#include <algorithm>
#include <cmath>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <unordered_map>
#include <utility>
#include <vector>

namespace py = pybind11;
using namespace py::literals;

double min(double x, double y) { return x < y ? x : y; }

namespace Type1 {
// 只算bin数目

struct Shelf {
  double w, h;
};

bool try_slot(std::unordered_map<double, std::vector<double>> &slots, double w,
              double h) {
  for (auto &s : slots[w]) {
    if (h <= s) {
      s -= h;
      return true;
    }
  }
  return false;
}

bool try_shelf(std::vector<Shelf> &shelves,
               std::unordered_map<double, std::vector<double>> &slots,
               double min_l, double w, double h) {
  std::pair<double, double> best_m = {-INFINITY, -INFINITY};
  Shelf *best_s = nullptr;
  for (auto &s : shelves) {
    if (h <= s.h) {
      auto w_ = s.w - w;
      if (w_ >= 0) {
        std::pair<double, double> m = {w_ * s.h, w_};
        if (m > best_m) {
          best_m = m;
          best_s = &s;
        }
      }
    }
  }
  if (!best_s) {
    return false;
  }
  best_s->w = best_m.second;
  auto h_ = best_s->h - h;
  if (h_ >= min_l && w >= min_l) {
    slots[w].push_back(h_);
  }
  return true;
}

int get_bin_out_fast(
    py::array_t<double, py::array::c_style | py::array::forcecast> _items,
    bool rotatable, double bin_width, double bin_height) {
  auto items = _items.unchecked<2>();
  auto n = items.shape(0);
  auto last = INFINITY;
  std::vector<double> min_ls(n);
  for (auto i = n - 1; i >= 0; --i) {
    min_ls.at(i) = last;
    last = min(last, min(items(i, 0), items(i, 1)));
  }
  std::vector<double> bin_hs;
  std::vector<Shelf> shelves;
  std::unordered_map<double, std::vector<double>> slots;
  for (int i = 0; i < n; ++i) {
    auto w = items(i, 0);
    auto h = items(i, 1);
    auto min_l = min_ls.at(i);
    if (try_slot(slots, w, h) || (rotatable && try_slot(slots, h, w))) {
      continue;
    }
    if (try_shelf(shelves, slots, min_l, w, h) ||
        (rotatable && try_shelf(shelves, slots, min_l, h, w))) {
      continue;
    }
    bool found = false;
    for (auto &bin_h : bin_hs) {
      if (h <= bin_h) {
        bin_h -= h;
        found = true;
        break;
      }
    }
    if (!found) {
      bin_hs.push_back(bin_height - h);
    }
    shelves.push_back({.w = bin_width - w, .h = h});
  }
  return bin_hs.size();
}
}; // namespace Type1

namespace Type2 {
// 额外记录每个item放入了哪个bin

struct Shelf {
  int id;
  double w, h;
};

struct Slot {
  int id;
  double h;
};

bool try_slot(std::vector<int> &item2bins,
              std::unordered_map<double, std::vector<Slot>> &slots, double w,
              double h) {
  for (auto &s : slots[w]) {
    if (h <= s.h) {
      s.h -= h;
      item2bins.push_back(s.id);
      return true;
    }
  }
  return false;
}

bool try_shelf(std::vector<int> &item2bins, std::vector<Shelf> &shelves,
               std::unordered_map<double, std::vector<Slot>> &slots,
               double min_l, double w, double h) {
  std::pair<double, double> best_m = {-INFINITY, -INFINITY};
  Shelf *best_s = nullptr;
  for (auto &s : shelves) {
    if (h <= s.h) {
      auto w_ = s.w - w;
      if (w_ >= 0) {
        std::pair<double, double> m = {w_ * s.h, w_};
        if (m > best_m) {
          best_m = m;
          best_s = &s;
        }
      }
    }
  }
  if (!best_s) {
    return false;
  }
  item2bins.push_back(best_s->id);
  best_s->w = best_m.second;
  auto h_ = best_s->h - h;
  if (h_ >= min_l && w >= min_l) {
    slots[w].push_back({.id = best_s->id, .h = h_});
  }
  return true;
}

std::pair<int, double> get_bin_out_fast(
    py::array_t<double, py::array::c_style | py::array::forcecast> _items,
    bool rotatable, double bin_width, double bin_height) {
  auto items = _items.unchecked<2>();
  auto n = items.shape(0);
  auto last = INFINITY;
  std::vector<double> min_ls(n);
  for (auto i = n - 1; i >= 0; --i) {
    min_ls.at(i) = last;
    last = min(last, min(items(i, 0), items(i, 1)));
  }
  std::vector<double> bin_hs;
  std::vector<Shelf> shelves;
  std::unordered_map<double, std::vector<Slot>> slots;
  std::vector<int> item2bins;
  for (int i = 0; i < n; ++i) {
    auto w = items(i, 0);
    auto h = items(i, 1);
    auto min_l = min_ls.at(i);
    if (try_slot(item2bins, slots, w, h) ||
        (rotatable && try_slot(item2bins, slots, h, w))) {
      continue;
    }
    if (try_shelf(item2bins, shelves, slots, min_l, w, h) ||
        (rotatable && try_shelf(item2bins, shelves, slots, min_l, h, w))) {
      continue;
    }
    int b_id = 0, size = bin_hs.size();
    for (; b_id < size; ++b_id) {
      if (h <= bin_hs.at(b_id)) {
        bin_hs.at(b_id) -= h;
        break;
      }
    }
    if (b_id == size) {
      bin_hs.push_back(bin_height - h);
    }
    shelves.push_back({.id = b_id, .w = bin_width - w, .h = h});
    item2bins.push_back(b_id);
  }
  auto bin_id = bin_hs.size() - 1;
  double area = 0;
  for (int i = 0; i < n; ++i) {
    if (item2bins.at(i) == bin_id) {
      area += items(i, 0) * items(i, 1);
    }
  }
  return {bin_id + 1, area / (bin_height * bin_width)};
}
std::pair<std::vector<int>, std::vector<double>> get_bin_out_fast_with_result(
    py::array_t<double, py::array::c_style | py::array::forcecast> _items,
    bool rotatable, double bin_width, double bin_height) {
  auto items = _items.unchecked<2>();
  auto n = items.shape(0);
  auto last = INFINITY;
  std::vector<double> min_ls(n);
  for (auto i = n - 1; i >= 0; --i) {
    min_ls.at(i) = last;
    last = min(last, min(items(i, 0), items(i, 1)));
  }
  std::vector<double> bin_hs;
  std::vector<Shelf> shelves;
  std::unordered_map<double, std::vector<Slot>> slots;
  std::vector<int> item2bins;
  item2bins.reserve(n);
  for (int i = 0; i < n; ++i) {
    auto w = items(i, 0);
    auto h = items(i, 1);
    auto min_l = min_ls.at(i);
    if (try_slot(item2bins, slots, w, h) ||
        (rotatable && try_slot(item2bins, slots, h, w))) {
      continue;
    }
    if (try_shelf(item2bins, shelves, slots, min_l, w, h) ||
        (rotatable && try_shelf(item2bins, shelves, slots, min_l, h, w))) {
      continue;
    }
    int b_id = 0, size = bin_hs.size();
    for (; b_id < size; ++b_id) {
      if (h <= bin_hs.at(b_id)) {
        bin_hs.at(b_id) -= h;
        break;
      }
    }
    if (b_id == size) {
      bin_hs.push_back(bin_height - h);
    }
    shelves.push_back({.id = b_id, .w = bin_width - w, .h = h});
    item2bins.push_back(b_id);
  }
  std::vector<double> usage(bin_hs.size(), 0);
  for (int i = 0; i < n; ++i) {
    usage.at(item2bins.at(i)) += items(i, 0) * items(i, 1);
  }
  auto a = bin_height * bin_width;
  for (auto &i : usage) {
    i /= a;
  }
  return {item2bins, usage};
}
}; // namespace Type2

namespace Type3 {
// 使用新的启发式方法选bin

struct Shelf {
  int id;
  double w, h;
};

struct Slot {
  int id;
  double h;
};

bool try_slot(std::vector<int> &item2bins,
              std::unordered_map<double, std::vector<Slot>> &slots, double w,
              double h) {
  for (auto &s : slots[w]) {
    if (h <= s.h) {
      s.h -= h;
      item2bins.push_back(s.id);
      return true;
    }
  }
  return false;
}

bool try_shelf(std::vector<int> &item2bins, std::vector<Shelf> &shelves,
               std::unordered_map<double, std::vector<Slot>> &slots,
               double min_l, double w, double h) {
  std::pair<double, double> best_m = {INFINITY, INFINITY};
  Shelf *best_s = nullptr;
  for (auto &s : shelves) {
    if (h <= s.h) {
      auto w_ = s.w - w;
      if (w_ >= 0) {
        std::pair<double, double> m = {w_, s.h - h};
        if (m < best_m) {
          best_m = m;
          best_s = &s;
        }
      }
    }
  }
  if (!best_s) {
    return false;
  }
  item2bins.push_back(best_s->id);
  best_s->w = best_m.first;
  auto h_ = best_m.second;
  if (h_ >= min_l && w >= min_l) {
    slots[w].push_back({.id = best_s->id, .h = h_});
  }
  return true;
}

std::pair<int, double> get_bin_out_fast(
    py::array_t<double, py::array::c_style | py::array::forcecast> _items,
    bool rotatable, double bin_width, double bin_height) {
  auto items = _items.unchecked<2>();
  auto n = items.shape(0);
  auto last = INFINITY;
  std::vector<double> min_ls(n);
  for (auto i = n - 1; i >= 0; --i) {
    min_ls.at(i) = last;
    last = min(last, min(items(i, 0), items(i, 1)));
  }
  std::vector<double> bin_hs;
  std::vector<Shelf> shelves;
  std::unordered_map<double, std::vector<Slot>> slots;
  std::vector<int> item2bins;
  for (int i = 0; i < n; ++i) {
    auto w = items(i, 0);
    auto h = items(i, 1);
    auto min_l = min_ls.at(i);
    if (try_slot(item2bins, slots, w, h) ||
        (rotatable && try_slot(item2bins, slots, h, w))) {
      continue;
    }
    if (try_shelf(item2bins, shelves, slots, min_l, w, h) ||
        (rotatable && try_shelf(item2bins, shelves, slots, min_l, h, w))) {
      continue;
    }
    int best_id = -1;
    double best_h = bin_height;
    for (int i = 0, j = bin_hs.size(); i < j; ++i) {
      auto h_ = bin_hs.at(i) - h;
      if (h_ >= 0 && h_ < best_h) {
        best_h = h_;
        best_id = i;
      }
    }
    if (best_id == -1) {
      best_id = bin_hs.size();
      bin_hs.push_back(bin_height - h);
    } else {
      bin_hs.at(best_id) = best_h;
    }
    shelves.push_back({.id = best_id, .w = bin_width - w, .h = h});
    item2bins.push_back(best_id);
  }
  auto bin_id = bin_hs.size() - 1;
  double area = 0;
  for (int i = 0; i < n; ++i) {
    if (item2bins.at(i) == bin_id) {
      area += items(i, 0) * items(i, 1);
    }
  }
  return {bin_id + 1, area / (bin_height * bin_width)};
}
std::pair<std::vector<int>, std::vector<double>> get_bin_out_fast_with_result(
    py::array_t<double, py::array::c_style | py::array::forcecast> _items,
    bool rotatable, double bin_width, double bin_height) {
  auto items = _items.unchecked<2>();
  auto n = items.shape(0);
  auto last = INFINITY;
  std::vector<double> min_ls(n);
  for (auto i = n - 1; i >= 0; --i) {
    min_ls.at(i) = last;
    last = min(last, min(items(i, 0), items(i, 1)));
  }
  std::vector<double> bin_hs;
  std::vector<Shelf> shelves;
  std::unordered_map<double, std::vector<Slot>> slots;
  std::vector<int> item2bins;
  item2bins.reserve(n);
  for (int i = 0; i < n; ++i) {
    auto w = items(i, 0);
    auto h = items(i, 1);
    auto min_l = min_ls.at(i);
    if (try_slot(item2bins, slots, w, h) ||
        (rotatable && try_slot(item2bins, slots, h, w))) {
      continue;
    }
    if (try_shelf(item2bins, shelves, slots, min_l, w, h) ||
        (rotatable && try_shelf(item2bins, shelves, slots, min_l, h, w))) {
      continue;
    }
    int best_id = -1;
    double best_h = bin_height;
    for (int i = 0, j = bin_hs.size(); i < j; ++i) {
      auto h_ = bin_hs.at(i) - h;
      if (h_ >= 0 && h_ < best_h) {
        best_h = h_;
        best_id = i;
      }
    }
    if (best_id == -1) {
      best_id = bin_hs.size();
      bin_hs.push_back(bin_height - h);
    } else {
      bin_hs.at(best_id) = best_h;
    }
    shelves.push_back({.id = best_id, .w = bin_width - w, .h = h});
    item2bins.push_back(best_id);
  }
  std::vector<double> usage(bin_hs.size(), 0);
  for (int i = 0; i < n; ++i) {
    usage.at(item2bins.at(i)) += items(i, 0) * items(i, 1);
  }
  auto a = bin_height * bin_width;
  for (auto &i : usage) {
    i /= a;
  }
  return {item2bins, usage};
}
}; // namespace Type3

namespace bp_mr {
struct Rect {
  double x, y, w, h;
};

struct Bin {
  int id;
  std::vector<Rect> rs;
};

using Score = std::pair<double, double>;

std::pair<std::vector<int>, std::vector<double>> bp_max_rect(
    py::array_t<double, py::array::c_style | py::array::forcecast> _items,
    bool rotatable, double bin_width, double bin_height) {
  auto items = _items.unchecked<2>();
  auto n = items.shape(0);
  std::vector<int> item2bins(n);
  std::vector<Bin> bins;
  std::vector<Rect> rs1, rs2;
  int bin_cnt = 0;
  for (int i = 0; i < n; ++i) {
    auto w = items(i, 0);
    auto h = items(i, 1);
    Bin *best_bin = nullptr;
    Rect *best_rect = nullptr;
    Score best_s = {INFINITY, INFINITY};
    for (auto &b : bins) {
      for (auto &r : b.rs) {
        if (r.w >= w && r.h >= h) {
          Score s = {r.w - w, h - r.h};
          if (s < best_s) {
            best_s = s;
            best_rect = &r;
            best_bin = &b;
          }
        }
      }
    }
    if (best_rect) {
      auto x = best_rect->x, y = best_rect->y;
      item2bins.at(i) = best_bin->id;
      if (w == best_rect->w && h == best_rect->h) {
        *best_rect = best_bin->rs.back();
        best_bin->rs.pop_back();
      }
      auto r = x + w, t = y + h;
      for (auto &rr : best_bin->rs) {
        auto x_ = rr.x, y_ = rr.y, w_ = rr.w, h_ = rr.h;
        auto r_ = x_ + w_, t_ = y_ + h_;
        if (std::max(r, r_) - std::min(x, x_) >= w + w_ ||
            std::max(t, t_) - std::min(y, y_) >= h + h_) {
          rs1.push_back(rr);
          continue;
        }
        if (x > x_)
          rs1.push_back({x_, y_, x - x_, h_});
        if (y > y_)
          rs1.push_back({x_, y_, w_, y - y_});
        if (r < r_)
          rs1.push_back({r, y_, r_ - r, h_});
        if (t < t_)
          rs1.push_back({x_, t, w_, t_ - t});
      }
      std::sort(rs1.begin(), rs1.end(), [](const Rect &a, const Rect &b) {
        return a.w * a.h > b.w * b.h;
      });
      for (auto &i : rs1) {
        bool covered = false;
        for (auto &j : rs2) {
          if (i.x >= j.x && i.y >= j.y && i.x + i.w <= j.x + j.w &&
              i.y + i.h <= j.y + j.h) {
            covered = true;
            break;
          }
        }
        if (!covered) {
          rs2.push_back(i);
        }
      }
      best_bin->rs.swap(rs2);
      rs1.clear();
      rs2.clear();
    } else {
      item2bins.at(i) = bin_cnt;
      Bin b;
      b.id = bin_cnt;
      if (w < bin_width) {
        b.rs.push_back({w, 0, bin_width - w, bin_height});
      }
      if (h < bin_height) {
        b.rs.push_back({0, h, bin_width, bin_height - h});
      }
      if (!b.rs.empty()) {
        bins.push_back(std::move(b));
      }
      ++bin_cnt;
    }
  }
  std::vector<double> usage(bins.size(), 0);
  for (int i = 0; i < n; ++i) {
    usage.at(item2bins.at(i)) += items(i, 0) * items(i, 1);
  }
  auto a = bin_height * bin_width;
  for (auto &i : usage) {
    i /= a;
  }
  return {item2bins, usage};
}
} // namespace bp_mr

PYBIND11_MODULE(pybp, m) {
  m.doc() = R"pbdoc(
    2D Bin Packing
    -----------------------

    .. currentmodule:: pybp

    .. autosummary::
        :toctree: _generate

        get_bin_out_fast_1
        get_bin_out_fast_2
    )pbdoc";
  m.def("get_bin_out_fast_1", Type1::get_bin_out_fast,
        R"pbdoc(
        Do bin-packing algorithm on the given input

        items: float np.array, Nx2

        rotatable: bool, whether the items are rotatable

        bin_width: float, width of the bins

        bin_height: float, height of the bins

        return: # of used bins
        )pbdoc",
        "items"_a, "rotatable"_a, "bin_width"_a = 2440, "bin_height"_a = 1220);
  m.def("get_bin_out_fast_2", Type2::get_bin_out_fast,
        R"pbdoc(
        Do bin-packing algorithm on the given input

        items: float np.array, Nx2

        rotatable: bool, whether the items are rotatable

        bin_width: float, width of the bins

        bin_height: float, height of the bins

        return: (# of used bins, % of usage of the last bin)
        )pbdoc",
        "items"_a, "rotatable"_a, "bin_width"_a = 2440, "bin_height"_a = 1220);
  m.def("get_bin_out_fast_3", Type2::get_bin_out_fast_with_result,
        R"pbdoc(
        Do bin-packing algorithm on the given input

        items: float np.array, Nx2

        rotatable: bool, whether the items are rotatable

        bin_width: float, width of the bins

        bin_height: float, height of the bins

        return: (bin_id of each item, usage of each bin)
        )pbdoc",
        "items"_a, "rotatable"_a, "bin_width"_a = 2440, "bin_height"_a = 1220);
  m.def("get_bin_out_fast_4", Type3::get_bin_out_fast_with_result,
        R"pbdoc(
        Do bin-packing algorithm on the given input

        items: float np.array, Nx2

        rotatable: bool, whether the items are rotatable

        bin_width: float, width of the bins

        bin_height: float, height of the bins

        return: (bin_id of each item, usage of each bin)
        )pbdoc",
        "items"_a, "rotatable"_a, "bin_width"_a = 2440, "bin_height"_a = 1220);
  m.def("bp_max_rect", bp_mr::bp_max_rect,
        R"pbdoc(
        Do MaxRect bin-packing algorithm on the given input

        items: float np.array, Nx2

        rotatable: bool, whether the items are rotatable (IGNORED)

        bin_width: float, width of the bins

        bin_height: float, height of the bins

        return: (bin_id of each item, usage of each bin)
        )pbdoc",
        "items"_a, "rotatable"_a = false, "bin_width"_a = 2440,
        "bin_height"_a = 1220);
#ifdef VERSION_INFO
#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)
  m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
  m.attr("__version__") = "dev";
#endif
}
