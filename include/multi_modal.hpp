#pragma once

#include <type_traits>
#include <numeric>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>
#include <algorithm>
#include <tuple>
#include <map>

using std::pair;

const double sqrt2 = 1.414213562373095;
const double pi = 3.141592653589793;
const double sqrtpi = 1.772453850905516;
const double eps = 1e-10;


template<typename X> struct norm {
  double operator()(X const & a) const {
    return std::abs(a);
  }
};
template<typename N> struct norm<std::vector<N>> {
  double operator()(std::vector<N> const & a) const {
    std::vector<N> sq = a;
    std::transform(a.begin(), a.end(), sq.begin(), [](N const & n) { return n * n; });
    return std::sqrt(std::accumulate(sq.begin(), sq.end(), 0.0));
  }
};
template<typename X, typename T> X scale(X const & x, T const & factor) {
    return x * factor;
}
template<typename U, typename T> std::vector<U> scale(std::vector<U> const & x, T const & factor) {
  std::vector<U> ret = x;
  std::transform(x.begin(), x.end(), ret.begin(), [&factor](auto const & v) {
    return v * factor;
  });
  return ret;
};

namespace std {
  template<typename V> struct plus<std::vector<V>> {
    std::vector<V> operator()(std::vector<V> const & a, std::vector<V> const & b) const {
      std::vector<V> ret = a;
      std::transform(a.begin(), a.end(), b.begin(), ret.begin(), std::plus<>());
      return ret;
    }
  };
  template<typename V> struct minus<std::vector<V>> {
    std::vector<V> operator()(std::vector<V> const & a, std::vector<V> const & b) const {
      std::vector<V> ret = a;
      std::transform(a.begin(), a.end(), b.begin(), ret.begin(), std::minus<>());
      return ret;
    }
  };
}

template<typename X> struct distribution {
  X mean;
  double m2;
  unsigned long count;

  distribution(X const & p)
    : mean(p), m2(0.), count(1)
  { }
  distribution()
    : mean(), m2(0.), count(0)
  { }

  static distribution<X> from_standard_deviation(X const & mean, double standard_deviation, unsigned long count) {
    distribution<X> ret(mean);
    ret.m2 = standard_deviation * standard_deviation * count;
    ret.count = count;
    return ret;
  }

  double variance() const {
    return m2 / (double)count;
  }
  double standard_deviation() const {
    return std::sqrt(variance());
  }
  double likelihood(X const & x) const {
    static std::minus<X> minus;
    static norm<X> norm;

    if (count == 1) return 0.;
    return 1. - std::erf(norm(minus(mean, x)) / standard_deviation() / sqrt2);
  }

  double density(double from_mean) const {
    if (m2 > 0)
      return std::exp(-from_mean * from_mean / variance() / 2.) / sqrt2 / sqrtpi / standard_deviation();
    
    if (from_mean == 0.) return 1.0;
    return 0.0;
  }
};

// template<typename T> struct distribution<std::vector<T>> {
//   std::vector<T> mean;
//   std::vector<double> m2;
//   unsigned long count;

//   distribution(std::vector<T> const & p) 
//     : mean(p), m2(p.size()), count(1) 
//   { 
//     std::fill(m2.begin(), m2.end(), 0.);
//   }

//   static distribution<std::vector<T>> from_standard_deviation(std::vector<T> const & mean, std::vector<double> const & standard_deviation, unsigned long count) {
//     std::vector<double> m2(standard_deviation);
//     std::transform(standard_deviation.begin(), standard_deviation.end(), m2.begin(), [count](auto const & stddev) -> double {
//       return stddev * stddev * (double)count;
//     });

//     return distribution<std::vector<T>>{mean, m2, count};
//   }

//   std::vector<double> variance() const {
//     std::vector<double> variance(m2);
//     std::transform(m2.begin(), m2.end(), variance.begin(), [&](auto const & m) {
//       return m / (double)count;
//     });
//     return variance;
//   }

//   double variance_determinant() const {
//     std::vector<double> var = variance();
//     return std::accumulate(var.begin(), var.end(), 1., std::multiplies<double>{});
//   }

//   double standard_deviation() const {
//     return std::sqrt(variance_determinant());
//   }

//   // double likelihood(std::vector<T> const & x) {
//   //   std::vector<double> variance(m2);
//   //   std::transform(m2.begin(), m2.end(), variance.begin(), [count](auto const & m) {
//   //     return m / (double)count;
//   //   });

//   //   std::vector<T> diff(mean);
//   //   std::transform(mean.begin(), mean.end(), x.begin(), diff.begin(), std::minus<double>{});

//   //   std::vector<double> scaled(m2);
//   //   std::transform(diff.begin(), diff.end(), variance.begin(), scaled.begin(), [count](auto const & diff, double v) {
//   //     return diff * diff / v;
//   //   });
//   //   double exp = std::accumulate(scaled.begin(), scaled.end(), 0.);

//   //   double det = std::accumulate(variance.begin(), variance.end(), 1., std::multiplies<double>{});

//   //   double g = std::exp(-exp / 2.) / 
//   // }

//   double density(std::vector<T> const & x) const {
//     auto var = variance();
//     auto det = variance_determinant();

//     if (det < eps)
//       return 0.;

//     std::vector<T> val(mean);
//     std::transform(mean.begin(), mean.end(), x.begin(), val.begin(), std::minus<double>{});
//     std::transform(val.begin(), val.end(), var.begin(), val.begin(), [](auto const & v, auto const & s) -> double {
//       return v * v / s;
//     });
//     double exponent = std::accumulate(val.begin(), val.end(), 0., std::plus<double>{});

//     return std::exp(-exponent / 2.) / std::pow(sqrtpi * sqrt2, x.size()) / std::sqrt(det);
//   }

//   std::vector<T> zero() {
//     std::vector<T> zed(mean.size());
//     std::fill(zed.begin(), zed.end(), 0.);
//     return zed;
//   }
// };

template<typename T>
std::ostream& operator<<(std::ostream & os, distribution<T> const & dist) {
  return os << "{ " << dist.mean << " / " << dist.standard_deviation() << " # " << dist.count << " }";
}

// L2 norm
template<typename X> 
double mixture_error(distribution<X> const & a, distribution<X> const & b) {
  static norm<X> norm;
  static std::minus<X> minus;

  distribution<X> c = mix(a, b);

  double mua = norm(minus(a.mean, c.mean));
  double mub = norm(minus(b.mean, c.mean));

  double alpha = (double)a.count / (double)c.count;

  double vara = a.variance();
  double varb = b.variance();
  double stda = a.standard_deviation();
  double stdb = b.standard_deviation();

  double ret = 0.;

  if (vara == 0. || varb == 0.) {
    return std::numeric_limits<double>::max();
  }

  ret += alpha * alpha / stda +
         (1. - alpha) * (1. - alpha) / stdb +
         1. / c.standard_deviation();
  ret /= sqrt2;

  ret +=   2. * alpha * (1. - alpha) / std::sqrt(vara + varb)
         - 2. * alpha / std::sqrt(vara + c.variance())
         - 2. * (1. - alpha) / std::sqrt(varb + c.variance());
  
  ret /= sqrt2 / sqrtpi;

  return ret;
}

// total variational distance sup | (A(+)B) - C | 
template<typename X>
double mixture_error2(distribution<X> const & a, distribution<X> const & b) {
  static norm<X> norm;
  static std::minus<X> minus;

  distribution<X> c = mix(a, b);

  double mac = norm(minus(a.mean, c.mean));
  double mbc = norm(minus(b.mean, c.mean));
  double mab = norm(minus(a.mean, b.mean));
  double alpha = (double)a.count / (double)c.count;

  double ea = std::abs(alpha * a.density(0)   + (1. - alpha) * b.density(mab) - c.density(mac));
  double eb = std::abs(alpha * a.density(mab) + (1. - alpha) * b.density(0)   - c.density(mbc));
  double ec = std::abs(alpha * a.density(mac) + (1. - alpha) * b.density(mbc) - c.density(0));

  // std::cout << "density: " << a.density(0) << "\n";

  // std::cout << "ac, bc, ab, alpha, ea, eb, ec: " << mac << " " << mbc << " " << mab << " " << alpha << " " << ea << " " << eb << " " << ec << " " << std::endl;

  return ea * a.standard_deviation() + eb * b.standard_deviation() + ec * c.standard_deviation();
}

// template<typename T>
// double mixture_error(distribution<std::vector<T>> const & a, distribution<std::vector<T>> const & b) {
//   distribution<std::vector<T>> c = mix(a, b);

//   double alpha = (double)a.count / (double)c.count;

//   double scale = std::pow(sqrtpi * sqrt2, a.mean.size());
//   double astd = std::sqrt(a.variance_determinant()),
//          bstd = std::sqrt(b.variance_determinant()),
//          cstd = std::sqrt(c.variance_determinant());

//   double ea = std::abs(alpha / scale / astd + (1. - alpha) * b.density(a.mean) - c.density(a.mean)),
//          eb = std::abs(alpha * a.density(b.mean) + (1. - alpha) / scale / bstd - c.density(b.mean)),
//          ec = std::abs(alpha * a.density(c.mean) + (1. - alpha) * b.density(c.mean) - 1. / scale / cstd);

//   return ea * astd + eb * bstd + ec * cstd;
// }


template<typename X>
distribution<X> mix(distribution<X> const & a, distribution<X> const & b) {
  static norm<X> norm;
  static std::plus<X> plus;
  static std::minus<X> minus;

  unsigned long count = a.count + b.count;

  // weight factor
  double a_left = (double)a.count / (double)count;
  double a_right = 1. - a_left;

  X mean = plus(
    scale(a.mean, a_left),
    scale(b.mean, a_right)
  );

  // distance between means
  double left_distance  = norm(minus(mean, a.mean));
  double right_distance = norm(minus(mean, b.mean));

  double variance = a_left  * (a.variance() + left_distance  * left_distance) +
                    a_right * (b.variance() + right_distance * right_distance);

  distribution<X> ret(mean);
  ret.m2 = variance * (double)count;
  ret.count = count;
  return ret;
}

// template<typename T>
// distribution<std::vector<T>> mix(distribution<std::vector<T>> const & a, distribution<std::vector<T>> const & b) {
//   static norm<std::vector<T>> norm;
//   static std::plus<std::vector<T>> plus;
//   static std::minus<std::vector<T>> minus;

//   unsigned long count = a.count + b.count;
//   double alpha = (double)a.count / (double)count;

//   std::vector<T> mean = plus(
//     scale(a.mean, alpha),
//     scale(b.mean, 1. - alpha)
//   );

//   std::vector<double> m2(a.m2);

//   for (int i = 0; i < m2.size(); i++) {
//     double a_dist = a.mean[i] - mean[i];
//     double b_dist = b.mean[i] - mean[i];
//     m2[i] = a.m2[i] + (double)a.count * a_dist * a_dist + 
//             b.m2[i] + (double)b.count * b_dist * b_dist;
//   }

//   distribution<std::vector<T>> ret(mean);
//   ret.m2 = m2;
//   ret.count = count;
//   return ret;
// }

template<typename X>
distribution<X> unmix(distribution<X> const & c, distribution<X> const & b) {
  static norm<X> norm;
  static std::plus<X> plus;
  static std::minus<X> minus;

  double left_factor = (double)c.count / (double)(c.count - b.count);
  double right_factor = (double)b.count / (double)(c.count - b.count);

  X mean = plus(
    scale(c.mean, -left_factor),
    scale(b.mean, right_factor)
  );

  // distance between means
  double left_distance  = norm(minus(mean, c.mean));
  double right_distance = norm(minus(mean, b.mean));

  double variance = left_factor  * (c.variance() + left_distance  * left_distance) -
                    right_factor * (b.variance() + right_distance * right_distance);

  unsigned long count = c.count - b.count;

  distribution<X> ret(mean);
  ret.m2 = variance * (double)count;
  ret.count = count;
  return ret;
}


template<typename X> class multi_modal {
  struct node {
    distribution<X> dist;
    double error;
    node * left, * right;
    unsigned long id;
  };

  node * root;
  unsigned long maximum_nodes;
  unsigned long count;
  unsigned long next_id;
private:
  void insert_helper(node * n, distribution<X> const & dist) {
    static norm<X> norm;
    static std::minus<X> minus;

    if (n->left == nullptr) {
      if (count < maximum_nodes) {
        n->left = new node(*n);
        n->left->id = next_id++;
        n->right = new node{dist, 0, nullptr, nullptr, next_id++};
        n->dist = mix(n->right->dist, n->left->dist);
        // this could be more expensive
        n->error = mixture_error(n->right->dist, n->left->dist);

        count++;
      } else {
        n->dist = mix(n->dist, dist);
      }
      return;
    }

    node ** op, ** other, ** adjust;

    auto left = norm(minus(n->left->dist.mean, dist.mean));
    auto right = norm(minus(n->right->dist.mean, dist.mean));

    if (left < right) {
      op = &(n->left);
      other = &(n->right);
    } else {
      op = &(n->right);
      other = &(n->left);
    }

    insert_helper(*op, dist);
    n->dist = mix(n->left->dist, n->right->dist);
    // this could be more expensive
    n->error = mixture_error(n->right->dist, n->left->dist);

    if (n->dist.variance() < (*op)->dist.variance()) {
      // handle the case where *op is a leaf
      if ((*op)->left == nullptr) {

      } else {
        if ((*op)->left->dist.variance() < (*op)->right->dist.variance()) {
          adjust = &((*op)->left);
        } else {
          adjust = &((*op)->right);
        }

        std::swap(*other, *adjust);

        (*op)->dist = mix((*op)->right->dist, (*op)->left->dist);
        (*op)->error = mixture_error((*op)->right->dist, (*op)->left->dist);
        n->dist = mix(n->right->dist, n->left->dist);
        n->error = mixture_error(n->left->dist, n->right->dist);
      }
    }
  }

  void delete_helper(node * n) {
    std::vector<node*> stack;
    stack.push_back(n);
    while(!stack.empty()) {
      node * cur = stack.back();
      stack.pop_back();
      if (cur->left != nullptr) {
        stack.push_back(cur->left);
        cur->left = nullptr;
      }
      if (cur->right != nullptr) {
        stack.push_back(cur->right);
        cur->right = nullptr;
      }
      delete cur;
      count--;
    }
  }

  template<typename Visitor> void visit_nodes(Visitor v) const {
    std::vector<std::pair<node*,unsigned long>> stack;
    stack.push_back({root,0});
    node * cur;
    unsigned long depth;
    while(!stack.empty()) {
      std::tie(cur, depth) = stack.back();
      stack.pop_back();

      if (v(cur, depth)) {
        if (cur->left != nullptr) stack.push_back({cur->left, depth+1});
        if (cur->right != nullptr) stack.push_back({cur->right, depth+1});
      }
    }
  }

  bool extract_peaks_helper(std::vector<pair<unsigned long, distribution<X>>> & peaks, node * cur, node * parent) const {
    if (parent != nullptr && cur->error < parent->error) {
      peaks.push_back({cur->id, cur->dist});
      return true;
    }

    bool left = false, right = false;
    if (cur->left != nullptr) {
      left = extract_peaks_helper(peaks, cur->left, cur);
      right = extract_peaks_helper(peaks, cur->right, cur);

      // if (left && !right) {
      //   peaks.push_back(cur->right->dist);
      // } else if (!left && right) {
      //   peaks.push_back(cur->left->dist);
      // }
    }

    // add the root if we haven't added anything else
    if (parent == nullptr && peaks.size() == 0) {
      peaks.push_back({cur->id, cur->dist});
    }

    return left || right;
  }

  bool find_peak_helper(X const & x, distribution<X> & dist, node * n) const {
    static norm<X> norm;
    static std::minus<X> minus;

    bool found = false;
    if (n->left == nullptr) return false;

    auto left = norm(minus(n->left->dist.mean, dist.mean));
    auto right = norm(minus(n->right->dist.mean, dist.mean));

    node * chosen = n->left, * other = n->right;
    if(right < left) {
      chosen = n->right;
      other = n->left;
    }

    if (chosen->error < n->error) {
      dist = chosen->dist;
      return true;
    } 

    bool ret = find_peak_helper(x, dist, chosen);

    if (!ret && other->error < n->error) {
      dist = chosen->dist;
      return true;
    }

    return ret;
  }

public:
  std::vector<pair<unsigned long, distribution<X>>> extract_peaks() const {
    std::vector<pair<unsigned long, distribution<X>>> ret;
    extract_peaks_helper(ret, root, nullptr);
    return ret;
  }

  template<typename Visitor> void visit(Visitor v) const {
    visit_nodes([&v](node * n, unsigned long depth) -> bool {
      return v(n->dist, depth);
    });
  }
  template<typename Visitor> void visit_children(Visitor v) const {
    visit_nodes([&v](node * n, unsigned long depth) -> bool {
      if (n->left == nullptr) return false;

      return v(n->dist, n->left->dist, n->right->dist, depth);
    });
  }

  void insert(X const & x) {
    distribution<X> dist(x);

    if (root == nullptr) {
      root = new node{dist, 0, nullptr, nullptr};
      return;
    }
    insert_helper(root, dist);
  }

  distribution<X> find_peak(X const & x) {
    distribution<X> ret(x);

    if (root == nullptr) {
      return ret;
    }

    find_peak_helper(ret, root, 0.);
    return ret;
  }

  unsigned long get_count() const { return count; }

  multi_modal(unsigned long max) 
    : root(nullptr), maximum_nodes(max), count(0), next_id(0)
  {}

  multi_modal() 
    : multi_modal(std::numeric_limits<unsigned long>::max()) 
  {}

  ~multi_modal() { 
    if (root == nullptr) return;

    delete_helper(root); 
    root = nullptr;
  }
};