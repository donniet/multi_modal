#pragma once

#include <type_traits>
#include <numeric>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

const double sqrt2 = 1.414213562373095;
const double pi = 3.141592653589793;
const double sqrtpi = 1.772453850905516;


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
  size_t count;

  static distribution<X> from_standard_deviation(X const & mean, double standard_deviation, size_t count) {
    return distribution<X>{mean, standard_deviation * standard_deviation * count, count};
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
    
    if (from_mean == mean) return 1.0;
    return 0.0;
  }
};

template<typename T>
std::ostream& operator<<(std::ostream & os, distribution<T> const & dist) {
  return os << "{ " << dist.mean << " / " << dist.standard_deviation() << " # " << dist.count << " }";
}

// double mixture_constant(double mux, double varx, double muy, double vary) {
//   double ret = 0.;
//   ret = (mux * vary + muy * varx) / (varx + vary) / 2. / varx / vary;
//   ret = (ret * ret);
//   ret += (mux * mux * vary + muy * muy * varx) / 2. / varx / vary;
//   return std::exp(-ret);
// }

double mixture_constant(double mux, double varx, double muy, double vary) {
  double ret = mux - muy;
  ret *= ret;

  std::cout << "\nmu2: " << ret << "\n";

  ret *= varx * vary;
  ret /= (varx + vary);
  return std::exp(ret);
}

// L2 norm
template<typename X> 
double mixture_error2(distribution<X> const & a, distribution<X> const & b) {
  static norm<X> norm;
  static std::minus<X> minus;

  distribution<X> c = mix(a, b);

  double mua = norm(minus(a.mean, c.mean));
  double mub = norm(minus(b.mean, c.mean));

  // std::cout << "\nmua: " << mua << " mub: " << mub << "\n";  

  double Cab = mixture_constant(mua, a.variance(), mub, b.variance());
  double Cac = mixture_constant(mua, a.variance(), 0, c.variance());
  double Cbc = mixture_constant(mub, b.variance(), 0, c.variance());

  // std::cout << "Cab " << Cab << " Cac " << Cac << " Cbc " << Cbc << "\n";

  double alpha = (double)a.count / (double)c.count;

  double ret = 0.;
  ret += alpha * alpha * a.standard_deviation();
  ret += (1. - alpha) * (1. - alpha) * b.standard_deviation();
  ret += c.standard_deviation();
  ret = ret / 2. / sqrtpi;

  double ret2 = 0.;
  ret2 += Cab * 2. * alpha * (1. - alpha) * a.standard_deviation() * b.standard_deviation() / std::sqrt(a.variance() + b.variance());
  ret2 -= Cac * 2 * alpha * a.standard_deviation() * c.standard_deviation() / std::sqrt(a.variance() + c.variance());
  ret2 -= Cbc * 2 * (1. - alpha) * b.standard_deviation() * c.standard_deviation() / std::sqrt(b.variance() + c.variance());
  ret2 = ret2 / sqrt2 / sqrtpi;

  return ret + ret2;
}

// total variational distance sup | (A(+)B) - C | 
template<typename X>
double mixture_error(distribution<X> const & a, distribution<X> const & b) {
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


template<typename X>
distribution<X> mix(distribution<X> const & a, distribution<X> const & b) {
  static norm<X> norm;
  static std::plus<X> plus;
  static std::minus<X> minus;

  size_t count = a.count + b.count;

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

  return distribution<X>{mean, variance * (double)count, count};
}

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

  size_t count = c.count - b.count;

  return distribution<X>{mean, variance * (double)count, count};
}


template<typename X> class multi_modal {
  struct node {
    distribution<X> dist;
    double error;
    node * left, * right;
  };

  node * root;
  size_t maximum_nodes;
  size_t count;
private:
  void insert_helper(node * n, distribution<X> const & dist) {
    static norm<X> norm;
    static std::minus<X> minus;

    if (n->left == nullptr) {
      if (count < maximum_nodes) {
        n->left = new node(*n);
        n->right = new node{dist, 0, nullptr, nullptr};
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
      if (cur->left != nullptr) {
        stack.push_back(cur->left);
        cur->left = nullptr;
      } else if (cur->right != nullptr) {
        stack.push_back(cur->right);
        cur->right = nullptr;
      } else {
        stack.pop_back();
        delete cur;
        count--;
      }
    }
  }

  template<typename Visitor> void visit_nodes(Visitor v) const {
    std::vector<std::pair<node*,size_t>> stack;
    stack.push_back({root,0});
    node * cur;
    size_t depth;
    while(!stack.empty()) {
      std::tie(cur, depth) = stack.back();
      stack.pop_back();

      if (v(cur, depth)) {
        if (cur->left != nullptr) stack.push_back({cur->left, depth+1});
        if (cur->right != nullptr) stack.push_back({cur->right, depth+1});
      }
    }
  }

public:
  std::vector<distribution<X>> extract_peaks() const {
    std::vector<std::pair<node*,node*>> stack;
    stack.push_back({root,nullptr});
    node * cur, * parent;

    std::vector<distribution<X>> ret;

    while(!stack.empty()) {
      std::tie(cur, parent) = stack.back();
      stack.pop_back();

      if (cur->left == nullptr) continue;
      
      if (parent != nullptr && cur->error < parent->error) {
        ret.push_back(cur->dist);
        continue;
      } 

      stack.push_back({cur->left, cur});
      stack.push_back({cur->right, cur});
    }

    return ret;
  }

  template<typename Visitor> void visit(Visitor v) const {
    visit_nodes([&v](node * n, size_t depth) -> bool {
      return v(n->dist, depth);
    });
  }
  template<typename Visitor> void visit_children(Visitor v) const {
    visit_nodes([&v](node * n, size_t depth) -> bool {
      if (n->left == nullptr) return false;

      return v(n->dist, n->left->dist, n->right->dist, depth);
    });
  }

  void insert(X const & x) {
    distribution<X> dist{x, 0, 1};
    if (root == nullptr) {
      root = new node{dist, 0, nullptr, nullptr};
      return;
    }
    insert_helper(root, dist);
  }

  size_t get_count() const { return count; }

  multi_modal(size_t max) 
    : root(nullptr), maximum_nodes(max), count(0) 
  {}

  multi_modal() 
    : multi_modal(std::numeric_limits<size_t>::max()) 
  {}

  ~multi_modal() { 
    if (root == nullptr) return;

    delete_helper(root); 
    root = nullptr;
  }
};