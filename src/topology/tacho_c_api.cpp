/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

// Native skew factorization using Tacho symbolic analysis; no MATLAB API.
#include <Kokkos_Core.hpp>
#include <Tacho_Blas_External.hpp>
#include <Tacho_GraphTools_Metis.hpp>
#include <Tacho_SymbolicTools.hpp>
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/max_cardinality_matching.hpp>
#include <cmath>
#include <complex>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace {
// Numerical buffers are bounded separately from the sparse input and ordering.
struct PfaffianMemoryLimit : std::runtime_error {
  PfaffianMemoryLimit()
      : std::runtime_error("Pfaffian numerical memory limit exceeded") {}
};

// Real skew-specific multifrontal factorization with delayed 2x2 pivots.
// Symbolic analysis and BLAS are supplied by Tacho; no pivot perturbations.
class SkewLDL {
  struct Panel {
    std::vector<int> ids;
    std::vector<double> lower, pivots;
    int eliminated = 0;
  };
  struct Update {
    std::vector<int> ids;
    std::vector<double> values;
    int delayed = 0;
  };
  int order_;
  double tolerance_, threshold_;
  long long live_ = 0;
  const long long budget_;
  long long &peak_bytes_;
  std::vector<Panel> panels_;
  std::vector<int> elimination_;
  int sign_ = 1;

  void allocate(long long count) {
    if (count < 0 || count > budget_ - live_)
      throw PfaffianMemoryLimit();
    live_ += count;
    peak_bytes_ =
        std::max(peak_bytes_, live_ * static_cast<long long>(sizeof(double)));
  }
  void release(long long count) { live_ -= count; }

  static void exchange(std::vector<double> &a, std::vector<int> &ids, int n,
                       int i, int j) {
    if (i == j)
      return;
    for (int k = 0; k < n; ++k)
      std::swap(a[i + size_t(k) * n], a[j + size_t(k) * n]);
    for (int k = 0; k < n; ++k)
      std::swap(a[k + size_t(i) * n], a[k + size_t(j) * n]);
    std::swap(ids[i], ids[j]);
  }

public:
  long long panel_entries = 0;

  SkewLDL(int n, double tolerance, long long max_bytes, long long &peak_bytes,
          double threshold = 0.01)
      : order_(n), tolerance_(tolerance), threshold_(threshold),
        budget_(max_bytes / sizeof(double)), peak_bytes_(peak_bytes) {}

  void reserve_rhs(int nrhs) { allocate(2LL * order_ * nrhs); }

  bool factor(const std::vector<std::map<int, double>> &rows,
              const Tacho::SymbolicTools &symbolic) {
    const auto order = symbolic.PermVector();
    const auto bounds = symbolic.Supernodes();
    const auto parents = symbolic.SupernodesTreeParent();
    const auto ptr = symbolic.gidSuperPanelPtr();
    const auto indices = symbolic.gidSuperPanelColIdx();
    const int ns = symbolic.NumSupernodes();
    std::vector<int> inverse(order_);
    for (int i = 0; i < order_; ++i)
      inverse[order(i)] = i;
    std::vector<std::vector<int>> children(ns);
    for (int s = 0; s < ns; ++s) {
      if (parents(s) >= 0) {
        if (parents(s) <= s || parents(s) >= ns)
          throw std::runtime_error("Expected postordered Tacho supernodes");
        children[parents(s)].push_back(s);
      }
    }
    std::vector<Update> updates(ns);
    std::vector<int> local(order_, -1);
    for (int s = 0; s < ns; ++s) {
      std::set<int> eligible, border;
      for (int i = bounds(s); i < bounds(s + 1); ++i)
        eligible.insert(i);
      for (int child : children[s])
        for (int i = 0; i < updates[child].delayed; ++i)
          eligible.insert(updates[child].ids[i]);
      for (auto k = ptr(s); k < ptr(s + 1); ++k)
        if (!eligible.count(indices(k)))
          border.insert(indices(k));
      for (int child : children[s])
        for (int i : updates[child].ids)
          if (!eligible.count(i))
            border.insert(i);
      if (parents(s) < 0 && !border.empty())
        throw std::runtime_error("Root contains unowned separator variables");
      Panel panel;
      panel.ids.assign(eligible.begin(), eligible.end());
      panel.ids.insert(panel.ids.end(), border.begin(), border.end());
      const int n = panel.ids.size(), initial_eligible = eligible.size();
      const long long square = static_cast<long long>(n) * n;
      allocate(square + initial_eligible / 2);
      panel.pivots.reserve(initial_eligible / 2);
      std::vector<double> a(square, 0.0);
      for (int i = 0; i < n; ++i)
        local[panel.ids[i]] = i;
      // Own each original upper entry exactly once, at its earlier endpoint.
      for (int i = bounds(s); i < bounds(s + 1); ++i) {
        for (const auto &[column, value] : rows[order(i)]) {
          const int j = inverse[column];
          if (j <= i)
            continue;
          const int row = local[i], col = local[j];
          if (row < 0 || col < 0)
            throw std::runtime_error(
                "Original entry missing from symbolic front");
          a[row + size_t(col) * n] += value;
          a[col + size_t(row) * n] -= value;
        }
      }
      for (int child : children[s]) {
        auto &update = updates[child];
        const int nu = update.ids.size();
        for (int j = 0; j < nu; ++j)
          for (int i = 0; i < nu; ++i) {
            const int row = local[update.ids[i]], col = local[update.ids[j]];
            if (row < 0 || col < 0)
              throw std::runtime_error("Missing child variable");
            a[row + size_t(col) * n] += update.values[i + size_t(j) * nu];
          }
        release(update.values.size());
        update = Update{};
      }
      for (int i : panel.ids)
        local[i] = -1;
      constexpr int block_size = 64;
      const long long scratch = static_cast<long long>(n) * (block_size + 2);
      allocate(scratch);
      std::vector<double> original_columns(size_t(n) * block_size), first(n),
          second(n);
      int active = initial_eligible, k = 0, panel_start = 0;
      auto swap_pending = [&](int i, int j) {
        exchange(a, panel.ids, n, i, j);
        for (int t = 0; t < k - panel_start; ++t)
          std::swap(original_columns[i + size_t(t) * n],
                    original_columns[j + size_t(t) * n]);
      };
      auto updated_column = [&](int column, std::vector<double> &values) {
        for (int i = k; i < n; ++i)
          values[i] = a[i + size_t(column) * n];
        const int pending = k - panel_start;
        if (pending)
          Tacho::Blas<double>::gemv('N', n - k, pending, -1.0,
                                    original_columns.data() + k, n,
                                    a.data() + column + size_t(panel_start) * n,
                                    n, 1.0, values.data() + k, 1);
        values[column] = 0.0;
      };
      auto flush = [&]() {
        const int pending = k - panel_start, remaining = n - k;
        if (pending && remaining) {
          Tacho::Blas<double>::gemm('N', 'T', remaining, remaining, pending,
                                    -1.0, original_columns.data() + k, n,
                                    a.data() + k + size_t(panel_start) * n, n,
                                    1.0, a.data() + k + size_t(k) * n, n);
          for (int j = k; j < n; ++j) {
            a[j + size_t(j) * n] = 0.0;
            for (int i = j + 1; i < n; ++i)
              a[j + size_t(i) * n] = -a[i + size_t(j) * n];
          }
        }
        panel_start = k;
      };
      while (k + 1 < active) {
        updated_column(k, first);
        int partner = k + 1;
        for (int j = k + 2; j < active; ++j)
          if (std::abs(first[j]) > std::abs(first[partner]))
            partner = j;
        updated_column(partner, second);
        const double candidate = std::abs(first[partner]);
        double column_max = 0.0;
        for (int i = k; i < n; ++i) {
          const double x = first[i], y = second[i];
          if (!std::isfinite(x) || !std::isfinite(y))
            return false;
          column_max = std::max({column_max, std::abs(x), std::abs(y)});
        }
        // Include the separator in the threshold test. A good local pivot can
        // otherwise have arbitrarily large multipliers outside this supernode.
        if (candidate <= tolerance_ || candidate < threshold_ * column_max) {
          if (parents(s) < 0) {
            if (candidate <= tolerance_)
              return false;
            // Rook search at the root: no ancestor can accept a delayed row.
            swap_pending(k, partner);
            continue;
          }
          swap_pending(k, active - 1);
          --active;
          continue;
        }
        swap_pending(k + 1, partner);
        std::swap(first[k + 1], first[partner]);
        std::swap(second[k + 1], second[partner]);
        const double pivot = -first[k + 1];
        sign_ *= pivot > 0 ? 1 : -1;
        panel.pivots.push_back(pivot);
        const int start = k + 2, remaining = n - start;
        for (int i = 0; i < remaining; ++i) {
          const double x = first[start + i], y = second[start + i];
          original_columns[start + i + size_t(k - panel_start) * n] = x;
          original_columns[start + i + size_t(k - panel_start + 1) * n] = y;
          a[start + i + size_t(k) * n] = y / pivot;
          a[start + i + size_t(k + 1) * n] = -x / pivot;
        }
        a[k + size_t(k) * n] = a[k + 1 + size_t(k + 1) * n] = 1.0;
        a[k + 1 + size_t(k) * n] = a[k + size_t(k + 1) * n] = 0.0;
        k += 2;
        if (k - panel_start == block_size)
          flush();
      }
      flush();
      panel.eliminated = k;
      const int rest = n - k, delayed = initial_eligible - k;
      if (parents(s) < 0 && rest)
        return false;
      for (int i = 0; i < k; ++i)
        elimination_.push_back(order(panel.ids[i]));
      auto &update = updates[s];
      update.ids.assign(panel.ids.begin() + k, panel.ids.end());
      update.delayed = delayed;
      if (k == 0) {
        release(initial_eligible / 2);
        update.values = std::move(a);
      } else if (rest == 0) {
        panel.lower = std::move(a);
        panel_entries += square;
        panels_.push_back(std::move(panel));
      } else {
        const long long lower_size = static_cast<long long>(n) * k;
        const long long update_size = static_cast<long long>(rest) * rest;
        allocate(lower_size + update_size);
        panel.lower.assign(a.begin(), a.begin() + lower_size);
        update.values.resize(update_size);
        for (int j = 0; j < rest; ++j)
          for (int i = 0; i < rest; ++i)
            update.values[i + size_t(j) * rest] = a[k + i + size_t(k + j) * n];
        release(square);
        panel_entries += lower_size;
        panels_.push_back(std::move(panel));
      }
      release(scratch);
    }
    if (elimination_.size() != size_t(order_))
      return false;
    std::vector<bool> visited(order_, false);
    for (int i = 0; i < order_; ++i) {
      if (visited[i])
        continue;
      visited[i] = true;
      for (int j = elimination_[i]; j != i; j = elimination_[j]) {
        visited[j] = true;
        sign_ = -sign_;
      }
    }
    return true;
  }

  int sign() const { return sign_; }

  std::vector<double> solve(const std::vector<double> &rhs,
                            const Tacho::SymbolicTools &symbolic, int nrhs) {
    allocate(rhs.size());
    const auto order = symbolic.PermVector();
    std::vector<double> x(rhs.size());
    for (int p = 0; p < nrhs; ++p)
      for (int i = 0; i < order_; ++i)
        x[i + size_t(p) * order_] = rhs[order(i) + size_t(p) * order_];
    for (const auto &panel : panels_) {
      const int n = panel.ids.size(), k = panel.eliminated, rest = n - k;
      allocate(static_cast<long long>(n) * nrhs);
      std::vector<double> front(size_t(k) * nrhs);
      for (int p = 0; p < nrhs; ++p)
        for (int i = 0; i < k; ++i)
          front[i + size_t(p) * k] = x[panel.ids[i] + size_t(p) * order_];
      Tacho::Blas<double>::trsm('L', 'L', 'N', 'U', k, nrhs, 1.0,
                                panel.lower.data(), n, front.data(), k);
      if (rest) {
        std::vector<double> changed(size_t(rest) * nrhs);
        Tacho::Blas<double>::gemm('N', 'N', rest, nrhs, k, 1.0,
                                  panel.lower.data() + k, n, front.data(), k,
                                  0.0, changed.data(), rest);
        for (int p = 0; p < nrhs; ++p)
          for (int i = 0; i < rest; ++i)
            x[panel.ids[k + i] + size_t(p) * order_] -=
                changed[i + size_t(p) * rest];
      }
      for (int p = 0; p < nrhs; ++p)
        for (int i = 0; i < k; i += 2) {
          x[panel.ids[i] + size_t(p) * order_] =
              -front[i + 1 + size_t(p) * k] / panel.pivots[i / 2];
          x[panel.ids[i + 1] + size_t(p) * order_] =
              front[i + size_t(p) * k] / panel.pivots[i / 2];
        }
      release(static_cast<long long>(n) * nrhs);
    }
    for (auto it = panels_.rbegin(); it != panels_.rend(); ++it) {
      const auto &panel = *it;
      const int n = panel.ids.size(), k = panel.eliminated, rest = n - k;
      allocate(static_cast<long long>(n) * nrhs);
      std::vector<double> front(size_t(k) * nrhs), tail(size_t(rest) * nrhs);
      for (int p = 0; p < nrhs; ++p) {
        for (int i = 0; i < k; ++i)
          front[i + size_t(p) * k] = x[panel.ids[i] + size_t(p) * order_];
        for (int i = 0; i < rest; ++i)
          tail[i + size_t(p) * rest] = x[panel.ids[k + i] + size_t(p) * order_];
      }
      if (rest)
        Tacho::Blas<double>::gemm('T', 'N', k, nrhs, rest, -1.0,
                                  panel.lower.data() + k, n, tail.data(), rest,
                                  1.0, front.data(), k);
      Tacho::Blas<double>::trsm('L', 'L', 'T', 'U', k, nrhs, 1.0,
                                panel.lower.data(), n, front.data(), k);
      for (int p = 0; p < nrhs; ++p)
        for (int i = 0; i < k; ++i)
          x[panel.ids[i] + size_t(p) * order_] = front[i + size_t(p) * k];
      release(static_cast<long long>(n) * nrhs);
    }
    std::vector<double> solution(rhs.size());
    for (int p = 0; p < nrhs; ++p)
      for (int i = 0; i < order_; ++i)
        solution[order(i) + size_t(p) * order_] = x[i + size_t(p) * order_];
    release(rhs.size());
    return solution;
  }
};

using Complex = std::complex<double>;
using Matrix = std::vector<std::map<int, Complex>>;
using Graph =
    boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>;
using Solver = Tacho::SymbolicTools;
constexpr long long default_memory = 4096LL * 1024 * 1024;

struct Runtime {
  bool owner = !Kokkos::is_initialized();
  Runtime() {
    if (Kokkos::is_finalized())
      throw std::runtime_error("Kokkos already finalized");
    if (owner)
      Kokkos::initialize();
  }
  ~Runtime() {
    if (owner && Kokkos::is_initialized() && !Kokkos::is_finalized())
      Kokkos::finalize();
  }
};

Complex entry(const Matrix &a, int i, int j) {
  auto p = a[i].find(j);
  return p == a[i].end() ? Complex{} : p->second;
}

double skew_entry(const Matrix &a, int i, int j) {
  const double x = entry(a, i, j).real(), y = entry(a, j, i).real();
  return x == -y ? x : 0.5 * x - 0.5 * y;
}

int parity(const std::vector<int> &permutation) {
  std::vector<bool> visited(permutation.size(), false);
  int sign = 1;
  for (size_t i = 0; i < permutation.size(); ++i) {
    if (visited[i])
      continue;
    visited[i] = true;
    for (size_t j = permutation[i]; j != i; j = permutation[j]) {
      visited[j] = true;
      sign = -sign;
    }
  }
  return sign;
}

// Q=(I-iC)/sqrt(2), C a signed involution. At most four sparse insertions per
// entry.
Matrix transform(const Matrix &a, bool multiply_i) {
  const int n = a.size(), quarter = n / 4;
  Matrix result(n);
  auto partner = [quarter](int i) {
    return (3 - i / quarter) * quarter + i % quarter;
  };
  auto sign = [quarter](int i) {
    return i / quarter == 0 || i / quarter == 3 ? 1.0 : -1.0;
  };
  const Complex factor = multiply_i ? Complex(0, 0.5) : Complex(0.5, 0);
  for (int i = 0; i < n; ++i) {
    for (const auto &[j, value] : a[i]) {
      const Complex v = factor * value;
      result[i][j] += v;
      result[partner(i)][partner(j)] += sign(i) * sign(j) * v;
      result[partner(i)][j] += Complex(0, sign(i)) * v;
      result[i][partner(j)] -= Complex(0, sign(j)) * v;
    }
  }
  return result;
}

// Factor each connected component independently. Pairing uses a general-graph
// perfect matching, not the bipartite matching assumed by Tacho's MATLAB
// driver.
int pfaffian(const Matrix &input, double tolerance, double &residual,
             long long &fill, long long max_memory, long long &peak_memory) {
  static Runtime runtime;
  const int n = input.size();
  if (n < 2 || n % 2)
    throw std::invalid_argument("Invalid skew order");
  Graph graph(n);
  double scale = 0;
  for (int i = 0; i < n; ++i) {
    for (const auto &[j, v] : input[i]) {
      if (!std::isfinite(v.real()) || !std::isfinite(v.imag()))
        throw std::invalid_argument("Nonfinite skew matrix");
      scale = std::max(scale, std::abs(v));
      // Tolerated one-sided roundoff must not give an asymmetric graph.
      if (i != j && v.real() != 0 && skew_entry(input, i, j) != 0 &&
          (j > i || entry(input, j, i).real() == 0))
        boost::add_edge(i, j, graph);
    }
  }
  if (scale == 0)
    return 0;
  for (int i = 0; i < n; ++i)
    for (const auto &[j, v] : input[i])
      if (std::abs(v.imag()) > tolerance * scale ||
          std::abs(v + entry(input, j, i)) > tolerance * scale)
        throw std::invalid_argument("Not real skew symmetric");

  std::vector<int> component(n);
  const int count = boost::connected_components(graph, component.data());
  std::vector<std::vector<int>> members(count);
  for (int i = 0; i < n; ++i)
    members[component[i]].push_back(i);
  const auto unmatched = boost::graph_traits<Graph>::null_vertex();
  std::vector<Graph::vertex_descriptor> mate(n, unmatched);
  std::vector<std::tuple<double, int, int>> edges;
  for (int i = 0; i < n; ++i)
    for (const auto &[j, v] : input[i])
      if (i != j && v.real() != 0 && skew_entry(input, i, j) != 0 &&
          (j > i || entry(input, j, i).real() == 0))
        edges.emplace_back(-std::abs(skew_entry(input, i, j)), i, j);
  std::sort(edges.begin(), edges.end());
  // Start with large entries: a purely structural matching can give a singular
  // intermediate Schur complement even when the complete matrix is well gapped.
  for (const auto &[weight, i, j] : edges)
    if (mate[i] == unmatched && mate[j] == unmatched) {
      mate[i] = j;
      mate[j] = i;
    }
  auto vertex_index = boost::get(boost::vertex_index, graph);
  boost::edmonds_augmenting_path_finder<Graph, Graph::vertex_descriptor *,
                                        decltype(vertex_index)>
      augmentor(graph, mate.data(), vertex_index);
  while (augmentor.augment_matching()) {
  }
  augmentor.get_current_matching(mate.data());
  if (2 * boost::matching_size(graph, mate.data()) != static_cast<size_t>(n))
    return 0;

  std::vector<int> permutation, inverse(n);
  permutation.reserve(n);
  for (const auto &vertices : members)
    for (int i : vertices)
      if (i < static_cast<int>(mate[i])) {
        permutation.push_back(i);
        permutation.push_back(mate[i]);
      }
  for (int i = 0; i < n; ++i)
    inverse[permutation[i]] = i;
  int result = parity(permutation), begin = 0;
  for (const auto &vertices : members) {
    const int m = vertices.size();
    if (m % 2)
      return 0;
    // Exact sparse 2x2 component; also avoids Tacho's unsupported one-node
    // graph.
    if (m == 2) {
      const double pivot =
          skew_entry(input, permutation[begin], permutation[begin + 1]);
      if (std::abs(pivot) <= tolerance * scale)
        return 0;
      result *= pivot > 0 ? 1 : -1;
      begin += m;
      continue;
    }
    std::vector<std::map<int, double>> rows(m);
    size_t nnz = 0;
    for (int i = 0; i < m; ++i)
      rows[i][i] = 0;
    for (int i = 0; i < m; ++i) {
      for (const auto &[j, v] : input[permutation[begin + i]]) {
        const double value = skew_entry(input, permutation[begin + i], j);
        if (value != 0) {
          const int column = inverse[j] - begin;
          if (column < 0 || column >= m)
            throw std::invalid_argument("Asymmetric component pattern");
          rows[i][column] = value / scale;
          rows[column][i] = -rows[i][column];
        }
      }
    }
    for (int i = 0; i < m; ++i)
      nnz += rows[i].size();
    if (nnz > static_cast<size_t>(std::numeric_limits<int>::max()))
      throw std::overflow_error("Tacho index overflow");
    Solver::size_type_array ap("rowptr", m + 1);
    Solver::ordinal_type_array aj("columns", nnz);
    size_t k = 0;
    for (int i = 0; i < m; ++i) {
      ap(i) = k;
      for (const auto &[j, v] : rows[i])
        aj(k++) = j;
    }
    ap(m) = k;
    Tacho::Graph pattern(m, nnz, ap, aj);
    Tacho::GraphTools_Metis ordering(pattern);
    ordering.reorder(0);
    Solver symbolic(m, ap, aj, ordering.PermVector(), ordering.InvPermVector());
    symbolic.symbolicFactorize(0);
    SkewLDL factors(m, tolerance, max_memory, peak_memory);
    if (!factors.factor(rows, symbolic))
      return 0;
    fill += factors.panel_entries;

    // Independent right-hand sides also probe nullspaces. The independent
    // MUMPS gap test is still necessary; these are not a global error bound.
    constexpr int probes = 4;
    factors.reserve_rhs(probes);
    std::vector<double> rhs(size_t(m) * probes);
    for (int p = 0; p < probes; ++p)
      for (int i = 0; i < m; ++i)
        rhs[i + size_t(p) * m] =
            std::sin((i + 1.0) * (p + 1.0)) + std::cos((i + 1.0) / (p + 2.0));
    const auto solution = factors.solve(rhs, symbolic, probes);
    double error = 0.0;
    for (int p = 0; p < probes; ++p) {
      double numerator = 0.0, denominator = 0.0;
      for (int i = 0; i < m; ++i) {
        double ax = 0.0;
        for (const auto &[j, value] : rows[i])
          ax += value * solution[j + size_t(p) * m];
        const double r = ax - rhs[i + size_t(p) * m];
        if (!std::isfinite(r))
          return 0;
        numerator += r * r;
        denominator += rhs[i + size_t(p) * m] * rhs[i + size_t(p) * m];
      }
      error = std::max(error, std::sqrt(numerator / denominator));
    }
    residual = std::max(residual, error);
    if (!std::isfinite(error) || error > tolerance)
      return 0;
    const int sign = factors.sign();
    result *= sign;
    begin += m;
  }
  return result;
}
} // namespace

extern "C" int
cp2k_tacho_pfaffian_limited(int n, int nnz, const int *rows, const int *cols,
                            const double *values, double tolerance, int *sign,
                            double *residual, long long *fill,
                            long long max_memory, long long *peak_memory) {
  *sign = 0;
  *peak_memory = 0;
  *residual = 0;
  *fill = 0;
  try {
    if (n < 2 || n % 2 || nnz < 0 || max_memory <= 0 || !(tolerance > 0) ||
        !std::isfinite(tolerance))
      return -1;
    Matrix a(n);
    for (int k = 0; k < nnz; ++k) {
      if (rows[k] < 1 || rows[k] > n || cols[k] < 1 || cols[k] > n)
        return -1;
      a[rows[k] - 1][cols[k] - 1] += values[k];
    }
    *sign = pfaffian(a, tolerance, *residual, *fill, max_memory, *peak_memory);
    return *sign == 0 ? 1 : 0;
  } catch (const std::invalid_argument &) {
    return -1;
  } catch (const PfaffianMemoryLimit &) {
    return -6;
  } catch (const std::bad_alloc &) {
    return -6;
  } catch (...) {
    return -4;
  }
}

extern "C" int cp2k_tacho_z2_limited(int order, int nnz, const int *rows,
                                     const int *cols, const double *a,
                                     const double *b, double tolerance,
                                     double metric_ratio, double gap,
                                     int *index, double *residual,
                                     long long *fill, long long max_memory,
                                     long long *peak_memory) {
  *index = 0;
  *peak_memory = 0;
  *residual = 0;
  *fill = 0;
  try {
    if (order < 8 || order % 8 || nnz < 0 || max_memory <= 0 ||
        !std::isfinite(tolerance) || !std::isfinite(metric_ratio) ||
        !std::isfinite(gap) || tolerance <= 0 || tolerance >= 1 ||
        metric_ratio <= 0 || metric_ratio > 1 || gap <= 0)
      return -1;
    const int n = order / 2;
    Matrix localizer(n), metric(n);
    std::vector<double> bnorm(order, 0);
    for (int k = 0; k < nnz; ++k) {
      const int i = rows[k] - 1, j = cols[k] - 1;
      if (j < 0 || i < j || i >= order || !std::isfinite(a[k]) ||
          !std::isfinite(b[k]))
        return -1;
      bnorm[i] += std::abs(b[k]);
      if (i != j)
        bnorm[j] += std::abs(b[k]);
      if (i < n) {
        localizer[i][j] += a[k];
        metric[i][j] += b[k];
        if (i != j) {
          localizer[j][i] += a[k];
          metric[j][i] += b[k];
        }
      } else if (j < n) {
        localizer[i - n][j] += Complex(0, a[k]);
        metric[i - n][j] += Complex(0, b[k]);
      }
    }
    Matrix skew = transform(localizer, true),
           real_metric = transform(metric, false);
    double scale = 0, structure_error = 0, metric_scale = 0, metric_error = 0;
    std::vector<double> discarded_rows(n, 0), discarded_cols(n, 0);
    for (int i = 0; i < n; ++i) {
      for (auto &[j, v] : skew[i]) {
        if (!std::isfinite(v.real()) || !std::isfinite(v.imag()))
          return -1;
        scale = std::max(scale, std::abs(v));
        const double projected = skew_entry(skew, i, j);
        const double defect = std::abs(v - projected);
        discarded_rows[i] += defect;
        discarded_cols[j] += defect;
        if (!skew[j].count(i)) {
          discarded_rows[j] += std::abs(projected);
          discarded_cols[i] += std::abs(projected);
        }
        v = v.real();
      }
      for (const auto &[j, v] : real_metric[i]) {
        if (!std::isfinite(v.real()) || !std::isfinite(v.imag()))
          return -1;
        metric_scale = std::max(metric_scale, std::abs(v));
        metric_error = std::max(metric_error, std::abs(v.imag()));
      }
    }
    structure_error = std::max(
        *std::max_element(discarded_rows.begin(), discarded_rows.end()),
        *std::max_element(discarded_cols.begin(), discarded_cols.end()));
    if (structure_error > tolerance * scale ||
        metric_error > tolerance * metric_scale)
      return -3;
    const double metric_min =
        metric_ratio * *std::max_element(bnorm.begin(), bnorm.end());
    if (structure_error >= 0.25 * gap * metric_min)
      return 1;
    // Positive TR-compatible S deforms continuously to I without a Pfaffian
    // zero. In this fixed component ordering the atomic reference has sign
    // (-1)^(n/4).
    const int reference = (n / 4) % 2 ? -1 : 1;
    const int sign =
        pfaffian(skew, tolerance, *residual, *fill, max_memory, *peak_memory);
    if (!sign)
      return 1;
    *index = (1 - reference * sign) / 2;
    return 0;
  } catch (const PfaffianMemoryLimit &) {
    return -6;
  } catch (const std::bad_alloc &) {
    return -6;
  } catch (...) {
    return -4;
  }
}

// Preserve the original C entry points for external validation clients.
extern "C" int cp2k_tacho_pfaffian(int n, int nnz, const int *rows,
                                   const int *cols, const double *values,
                                   double tolerance, int *sign,
                                   double *residual, long long *fill) {
  long long peak_memory = 0;
  return cp2k_tacho_pfaffian_limited(n, nnz, rows, cols, values, tolerance,
                                     sign, residual, fill, default_memory,
                                     &peak_memory);
}

extern "C" int cp2k_tacho_z2(int order, int nnz, const int *rows,
                             const int *cols, const double *a, const double *b,
                             double tolerance, double metric_ratio, double gap,
                             int *index, double *residual, long long *fill) {
  long long peak_memory = 0;
  return cp2k_tacho_z2_limited(order, nnz, rows, cols, a, b, tolerance,
                               metric_ratio, gap, index, residual, fill,
                               default_memory, &peak_memory);
}
