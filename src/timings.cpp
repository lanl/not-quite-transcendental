//======================================================================
// © 2022-2025. Triad National Security, LLC. All rights reserved.  This
// program was produced under U.S. Government contract 89233218CNA000001
// for Los Alamos National Laboratory (LANL), which is operated by Triad
// National Security, LLC for the U.S.  Department of Energy/National
// Nuclear Security Administration. All rights in the program are
// reserved by Triad National Security, LLC, and the U.S. Department of
// Energy/National Nuclear Security Administration. The Government is
// granted for itself and others acting on its behalf a nonexclusive,
// paid-up, irrevocable worldwide license in this material to reproduce,
// prepare derivative works, distribute copies to the public, perform
// publicly and display publicly, and to permit others to do so.
//======================================================================

// C headers
#include <cmath>
#include <cstdio>
#include <cstdlib>

// C++ headers
#include <chrono>
#include <iostream>
#include <limits>

#include "Kokkos_Core.hpp"

// Local headers
#include "logs.hpp"

using duration_t = std::chrono::nanoseconds;

constexpr int NTRIALS = 10;
template <typename Operation, typename Duration = duration_t>
auto TimeOperation(const Operation &op) {
  Kokkos::fence();
  auto start = std::chrono::high_resolution_clock::now();
  //for (int t = 0; t < NTRIALS; ++t) {
    op();
    op();
    op();
    op();
    op();
    op();
    op();
    op();
    op();
    op();
  //  }
  Kokkos::fence();
  auto stop = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<Duration>(stop - start);
}

template <typename Kernel, typename Duration = duration_t>
auto TimeKernel(const char *name, const int ntimes, const Kernel &k) {
  return TimeOperation([=]() { Kokkos::parallel_for(name, ntimes, k); });
}

KOKKOS_INLINE_FUNCTION
double get_dx(const double min, const double max, const int N) {
  return (max - min) / (N - 1);
}

KOKKOS_INLINE_FUNCTION
double get_diff(const double a, const double b) {
  return 2. * std::abs(a - b) /
         (std::abs(a) + std::abs(b) + std::numeric_limits<double>::epsilon());
}

int main(int argc, char *argv[]) {
  Kokkos::initialize(argc, argv);
  {
    if (argc < 2) {
      std::cerr << "Usage: " << argv[0] << " npoints" << std::endl;
      std::exit(1);
    }

    const int NPOINTS = std::atoi(argv[1]);
    double min, max, dx;
    Kokkos::View<double *> lx("lx", NPOINTS);
    Kokkos::View<double *> ex("ex", NPOINTS);
    Kokkos::View<double *> fastex("fastex", NPOINTS);
    Kokkos::View<double *> scratch("scratch", NPOINTS);

    // exp10
    min = -10;
    max = 10;
    dx = get_dx(min, max, NPOINTS);
    Kokkos::parallel_for(
        "setup", NPOINTS, KOKKOS_LAMBDA(const int i) { lx[i] = min + dx * i; });
    auto dt_std_exp10 = TimeKernel(
        "std::exp10", NPOINTS, KOKKOS_LAMBDA(const int i) {
	  ex[i] = std::pow(10, lx[i]);
        });
    auto dt_fast_exp10 = TimeKernel(
        "fast::exp10", NPOINTS, KOKKOS_LAMBDA(const int i) {
	  fastex[i] = FastMath::pow10(lx[i]);
        });

    // log10
    auto dt_std_log10 = TimeKernel(
        "std::log10", NPOINTS, KOKKOS_LAMBDA(const int i) {
		  scratch[i] = std::log10(ex[i]);
        });
    double std_diff = 0;
    Kokkos::parallel_reduce(
        "check diff", NPOINTS,
        KOKKOS_LAMBDA(const int i, double &d) { d += get_diff(scratch[i], lx[i]); },
        std_diff);

    auto dt_fast_log10 = TimeKernel(
        "fast::log10", NPOINTS, KOKKOS_LAMBDA(const int i) {
	  scratch[i] = FastMath::log10(fastex[i]);
        });
    double fast_diff = 0;
    Kokkos::parallel_reduce(
        "check diff", NPOINTS,
        KOKKOS_LAMBDA(const int i, double &d) { d += get_diff(scratch[i], lx[i]); },
        fast_diff);

    printf("Invertibility of log10 (smaller is better):\n"
           "\tstdlib: %.14e\n"
           "\tfast:   %.14e\n",
           std_diff / static_cast<double>(NPOINTS),
           fast_diff / static_cast<double>(NPOINTS));

    /*
    // tanh
    auto dt_std_tanh = TimeKernel(
        "std::tanh", NPOINTS, KOKKOS_LAMBDA(const int i) {
	  scratch[i] = std::tanh(lx[i]);
        });
    auto dt_fast_tanh = TimeKernel(
        "fast::tanh", NPOINTS, KOKKOS_LAMBDA(const int i) {
	  scratch[i] = FastMath::tanh(lx[i]);
        });
    */

    // Report
    printf("Timings for pow10:\n"
           "\tstdlib: %.14e (ns/point)\n"
           "\tfast:   %.14e (ns/point)\n"
           "\tratio:  %.14e\n"
           "Timings for log10:\n"
           "\tstdlib: %.14e (ns/point)\n"
           "\tfast:   %.14e (ns/point)\n"
           "\tratio:  %.14e\n",
           // "Timings for tanh:\n"
           // "\tstdlib: %.14e (ns/point)\n"
           // "\tfast:   %.14e (ns/point)\n"
           // "\tratio:  %.14e\n",
           dt_std_exp10.count() / static_cast<double>(NPOINTS * NTRIALS),
           dt_fast_exp10.count() / static_cast<double>(NPOINTS * NTRIALS),
           dt_std_exp10.count() / static_cast<double>(dt_fast_exp10.count()),
           dt_std_log10.count() / static_cast<double>(NPOINTS * NTRIALS),
           dt_fast_log10.count() / static_cast<double>(NPOINTS * NTRIALS),
           dt_std_log10.count() / static_cast<double>(dt_fast_log10.count())
           );
           // dt_std_tanh.count() / static_cast<double>(NPOINTS * NTRIALS),
           // dt_fast_tanh.count() / static_cast<double>(NPOINTS * NTRIALS),
           // dt_std_tanh.count() / static_cast<double>(dt_fast_tanh.count()));
  }
  Kokkos::finalize();

  return 0;
}
