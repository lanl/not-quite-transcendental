//======================================================================
// © 2022. Triad National Security, LLC. All rights reserved.  This
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

// Code taken from singularity-eos
// https://github.com/lanl/singularity-eos

#pragma once

#include <cassert>
#include <cmath>

#include "Kokkos_Core.hpp"

namespace FastMath {

template<typename T>
KOKKOS_FORCEINLINE_FUNCTION
auto sgn(const T &x) {
  return (T(0) < x) - (x < T(0));
}

KOKKOS_FORCEINLINE_FUNCTION
auto as_int(double f) {
  return *reinterpret_cast<long long int*>(&f);
}
KOKKOS_FORCEINLINE_FUNCTION
auto as_double(long long int i) {
  return *reinterpret_cast<double*>(&i);
}

// Reference implementations, however the integer cast implementation
// below is probably faster.
/*
KOKKOS_FORCEINLINE_FUNCTION
double lg(const double x) {
  int n;
  assert(x > 0 && "log divergent for x <= 0");
  const double y = frexp(x, &n);
  return 2 * (y - 1) + n;
}


KOKKOS_FORCEINLINE_FUNCTION
double pow2(const double x) {
  const int flr = std::floor(x);
  const double remainder = x - flr;
  const double mantissa = 0.5 * (remainder + 1);
  const double exponent = flr + 1;
  return ldexp(mantissa, exponent);
}
*/

KOKKOS_FORCEINLINE_FUNCTION
double lg(const double x) {
  // Magic numbers constexpr because C++ doesn't constexpr reinterpret casts
  // these are floating point numbers as reinterpreted as integers.
  // as_int(1.0)
  constexpr long long int one_as_int = 4607182418800017408;
  // 1./static_cast<double>(as_int(2.0) - as_int(1.0))
  constexpr double scale_down = 2.22044604925031e-16;
  return static_cast<double>(as_int(x) - one_as_int) * scale_down;
}

KOKKOS_FORCEINLINE_FUNCTION
double pow2(const double x) {
  // Magic numbers constexpr because C++ doesn't constexpr reinterpret casts
  // these are floating point numbers as reinterpreted as integers.
  // as_int(1.0)
  constexpr long long int one_as_int = 4607182418800017408;
  // as_int(2.0) - as_int(1.0)
  constexpr double scale_up = 4503599627370496;
  return as_double(static_cast<long long int>(x*scale_up) + one_as_int);
}

KOKKOS_FORCEINLINE_FUNCTION
double ln(const double x) {
  constexpr double ILOG2E = 0.6931471805599453;
  return ILOG2E * lg(x);
}

KOKKOS_FORCEINLINE_FUNCTION
double exp(const double x) {
  constexpr double LOG2E = 1.4426950408889634;
  return pow2(LOG2E * x);
}

KOKKOS_FORCEINLINE_FUNCTION
double log10(const double x) {
  constexpr double LOG2OLOG10 = 0.301029995663981195;
  return LOG2OLOG10 * lg(x);
}

KOKKOS_FORCEINLINE_FUNCTION
double pow10(const double x) {
  constexpr double LOG10OLOG2 = 3.321928094887362626;
  return pow2(LOG10OLOG2 * x);
}

KOKKOS_FORCEINLINE_FUNCTION
double tanh(const double x) {
  const double expx = exp(2 * x);
  return (expx - 1) / (expx + 1);
}

KOKKOS_FORCEINLINE_FUNCTION
double sinh(const double x) {
  constexpr double IE = 1.0/M_E;
  constexpr double LG2 = 1.4426950408889634074;
  const double a2x = 2*std::abs(x);
  const bool mask = (a2x < M_E);
  return mask * a2x*IE + !mask*LG2*sgn(x)*lg(a2x);
}

KOKKOS_FORCEINLINE_FUNCTION
double asinh(const double x) {
  constexpr double LG2 = 1.4426950408889634074;
  const double ax = std::abs(x);
  const bool mask = (ax < 1.0);
  return mask * 0.5*M_E*x + !mask*0.5*sgn(x)*pow2(LG2*ax);
}

} // namespace FastMath
