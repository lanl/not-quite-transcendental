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
#include <cstdint>
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
  return *reinterpret_cast<std::int64_t*>(&f);
}
KOKKOS_FORCEINLINE_FUNCTION
auto as_double(std::int64_t i) {
  return *reinterpret_cast<double*>(&i);
}

// First order interpolation based NQTs
// ----------------------------------------------------------------------
// Reference implementations, however the integer cast implementation
// below is probably faster.
KOKKOS_FORCEINLINE_FUNCTION
double lg_o1_portable(const double x) {
  int e;
  assert(x > 0 && "log divergent for x <= 0");
  const double m = frexp(x, &e);
  return 2 * (m - 1) + e;
}


KOKKOS_FORCEINLINE_FUNCTION
double pow2_o1_portable(const double x) {
  const int flr = std::floor(x);
  const double remainder = x - flr;
  const double mantissa = 0.5 * (remainder + 1);
  const double exponent = flr + 1;
  return ldexp(mantissa, exponent);
}

// Integer aliased versions
KOKKOS_FORCEINLINE_FUNCTION
double lg_o1_aliased(const double x) {
  // Magic numbers constexpr because C++ doesn't constexpr reinterpret casts
  // these are floating point numbers as reinterpreted as integers.
  // as_int(1.0)
  constexpr std::int64_t one_as_int = 4607182418800017408;
  // 1./static_cast<double>(as_int(2.0) - as_int(1.0))
  constexpr double scale_down = 2.22044604925031e-16;
  return static_cast<double>(as_int(x) - one_as_int) * scale_down;
}

KOKKOS_FORCEINLINE_FUNCTION
double pow2_o1_aliased(const double x) {
  // Magic numbers constexpr because C++ doesn't constexpr reinterpret casts
  // these are floating point numbers as reinterpreted as integers.
  // as_int(1.0)
  constexpr std::int64_t one_as_int = 4607182418800017408;
  // as_int(2.0) - as_int(1.0)
  constexpr double scale_up = 4503599627370496;
  return as_double(static_cast<std::int64_t>(x*scale_up) + one_as_int);
}
// ----------------------------------------------------------------------

// Second-order interpolation based NQTs
// These implementations are due to Peter Hammond
// ----------------------------------------------------------------------
// Portable versions that use frexp/ldexp rather than integer aliasing
KOKKOS_FORCEINLINE_FUNCTION
double lg_o2_portable(const double x) {
  constexpr double four_thirds = 4./3.;
  int e;
  assert(x > 0 && "log divergent for x <= 0");
  const double m = frexp(x, &e);
  return e - four_thirds*(m - 2)*(m - 1);
}

// This version uses the exact formula 
KOKKOS_FORCEINLINE_FUNCTION
double pow2_o2_portable(const double x) {
  // log2(mantissa). should go between -1 and 0
  const int flr = std::floor(x);
  const double lm = x - flr - 1; 
  const double mantissa = 0.5*(3 - std::sqrt(1 - 3*lm));
  const double exponent = flr + 1;
  return ldexp(mantissa, exponent);
}

// Integer aliased/bithacked versions
KOKKOS_FORCEINLINE_FUNCTION
double lg_o2_aliased(const double x) {
  // as_int(1.0) == 2^62 - 2^52
  constexpr std::int64_t one_as_int = 4607182418800017408;
  // 1/(as_int(2.0) - as_int(1.0)) == 2^-52
  constexpr double scale_down = 2.220446049250313e-16;
  // 2^52 - 1
  constexpr std::int64_t mantissa_mask = 4503599627370495;
  // 2^26 - 1
  constexpr std::int64_t low_mask = 67108863;
  
  const std::int64_t x_as_int = as_int(x) - one_as_int;
  const std::int64_t frac_as_int = x_as_int & mantissa_mask;
  const std::int64_t frac_high = frac_as_int>>26;
  const std::int64_t frac_low  = frac_as_int & low_mask;
  const std::int64_t frac_squared = frac_high*frac_high + ((frac_high*frac_low)>>25);
  
  return static_cast<double>(x_as_int +
                             ((frac_as_int - frac_squared)/3)) * scale_down;
}

KOKKOS_FORCEINLINE_FUNCTION
double pow2_o2_aliased(const double x) {
  // as_int(1.0) == 2^62 - 2^52
  constexpr std::int64_t one_as_int = 4607182418800017408;
  // as_int(2.0) - as_int(1.0) == 2^52
  constexpr double scale_up = 4503599627370496;
  constexpr std::int64_t mantissa_mask = 4503599627370495; // 2^52 - 1
  constexpr std::int64_t a = 9007199254740992; // 2 * 2^52
  constexpr double b = 67108864; // 2^26
  constexpr std::int64_t c = 18014398509481984; // 4 * 2^52
  
  const std::int64_t x_as_int = static_cast<std::int64_t>(x*scale_up);
  const std::int64_t frac_as_int = x_as_int & mantissa_mask;
  const std::int64_t frac_sqrt = static_cast<std::int64_t>(
                                         b*std::sqrt(static_cast<double>(c-3*frac_as_int)));
  
  return as_double(x_as_int + a - frac_sqrt - frac_as_int + one_as_int);
}
// ----------------------------------------------------------------------

KOKKOS_FORCEINLINE_FUNCTION
double lg(const double x) {
#ifdef NQT_ORDER_1
#ifdef NQT_PORTABLE
  return lg_o1_portable(x);
#else
  return lg_o1_aliased(x);
#endif // PORTABLE
#else
#ifdef NQT_PORTABLE
  return lg_o2_portable(x);
#else
  return lg_o2_aliased(x);
#endif // PORTABLE
#endif // NQT_ORDER
}

KOKKOS_FORCEINLINE_FUNCTION
double pow2(const double x) {
#ifdef NQT_ORDER_1
#ifdef NQT_PORTABLE
  return pow2_o1_portable(x);
#else
  return pow2_o1_aliased(x);
#endif // PORTABLE
#else
#ifdef NQT_PORTABLE
  return pow2_o2_portable(x);
#else
  return pow2_o2_aliased(x);
#endif // PORTABLE
#endif // NQT_ORDER
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
double sinh(const double x) {
  constexpr double IE = 1.0/M_E;
  constexpr double LG2 = 1.4426950408889634074;
  const double a2x = 2*std::abs(x);
  const double mask = (a2x < M_E); // to make expr below single type
  return mask * a2x*IE + (1.-mask)*LG2*sgn(x)*lg(a2x);
}

KOKKOS_FORCEINLINE_FUNCTION
double asinh(const double x) {
  constexpr double LG2 = 1.4426950408889634074;
  const double ax = std::abs(x);
  const double mask = (ax < 1.0);
  return mask * 0.5*M_E*x + (1.-mask)*0.5*sgn(x)*pow2(LG2*ax);
}

} // namespace FastMath
