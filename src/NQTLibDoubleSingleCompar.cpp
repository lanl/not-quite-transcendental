#include <cmath>
#include <cassert>
#include <cstdlib> 
#include <iostream>
#include <random>
#include <chrono>

template<typename T>
auto sgn(const T &x) {
  return (T(0) < x) - (x < T(0));
}

// Double Precision
int64_t as_int(double f) {
  return *reinterpret_cast<int64_t*>(&f);
}

double as_real(int64_t i) {
  return *reinterpret_cast<double*>(&i);
}

// Single Precision
int32_t as_int(float f) {
  return *reinterpret_cast<int32_t*>(&f);
}

float as_real(int32_t i) {
  return *reinterpret_cast<float*>(&i);
}

// Double Precision
double NQT_log_O1(const double x) {
  // as_int(1.0) == 2^30 - 2^23
  constexpr int64_t one_as_int = 4607182418800017408;
  // 1./static_cast<double>(as_int(2.0) - as_int(1.0))
  constexpr double scale_down = 2.220446049250313e-16;
  return static_cast<double>(as_int(x) - one_as_int) * scale_down;
}

double NQT_exp_O1(const double x) {
  // as_int(1.0) == 2^30 - 2^23
  constexpr int64_t one_as_int = 4607182418800017408;
  // 1./static_cast<double>(as_int(2.0) - as_int(1.0))
  constexpr double scale_up = 4503599627370496;
  return as_real(static_cast<int64_t>(x*scale_up) + one_as_int);
}

double NQT_log_O2(const double x) {
      // as_int(1.0) == 2^62 - 2^52
      constexpr int64_t one_as_int = 4607182418800017408;
      // 1/(as_int(2.0) - as_int(1.0)) == 2^-52
      constexpr double scale_down = 2.220446049250313e-16;
      // 2^52 - 1
      constexpr int64_t mantissa_mask = 4503599627370495;
      // 2^26 - 1
      constexpr int64_t low_mask = 67108863;

      const int64_t x_as_int = as_int(x) - one_as_int;
      const int64_t frac_as_int = x_as_int & mantissa_mask;
      const int64_t frac_high = frac_as_int>>26;
      const int64_t frac_low  = frac_as_int & low_mask;
      const int64_t frac_squared = frac_high*frac_high + ((frac_high*frac_low)>>25);
      
      return static_cast<double>(x_as_int + ((frac_as_int - frac_squared)/3)) * scale_down;
  }

double NQT_exp_O2(const double x) {
      // as_int(1.0) == 2^62 - 2^52
      constexpr int64_t one_as_int = 4607182418800017408;
      // as_int(2.0) - as_int(1.0) == 2^52
      constexpr double scale_up = 4503599627370496.;
      constexpr int64_t mantissa_mask = 4503599627370495; // 2^52 - 1
      constexpr int64_t a = 9007199254740992; // 2 * 2^52
      constexpr int64_t c = 18014398509481984; // 4 * 2^52

      const int64_t x_as_int = static_cast<int64_t>(x*scale_up);
      const int64_t frac_as_int = x_as_int & mantissa_mask;
      const int64_t frac_sqrt = static_cast<int64_t>(std::sqrt(scale_up*static_cast<double>(c-3*frac_as_int)));

      return as_real(x_as_int + a - frac_sqrt - frac_as_int + one_as_int);
  }

// Single Precision
float NQT_log_O1(const float x) {
  // as_int(1.0) == 2^30 - 2^23
  constexpr int32_t one_as_int = 1065353216;
  // 1./static_cast<double>(as_int(2.0) - as_int(1.0))
  // constexpr float scale_down = 1.1920928955078125e-07f;
  constexpr float scale_down = 1.f/8388608.f;
  return static_cast<float>(as_int(x) - one_as_int) * scale_down;
}

float NQT_exp_O1(const float x) {
  // as_int(1.0) == 2^30 - 2^23
  constexpr int32_t one_as_int = 1065353216;
  // static_cast<double>(as_int(2.0) - as_int(1.0))
  constexpr float scale_up = 8388608;
  return as_real(static_cast<int32_t>(x*scale_up) + one_as_int);
}

float NQT_log_O2(const float x) {
  // as_int(1.0) == 2^30 - 2^23
  constexpr int32_t one_as_int = 1065353216;
  // 1/(as_int(2.0) - as_int(1.0)) == 2^-23
  constexpr float scale_down = 1.1920929e-07f;
  // 2^23 - 1
  constexpr int32_t mantissa_mask = 8388607;
  // 2^11 - 1
  constexpr int32_t low_mask = 2047;

  const int32_t x_as_int = as_int(x) - one_as_int;
  const int32_t frac_as_int = x_as_int & mantissa_mask;
  const int32_t frac_high = frac_as_int>>11;
  const int32_t frac_low  = frac_as_int & low_mask;
  const int32_t frac_squared = (frac_high*frac_high + ((frac_high*frac_low)>>10))>>1;
  
  return static_cast<float>(x_as_int + ((frac_as_int - frac_squared)/3)) * scale_down;
  }

float NQT_exp_O2(const float x) {
  // as_int(1.0) == 2^30 - 2^23
  constexpr int32_t one_as_int = 1065353216;
  // as_int(2.0) - as_int(1.0) == 2^23
  constexpr float scale_up = 8388608;
  constexpr int32_t mantissa_mask = 8388607; // 2^23 - 1
  constexpr int32_t a = 16777216; // 2 * 2^23
  constexpr int32_t c = 33554432; // 4 * 2^23

  const int32_t x_as_int = static_cast<int32_t>(x*scale_up);
  const int32_t frac_as_int = x_as_int & mantissa_mask;
  const int32_t frac_sqrt = static_cast<int32_t>(std::sqrt(scale_up*static_cast<float>(c-3*frac_as_int)));

  return as_real(x_as_int + a - frac_sqrt - frac_as_int + one_as_int);
  }


int main () {
  const int32_t samples = 100000000;

  float*  rands_float  = new float[samples];
  double* rands_double = new double[samples];
  float*  logs_float  = new float[samples];
  double* logs_double = new double[samples];

  for (int32_t i=0; i<samples; ++i) {
    double rand_double_01 = static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
    float  rand_float_01  = static_cast<float>(rand_double_01);

    rands_float[i]  = 1.0 + rand_float_01;
    rands_double[i] = 1.0 + rand_double_01;
    logs_float[i]   = 1.0 + rand_float_01;
    logs_double[i]  = 1.0 + rand_double_01;
  }


  std::cout << "Runtime (seconds):" << std::endl;

  std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();
  std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();

  // Time single log
  t_start = std::chrono::high_resolution_clock::now();
  for (int32_t i=0; i<samples; ++i) {
    logs_float[i] = std::log2(rands_float[i]);
  }
  t_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span_single_log = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start);
  std::cout << "single log STL: " << time_span_single_log.count() << std::endl;

  // Time single log O1
  t_start = std::chrono::high_resolution_clock::now();
  for (int32_t i=0; i<samples; ++i) {
    logs_float[i] = NQT_log_O1(rands_float[i]);
  }
  t_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span_single_log_O1 = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start);
  std::cout << "single log O1:  " << time_span_single_log_O1.count() << std::endl;

  // Time single log O2
  t_start = std::chrono::high_resolution_clock::now();
  for (int32_t i=0; i<samples; ++i) {
    logs_float[i] = NQT_log_O2(rands_float[i]);
  }
  t_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span_single_log_O2 = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start);
  std::cout << "single log O2:  " << time_span_single_log_O2.count() << std::endl;


  // Time double log
  t_start = std::chrono::high_resolution_clock::now();
  for (int32_t i=0; i<samples; ++i) {
    logs_double[i] = std::log2(rands_double[i]);
  }
  t_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span_double_log = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start);
  std::cout << "double log STL: " << time_span_double_log.count() << std::endl;

  // Time double log O1
  t_start = std::chrono::high_resolution_clock::now();
  for (int32_t i=0; i<samples; ++i) {
    logs_double[i] = NQT_log_O1(rands_double[i]);
  }
  t_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span_double_log_O1 = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start);
  std::cout << "double log O1:  " << time_span_double_log_O1.count() << std::endl;

  // Time double log O2
  t_start = std::chrono::high_resolution_clock::now();
  for (int32_t i=0; i<samples; ++i) {
    logs_double[i] = NQT_log_O2(rands_double[i]);
  }
  t_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span_double_log_O2 = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start);
  std::cout << "double log O2:  " << time_span_double_log_O2.count() << std::endl;



  // Time single exp
  t_start = std::chrono::high_resolution_clock::now();
  for (int32_t i=0; i<samples; ++i) {
    logs_float[i] = std::exp2(rands_float[i]);
  }
  t_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span_single_exp = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start);
  std::cout << "single exp STL: " << time_span_single_exp.count() << std::endl;

  // Time single log O1
  t_start = std::chrono::high_resolution_clock::now();
  for (int32_t i=0; i<samples; ++i) {
    logs_float[i] = NQT_exp_O1(rands_float[i]);
  }
  t_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span_single_exp_O1 = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start);
  std::cout << "single exp O1:  " << time_span_single_exp_O1.count() << std::endl;

  // Time single log O2
  t_start = std::chrono::high_resolution_clock::now();
  for (int32_t i=0; i<samples; ++i) {
    logs_float[i] = NQT_exp_O2(rands_float[i]);
  }
  t_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span_single_exp_O2 = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start);
  std::cout << "single exp O2:  " << time_span_single_exp_O2.count() << std::endl;


  // Time double log
  t_start = std::chrono::high_resolution_clock::now();
  for (int32_t i=0; i<samples; ++i) {
    logs_double[i] = std::exp2(rands_double[i]);
  }
  t_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span_double_exp = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start);
  std::cout << "double exp STL: " << time_span_double_exp.count() << std::endl;

  // Time double log O1
  t_start = std::chrono::high_resolution_clock::now();
  for (int32_t i=0; i<samples; ++i) {
    logs_double[i] = NQT_exp_O1(rands_double[i]);
  }
  t_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span_double_exp_O1 = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start);
  std::cout << "double exp O1:  " << time_span_double_exp_O1.count() << std::endl;

  // Time double log O2
  t_start = std::chrono::high_resolution_clock::now();
  for (int32_t i=0; i<samples; ++i) {
    logs_double[i] = NQT_exp_O2(rands_double[i]);
  }
  t_end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_span_double_exp_O2 = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start);
  std::cout << "double exp O2:  " << time_span_double_exp_O2.count() << std::endl;

  delete rands_float;
  delete rands_double;
  delete logs_float;
  delete logs_double;

  std::cout << std::endl;

  std::cout << "Speedups (NQT vs STL):" << std::endl;
  std::cout << "log, single, O1: " << 100.0*(time_span_single_log.count() / time_span_single_log_O1.count()) - 100.0 << "%" << std::endl; 
  std::cout << "exp, single, O1: " << 100.0*(time_span_single_exp.count() / time_span_single_exp_O1.count()) - 100.0 << "%" << std::endl; 

  std::cout << "log, single, O2: " << 100.0*(time_span_single_log.count() / time_span_single_log_O2.count()) - 100.0 << "%" << std::endl; 
  std::cout << "exp, single, O2: " << 100.0*(time_span_single_exp.count() / time_span_single_exp_O2.count()) - 100.0 << "%" << std::endl; 

  std::cout << "log, double, O1: " << 100.0*(time_span_double_log.count() / time_span_double_log_O1.count()) - 100.0 << "%" << std::endl; 
  std::cout << "exp, double, O1: " << 100.0*(time_span_double_exp.count() / time_span_double_exp_O1.count()) - 100.0 << "%" << std::endl; 

  std::cout << "log, double, O2: " << 100.0*(time_span_double_log.count() / time_span_double_log_O2.count()) - 100.0 << "%" << std::endl; 
  std::cout << "exp, double, O2: " << 100.0*(time_span_double_exp.count() / time_span_double_exp_O2.count()) - 100.0 << "%" << std::endl; 
  
  std::cout << std::endl;

  std::cout << "Speedups (single vs double):" << std::endl;
  std::cout << "log, STL: " << 100.0*(time_span_double_log.count() / time_span_single_log.count()) - 100.0 << "%" << std::endl; 
  std::cout << "exp, STL: " << 100.0*(time_span_double_exp.count() / time_span_single_exp.count()) - 100.0 << "%" << std::endl; 

  std::cout << "log, O1:  " << 100.0*(time_span_double_log_O1.count() / time_span_single_log_O1.count()) - 100.0 << "%" << std::endl; 
  std::cout << "exp, O1:  " << 100.0*(time_span_double_exp_O1.count() / time_span_single_exp_O1.count()) - 100.0 << "%" << std::endl; 

  std::cout << "log, O2:  " << 100.0*(time_span_double_log_O2.count() / time_span_single_log_O2.count()) - 100.0 << "%" << std::endl; 
  std::cout << "exp, O2:  " << 100.0*(time_span_double_exp_O2.count() / time_span_single_exp_O2.count()) - 100.0 << "%" << std::endl; 

  return 0;
}

