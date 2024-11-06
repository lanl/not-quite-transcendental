#------------------------------------------------------------------------------#
# © 2022. Triad National Security, LLC. All rights reserved.  This
# program was produced under U.S. Government contract 89233218CNA000001
# for Los Alamos National Laboratory (LANL), which is operated by Triad
# National Security, LLC for the U.S.  Department of Energy/National
# Nuclear Security Administration. All rights in the program are
# reserved by Triad National Security, LLC, and the U.S. Department of
# Energy/National Nuclear Security Administration. The Government is
# granted for itself and others acting on its behalf a nonexclusive,
# paid-up, irrevocable worldwide license in this material to reproduce,
# prepare derivative works, distribute copies to the public, perform
# publicly and display publicly, and to permit others to do so.
#--------------p----------------------------------------------------------------#

cmake_minimum_required(VERSION 3.14)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_EXPORT_COMPILE_COMMANDS On)

set(TRANSCENDENTAL_VERSION 2.0.0)
project(Transcendental VERSION ${TRANSCENDENTAL_VERSION})

option(USE_CUDA "Use Cuda with Kokkos" OFF)
option(NQT_ORDER_1 "Use first-order NQT rather than quadratic" OFF)
option(NQT_PORTABLE "Use portable frexp/ldexp NQT rather than bithacked" OFF)

# Don't allow in-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
  message(FATAL_ERROR
    "You cannot build in a source directory (or any directory with a CMakeLists.txt file). "
    "Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

set(Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION ON CACHE BOOL
  "Kokkos aggressive vectorization" FORCE)
if (USE_CUDA)
  set(Kokkos_ENABLE_CUDA ON CACHE BOOL
    "Cuda backend for Kokkos" FORCE)
  set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL
    "Enable cuda lambda" FORCE)
endif()
add_subdirectory(kokkos)

add_executable(timelogs
  logs.hpp
  timings.cpp
  )

target_include_directories(timelogs
  PRIVATE
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

target_link_libraries(timelogs
  PUBLIC
  Kokkos::kokkos)

if (NQT_ORDER_1)
  target_compile_definitions(timelogs PUBLIC NQT_ORDER_1)
endif()
if (NQT_PORTABLE)
  target_compile_definitions(timelogs PUBLIC NQT_PORTABLE)
endif()

if (USE_CUDA)
  target_compile_options(
    timelogs
    PUBLIC
    "--expt-relaxed-constexpr"
    )
endif()
# target_compile_options(timelogs PUBLIC "-qopt-report=1 -qopt-report-phase=vec")
