Not-Quite-Transcendental Functions
===

Reference Implementations of Not-Quite-Transcendental Functions

## Introduction

This is the companion repository for [ArXiv:2206.08957](https://arxiv.org/abs/2206.08957). 
Here we present
the reference impementations for not-quite-transendental functions as
well as machinery to time them. For completeness, we also include the
paper in TeX in the paper subdirectory. The logarithms are in
`src/logs.hpp`, where we include a portable version as well as a very
performant version based on integer aliasing.

## Dependencies

To enable performance portability, we utilize the Kokkos library,
included as a submodule. We use `cmake` as our build system.

## Building and running

To compile for CPU, this is enough:
```bash
mkdir bin
cd bin
cmake -DUSE_CUDA=OFF src
make -j
timelogs npoints
```
where here `npoints` is the number of points to use. Use a number large enough to saturate your target architecture. If you want to compile for GPU, change `-DUSE_CUDA` to `ON`.

## Building the paper

If you have `latexmk`, just type `make` in the `paper` directory.

# Report Number

The paper is authorized for unlimited release with LA-UR-22-25573.

## Authors

This code was developed by
- Jonah Miller
- Josh Dolence
- Daniel Holladay

## Copyright

© (or copyright) 2022-2025. Triad National Security, LLC. All rights
reserved.  This program was produced under U.S. Government contract
89233218CNA000001 for Los Alamos National Laboratory (LANL), which is
operated by Triad National Security, LLC for the U.S.  Department of
Energy/National Nuclear Security Administration. All rights in the
program are reserved by Triad National Security, LLC, and the
U.S. Department of Energy/National Nuclear Security
Administration. The Government is granted for itself and others acting
on its behalf a nonexclusive, paid-up, irrevocable worldwide license
in this material to reproduce, prepare derivative works, distribute
copies to the public, perform publicly and display publicly, and to
permit others to do so.

This program is open source under the BSD-3 License.  Redistribution
and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE
