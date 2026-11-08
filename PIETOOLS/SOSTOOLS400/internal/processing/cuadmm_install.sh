#!/usr/bin/env bash
# CUADMM_INSTALL.SH --- Build cuadmm_exe, the cuADMM executable that
# CUADMM_SEDUMI runs, on Ubuntu or Debian, including Ubuntu under WSL.
#
#   bash cuadmm_install.sh [DIR] [CMAKE_ARG ...]
#
# Description: installs the missing build packages (with sudo apt-get),
# clones cuADMM at commit fcadc30 with its submodule into DIR, applies
# cuadmm_exe.patch from the folder of this script, compiles cuadmm_exe for
# the compute capability of each GPU that nvidia-smi reports, and solves a
# 2 x 2 test SDP with it. A rerun resumes after a failed step.
#
# INPUT:
# DIR:        (optional) path, the folder for the cuADMM sources and build;
#             default ~/cuADMM, where CUADMM_SEDUMI looks for the
#             executable. It must be absent, empty, or made by this script.
# CMAKE_ARG:  (optional) further arguments for the CMake configure step,
#             e.g. -DFETCHCONTENT_SOURCE_DIR_GOOGLETEST=<folder> to use a
#             local googletest in place of the download.
#
# OUTPUT:
# DIR/build/cuadmm_exe, the executable. For DIR other than ~/cuADMM, set
# the environment variable CUADMM_EXE to its path, or pass the path to
# CUADMM_SEDUMI as pars.exe.
#
# NOTES:
# On Windows, run it in the default WSL distribution, which CUADMM_SEDUMI
# uses: in PowerShell, in the folder of this script,
#   wsl bash cuadmm_install.sh
# The GPU needs the NVIDIA driver; under WSL, the Windows driver serves.
# Downloads: cuADMM and its submodule psd_projection (GitHub), googletest
# (fetched by CMake while configuring), and the missing packages; on Ubuntu
# 26.04 the packages take about 2.3 GB, most of it the CUDA toolkit.
#
# For support, contact M. Peet, Arizona State University at mpeet@asu.edu

##################################################################
# Copyright (C)2026  M. Peet
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
##################################################################
#
# If you modify this code, document all changes carefully and include date
# authorship, and a brief description of modifications
#
# 11/08/26 - MMP - initial coding

set -euo pipefail
die() { echo "cuadmm_install.sh: $*" >&2; exit 1; }

# Sources: the cuADMM commit that cuadmm_exe.patch is made against; the
# patch sits next to this script
repo=https://github.com/ComputationalRobotics/cuADMM.git
commit=fcadc30fee799464b603419cc21eb05c4771cad2
patch="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/cuadmm_exe.patch"
[ -f "$patch" ] || die "no cuadmm_exe.patch next to this script ($patch)."

# DIR: the first argument, unless it is a CMake argument; the rest go to CMake
case "${1:-}" in
    ''|-*)  dir="$HOME/cuADMM" ;;
    *)      dir="$1"; shift ;;
esac
case "$dir" in /*) ;; *) dir="$PWD/$dir" ;; esac
# DIR must be absent, empty, or made by this script, which leaves the marker
# .git/cuadmm_install in it: no other tree is changed
if [ ! -e "$dir/.git/cuadmm_install" ] && [ -n "$(ls -A "$dir" 2>/dev/null)" ]; then
    die "$dir is neither an empty folder nor one made by this script; pass another folder as the first argument."
fi

# Platform: Linux or WSL with apt; the user whose home folder CUADMM_SEDUMI
# searches, so not under sudo; an NVIDIA GPU, whose compute capability
# nvidia-smi reports as e.g. 8.9 (read as 89)
[ "$(uname -s)" = Linux ] || die "run this script under Linux or WSL; on Windows, in PowerShell, in its folder: wsl bash cuadmm_install.sh"
command -v apt-get >/dev/null || die "no apt-get: on distributions other than Ubuntu and Debian, follow the build steps at the top of $patch."
[ -z "${SUDO_USER:-}" ] || die "run this script without sudo; it calls sudo itself to install packages."
command -v nvidia-smi >/dev/null || die "no nvidia-smi: install the NVIDIA driver (for WSL, the Windows driver)."
cc=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | tr -d '. \r' | sort -un) \
    || die "nvidia-smi reports no compute capability; is a GPU visible?"
[ -n "$cc" ] || die "nvidia-smi reports no GPU."

# Packages: install those missing. An nvcc from NVIDIA's own CUDA packages,
# in /usr/local/cuda/bin, is not on the default PATH, and serves as well
if ! command -v nvcc >/dev/null && [ -x /usr/local/cuda/bin/nvcc ]; then
    PATH="/usr/local/cuda/bin:$PATH"
fi
pkgs="git ca-certificates build-essential cmake libsuitesparse-dev libopenblas-dev liblapack-dev"
command -v nvcc >/dev/null || pkgs="$pkgs nvidia-cuda-toolkit"
missing=""
for p in $pkgs; do
    s=$(dpkg-query -W -f='${Status}' "$p" 2>/dev/null || true)
    case "$s" in *' installed') ;; *) missing="$missing $p" ;; esac
done
if [ -n "$missing" ]; then
    echo "Installing the missing packages:$missing"
    sudo apt-get update
    sudo apt-get install -y --no-install-recommends $missing
fi

# Architectures: native code and PTX for the compute capability of each GPU.
# For one that nvcc does not know (a GPU newer than the toolkit), PTX for the
# newest architecture nvcc knows, which the driver compiles at load
known=" $(nvcc --list-gpu-arch 2>/dev/null | sed -n 's/^compute_\([0-9]*\)$/\1/p' | sort -n | tr '\n' ' ' || true)"
flags=""
for a in $cc; do
    case "$known" in
        ' '|*" $a "*)
            flags="$flags -gencode arch=compute_$a,code=sm_$a -gencode arch=compute_$a,code=compute_$a" ;;
        *)
            b=""
            for k in $known; do [ "$k" -ge "$a" ] || b=$k; done
            [ -n "$b" ] || die "nvcc compiles for no architecture at or below compute capability $a."
            echo "nvcc does not know compute_$a: compiling PTX for compute_$b, which the driver compiles for the GPU at load."
            flags="$flags -gencode arch=compute_$b,code=compute_$b" ;;
    esac
done

# Sources: cloned into DIR on the first run, which leaves the marker; a rerun
# reuses the clone
if [ ! -e "$dir/.git/cuadmm_install" ]; then
    git clone "$repo" "$dir"
    touch "$dir/.git/cuadmm_install"
fi
cd "$dir"
dir=$(pwd)
git -c advice.detachedHead=false checkout -q "$commit"
git submodule update --init

# Patch: applied once. Carriage returns are stripped first: the patch has
# none, and a CRLF copy (a Windows checkout without .gitattributes) does not
# apply
if tr -d '\r' < "$patch" | git apply --reverse --check - 2>/dev/null; then
    echo "cuadmm_exe.patch is applied already."
else
    tr -d '\r' < "$patch" | git apply -
fi

# Build: Release, for the architectures above, with at most 8 parallel jobs,
# the number tested. Configuring fetches googletest, which only the test
# target uses
j=$(nproc);     [ "$j" -le 8 ] || j=8
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CUDA_FLAGS="${flags# }" "$@"
cmake --build build --target cuadmm_exe -j "$j"

# Test: min <C,X> s.t. X11 = X22 = 1, X 2 x 2 positive semidefinite, with
# C = [0 1; 1 0]; the optimal value is -2. The files are in cuADMM's input
# format, as CUADMM_SEDUMI writes them: svec coordinates (X11, sqrt(2)*X12,
# X22), 0-based COO lines "row column value"
t=$(mktemp -d)
printf 's 2\n'                      > "$t/blk.txt"
printf '2\n'                        > "$t/con_num.txt"
printf '0 0 1\n2 1 1\n'             > "$t/At.txt"
printf '0 0 1\n1 0 1\n'             > "$t/b.txt"
printf '1 0 1.4142135623730951\n'   > "$t/C.txt"
# Under WSL the WSL driver library goes first: the libcuda.so.1 installed
# with nvidia-cuda-toolkit otherwise shadows it, and every CUDA call fails
lp="${LD_LIBRARY_PATH:-}"
[ ! -d /usr/lib/wsl/lib ] || lp="/usr/lib/wsl/lib${lp:+:$lp}"
LD_LIBRARY_PATH="$lp" build/cuadmm_exe "$t/" 1e-6 10000 1 1 > "$t/log" 2>&1 || true
# Success: the line CUADMM_TIMING (the exit status is 0 even when CUDA fails)
# and the objective within 1e-4 of -2
obj=$(sed -n 's/^ *primal objective = *\([^ ]*\).*/\1/p' "$t/log" | tail -n 1)
if grep -q CUADMM_TIMING "$t/log" && awk -v o="$obj" 'BEGIN { exit !(o != "" && o+2 < 1e-4 && o+2 > -1e-4) }'; then
    rm -rf "$t"
else
    tail -n 20 "$t/log" >&2
    die "the test SDP failed; its log is $t/log."
fi

# Report: the executable, and where CUADMM_SEDUMI finds it
exe="$dir/build/cuadmm_exe"
echo "cuadmm_exe solved the test SDP (objective $obj; exact -2): $exe"
if [ "$dir" != "$HOME/cuADMM" ]; then
    echo "CUADMM_SEDUMI looks for ~/cuADMM/build/cuadmm_exe: set the environment variable CUADMM_EXE to $exe (in MATLAB, setenv('CUADMM_EXE','$exe')), or pass it as pars.exe."
fi
echo "In SOSTOOLS, options.solver = 'cuadmm' selects it in sossolve."
