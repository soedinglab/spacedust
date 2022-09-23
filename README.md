# Clustersearch
Discovery of conserved gene clusters in multiple genomes

## Installation

Compiling clustersearch from source has the advantage of system-specific optimizations, which should improve its performance. To compile clustersearch `git`, `g++` (4.9 or higher) and `cmake` (3.0 or higher) are required. Afterwards, the clustersearch binary will be located in the `build/bin` directory.

    git clone https://github.com/soedinglab/clustersearch.git
    cd clustersearch
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
    make -j
    make install
    export PATH=$(pwd)/clustersearch/bin/:$PATH

:exclamation: If you want to compile clustersearch on macOS, please install and use `gcc` from Homebrew. The default macOS `clang` compiler does not support OpenMP and clustersearch will not be able to run multithreaded. Adjust the `cmake` call above to:

    CC="$(brew --prefix)/bin/gcc-10" CXX="$(brew --prefix)/bin/g++-10" cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
