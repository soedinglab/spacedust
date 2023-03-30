# Spacedust: Discovery of conserved gene clusters in multiple genomes

Spacedust is a modular toolkit for identification of conserved gene clusters among multiple genomes based on homology and conservation of gene neighborhood. Spacedust adapts the fast and sensitive structure comparisons of [Foldseek](https://github.com/steineggerlab/foldseek) and homology search capabilities of [MMseqs2](https://github.com/soedinglab/MMseqs2). It introduces a novel approach of aggregating sets of homologous hits between pairs of genomes and identifies cluster of hits with conserved gene neighborhood between each using agglomerative hierarchical clustering algorithm. Spacedust is GPLv3-licensed open source software implemented in C++ and available for Linux and macOS. The software is designed to run efficiently on multiple cores.

<p align="center"><img src="https://github.com/soedinglab/spacedust/blob/master/.github/spacedust.png" height="350"/></p>

## Installation

Spacedust can be used by downloading a [statically compiled version](https://mmseqs.com/spacedust/) or [compiling from source](#compile-from-source). It requires a 64-bit system (check with `uname -a | grep x86_64`) with at least the SSE4.1 instruction set (check by executing `cat /proc/cpuinfo | grep sse4_1` on Linux or `sysctl -a | grep machdep.cpu.features | grep SSE4.1` on MacOS).

    # static Linux AVX2 build (check using: cat /proc/cpuinfo | grep avx2)
    wget https://mmseqs.com/spacedust/spacedust-linux-avx2.tar.gz; tar xvzf spacedust-linux-avx2.tar.gz; export PATH=$(pwd)/spacedust/bin/:$PATH
    # static Linux SSE4.1 build (check using: cat /proc/cpuinfo | grep sse4_1)
    wget https://mmseqs.com/spacedust/spacedust-linux-sse41.tar.gz; tar xvzf spacedust-linux-sse41.tar.gz; export PATH=$(pwd)/spacedust/bin/:$PATH
    # static macOS build (universal binary with SSE4.1/AVX2/M1 NEON)
    wget https://mmseqs.com/spacedust/spacedust-osx-universal.tar.gz; tar xvzf spacedust-osx-universal.tar.gz; export PATH=$(pwd)/spacedust/bin/:$PATH

Other precompiled binaries for ARM64, PPC64LE amd SSE2 are available at [https://mmseqs.com/spacedust](https://mmseqs.com/spacedust).

## Input

The starting point is either (1) FASTA (`.fna` or `.fna.gz`) files of nulceotide genome or contig sequences coupled with protein-coding sequences (**CDS**) annotated in GFF3 (`.gff3`) format (typically by running a gene prediction program like [Prodigal](https://github.com/hyattpd/Prodigal)). **Note**: this only works with (prokaryotic) genomes without intron/exon structures or (2) FASTA (`.faa` or `.faa.gz`) files of protein sequences with header format from [Prodigal](https://github.com/hyattpd/Prodigal).
 <!-- for prokaryotes or [metaeuk](https://github.com/soedinglab/metaeuk) for eukaryotes.  -->
Input genomes are supplied as separate FASTA and GFF3 files (one genome per file). 

## Dependencies

To enable structure comparisons, spacedust requires the installation of [Foldseek](https://github.com/steineggerlab/foldseek) in the working directory.
 <!-- (the binary file `/spacedust/build/bin/foldseek` should exist in the working directory). -->

## Running Spacedust

### Main Modules

<!-- * `easy-clustersearch`     search for conserved gene clusters between genomes (fasta/db) -->
* `createsetdb`       create sequence database from FASTA input
* `clustersearch`   search for conserved gene clusters between genomes (using PSI-BLAST like iterative searches)
* `aa2foldseek` convert a sequence database to structure database by mapping to Uniprot/Alphafold

### Important parameters

    # createsetdb
    --gff-type                             Type of feature in the GFF file to filter by (default: "", all features)
    --gff-dir                              Path to gff dir file

    # clustersearch
    --search-mode                          0: sequence search with MMseqs2, 1: structure comparison with Foldseek (default:0)
    --num-iterations                       Number of iterative profile search iterations (default:1)
    --profile-cluster-search               Perform profile(target)-sequence searches
    --filter-self-match                    Remove hits between the same set
    --max-gene-gap                         Maximum number of genes allowed between two clusters to merge (default:3)
    --cluster-size                         Minimum number of genes to define cluster (default:2)

### Quick start
<!-- The `easy-clustersearch` workflow combines the clustersearch modules into a single step: createsetdb, aa2foldseek and (iterative)clustersearch.

    spacedust easy-clustersearch examples/*.fna targetSetDB clusterResult tmpFolder --gff-dir gffDir.txt --gff-type CDS
    spacedust easy-clustersearch examples/*.faa targetSetDB clusterResult tmpFolder -->

### Creating databases

To start, you need to create a database of the input genomes `setDB`. Before search, query or target sequences contained in FASTA files need to be converted to database format by calling `createsetdb`. This command first creates a sequence DB with mapped strand and genomic coordinates extracted from GFF3 or header, and finally generates associated metadata. For nucleotide input, the respective GFF3 files should be given using argument `--gff-dir` as a list of paths to the GFF3 files.

    spacedust createsetdb genome1.fna [...genomeN.fna] setDB tmpFolder --gff-dir gffDir.txt --gff-type CDS
    spacedust createsetdb genome1.faa [...genomeN.faa] setDB tmpFolder

To enable protein structure search with Foldseek, the protein sequences are mapped to Foldseek structure sequence DB like AlphaFoldDB. This requires pre-downloading the reference FoldseekDB.

    # Download reference FoldseekDB
    path/to/foldseek databases Alphafold/UniProt refFoldseekDB tmpFolder

    # Convert to structure sequence DB
    spacedust aa2foldseek setDB refFoldseekDB tmpFolder

### Spacedust (with MMseqs2 and/or Foldseek)

Spacedust will first conduct an all-against-all homology search/structure comparison between two sets of protein-coding genes derived from multiple genomes, and then find clusters of homologous hits based on conservation of gene neighborhood. Structure comparison with Foldseek is invoked by `--search-mode 1`. The sequences which could be mapped to a structure by `aa2foldseek` will be searched with Foldseek, and the rest will be searched with MMseqs2. For a more sensitive search, iterative searches in MMseqs2 and Foldseek can be done by setting `--num-iterations`.

    # Search querySetDB against targetSetDB (using MMseqs)
    spacedust clustersearch querySetDB targetSetDB result.tsv tmpFolder

    # Search querySetDB against targetSetDB turned into profile
    spacedust clustersearch querySetDB targetSetDB result.tsv tmpFolder --profile-cluster-search

    # Iterative cluster search (like PSI-BLAST) with 2 iterations
    spacedust clustersearch querySetDB targetSetDB result.tsv tmpFolder --num-iterations 2

    # Search querySetDB against targetSetDB (using Foldseek and MMseqs)
    spacedust clustersearch querySetDB targetSetDB result.tsv tmpFolder --search-mode 1

### The Spacedust output

Upon completion, Spacedust outputs a tab-separated text file (`.tsv`). Each reported cluster consist of one summary line followed by multiple lines, one line for each pairwise hit between the query and target genome.

    #clusterID  query_acc  target_acc   clusterMatchPvalue multihitPvalue  num_hits
    >queryID    targetID    bestHitPvalue   seqIdentity eVal    qStart  qEnd    qLen    tStart  tEnd    tLen    alnCigar

The summary line starts with `#`: a unique cluster ID, query genome accession, target genome accession, cluster match P-value (joint P-value of clustering and ordering), multihit P-value and number of hits in the cluster.

Each following line starts with `>` and describes an individual member hit of the cluster (i.e. one target sequence aligned to the query) in MMseqs2 alignment-result-like format, which has the following columns separated by tab characters: query protein ID, target protein ID, besthit P-value, sequence identity, pairwise E-value, query protein start, end and length, target protein start, end and length, alnCigar, which is a string that encodes the alignment in compressed format. The columns are further described in the [MMseqs2 wiki](https://github.com/soedinglab/MMseqs2/wiki#internal-alignment-format). The query and target protein IDs contain the accession, protein position index, start and end coordinates (end and start coordinates if on reverse strand), separated by `_`.

### Removing temporary files

During the workflow execution, spacedust will keep all intermediate outputs in `tmpFolder`, passing the `--remove-tmp-files` parameter will clear out the `tmpFolder` after workflows have finished.

## Compile from source

Compiling spacedust from source has the advantage of system-specific optimizations, which should improve its performance. To compile spacedust `git`, `g++` (4.9 or higher) and `cmake` (3.0 or higher) are required. Afterwards, the spacedust binary will be located in the `build/bin` directory.

    git clone https://github.com/soedinglab/spacedust.git
    cd spacedust
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
    make -j
    make install
    export PATH=$(pwd)/bin/:$PATH

:exclamation: If you want to compile spacedust on macOS, please install and use `gcc` from Homebrew. The default macOS `clang` compiler does not support OpenMP and spacedust will not be able to run multithreaded. Adjust the `cmake` call above to:

    CC="$(brew --prefix)/bin/gcc-10" CXX="$(brew --prefix)/bin/g++-10" cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
