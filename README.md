# Clustersearch: Discovery of conserved gene clusters in multiple genomes

Clustersearch is a modular toolkit for identification of conserved gene clusters among multiple genomes based on homology and conservation of gene neighborhood. Clustersearch adapts the fast and sensitive structure comparisons of [Foldseek](https://github.com/steineggerlab/foldseek) and homology search capabilities of [MMseqs2](https://github.com/soedinglab/MMseqs2). It introduces a novel approach of aggregating sets of homologous hits between pairs of genomes and identifies cluster of hits with conserved gene neighborhood between each using agglomerative hierarchical clustering algorithm. Clustersearch is GPLv3-licensed open source software implemented in C++ and available for Linux and macOS. The software is designed to run efficiently on multiple cores.

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

## Input

The starting point is either (1) FASTA (`.fna` or `.fna.gz`) files of nulceotide genome or contig sequences coupled with protein-coding sequences (**CDS**) annotated in GFF3 (`.gff3`) format (typically by running a gene prediction program like [Prodigal](https://github.com/hyattpd/Prodigal)). **Note**: this only works with (prokaryotic) genomes without intron/exon structures or (2) FASTA (`.faa` or `.faa.gz`) files of protein sequences with header format from [Prodigal](https://github.com/hyattpd/Prodigal)
 <!-- for prokaryotes or [metaeuk](https://github.com/soedinglab/metaeuk) for eukaryotes.  -->
Input genomes are supplied as separate FASTA and GFF3 files (one genome per file). 

## Dependencies

To enable structure comparisons, clustersearch requires the installation of [Foldseek](https://github.com/steineggerlab/foldseek) in the working directory.
 <!-- (the binary file `/foldseek/build/bin/foldseek` should exist in the working directory). -->

## Running Clustersearch

### Main Modules

<!-- * `easy-clustersearch`     search for conserved gene clusters between genomes (fasta/db) -->
* `createsetdb`       create sequence database from FASTA input
* `clustersearch`   search for conserved gene clusters between genomes (using PSI-BLAST like iterative searches)
* `aa2foldseek` convert a sequence database to structure database by mapping to Uniprot/Alphafold

### Important parameters

    # createsetdb
    --gff-type                             Type in the GFF file to filter by
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

    clustersearch easy-clustersearch examples/*.fna targetSetDB clusterResult tmpFolder --gff-dir gffDir.txt --gff-type CDS
    clustersearch easy-clustersearch examples/*.faa targetSetDB clusterResult tmpFolder -->

### Creating databases

To start, you need to create a database of the input genomes `setDB`. Before search, query or target sequences contained in FASTA files need to be converted to database format by calling `createsetdb`. This command first creates a sequence DB with mapped strand and genomic coordinates extracted from GFF3 or header, and finally generates associated metadata. For nucleotide input, the respective GFF3 files should be given using argument `--gff-dir` as a list of paths to the GFF3 files.

    clustersearch createsetdb genome1.fna [...genomeN.fna] setDB tmpFolder --gff-dir gffDir.txt --gff-type CDS
    clustersearch createsetdb genome1.faa [...genomeN.faa] setDB tmpFolder

To enable protein structure search with Foldseek, the protein sequences are mapped to Foldseek structure sequence DB like AlphaFoldDB. This requires pre-downloading the reference FoldseekDB.

    # Download reference FoldseekDB
    path/to/foldseek databases Alphafold/UniProt-NO-CA refFoldseekDB tmpFolder

    # Convert to structure sequence DB
    clustersearch aa2foldseek setDB refFoldseekDB outDB tmpFolder

### Clustersearch (with MMseqs2 or Foldseek)

Clustersearch will first conduct an all-against-all homology search/structure comparison between two sets of protein-coding genes derived from multiple genomes, and then find clusters of homologous hits based on conservation of gene neighborhood. Structure comparison with Foldseek is invoked by `--search-mode 1`. For a more sensitive search, iterative searches in MMseqs2 and Foldseek can be done by setting `--num-iterations`.

    # Search querySetDB against targetSetDB (using MMseqs)
    clustersearch clustersearch querySetDB targetSetDB resultDB tmpFolder

    # Search querySetDB against targetSetDB turned into profile
    clustersearch clustersearch querySetDB targetSetDB resultDB tmpFolder --profile-cluster-search

    # Iterative cluster search (like PSI-BLAST) with 2 iterations
    clustersearch clustersearch querySetDB targetSetDB resultDB tmpFolder --num-iterations 2

    # Search querySetDB against targetSetDB (using Foldseek)
    clustersearch clustersearch querySetDB targetSetDB resultDB tmpFolder --search-mode 1

### The Clustersearch output

Upon completion, clustersearch outputs two files: an **alignment result file** for clusters of hits (in MMseqs internal alignment format), and a **header file** (`_h`) describing each cluster. Each cluster entry is separated by a \0 byte. Each line in the alignment file describes an individual hit, i.e. one target sequence aligned to the query, which has the following columns separated by tab characters:

    queryID  targetID  bestHitPvalue  seqIdentity  eVal  qStart  qEnd  qLen  tStart  tEnd  tLen  alnCigar

The rest of the columns are described in the [MMseqs2 wiki](https://github.com/soedinglab/MMseqs2/wiki#internal-alignment-format).

Each line in the header file describes an individual cluster of hits. The columns of header file: query genome ID, target genome ID, cluster match P-value, multihit E-value and number of hits in the cluster.

### Removing temporary files

During the workflow execution, clustersearch will keep all intermediate outputs in `tmpFolder`, passing the `--remove-tmp-files` parameter will clear out the `tmpFolder` after workflows have finished.
