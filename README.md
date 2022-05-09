# Introduction

Scallop2 is a reference-based transcript assembler
specifically optimized for paired-/multiple-end RNA-seq data.
The development of Scallop2 has been based on the [Scallop](https://github.com/Kingsford-Group/scallop) assembler.
Scripts and documentation that evaluate Scallop2 is available at [scallop2-test](https://github.com/Shao-Group/scallop2-test).

# Installation

Scallop2 can be easily installed with conda: [![Anaconda-Server Badge](https://anaconda.org/bioconda/scallop2/badges/installer/conda.svg)](https://anaconda.org/bioconda/scallop2)

Install from source code: download the source code of Scallop2 from
[here](https://github.com/Shao-Group/scallop2/releases/download/v1.1.2/scallop2-1.1.2.tar.gz).
Scallop2 uses additional libraries of Boost and htslib. 
If they have not been installed in your system, you first
need to download and install them. You might also need to
export the runtime library path to certain environmental
variable (for example, `LD_LIBRARY_PATH`, for most linux distributions).
After install these dependencies, you then compile the source code of Scallop2.
If some of the above dependencies are not installed to the default system 
directories (for example, `/usr/local`, for most linux distributions),
their corresponding installing paths should be specified to `configure` of Scallop2.
The entire installation typically takes a few minutes to complete.

## Download Boost
If Boost has not been downloaded/installed, download Boost
[(license)](http://www.boost.org/LICENSE_1_0.txt) from (http://www.boost.org).
Uncompress it somewhere (compiling and installing are not necessary).

## Install htslib
If htslib has not been installed, download htslib 
[(license)](https://github.com/samtools/htslib/blob/develop/LICENSE)
from (http://www.htslib.org/) with version 1.5 or higher.
Note that htslib relies on zlib. So if zlib has not been installed in your system,
you need to install zlib first. To do so, download zlib
[(license)](https://zlib.net/zlib_license.html) at (https://zlib.net/).
Use the following commands to install zlib:
```
./configure
make
make install
```
After installing zlib, use the following commands to build htslib:
```
./configure --disable-bz2 --disable-lzma --disable-gcs --disable-s3 --enable-libcurl=no
make
make install
```
The default installation location of htslib is `/usr/lib`.
If you would install it to a different location, replace the above `configure` line with
the following (by adding `--prefix=/path/to/your/htslib` to the end):
```
./configure --disable-bz2 --disable-lzma --disable-gcs --disable-s3 --enable-libcurl=no --prefix=/path/to/your/htslib
```
In this case, you also need to export the runtime library path (note that there
is an additional `lib` following the installation path):
```
export LD_LIBRARY_PATH=/path/to/your/htslib/lib:$LD_LIBRARY_PATH
```

## Build Scallop2

Use the following to compile Scallop2:
```
./configure --with-htslib=/path/to/your/htslib --with-boost=/path/to/your/boost
make
```

If some of the dependencies are installed in the default system directory (for example, `/usr/lib`),
then the corresponding `--with-` option might not be necessary.
The executable file `scallop2` will appear at `src/scallop2`.


# Usage

The usage of `scallop2` is:
```
./scallop2 -i <input.bam> -o <output.gtf> [options]
```

The `input.bam` is the read alignment file generated by some RNA-seq aligner, (for example, STAR or HISAT2).
Make sure that it is sorted; otherwise run `samtools` to sort it:
```
samtools sort input.bam > input.sort.bam
```

The reconstructed transcripts shall be written as gtf format into `output.gtf`.

Scallop2 support the following parameters. Please refer
to additional explanations below the table.

 Parameters | Default Value | Description
 ------------------------- | ------------- | ----------
 --help  | | print usage of Scallop2 and exit
 --version | | print version of Scallop2 and exit
 --preview | | show the inferred `library_type` and exit
 --verbose | 1 | chosen from {0, 1, 2}
 -f/--transcript_fragments    | | file to which the assembled non-full-length transcripts will be written to
 --library_type               | empty | chosen from {empty, unstranded, first, second}
 --assemble_duplicates		  | 10 | the number of consensus runs of the decomposition
 --min_transcript_coverage    | 1.5 | the minimum coverage required to output a multi-exon transcript
 --min_single_exon_coverage   | 20 | the minimum coverage required to output a single-exon transcript
 --min_transcript_length_base      |150 | the minimum base length of a transcript
 --min_transcript_length_increase  | 50 | the minimum increased length of a transcript with each additional exon
 --min_mapping_quality        | 1 | ignore reads with mapping quality less than this value
 --max_num_cigar              | 1000 | ignore reads with CIGAR size larger than this value
 --min_bundle_gap             | 100 | the minimum distances required to start a new bundle
 --min_num_hits_in_bundle     | 5 | the minimum number of reads required in a bundle
 --min_flank_length           | 3 | the minimum match length required in each side for a spliced read
 --min_splice_bundary_hits    | 1 | the minimum number of spliced reads required to support a junction

1. For `--verbose`, 0: quiet; 1: one line for each splice graph; 2: details of graph decomposition.

2. `--library_type` is highly recommended to provide. The `unstranded`, `first`, and `second`
correspond to `fr-unstranded`, `fr-firststrand`, and `fr-secondstrand` used in standard Illumina
sequencing libraries. If none of them is given, i.e., it is `empty` by default, then Scallop2
will try to infer the `library_type` by itself (see `--preview`). Notice that such inference is based
on the `XS` tag stored in the input `bam` file. If the input `bam` file do not contain `XS` tag,
then it is essential to provide the `library_type` to Scallop2. You can try `--preview` to see
the inferred `library_type`.

3. `--min_transcript_coverage` is used to filter lowly expressed transcripts: Scallop2 will filter
out transcripts whose (predicted) raw counts (number of moleculars) is less than this number.

4. `--min_transcript_length_base` and `--min_transcript_length_increase` is combined to filter
short transcripts: the minimum length of a transcript is given by `--min_transcript_length_base`
\+ `--min_transcript_length_increase` * num-of-exons-in-this-transcript. Transcripts that are less
than this number will be filtered out.

# Examples 
The alignments of a set of 10 Illumina paired-end RNA-seq samples
are available at [doi:10.26208/8c06-w247](https://doi.org/10.26208/8c06-w247).
Scallop2 usually take a few minteus to an hour to assemble a typical aligned Illumina RNA-seq sample.

## Running Scallop2 on a small example
A small example of input data `example-input.bam` is available in the `example` directory.

Suppose we have installed Scallop2 following the steps in the `Installation` section, we have the executable file `scallop2` at `src/scallop2`.

Commands to enter `example` directory and run Scallop2 using `example-input.bam` as input:
```
cd ./example
../src/scallop2 -i example-input.bam -o example-output.gtf
```

An output file named `example-output.gtf` will appear in the `example` directory.
The output file stores the reconstructed transcripts assembled by Scallop2 in GTF format. 

The Gene transfer format (GTF) is a file format used to hold information about gene structure, detailed documentation about [GTF format](https://useast.ensembl.org/info/website/upload/gff.html) is available from Ensembl.

We can use [IGV](https://software.broadinstitute.org/software/igv/home) to visualize the input reads (BAM file) and the output assembled transcripts (GTF file).
A snapshot of IGV visualization of `example-input.bam` and `example-output.gtf` at gene loci chr1:46,303,377-46,321,251 is avaliable at `example/example-IGV-snapshot.png`.
