# fastx2bam
Convert fasta/fastq files to bam file (unmapped hopefully PacBio compatible)

# Usage

Convert FASTA/FASTQ files to BAM (unmapped).

## Command Line Usage

    fastx2bam [OPTIONS] <input.fastx> <output.bam>

## Positional Arguments

- **<input.fastx>**: Path to the input FASTA or FASTQ file (can be gzipped).
- **<output.bam>**: Path to the output BAM file.

## Options

- `--threads N`: Number of threads to use with samtools (default: 1).
- `--rename`: Rename query names to simple sequential numbers (e.g., 1, 2, ...).
- `--prefix STR`: Add a prefix to all query names (default: "" [empty string]).
- `--suffix STR`: Add a suffix to all query names (default: "" [empty string]).
- `--header STR`: Optional SAM header (file should be plain text) to use instead of the default.
- `--help`: Show this help message and exit.

## Examples

Convert a FASTQ file without renaming:
    fastx2bam input.fq output.bam

Rename query names in the conversion:
    fastx2bam --rename input.fq output.bam

Add a prefix and suffix to the query names:
    fastx2bam --rename --prefix R --suffix /ccs input.fq output.bam

Use 4 threads with samtools with custom header:
    fastx2bam --threads 4 input.fq output.bam --header header.sam

## Installation
```
make
make install INSTALLATION_DIR=/path/to/installation/dir
fastx2bam --version
```

## Dependencies
### Compilation
- gcc
- zlib

### Running
- samtools
