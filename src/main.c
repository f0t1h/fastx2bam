/* The MIT License

Copyright (c) 2025 Fatih Karaoglanoglu

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#define F2B_VERSION "0.1"

#ifdef _WIN32
#define SAMTOOLS_CHECK_CMD "where samtools > nul 2>&1"
#else
#define SAMTOOLS_CHECK_CMD "command -v samtools > /dev/null 2>&1"
#endif

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)
#define UNUSED(x) (void)(x)

#define QUALBUFFER 16384
#define MIN(a, b) ((a) < (b) ? (a) : (b))

void print_sam(FILE *stream, kseq_t *seq, char *mock_qual, int rename, int count, const char *prefix, const char *suffix) {
    if (rename) {
        fprintf(stream, "%s%d%s\t4\t*\t0\t255\t*\t*\t0\t%ld\t%.*s\t",
                prefix, count, suffix, seq->seq.l, (int)seq->seq.l, seq->seq.s);
    } else {
        fprintf(stream, "%s%s%s\t4\t*\t0\t255\t*\t*\t0\t%ld\t%.*s\t",
                prefix, seq->name.s, suffix, seq->seq.l, (int)seq->seq.l, seq->seq.s);
    }

    if (seq->qual.l > 0) {
        fprintf(stream, "%s\n", seq->qual.s);
    } else {
        int so_far;
        int remaining = seq->seq.l;
        while (remaining > 0) {
            fprintf(stream, "%.*s%n", MIN(remaining, QUALBUFFER), mock_qual, &so_far);
            remaining -= so_far;
        }
        fputc('\n', stream);
    }
}


void print_help(int argc, char **argv, FILE *stream){
    UNUSED(argc);
    // fprintf(stream, "Usage: %s [--threads N] [--rename] [--prefix P] [--suffix S] <input.fastx> <output.bam>\n", argv[0]); 
    fprintf(stream, "Usage:\n  %s [OPTIONS] <input.fastx> <output.bam>\n", argv[0]);
    fprintf(stream, "\nDescription:\n");
    fprintf(stream, "  Converts FASTQ/FASTA input to BAM using samtools view. Supports\n");
    fprintf(stream, "  renaming query names, adding prefixes/suffixes, and setting thread count.\n\n");

    fprintf(stream, "Positional arguments:\n");
    fprintf(stream, "  <input.fastx>     Path to input FASTA or FASTQ file (can be gzipped).\n");
    fprintf(stream, "  <output.bam>      Path to output BAM file.\n\n");

    fprintf(stream, "Options:\n");
    fprintf(stream, "  --threads N       Number of threads to use with samtools (default: 1).\n");
    fprintf(stream, "  --rename          Rename query names to simple sequential integers (1, 2, ...).\n");
    fprintf(stream, "  --prefix STR      Add a prefix to all query names (default: \"\").\n");
    fprintf(stream, "  --suffix STR      Add a suffix to all query names (default: \"\").\n");
    fprintf(stream, "  --header STR      Optional sam header (file should be plain text) to use instead of the default.\n");
    fprintf(stream, "  --help            Show this help message and exit.\n");
    fprintf(stream, "  --version         Version string.\n\n");
    
    fprintf(stream, "Examples:\n");
    fprintf(stream, "  %s input.fq output.bam\n", argv[0]);
    fprintf(stream, "  %s --rename input.fq output.bam\n", argv[0]);
    fprintf(stream, "  %s --rename --prefix R --suffix /ccs input.fq output.bam\n", argv[0]);
    fprintf(stream, "  %s --threads 4 input.fq output.bam\n", argv[0]);
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        print_help(argc, argv, stderr);
        return 1;
    }

    int num_threads = 1;
    int rename_flag = 0;
    const char *prefix = "";
    const char *suffix = "";
    const char *header = "";
    const char *input_path = NULL;
    const char *output_path = NULL;

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h")==0 ){
            print_help(argc, argv, stdout);
            return 0;
        }
        else if (strcmp(argv[i], "--version") == 0 ){
            puts(F2B_VERSION);
            return 0;
        }
        else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) {
            num_threads = atoi(argv[++i]);
            if (num_threads < 1) {
                fprintf(stderr, "Invalid thread count: %s\n", argv[i]);
                return 1;
            }
        } else if (strcmp(argv[i], "--rename") == 0) {
            rename_flag = 1;
        } else if (strcmp(argv[i], "--prefix") == 0 && i + 1 < argc) {
            prefix = argv[++i];
        } else if (strcmp(argv[i], "--suffix") == 0 && i + 1 < argc) {
            suffix = argv[++i];
        } else if (strcmp(argv[i], "--header") == 0 && i + 1 < argc) {
            header = argv[++i];
        } else if (!input_path) {
            input_path = argv[i];
        } else if (!output_path) {
            output_path = argv[i];
        } else {
            fprintf(stderr, "Unexpected argument: %s\n", argv[i]);
            return 1;
        }
    }

    if (!input_path || !output_path) {
        fprintf(stderr, "Error: missing input or output file.\n");
        fprintf(stderr, "Usage: %s [--threads N] [--rename] [--prefix P] [--suffix S] <input.fastx> <output.bam>\n", argv[0]);
        return 1;
    }

    if (system(SAMTOOLS_CHECK_CMD) != 0) {
        fprintf(stderr, "Error: 'samtools' not found in PATH.\n");
        return 1;
    }

    gzFile in = gzopen(input_path, "rb");
    if (!in) {
        fprintf(stderr, "Error opening input file: %s\n", input_path);
        return 1;
    }

    kseq_t *seq = kseq_init(in);

    char cmd[2048];
    snprintf(cmd, sizeof(cmd), "samtools view -@ %d -b /dev/stdin > \"%s\"", num_threads, output_path);
    FILE *pipe = popen(cmd, "w");
    if (!pipe) {
        fprintf(stderr, "Error opening output pipe to samtools\n");
        kseq_destroy(seq);
        gzclose(in);
        return 1;
    }

    char mock_qual[QUALBUFFER + 1];
    for (int i = 0; i < QUALBUFFER; ++i) mock_qual[i] = '@';
    mock_qual[QUALBUFFER] = '\0';

    if(strcmp(header, "") == 0) {
        const char *fake_hd_line = "@HD\tVN:1.5\tSO:unknown\tpb:3.0.1\n";
        const char *fake_rg_line =
            "@RG\tID:4f25f78c\tPL:PACBIO\t"
            "DS:READTYPE=CCS;BINDINGKIT=101-820-500;SEQUENCINGKIT=101-826-100;"
            "BASECALLERVERSION=5.0.0;FRAMERATEHZ=100.000000\tLB:SQlle "
            "zeBAM\tPU:m64187e_211217_130958\t"
            "PM:SEQUELII\tCM:S/P4.1-C2/5.0-8M\n";
        const char *fake_pg_line =
            "@PG\tID:ccs-6.0.0\tPN:ccs\tVN:6.0.0\tDS:Generate circular "
            "consensus sequences (ccs) from subreads.\tCL:ccs ...\n";
        fputs( fake_hd_line, pipe);
        fputs( fake_rg_line, pipe);
        fputs( fake_pg_line, pipe);
    }
    else{
        FILE *hf = fopen(header, "r");
        if (hf == NULL) {
            fprintf(stderr, "Error: Could not open header file %s\n", header);
            return -1;
        }
        char *line = NULL;
        size_t len = 0;
        ssize_t read;
        while ((read = getline(&line, &len, hf)) != -1) {
            fputs(line, pipe);
        }
        free(line);
        fclose(hf);
    }

    int read_count = 0;
    while (kseq_read(seq) >= 0) {
        ++read_count;
        print_sam(pipe, seq, mock_qual, rename_flag, read_count, prefix, suffix);
    }

    kseq_destroy(seq);
    gzclose(in);
    pclose(pipe);

    return 0;
}
