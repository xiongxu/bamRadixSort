/*  bam_sort.c -- sorting and merging.
gcc -g -O3 -Wall -c bam_radixSort.c -o bam_radixSort.o -I../htslib
gcc -g -O3 -Wall bam_radixSort.o radixSort.o sam_opts.o sam_utils.o ../htslib/libhts.a -o bam_radixSort -lz -lpthread

    Copyright (C) 2008-2016 Genome Research Ltd.
    Portions copyright (C) 2009-2012 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>
    Author: Martin Pollard <mp15@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <stdbool.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <unistd.h>
#include <getopt.h>
#include <assert.h>
#include "htslib/sam.h"
#include "sam_opts.h"
#include "samtools.h"
#include "radixSort.h"

const size_t SORT_DEFAULT_MEGS_PER_THREAD = 5;
#define SORT_LIST_SIZE 1<<20

static int change_SO(bam_hdr_t *h, const char *so)
{
    char *p, *q, *beg = NULL, *end = NULL, *newtext;
    if (h->l_text > 3) {
        if (strncmp(h->text, "@HD", 3) == 0) {
            if ((p = strchr(h->text, '\n')) == 0) return -1;
            *p = '\0';
            if ((q = strstr(h->text, "\tSO:")) != 0) {
                *p = '\n'; // change back
                if (strncmp(q + 4, so, p - q - 4) != 0) {
                    beg = q;
                    for (q += 4; *q != '\n' && *q != '\t'; ++q);
                    end = q;
                } else return 0; // no need to change
            } else beg = end = p, *p = '\n';
        }
    }
    if (beg == NULL) { // no @HD
        h->l_text += strlen(so) + 15;
        newtext = (char*)malloc(h->l_text + 1);
        sprintf(newtext, "@HD\tVN:1.3\tSO:%s\n", so);
        strcat(newtext, h->text);
    } else { // has @HD but different or no SO
        h->l_text = (beg - h->text) + (4 + strlen(so)) + (h->text + h->l_text - end);
        newtext = (char*)malloc(h->l_text + 1);
        strncpy(newtext, h->text, beg - h->text);
        sprintf(newtext + (beg - h->text), "\tSO:%s", so);
        strcat(newtext, end);
    }
    free(h->text);
    h->text = newtext;
    return 0;
}

int bam_Radix_sort(const char *fn,const char *fnout, const char *modeout,
                      size_t _max_mem,const htsFormat *in_fmt, const htsFormat *out_fmt)
{
    int ret = -1, i,j,k=0, n=0;

    samFile *fp = sam_open_format(fn, "r", in_fmt);
    if (fp == NULL) {
        print_error_errno("sort", "can't open \"%s\"", fn);
        return -2;
    }
    bam_hdr_t *header = sam_hdr_read(fp);
    if (header == NULL) {
        print_error("sort", "failed to read header from \"%s\"", fn);
        goto err;
    }
    change_SO(header, "coordinate");

    mapInfo_t * keyList = (mapInfo_t *)calloc(SORT_LIST_SIZE , sizeof(mapInfo_t));
    uint32_t keyListSize = SORT_LIST_SIZE;
    bam1_t *b=bam_init1() ;
    for (;;) {
        if ((ret = sam_read1(fp, header, b)) < 0) break;    
        keyList[n].index=n;
        keyList[n].chr = b->core.tid;
        keyList[n].pos = (uint32_t)((int32_t)b->core.pos+1)<<1 | bam_is_rev(b);    
        if(++n > keyListSize) {   //增加内存
            keyListSize += SORT_LIST_SIZE;
            keyList = (mapInfo_t *)realloc(keyList, keyListSize * sizeof(mapInfo_t));
            if(keyList == NULL) {
                perror("realloce keyList ");
                exit(1);
            }
        }
    }
    if (ret != -1) {
        print_error("sort", "truncated file. Aborting");
        ret = -1;
        goto err;
    }
    sam_close(fp);
    radixSort(n, keyList);

    samFile *fo = sam_open_format(fnout, modeout, out_fmt);
    if (fo == NULL) return -1;
    if (sam_hdr_write(fo, header) != 0) goto err;
    bam_hdr_destroy(header);
    header=NULL;
    uint8_t* flags = (uint8_t *)calloc(n , sizeof(uint8_t));
    uint32_t bufferReads = _max_mem,outReads=0;
    bam1_t **outPtrs = (bam1_t **)calloc(n , sizeof(bam1_t *)),
            **buf = (bam1_t **)calloc(bufferReads , sizeof(bam1_t *));
    for (i = 0; i < bufferReads; ++i) buf[i] = bam_init1();
    for (i = 0; i < (n / bufferReads) + 1; i++) {
        memset(flags,0,n * sizeof(uint8_t));
        for (j = 0; j < bufferReads && j + outReads < n; ++j) {
            flags[keyList[j + outReads].index]++;
        }
        fp = sam_open_format(fn, "r", in_fmt);
        if (fp == NULL) {
            print_error_errno("sort", "can't open \"%s\"", fn);
            return -2;
        }
        header = sam_hdr_read(fp);
        if (header == NULL) {
            print_error("sort", "failed to read header from \"%s\"", fn);
            goto err;
        }
        switch (fp->format.format) {
            case bam: {
                for (j=0,k=0;j<n;++j) {
                    if(flags[j]) {
                        ret = bam_read1(fp->fp.bgzf, buf[k]);
                        if (ret >= 0) {
                            if (buf[k]->core.tid  >= header->n_targets || buf[k]->core.tid  < -1 || 
                                buf[k]->core.mtid >= header->n_targets || buf[k]->core.mtid < -1){
                                ret = -3;
                                break;
                            }
                        }else{
                            break;
                        }
                        outPtrs[j] = buf[k++];
                    }else{
                        ret = bam_read1(fp->fp.bgzf, b);
                        if (ret >= 0) {
                            if (b->core.tid  >= header->n_targets || b->core.tid  < -1 || 
                                b->core.mtid >= header->n_targets || b->core.mtid < -1){
                                ret = -3;
                                break;
                            }
                        }else{
                            break;
                        }
                    }
                }
                break;
            }
            case cram: {
                for (j=0,k=0;j<n;++j) { 
                    if(flags[j]) {
                        ret = cram_get_bam_seq(fp->fp.cram, &buf[k]);
                        if (ret<0) {
                            ret = cram_eof(fp->fp.cram) ? -1 : -2;
                            break;
                        }
                        outPtrs[j] = buf[k++];
                    }else{
                        ret = cram_get_bam_seq(fp->fp.cram, &b);
                        if (ret<0) {
                            ret = cram_eof(fp->fp.cram) ? -1 : -2;
                            break;
                        }
                    }
                }
                break;
            }
            case sam: {
                for (j=0,k=0;j<n;++j) {
                    int ret;
err_recover:
                    if (fp->line.l == 0) {
                        ret = hts_getline(fp, KS_SEP_LINE, &fp->line);
                        if (ret < 0) break;
                    }
                    if(flags[j]) {
                        ret = sam_parse1(&fp->line, header, buf[k]);
                        fp->line.l = 0;
                        if (ret < 0) {
                            if (hts_verbose >= 1)
                                fprintf(stderr, "[W::%s] parse error at line %lld\n", __func__, (long long)fp->lineno);
                            if (header->ignore_sam_err) goto err_recover;
                            break;
                        }
                        outPtrs[j] = buf[k++];
                    }else{
                        ret = sam_parse1(&fp->line, header, b);
                        fp->line.l = 0;
                        if (ret < 0) {
                            if (hts_verbose >= 1)
                                fprintf(stderr, "[W::%s] parse error at line %lld\n", __func__, (long long)fp->lineno);
                            if (header->ignore_sam_err) goto err_recover;
                            break;
                        }
                    }
                }
            }
            default:
                abort();
        }
        for (j = 0; j < bufferReads && j + outReads < n; j++) {
            if (ret = sam_write1(fo, header, outPtrs[keyList[j + outReads].index] ) < 0) goto err;
        }
        outReads += j;
        bam_hdr_destroy(header);
        header=NULL;
        sam_close(fp);
    }
    free(flags);
    for (i = 0; i < bufferReads; ++i) bam_destroy1(buf[i]);
    free(buf);
    free(outPtrs);
    ret = 0;
 err:
    // free
    bam_destroy1(b);
    if (sam_close(fo) < 0) return -1;
    return ret;
}

static void sort_usage(FILE *fp)
{
    fprintf(fp,
"Usage: samtools sort [options...] [in.bam]\n"
"Options:\n"
"  -l INT     Set compression level, from 0 (uncompressed) to 9 (best)\n"
"  -m INT     Set maximum memory per thread; suffix K/M/G recognized [768M]\n"
"  -n         Sort by read name\n"
"  -o FILE    Write final output to FILE rather than standard output\n"
"  -T PREFIX  Write temporary files to PREFIX.nnnn.bam\n");
    sam_global_opt_help(fp, "-.O..@");
}

int main(int argc, char *argv[])
{
    size_t max_mem = SORT_DEFAULT_MEGS_PER_THREAD << 20;
    int c, nargs, ret, level = -1;
    char *fnout = "-", modeout[12];
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;

    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 'O', 0, 0, '@'),
        { "threads", required_argument, NULL, '@' },
        { NULL, 0, NULL, 0 }
    };

    while ((c = getopt_long(argc, argv, "l:m:o:O:T:@:", lopts, NULL)) >= 0) {
        switch (c) {
        case 'o': fnout = optarg; break;
        case 'm': {
                char *q;
                max_mem = strtol(optarg, &q, 0);
                if (*q == 'k' || *q == 'K') max_mem <<= 10;
                else if (*q == 'm' || *q == 'M') max_mem <<= 20;
                else if (*q == 'g' || *q == 'G') max_mem <<= 30;
                break;
            }
        case 'l': level = atoi(optarg); break;

        default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
                  /* else fall-through */
        case '?': sort_usage(stderr); ret = EXIT_FAILURE; goto sort_end;
        }
    }

    nargs = argc - optind;
    if (nargs == 0 && isatty(STDIN_FILENO)) {
        sort_usage(stdout);
        ret = EXIT_SUCCESS;
        goto sort_end;
    }
    else if (nargs >= 2) {
        if (nargs == 2) fprintf(stderr, "[bam_sort] Use -T PREFIX / -o FILE to specify temporary and final output files\n");
        sort_usage(stderr);
        ret = EXIT_FAILURE;
        goto sort_end;
    }
    strcpy(modeout, "wb");
    sam_open_mode(modeout+1, fnout, NULL);
    if (level >= 0) sprintf(strchr(modeout, '\0'), "%d", level < 9? level : 9);

    ret = bam_Radix_sort((nargs > 0)? argv[optind] : "-",fnout, modeout, max_mem, &ga.in, &ga.out);
    if (ret >= 0) ret = EXIT_SUCCESS;

sort_end:
    sam_global_args_free(&ga);

    return ret;
}
