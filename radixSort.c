//  radixSort.c
//  infiniAnalyzer
//
//  Created by xubl on 01/05/2017.
//  Modified by xux on 05/06/2017.
//  Copyright © 2017 infiniGenomics. All rights reserved.

#include "radixSort.h"
int radixSort(int n, mapInfo_t * sourceList) {
    mapInfo_t * dstList = (mapInfo_t *)calloc(n, sizeof(mapInfo_t));
    if(dstList == NULL) perror("calloc on radixSort function");
    uint64_t index[_RADIX_] = {0}, i ,sum ,counter[_RADIX_] = {0};
#define COUNT(__mask, __bits_n, __src, __dst, __attr)                                                       \
    for (i = 0; i < n; i++) counter[((__src + i)->__attr >> (__bits_n)) & __mask]++;                        \
    sum = 0;                                                                                                \
    for (i = 0; i < _RADIX_; i++) {                                                                         \
        index[i] = sum;                                                                                     \
        sum += counter[i];                                                                                  \
    }                                                                                                       \
    for (i = 0; i < n; i++) {                                                                               \
        *(__dst + index[((__src + i)->__attr >> (__bits_n)) & __mask]++) = *(__src + i);                    \
    }                                                                                                       \
    for(i = 0; i < _RADIX_; i++) {                                                                          \
        counter[i] = 0;                                                                                     \
        index[i] = 0;                                                                                       \
    }
    COUNT(_NORMAL_MASK_, 0, sourceList, dstList, pos)               /* count 0-9 bits */
    COUNT(_NORMAL_MASK_, _BITS_N_, dstList, sourceList, pos)        /* count 10-19 bits */
    COUNT(_NORMAL_MASK_, _BITS_N_ * 2, sourceList, dstList, pos)    /* count 20-29 bits */
    COUNT(_LAST_MASK_, 0, dstList, sourceList, chr)                 /* count chromosome */
    free(dstList);
    return 0;
}
