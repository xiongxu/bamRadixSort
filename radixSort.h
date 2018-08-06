//
//  radixSort.h
//  infiniAnalyzer
//
//  Created by xubl on 01/05/2017.
//  Modified by xux on 05/06/2017.
//  Copyright © 2017 infiniGenomics. All rights reserved.
//

#ifndef radixSort_h
#define radixSort_h

#define _RADIX_ 1024
#define _BITS_N_ 10    // 2^10 = 1024
#define _NORMAL_MASK_ 0x3FF // 0x3FF = 1024-1 perl -e 'printf("%b\n",0x3FF);' 1111111111
#define _LAST_MASK_ 0xFF // 3 bits 111

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include "string.h"

typedef struct _mapInfo_t_ {
    uint32_t chr;
    uint32_t index;
    uint32_t pos;
}mapInfo_t;

int radixSort(int n, mapInfo_t * sourceList);

#endif /* radixSort_h */
