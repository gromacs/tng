/* This code is part of the tng compression routines.
 *
 * Written by Daniel Spangberg
 * Copyright (c) 2010, 2013, The GROMACS development team.
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 */


#ifndef WARNMALLOC_H
#define WARNMALLOC_H

#include "compression/tng_compress.h"

void DECLSPECDLLEXPORT *Ptngc_warnmalloc_x(size_t size, char *file, int line);

#define warnmalloc(size) Ptngc_warnmalloc_x(size,__FILE__,__LINE__)

void DECLSPECDLLEXPORT *Ptngc_warnrealloc_x(void *old, size_t size, char *file, int line);

#define warnrealloc(old,size) Ptngc_warnrealloc_x(old,size,__FILE__,__LINE__)


#endif
