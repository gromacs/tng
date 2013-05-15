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


/* Here presence of 64 bit integers are checked for. If on windows
   define USE_WINDOWS if it does not get defined automatically. 
   If on UNIX, define the sizes of SIZEOF_INT, SIZEOF_LONG, SIZEOF_LONG_LONG
   If none of these symbols are defined, stdint.h is included (if on C99
   or USE_STDINT is defined) and 64 bit integers are assumed to be
   present. If none of this is fulfilled, 64 bit integers are not
   available. */

#ifndef MY64BIT_H
#define MY64BIT_H

#if HAVE_CONFIG_H
#include <config.h>
#endif

#if HAVE_STDINT_H
#include <stdint.h>
#endif

/* The USE_WINDOWS symbol should be automatically defined in tng_compress.h */
#include "tng_compress.h"

#ifdef USE_WINDOWS
typedef __int64 my_int64_t;
typedef unsigned __int64 my_uint64_t;
#define HAVE64BIT
#else  /* USE_WINDOWS */
/* UNIX. Use config.h */
#if SIZEOF_INT >= 8 
typedef int my_int64_t;
typedef unsigned int my_uint64_t;
#define HAVE64BIT
#else /* SIZEOF_INT */
#if SIZEOF_LONG >= 8 
typedef long my_int64_t;
typedef unsigned long my_uint64_t;
#define HAVE64BIT
#else /* SIZEOF_LONG */
#if SIZEOF_LONG_LONG >= 8 
typedef long long my_int64_t;
typedef unsigned long long my_uint64_t;
#define HAVE64BIT
#else /* SIZEOF_LONG_LONG */
#if HAVE_STDINT_H
typedef int64_t my_int64_t;
typedef uint64_t my_uint64_t;
#define HAVE64BIT
#endif /* STDINT_H */
#endif /* SIZEOF_LONG_LONG */
#endif /* SIZEOF_LONG */
#endif /* SIZEOF_INT */
#endif  /* USE_WINDOWS */

#endif
