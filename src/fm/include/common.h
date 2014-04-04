/*
 * common.h: this file is part of the FM project.
 * 
 * FM, a fast and optimized C implementation of Fourier-Motzkin 
 * projection algorithm.
 * 
 * Copyright (C) 2006-2008 Louis-Noel Pouchet
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 3
 * of the License, or (at your option) any later version.
 * 
 * The complete GNU Lesser General Public Licence Notice can be found
 *  as the `COPYING.LESSER' file in the root directory.
 * 
 * Author:
 * Louis-Noel Pouchet <Louis-Noel.Pouchet@inria.fr>
 * 
 */
#ifndef FM_COMMON_H
# define FM_COMMON_H

# if HAVE_CONFIG_H
#  include "config.h"
# endif

#include <ctype.h>

# include <stdio.h>
# include <sys/types.h>

# if STDC_HEADERS
#  include <stdlib.h>
#  include <string.h>
# elif HAVE_STRINGS_H
#  include <strings.h>
# endif /*STDC_HEADERS*/

# if HAVE_UNISTD_H
#  include <unistd.h>
# endif


# ifdef __cplusplus
#  define BEGIN_C_DECLS         extern "C" {
#  define END_C_DECLS           }
# else
#  define BEGIN_C_DECLS
#  define END_C_DECLS
# endif


# ifdef __GNUC__
#  ifndef const
#    define const       __const
#  endif
#  ifndef signed
#    define signed      __signed
#  endif
#  ifndef volatile
#    define volatile    __volatile
#  endif
# else
#  ifdef __STDC__
#    undef  signed
#    define signed
#    undef  volatile
#    define volatile
#  endif
# endif

# ifdef __STDC__
#  define STR(x)          #x
#  define CONC(x, y)      x##y
# else
#  define STR(x)          "x"
#  define CONC(x, y)      x/**/y
# endif


# ifndef EXIT_SUCCESS
#  define EXIT_SUCCESS  0
#  define EXIT_FAILURE  1
# endif


# define XCALLOC(type, num)                                  \
        ((type *) xcalloc ((num), sizeof(type)))
# define XMALLOC(type, num)                                  \
        ((type *) xmalloc ((num) * sizeof(type)))
# define XREALLOC(type, p, num)                              \
        ((type *) xrealloc ((p), (num) * sizeof(type)))
# define XFREE(stale)                            do {        \
        if (stale) { free (stale);  stale = 0; }            \
                                                } while (0)

# define FM_STATUS_ERROR 1
# define FM_STATUS_OK    0

#define FM_SIZE_BUF 5000


BEGIN_C_DECLS

extern void *xcalloc    (size_t num, size_t size);
extern void *xmalloc    (size_t num);
extern void *xrealloc   (void *p, size_t num);
extern char *xstrdup    (const char *string);
extern char *xstrerror  (int errnum);

typedef long int z_type_t;
# define Z_STRING_MODIFIER "%ld"


END_C_DECLS

#endif // FM_COMMON_H
