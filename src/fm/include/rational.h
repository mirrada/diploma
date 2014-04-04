/*
 * rational.h: this file is part of the FM project.
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
#ifndef FM_RATIONAL_H
# define FM_RATIONAL_H

# include <assert.h>
# include <common.h>
# include <macros.h>

BEGIN_C_DECLS

struct s_fm_rational
{
  z_type_t        num;
  z_type_t        denum;
};

typedef struct s_fm_rational s_fm_rational_t;

extern
inline
z_type_t
fm_z_gcd(z_type_t a, z_type_t b);

extern
inline
z_type_t
fm_z_lcm(z_type_t a, z_type_t b);


extern
inline
s_fm_rational_t*
fm_rational_alloc ();

extern
inline
void
fm_rational_init (s_fm_rational_t* r);

extern
inline
void
fm_rational_print (FILE* stream, s_fm_rational_t* r);

extern
inline
void
fm_rational_free (s_fm_rational_t* r);

extern
inline
int
fm_rational_assign (s_fm_rational_t* r, z_type_t num, z_type_t denum);

extern
inline
int
fm_rational_copy (s_fm_rational_t* r, s_fm_rational_t* s);

extern
inline
int
fm_rational_assign_int (s_fm_rational_t* r, z_type_t num);

extern
inline
int
fm_rational_cmp (s_fm_rational_t* r1, s_fm_rational_t* r2);

extern
inline
void
fm_rational_add (s_fm_rational_t*r,  s_fm_rational_t* r1, s_fm_rational_t* r2);

extern
inline
void
fm_rational_sub (s_fm_rational_t*r,  s_fm_rational_t* r1, s_fm_rational_t* r2);

extern
inline
void
fm_rational_mul (s_fm_rational_t*r,  s_fm_rational_t* r1, s_fm_rational_t* r2);

extern
inline
void
fm_rational_div (s_fm_rational_t*r,  s_fm_rational_t* r1, s_fm_rational_t* r2);

extern
inline
void
fm_rational_opp (s_fm_rational_t*r,  s_fm_rational_t* r1);

extern
inline
int
fm_rational_equal (s_fm_rational_t* r1, s_fm_rational_t* r2);


END_C_DECLS


#endif // FM_RATIONAL_H
