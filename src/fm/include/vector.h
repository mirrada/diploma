/*
 * vector.h: this file is part of the FM project.
 *
 * FM, a fast and optimized C implementation of Fourier-Motzkin
 * projection algorithm.
 *
 * Copyright (C) 2006,2007,2008 Louis-Noel Pouchet
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
#ifndef FM_VECTOR_H
# define FM_VECTOR_H

# include <common.h>
# include <rational.h>

BEGIN_C_DECLS

struct s_fm_vector
{
  size_t                size;
  s_fm_rational_t*        vector;
  z_type_t                key;
};

typedef struct s_fm_vector s_fm_vector_t;


extern
void
fm_vector_compute_key(z_type_t* key, s_fm_vector_t* v);

extern
inline
s_fm_vector_t*
fm_vector_alloc (size_t size);

extern
inline
s_fm_vector_t*
fm_vector_dup (s_fm_vector_t* v1);

extern
inline
int
fm_vector_init (s_fm_vector_t* v, size_t size);

extern
inline
void
fm_vector_free (s_fm_vector_t* v);

extern
inline
void
fm_vector_read (FILE* stream, s_fm_vector_t* v, unsigned size);

extern
inline
void
fm_vector_print (FILE* stream, s_fm_vector_t* v);

extern
inline
void
fm_vector_set_ineq (s_fm_vector_t* v);

extern
inline
void
fm_vector_set_eq (s_fm_vector_t* v);

extern
inline
int
fm_vector_assign (s_fm_vector_t* v, s_fm_vector_t* v1);

extern
inline
int
fm_vector_assign_at (s_fm_vector_t* v, s_fm_vector_t* v1, unsigned idx);

extern
inline
int
fm_vector_assign_idx (s_fm_vector_t* v, s_fm_rational_t* r, unsigned idx);

extern
inline
int
fm_vector_expand (s_fm_vector_t* v, s_fm_vector_t* v1);

extern
inline
int
fm_vector_expand_at (s_fm_vector_t* v, s_fm_vector_t* v1, unsigned idx);

extern
inline
int
fm_vector_shrink (s_fm_vector_t* v, s_fm_vector_t* v1, unsigned idx);

extern
inline
int
fm_vector_assign_int_idx (s_fm_vector_t* v, z_type_t i, unsigned idx);

extern
inline
int
fm_vector_is_null (s_fm_vector_t* v);

extern
inline
int
fm_vector_is_empty (s_fm_vector_t* v);

extern
inline
int
fm_vector_is_valid (s_fm_vector_t* v);

extern
inline
int
fm_vector_is_scalar_cst (s_fm_vector_t* v);

extern
inline
int
fm_vector_opp (s_fm_vector_t* v, s_fm_vector_t* v1);

extern
inline
int
fm_vector_add (s_fm_vector_t* v, s_fm_vector_t* v1, s_fm_vector_t* v2);

extern
inline
int
fm_vector_sub (s_fm_vector_t* v, s_fm_vector_t* v1, s_fm_vector_t* v2);

extern
inline
int
fm_vector_normalize_idx (s_fm_vector_t* v, s_fm_vector_t* v1, unsigned idx);

extern
inline
int
fm_vector_to_z (s_fm_vector_t* v, s_fm_vector_t* v1);

extern
inline
int
fm_vector_resize (s_fm_vector_t* v, s_fm_vector_t* v1);

extern
inline
int
fm_vector_equal (s_fm_vector_t* v1, s_fm_vector_t* v2);

extern
inline
int
fm_vector_do_subsume (s_fm_vector_t* v1, s_fm_vector_t* v2);


extern
s_fm_vector_t*
fm_vector_read_str (char* stream);


END_C_DECLS


#endif // FM_VECTOR_H
