/*
 * system.h: this file is part of the FM project.
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
#ifndef FM_SYSTEM_H
# define FM_SYSTEM_H

# include <stdio.h>
# include <vector.h>

// #define FM_SOLVER_REDREC_DESCENDANT        32
// #define FM_SOLVER_REDREC_IRIGOIN        64
#define FM_SYSTEM_REDREC_DESCENDANT        32 // Beware to match the values in
#define FM_SYSTEM_REDREC_IRIGOIN        64 // solver.h (convenience).

#define FM_SYSTEM_ALLOC_BUFFER_SIZE        64

BEGIN_C_DECLS

struct s_fm_system
{
  s_fm_vector_t** lines;
  unsigned nb_lines;
  unsigned nb_cols;
  // Internal service field.
  unsigned allocated;
};

typedef struct s_fm_system s_fm_system_t;


extern
s_fm_system_t*
fm_system_alloc (size_t nb_lines, size_t nb_cols);

extern
s_fm_system_t*
fm_system_dup (s_fm_system_t* s);

extern
void
fm_system_free (s_fm_system_t* s);

extern
void
fm_system_normalize_ineq (s_fm_system_t* s);

extern
s_fm_system_t*
fm_system_read (FILE* stream);

extern
int
fm_system_print (FILE* stream, s_fm_system_t* s);

extern
int
fm_system_sort (s_fm_system_t* s, unsigned idx, unsigned* n1, unsigned* n2);

extern
int
fm_system_remove_line (s_fm_system_t* s, unsigned idx);

extern
int
fm_system_add_line_at (s_fm_system_t* s, s_fm_vector_t* v, unsigned idx);

extern
int
fm_system_add_line (s_fm_system_t* s, s_fm_vector_t* v);

extern
int
fm_system_add_column (s_fm_system_t* s, int pos);

extern
s_fm_system_t*
fm_system_to_z (s_fm_system_t* s);

extern
s_fm_system_t*
fm_system_swap_column (s_fm_system_t* s, unsigned col_src, unsigned col_dst);

extern
s_fm_system_t*
fm_system_split (s_fm_system_t* s, unsigned col);

extern
int
fm_system_remove_duplicate (s_fm_system_t* s);

extern
void
fm_system_equalities_sort (s_fm_system_t* s);

extern
int
fm_system_equalities_find (s_fm_system_t* s);

extern
int
fm_system_point_included (s_fm_system_t* s, s_fm_vector_t* v);

extern
int**
fm_system_getconnected (s_fm_system_t* system);


/* Must have piplib on the system. */
# ifdef HAVE_LIBPIPLIB

// Forward declaration.
extern
int
fm_piptools_check_int (s_fm_system_t* sys);


extern
s_fm_system_t*
fm_system_simplify (s_fm_system_t* system,
                    int simplify_mode);
# endif


END_C_DECLS


#endif // FM_SYSTEM_H
