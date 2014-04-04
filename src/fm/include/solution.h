/*
 * solution.h: this file is part of the FM project.
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
#ifndef FM_SOLUTION_H
# define FM_SOLUTION_H

# include <list.h>
# include <vector.h>
# include <system.h>
# include <options.h>

# define FM_SOLUTION_PRINT_POINT 1

BEGIN_C_DECLS

struct s_fm_ball
{
  s_fm_list_t* positive;
  s_fm_list_t* negative;
};

typedef struct s_fm_ball s_fm_ball_t;


struct s_fm_solution
{
  s_fm_ball_t* solution;
  unsigned size;
};

typedef struct s_fm_solution s_fm_solution_t;

typedef void (*point_fun_t)(s_fm_solution_t*, s_fm_vector_t*, int, void*);

extern
s_fm_solution_t*
fm_solution_alloc (size_t size);

extern
void
fm_solution_free (s_fm_solution_t* s);

extern
s_fm_solution_t*
fm_solution_dup (s_fm_solution_t* s);

extern
void
fm_solution_print (FILE* stream, s_fm_solution_t* s);


extern
s_fm_system_t*
fm_solution_to_system (s_fm_solution_t* s);


extern
s_fm_system_t*
fm_solution_to_system_at (s_fm_solution_t* s, int idx);


extern
int
fm_solution_add_line_at (s_fm_solution_t* s, s_fm_vector_t* v, unsigned idx);

extern
int
fm_solution_add_unique_line_at (s_fm_solution_t* s,
                                s_fm_vector_t* v,
                                unsigned idx);

extern
s_fm_solution_t*
fm_system_to_solution (s_fm_system_t* s);


extern
s_fm_system_t*
fm_system_reduce (s_fm_system_t* in, s_fm_solution_t* redeq);


extern
void
fm_system_subst_in_vector (s_fm_vector_t* v,
                           s_fm_vector_t* v1,
                           s_fm_vector_t* pattern,
                           s_fm_solution_t* redeq,
                           int iv1);

extern
int
fm_solution_equalities_find (s_fm_solution_t* s);

extern
void
fm_solution_cut (s_fm_solution_t* s, int dim);

extern
int
fm_solution_point_included (s_fm_solution_t* s, s_fm_vector_t* v);

extern
unsigned long long int
fm_solution_count (s_fm_solution_t* sol, int bound_limit, int mask);

extern
void
fm_solution_traverse (s_fm_solution_t* sol, int bound_limit, int mask,
                      point_fun_t f, void* data);

/* Must have piplib on the system. */
# ifdef HAVE_LIBPIPLIB

extern
s_fm_solution_t*
fm_solution_simplify (s_fm_solution_t* s,
                      int simplify_mode);

# endif

END_C_DECLS


#endif // FM_SOLUTION_H
