/*
 * compsol.h: this file is part of the FM project.
 *
 * FM, a fast and optimized C implementation of Fourier-Motzkin
 * projection algorithm.
 *
 * Copyright (C) 2007,2008 Louis-Noel Pouchet
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
#ifndef FM_COMPSOL_H
# define FM_COMPSOL_H

# include <stdio.h>
# include <vector.h>
# include <system.h>
# include <solution.h>
# include <solver.h>


BEGIN_C_DECLS

/**
 * A compacted solution, composed of a convex polyhedron and a set of
 * equalities for linearly dependent variables of the system.
 *
 *
 */
struct s_fm_compsol
{
  s_fm_solution_t* redeq;        /// List of equations to build dependent
                                /// variables
  s_fm_solution_t* redfree;        /// List of constraints for free variables
  s_fm_solution_t* poly;        /// Compacted (convex) polyhedron
  unsigned size;                /// Size of the expanded polyhedron
  unsigned nb_reduc;                /// Number of equality reduced variables
  unsigned nb_free;                /// Number of free variables
  unsigned empty;                /// Boolean for emptiness of the polyhedron
};

typedef struct s_fm_compsol s_fm_compsol_t;

extern
s_fm_compsol_t*
fm_compsol_alloc (unsigned size);

extern
s_fm_compsol_t*
fm_compsol_dup (s_fm_compsol_t* s);

extern
void
fm_compsol_free (s_fm_compsol_t* s);

extern
s_fm_compsol_t*
fm_compsol_init_sol (s_fm_solution_t* in);

extern
s_fm_compsol_t*
fm_compsol_init_sol_free (s_fm_solution_t* in, s_fm_solution_t* sfree);

extern
s_fm_compsol_t*
fm_compsol_init_sys (s_fm_system_t* in);

extern
void
fm_compsol_add_unique (s_fm_compsol_t* s, s_fm_vector_t* v);

extern
s_fm_solution_t*
fm_compsol_expand (s_fm_compsol_t* s);

extern
s_fm_solution_t*
fm_compsol_redundancy_elimination (s_fm_solution_t* s);

extern
s_fm_system_t*
fm_compsol_redundancy_elimination_sys (s_fm_system_t* s);

END_C_DECLS


#endif // FM_COMPSOL_H
