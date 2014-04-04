/*
 * solver.h: this file is part of the FM project.
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
#ifndef FM_SOLVER_H
# define FM_SOLVER_H

# include <system.h>
# include <solution.h>
# include <options.h>


BEGIN_C_DECLS

extern
s_fm_solution_t*
fm_solver (s_fm_system_t* system, int solver_type);

extern
s_fm_solution_t*
fm_solver_solution_at (s_fm_system_t* system, int solver_type, unsigned last);

extern
s_fm_solution_t*
fm_solver_solution_to (s_fm_system_t* system, int solver_type, unsigned to);

extern
s_fm_rational_t**
fm_solver_minlexico(s_fm_solution_t* sol, z_type_t min, int is_integral);

extern
s_fm_rational_t**
fm_solver_maxlexico(s_fm_solution_t* sol, z_type_t max, int is_integral);

extern
void
fm_solver_compute_min (s_fm_rational_t** lb,
                       s_fm_list_t* l,
                       s_fm_vector_t* vect,
                       unsigned idx,
                       int is_int);

extern
void
fm_solver_compute_max (s_fm_rational_t** Ub,
                       s_fm_list_t* l,
                       s_fm_vector_t* vect,
                       unsigned idx,
                       int is_int);


s_fm_solution_t*
fm_solver_linind (s_fm_system_t* A);

int
fm_solver_gauss (s_fm_system_t* A);


END_C_DECLS


#endif // FM_SOLVER_H
