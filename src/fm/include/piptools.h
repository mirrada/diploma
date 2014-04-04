/*
 * piptools.h: this file is part of the FM project.
 *
 * FM, a fast and optimized C implementation of Fourier-Motzkin
 * projection algorithm.
 *
 * Copyright (C) 2007-2008 Louis-Noel Pouchet
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
#ifndef FM_PIPTOOLS_H
# define FM_PIPTOOLS_H

# include <assert.h>
# include <common.h>

/* Must have piplib on the system. */

// # ifdef HAVE_LIBPIPLIB

#  include <piplib/piplib64.h>

#  include <system.h>
#  include <solution.h>


#  define FM_PIPTOOLS_RAT 0
#  define FM_PIPTOOLS_INT 1


BEGIN_C_DECLS

extern
int
fm_piptools_check_rat (s_fm_system_t* sys);

extern
int
fm_piptools_check_int (s_fm_system_t* sys);

extern
int
fm_piptools_check (s_fm_system_t* system, int mode);

extern
PipQuast*
fm_piptools_pip (s_fm_system_t* sys, s_fm_system_t* context, int mode);

extern
int
fm_piptools_check_sol (s_fm_solution_t* sol, int mode);

extern
int
fm_piptools_check_sol_msg (char* msg,
                           FILE* stream,
                           s_fm_solution_t* sol,
                           int mode);

extern
int
fm_piptools_pipmatrix_equal (PipMatrix* a,
                             PipMatrix* b);

extern
PipMatrix*
fm_piptools_st_to_pipmatrix (s_fm_system_t* sys);

extern
s_fm_system_t*
fm_piptools_pm_to_system (PipMatrix* m);


END_C_DECLS


// # endif // HAVE_LIBPIPLIB

#endif // FM_PIPTOOLS_H
