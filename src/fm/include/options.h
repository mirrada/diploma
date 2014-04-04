/*
 * options.h: this file is part of the FM project.
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
#ifndef FM_DEFINES_H
# define FM_DEFINES_H


/*
 * Solver options.
 */
# define FM_SOLVER_DEFAULT                0
# define FM_SOLVER_FAST                        1
# define FM_SOLVER_MEMSIZE                2
# define FM_SOLVER_NORMALIZE_EQ                4
# define FM_SOLVER_DEBUG                       8
# define FM_SOLVER_SIMPLIFY                16
# define FM_SOLVER_REDREC_DESCENDANT        32
# define FM_SOLVER_REDREC_IRIGOIN        64
# define FM_SOLVER_AUTO_SIMPLIFY        128

# define FM_SOLVER_VERBOSE                256

# define FM_SOLVER_AUTOSIMP_THRESOLD        2000
# define FM_SOLVER_ERRPTR                0x1

/*
 * Lexico-minimum computation options.
 */
# define FM_MINLEXICO_RAT 0
# define FM_MINLEXICO_INT 1


#endif // FM_DEFINES_H
