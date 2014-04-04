/*
 * macros.h: this file is part of the FM project.
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


#ifndef FM_MACROS_H
# define FM_MACROS_H

# include <stdlib.h>
# include <common.h>




# define Z_INIT(v1) ( v1 = 0 )
# define Z_CLEAR(v1)
# define Z_ASSIGN(v, v1) ( v = v1 )
# define Z_ASSIGN_SI(v, val) ( v = val )


# define Z_CMP_SI(v1, op, v2) (v1 op v2)
# define Z_CMP(v1, op, v2) (v1 op v2)


# define Z_PRINT(out, v) fprintf(out, Z_STRING_MODIFIER, v)


# define Z_INC(v, v1) ((v = v1 + 1))
# define Z_INC_(v1) ((v1++))

# define Z_DEC(v, v1) ((v = v1 - 1))
# define Z_DEC_(v1) ((v1--))

# define Z_ABS(v, v1) ( v = v1 > 0 ? v1 : -v1 )
# define Z_ABS_(v1) (v1 > 0 ? v1 : -v1)

# define Z_OPP(v, v1) ((v = -v1))
# define Z_OPP_(v1) (-v1)


# define Z_GCD(v, v1, v2) ( v = fm_z_gcd(v1, v2) )
# define Z_GCD_(v1, v2) (fm_z_gcd(v1, v2))

# define Z_LCM(v, v1, v2) ( v = fm_z_lcm(v1, v2) )
# define Z_LCM_(v1, v2) (fm_z_lcm(v1, v2))


# define Z_ADD(v, v1, v2) ( v = v1 + v2 )
# define Z_ADD_(v1, v2) (v1 + v2)

# define Z_SUB(v, v1, v2) ( v = v1 - v2 )
# define Z_SUB_(v1, v2) (v1 - v2)

# define Z_MUL(v, v1, v2) ( v = v1 * v2 )
# define Z_MUL_(v1, v2) (v1 * v2)

# define Z_DIV(v, v1, v2) ( v = v1 / v2 )
# define Z_DIV_(v1, v2) (v1 / v2)

# define Z_MOD(v, v1, v2) ( v = v1 % v2 )
# define Z_MOD_(v1, v2) (v1 % v2)

#endif // FM_MACROS_H
