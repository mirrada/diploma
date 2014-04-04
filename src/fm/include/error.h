/*
 * error.h: this file is part of the FM project.
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
#ifndef FM_ERROR_H
# define FM_ERROR_H 

#include <common.h>

BEGIN_C_DECLS

extern const char *program_name;
extern void        set_program_name (const char *argv0);

extern void        fm_warning (const char *message);
extern void        fm_error (const char *message);
extern void        fm_fatal (const char *message);

END_C_DECLS


#endif // FM_ERROR_H
