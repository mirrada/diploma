/*
 * list.h: this file is part of the FM project.
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
#ifndef FM_LIST_H
# define FM_LIST_H

# include <assert.h>
# include <common.h>

BEGIN_C_DECLS

struct s_fm_list
{
  void*                        data;
  struct s_fm_list*        next;
};

typedef struct s_fm_list s_fm_list_t;
typedef void (*free_fun_t) (void*);
typedef int (*cmp_fun_t) (void*, void*);

extern s_fm_list_t*        fm_list_new(void *data);
extern void                fm_list_dummy_free (void*);
extern void                fm_list_free (s_fm_list_t* l, free_fun_t f);
extern s_fm_list_t*        fm_list_add_head (s_fm_list_t** head, void* data);
extern s_fm_list_t*        fm_list_add_head_unique (s_fm_list_t** head,
                                                 void* data,
                                                 cmp_fun_t f);
extern s_fm_list_t*        fm_list_cons(s_fm_list_t *head, s_fm_list_t *tail);
extern s_fm_list_t*        s_fm_list_tail(s_fm_list_t *head);
extern size_t                fm_list_length(s_fm_list_t *head);

extern int                fm_list_remove (s_fm_list_t** head, void* data);

END_C_DECLS


#endif // FM_LIST_H
