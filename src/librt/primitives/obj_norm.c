/*                    O B J _ N O R M . C
 * BRL-CAD
 *
 * Copyright (c) 2010-2012 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */

#include "common.h"

#include "raytrace.h"


int
rt_obj_norm(struct hit *hitp, struct soltab *stp, struct xray *rp)
{
    int id;
    const struct rt_functab *ft;

    if (!hitp || !stp)
	return -1;

    RT_CK_SOLTAB(stp);
    RT_CK_HIT(hitp);
    if (rp) RT_CK_RAY(rp);

    id = stp->st_id;
    if (id < 0)
	return -2;

    ft = &rt_functab[id];
    if (!ft)
	return -3;
    if (!ft->ft_norm)
	return -4;

    ft->ft_norm(hitp, stp, rp);
    return 0;
}


/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
