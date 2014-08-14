/*                          R E C T . C
 * BRL-CAD
 *
 * Copyright (c) 1998-2014 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file libdm/rect.c
 *
 * Rubber band rectangle.
 *
 */

#include "common.h"

#include <math.h>
#include <stdio.h>

#include "bio.h"
#include "bu.h"
#include "vmath.h"
#include "dm.h"
#include "dm_private.h"

void
dm_draw_rect(dm *dmp, struct ged_rect_state *grsp)
{
    if (ZERO(grsp->grs_width) &&
	ZERO(grsp->grs_height))
	return;

    /* draw rectangle */
    dm_set_fg(dmp,
		   (unsigned char)grsp->grs_color[0],
		   (unsigned char)grsp->grs_color[1],
		   (unsigned char)grsp->grs_color[2], 1, 1.0);
    dm_set_line_attr(dmp, grsp->grs_line_width, grsp->grs_line_style);

    dm_draw_line_2d(dmp,
		    grsp->grs_x,
		    grsp->grs_y * dmp->dm_aspect,
		    grsp->grs_x,
		    (grsp->grs_y + grsp->grs_height) * dmp->dm_aspect);
    dm_draw_line_2d(dmp,
		    grsp->grs_x,
		    (grsp->grs_y + grsp->grs_height) * dmp->dm_aspect,
		    grsp->grs_x + grsp->grs_width,
		    (grsp->grs_y + grsp->grs_height) * dmp->dm_aspect);
    dm_draw_line_2d(dmp,
		    grsp->grs_x + grsp->grs_width,
		    (grsp->grs_y + grsp->grs_height) * dmp->dm_aspect,
		    grsp->grs_x + grsp->grs_width,
		    grsp->grs_y * dmp->dm_aspect);
    dm_draw_line_2d(dmp,
		    grsp->grs_x + grsp->grs_width,
		    grsp->grs_y * dmp->dm_aspect,
		    grsp->grs_x,
		    grsp->grs_y * dmp->dm_aspect);
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
