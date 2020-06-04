/*                         S C R E E N G R A B . C
 * BRL-CAD
 *
 * Copyright (c) 2008-2020 United States Government as represented by
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
/** @file libged/screengrab.c
 *
 * The screengrab command.
 *
 */

#include "common.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "icv.h"
#include "dm.h"

#include "./ged_private.h"

static int
image_mime(struct bu_vls *msg, size_t argc, const char **argv, void *set_mime)
{
    int type_int;
    bu_mime_image_t type = BU_MIME_IMAGE_UNKNOWN;
    bu_mime_image_t *set_type = (bu_mime_image_t *)set_mime;

    BU_OPT_CHECK_ARGV0(msg, argc, argv, "mime format");

    type_int = bu_file_mime(argv[0], BU_MIME_IMAGE);
    type = (type_int < 0) ? BU_MIME_IMAGE_UNKNOWN : (bu_mime_image_t)type_int;
    if (type == BU_MIME_IMAGE_UNKNOWN) {
        if (msg) bu_vls_sprintf(msg, "Error - unknown geometry file type: %s \n", argv[0]);
        return -1;
    }
    if (set_type) (*set_type) = type;
    return 1;
}

int
ged_screen_grab(struct ged *gedp, int argc, const char *argv[])
{

    int i;
    int print_help = 0;
    int width = 0;
    int height = 0;
    int scr_xoff = 0;
    int scr_yoff = 0;
    int bytes_per_pixel = 0;
    int bytes_per_line = 0;
    int grab_fb = 0;
    unsigned char **rows = NULL;
    unsigned char *idata = NULL;
    struct icv_image *bif = NULL;	/**< icv image container for saving images */
    struct dm *dmp = NULL;
    struct fb *fbp = NULL;
    bu_mime_image_t type = BU_MIME_IMAGE_AUTO;
    static char usage[] = "Usage: screengrab [-h] [-F] [-X scr_xoff] [-Y scr_yoff] [-w scr_width] [-n scr_height] [--format fmt] [file.img]\n";

    struct bu_opt_desc d[8];
    BU_OPT(d[0], "h", "help",           "",            NULL,      &print_help,       "Print help and exit");
    BU_OPT(d[1], "F", "fb",             "",     NULL,             &grab_fb,          "screengrab framebuffer instead of scene display");
    BU_OPT(d[2], "X", "scr_xoff",       "#",    &bu_opt_int,      &scr_xoff,         "X offset");
    BU_OPT(d[3], "Y", "scr_yoff",       "#",    &bu_opt_int,      &scr_yoff,         "Y offset");
    BU_OPT(d[4], "W", "scr_maxwidth",   "#",    &bu_opt_int,      &width,            "width of image to grab");
    BU_OPT(d[5], "N", "scr_maxheight",  "#",    &bu_opt_int,      &height,           "height of image to grab");
    BU_OPT(d[6], "",  "format",         "fmt",  &image_mime,      &type,             "output image file format");
    BU_OPT_NULL(d[7]);

    GED_CHECK_DATABASE_OPEN(gedp, GED_ERROR);
    GED_CHECK_VIEW(gedp, GED_ERROR);
    GED_CHECK_DRAWABLE(gedp, GED_ERROR);
    GED_CHECK_ARGC_GT_0(gedp, argc, GED_ERROR);

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (!gedp->ged_dmp) {
	bu_vls_printf(gedp->ged_result_str, ": no display manager currently active");
	return GED_ERROR;
    }

    dmp = (struct dm *)gedp->ged_dmp;

    /* must be wanting help */
    if (argc == 1) {
	_ged_cmd_help(gedp, usage, d);
	return GED_HELP;
    }

    argc-=(argc>0); argv+=(argc>0); /* done with command name argv[0] */

    int opt_ret = bu_opt_parse(NULL, argc, argv, d);

    if (print_help) {
	_ged_cmd_help(gedp, usage, d);
	return GED_HELP;
    }

    argc = opt_ret;

    if (grab_fb) {
	fbp = dm_get_fb(dmp);
	if (!fbp) {
	    bu_vls_printf(gedp->ged_result_str, ": display manager does not have a framebuffer");
	    return GED_ERROR;
	}
    }

    /* must be wanting help */
    if (!argc) {
	_ged_cmd_help(gedp, usage, d);
	return GED_HELP;
    }

    bytes_per_pixel = 3;
    bytes_per_line = width * bytes_per_pixel;

    /* create image file */
    if (!grab_fb) {

	dm_get_display_image(dmp, &idata);
	if (!idata) {
	    bu_vls_printf(gedp->ged_result_str, "%s: display manager did not return image data.", argv[1]);
	    return GED_ERROR;
	}
	bif = icv_create(dm_get_width(dmp), dm_get_height(dmp), ICV_COLOR_SPACE_RGB);
	if (bif == NULL) {
	    bu_vls_printf(gedp->ged_result_str, ": could not create icv_image write structure.");
	    return GED_ERROR;
	}
	rows = (unsigned char **)bu_calloc(dm_get_height(dmp), sizeof(unsigned char *), "rows");
	for (i = 0; i < dm_get_height(dmp); ++i) {
	    rows[i] = (unsigned char *)(idata + ((dm_get_height(dmp)-i-1)*bytes_per_line));
	    /* TODO : Add double type data to maintain resolution */
	    icv_writeline(bif, i, rows[i], ICV_DATA_UCHAR);
	}
	bu_free(rows, "rows");
	bu_free(idata, "image data");

    } else {
	bif = fb_write_icv(fbp, 0, 0, fb_getwidth(fbp), fb_getheight(fbp));
	if (bif == NULL) {
	    bu_vls_printf(gedp->ged_result_str, ": could not create icv_image from framebuffer.");
	    return GED_ERROR;
	}
    }

    /* Crop the image, if the settings indicate we need to */
    if (scr_xoff || scr_yoff || width || height) {
        width = (width) ? width : (int)bif->width - scr_xoff;
        height = (height) ? height : (int)bif->height - scr_yoff;
        icv_rect(bif, scr_xoff, scr_yoff, width, height);
    }

    icv_write(bif, argv[0], type);
    icv_destroy(bif);

    return GED_OK;
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
