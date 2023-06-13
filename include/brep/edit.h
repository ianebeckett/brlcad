/*                            E D I T . H
 * BRL-CAD
 *
 * Copyright (c) 2007-2023 United States Government as represented by
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
/** @addtogroup brep_edit
 * @brief
 * Implementation of edit support for brep.
 */
#ifndef BREP_EDIT_H
#define BREP_EDIT_H

#include "common.h"
#include "brep/defines.h"

/** @{ */
/** @file brep/edit.h */

/// function below are C interface
__BEGIN_DECLS

/**
 * create an empty ON_Brep.
 */
BREP_EXPORT extern void *
brep_create();

__END_DECLS

/// function below are C++ interface
#ifdef __cplusplus
extern "C++" {
#include <vector>

/**
 * create an nurbs curve using a template
 * position can be specified if argc == 3
 */
BREP_EXPORT extern ON_NurbsCurve *
brep_make_curve(int argc, const char **argv);

/**
 * create an nurbs curve given detailed information
 */
BREP_EXPORT extern ON_NurbsCurve *
brep_in_curve(int argc, const char **argv);

/**
 * move control vertex of a curve
 */
BREP_EXPORT extern bool
brep_curve_move_cv(ON_Brep* brep, int curve_id, int cv_id, ON_4dPoint point);

/**
 * Reverse parameterizatrion by negating all knots
 * and reversing the order of the control vertices.
 */
BREP_EXPORT extern bool
brep_curve_reverse(ON_Brep* brep, int curve_id);

/**
 * insert knots into a curve
 */
BREP_EXPORT extern bool
brep_curve_insert_knot(ON_Brep* brep, int curve_id, double knot, int multiplicity);

/**
 * trim a curve using a parameter range
 */
BREP_EXPORT extern bool
brep_curve_trim(ON_Brep* brep, int curve_id, double t0, double t1);

/**
 * create an nurbs curve using a template
 * position can be specified if argc == 3
 * return id of the surface
 */
BREP_EXPORT extern int
brep_make_surface(ON_Brep* brep, std::vector<double> position);

} /* extern C++ */
#endif

#endif  /* BREP_EDIT_H */
/** @} */
/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
