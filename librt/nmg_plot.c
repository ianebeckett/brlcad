/* XXX Move to raytrace.h */
#define RT_CK_LIST(_hd)	RT_CKMAG(_hd, RT_LIST_HEAD_MAGIC, "struct rt_list");
/* Near definition of RT_GET_VLIST(), note: */
/*
 *  Applications that are going to use RT_ADD_VLIST are required to
 *  first execute this macro:	RT_LIST_INIT( &rt_g.rtg_vlfree );
 */

/*
 *			N M G _ P L O T . C
 *
 *  This file contains routines that create VLISTs and UNIX-plot files.
 *  Some routines are essential to the MGED interface, some are
 *  more for diagnostic and visualization purposes.
 *
 *  There are several distinct families -
 *	nmg_ENTITY_to_vlist	Wireframes & polgyons.  For MGED "ev".
 *	nmg_pl_ENTITY		Fancy edgeuse drawing, to plot file.
 *	nmg_vlblock_ENTITY	Fancy edgeuse drawing, into vlblocks.
 *	show_broken_ENTITY	Graphical display of classifier results.
 *	...as well as assorted wrappers for debugging use.
 *
 *  In the interest of having only a single way of creating the fancy
 *  drawings, the code is migrating to creating everything first as
 *  VLBLOCKS, and converting that to UNIX-plot files or other formats
 *  as appropriate.
 *
 *  Authors -
 *	Lee A. Butler
 *	Michael John Muuss
 *  
 *  Source -
 *	The U. S. Army Research Laboratory
 *	Aberdeen Proving Ground, Maryland  21005-5068  USA
 *  
 *  Distribution Notice -
 *	Re-distribution of this software is restricted, as described in
 *	your "Statement of Terms and Conditions for the Release of
 *	The BRL-CAD Pacakge" agreement.
 *
 *  Copyright Notice -
 *	This software is Copyright (C) 1993 by the United States Army
 *	in all countries except the USA.  All rights reserved.
 */
#ifndef lint
static char RCSid[] = "@(#)$Header$ (ARL)";
#endif

#include "conf.h"
#include <stdio.h>
#include <fcntl.h>
#include <signal.h>
#include "machine.h"
#include "vmath.h"
#include "nmg.h"
#include "raytrace.h"
#include "nurb.h"
#include "plot3.h"
#include "rtstring.h"

#define US_DELAY	10	/* Additional delay between frames */

void		(*nmg_plot_anim_upcall)();	/* For I/F with MGED */
void		(*nmg_vlblock_anim_upcall)();	/* For I/F with MGED */
void		(*nmg_mged_debug_display_hack)();
double nmg_eue_dist = 0.05;

/************************************************************************
 *									*
 *		NMG to VLIST routines, for MGED "ev" command.		*
 * XXX should take a flags array, to ensure each item done only once!   *
 *									*
 ************************************************************************/

/*
 *			N M G _ V U _ T O _ V L I S T
 *
 *  Plot a single vertexuse
 */
void
nmg_vu_to_vlist( vhead, vu )
struct rt_list		*vhead;
CONST struct vertexuse	*vu;
{
	struct vertex	*v;
	register struct vertex_g *vg;

	NMG_CK_VERTEXUSE(vu);
	v = vu->v_p;
	NMG_CK_VERTEX(v);
	vg = v->vg_p;
	if( vg )  {
		/* Only thing in this shell is a point */
		NMG_CK_VERTEX_G(vg);
		RT_ADD_VLIST( vhead, vg->coord, RT_VLIST_LINE_MOVE );
		RT_ADD_VLIST( vhead, vg->coord, RT_VLIST_LINE_DRAW );
	}
}

/*
 *			N M G _ E U _ T O _ V L I S T
 *
 *  Plot a list of edgeuses.  The last edge is joined back to the first.
 */
void
nmg_eu_to_vlist( vhead, eu_hd )
struct rt_list		*vhead;
CONST struct rt_list	*eu_hd;
{
	struct edgeuse		*eu;
	struct edgeuse		*eumate;
	struct vertexuse	*vu;
	struct vertexuse	*vumate;
	register struct vertex_g *vg;
	register struct vertex_g *vgmate;

	/* Consider all the edges in the wire edge list */
	for( RT_LIST_FOR( eu, edgeuse, eu_hd ) )  {
		/* This wire edge runs from vertex to mate's vertex */
		NMG_CK_EDGEUSE(eu);
		vu = eu->vu_p;
		NMG_CK_VERTEXUSE(vu);
		NMG_CK_VERTEX(vu->v_p);
		vg = vu->v_p->vg_p;

		eumate = eu->eumate_p;
		NMG_CK_EDGEUSE(eumate);
		vumate = eumate->vu_p;
		NMG_CK_VERTEXUSE(vumate);
		NMG_CK_VERTEX(vumate->v_p);
		vgmate = vumate->v_p->vg_p;

		if( !vg || !vgmate ) {
			rt_log("nmg_eu_to_vlist() no vg or mate?\n");
			continue;
		}
		NMG_CK_VERTEX_G(vg);
		NMG_CK_VERTEX_G(vgmate);

		RT_ADD_VLIST( vhead, vg->coord, RT_VLIST_LINE_MOVE );
		RT_ADD_VLIST( vhead, vgmate->coord, RT_VLIST_LINE_DRAW );
	}
}

/*
 *			N M G _ L U _ T O _ V L I S T
 *
 *  Plot a single loopuse into a rt_vlist chain headed by vhead.
 */
void
nmg_lu_to_vlist( vhead, lu, poly_markers, normal )
struct rt_list		*vhead;
CONST struct loopuse	*lu;
int			poly_markers;
CONST vectp_t		normal;
{
	CONST struct edgeuse		*eu;
	CONST struct vertexuse		*vu;
	CONST struct vertex		*v;
	register CONST struct vertex_g	*vg;
	CONST struct vertex_g		*first_vg;
	int		isfirst;
	point_t		centroid;
	int		npoints;

	RT_CK_LIST(vhead);

	NMG_CK_LOOPUSE(lu);
	if( RT_LIST_FIRST_MAGIC(&lu->down_hd)==NMG_VERTEXUSE_MAGIC )  {
		/* Process a loop of a single vertex */
		vu = RT_LIST_FIRST(vertexuse, &lu->down_hd);
		nmg_vu_to_vlist( vhead, vu );
		return;
	}

	/* Consider all the edges in the loop */
	isfirst = 1;
	first_vg = (struct vertex_g *)0;
	npoints = 0;
	VSETALL( centroid, 0 );
	for( RT_LIST_FOR( eu, edgeuse, &lu->down_hd ) )  {
		/* Consider this edge */
		NMG_CK_EDGEUSE(eu);
		vu = eu->vu_p;
		NMG_CK_VERTEXUSE(vu);
		v = vu->v_p;
		NMG_CK_VERTEX(v);
		vg = v->vg_p;
		if( !vg ) {
			continue;
		}
		NMG_CK_VERTEX_G(vg);
		VADD2( centroid, centroid, vg->coord );
		npoints++;
		if (isfirst) {
			if( poly_markers) {
				/* Insert a "start polygon, normal" marker */
				RT_ADD_VLIST( vhead, normal, RT_VLIST_POLY_START );
				RT_ADD_VLIST( vhead, vg->coord, RT_VLIST_POLY_MOVE );
			} else {
				/* move */
				RT_ADD_VLIST( vhead, vg->coord, RT_VLIST_LINE_MOVE );
			}
			isfirst = 0;
			first_vg = vg;
		} else {
			if( poly_markers) {
				RT_ADD_VLIST( vhead, vg->coord, RT_VLIST_POLY_DRAW );
			} else {
				/* Draw */
				RT_ADD_VLIST( vhead, vg->coord, RT_VLIST_LINE_DRAW );
			}
		}
	}

	/* Draw back to the first vertex used */
	if( !isfirst && first_vg )  {
		if( poly_markers )  {
			/* Draw, end polygon */
			RT_ADD_VLIST( vhead, first_vg->coord, RT_VLIST_POLY_END );
		} else {
			/* Draw */
			RT_ADD_VLIST( vhead, first_vg->coord, RT_VLIST_LINE_DRAW );
		}
	}
	if( poly_markers > 1 && npoints > 2 )  {
		/* Draw surface normal as a little vector */
		double	f;
		vect_t	tocent;
		point_t	tip;
		f = 1.0 / npoints;
		VSCALE( centroid, centroid, f );
		RT_ADD_VLIST( vhead, centroid, RT_VLIST_LINE_MOVE );
		VSUB2( tocent, first_vg->coord, centroid );
		f = MAGNITUDE( tocent ) * 0.5;
		VJOIN1( tip, centroid, f, normal );
		RT_ADD_VLIST( vhead, tip, RT_VLIST_LINE_DRAW );

		/* For any vertexuse attributes with normals, draw them too */
		for( RT_LIST_FOR( eu, edgeuse, &lu->down_hd ) )  {
			struct vertexuse_a_plane	*vua;
			/* Consider this edge */
			vu = eu->vu_p;
			if( !vu->a.magic_p || *vu->a.magic_p != NMG_VERTEXUSE_A_PLANE_MAGIC )  continue;
			vua = vu->a.plane_p;
			v = vu->v_p;
			vg = v->vg_p;
			if( !vg )  continue;
			NMG_CK_VERTEX_G(vg);
			RT_ADD_VLIST( vhead, vg->coord, RT_VLIST_LINE_MOVE );
			VJOIN1( tip, vg->coord, f, vua->N );
			RT_ADD_VLIST( vhead, tip, RT_VLIST_LINE_DRAW );
		}
	}
}

/*
 *			N M G _ S N U R B _ F U _ T O _ V L I S T
 */
void
nmg_snurb_fu_to_vlist( vhead, fu, poly_markers )
struct rt_list		*vhead;
CONST struct faceuse	*fu;
int			poly_markers;
{
	struct face_g_snurb	*fg;

	RT_CK_LIST(vhead);

	NMG_CK_FACEUSE(fu);
	NMG_CK_FACE(fu->f_p);
	fg = fu->f_p->g.snurb_p;
	NMG_CK_FACE_G_SNURB(fg);

	/* First step, draw the surface. */
	/* XXX For now, draw the whole surface, not just the interior */
	nmg_snurb_to_vlist( vhead, fg, 10 );

	/* Second step, draw the trimming curves */
}

/*
 *			N M G _ S _ T O _ V L I S T
 *
 *  Plot the entire contents of a shell.
 *
 *  poly_markers =
 *	0 for vectors
 *	1 for polygons
 *	2 for polygons and surface normals drawn with vectors
 */
void
nmg_s_to_vlist( vhead, s, poly_markers )
struct rt_list		*vhead;
CONST struct shell	*s;
int			poly_markers;
{
	struct faceuse	*fu;
	struct face_g_plane	*fg;
	register struct loopuse	*lu;
	vect_t		normal;

	NMG_CK_SHELL(s);

	/* faces */
	for( RT_LIST_FOR( fu, faceuse, &s->fu_hd ) )  {
		vect_t		n;

		/* Consider this face */
		NMG_CK_FACEUSE(fu);
		if (fu->orientation != OT_SAME)  continue;
		NMG_CK_FACE(fu->f_p);

		if( fu->f_p->g.magic_p && *fu->f_p->g.magic_p == NMG_FACE_G_SNURB_MAGIC )  {
			nmg_snurb_fu_to_vlist( vhead, fu, poly_markers );
			continue;
		}

		/* Handle planar faces directly */
		fg = fu->f_p->g.plane_p;
		NMG_CK_FACE_G_PLANE(fg);
		NMG_GET_FU_NORMAL( n, fu );
		for( RT_LIST_FOR( lu, loopuse, &fu->lu_hd ) )  {
		   	nmg_lu_to_vlist( vhead, lu, poly_markers, n );
		}
	}

	/* wire loops.  poly_markers=0 so wires are always drawn as vectors */
	VSETALL(normal, 0);
	for( RT_LIST_FOR( lu, loopuse, &s->lu_hd ) )  {
		nmg_lu_to_vlist( vhead, lu, 0, normal );
	}

	/* wire edges */
	nmg_eu_to_vlist( vhead, &s->eu_hd );

	/* single vertices */
	if (s->vu_p)  {
		nmg_vu_to_vlist( vhead, s->vu_p );
	}
}

/*
 *			N M G _ R _ T O _ V L I S T
 */
void
nmg_r_to_vlist( vhead, r, poly_markers )
struct rt_list		*vhead;
CONST struct nmgregion	*r;
int			poly_markers;
{
	register struct shell	*s;

	NMG_CK_REGION( r );
	for( RT_LIST_FOR( s, shell, &r->s_hd ) )  {
		nmg_s_to_vlist( vhead, s, poly_markers );
	}
}

/*
 *			N M G _ M _ T O _ V L I S T
 *
 */
void
nmg_m_to_vlist( vhead, m, poly_markers )
struct rt_list	*vhead;
struct model	*m;
int		poly_markers;
{
	register struct nmgregion	*r;

	NMG_CK_MODEL( m );
	for( RT_LIST_FOR( r, nmgregion, &m->r_hd ) )  {
		NMG_CK_REGION( r );
		nmg_r_to_vlist( vhead, r, poly_markers );
	}
}
/************************************************************************
 *									*
 *		Routines to lay out the fancy edgeuse drawings		*
 *									*
 ************************************************************************/

#define LEE_DIVIDE_TOL	(1.0e-5)	/* sloppy tolerance */


/*
 *			N M G _ O F F S E T _ E U _ V E R T
 *
 *	Given an edgeuse, find an offset for its vertexuse which will place
 *	it "above" and "inside" the area of the face.
 *
 *  The point will be offset inwards along the edge
 *  slightly, to avoid obscuring the vertex, and will be offset off the
 *  face (in the direction of the face normal) slightly, to avoid
 *  obscuring the edge itself.
 *	
 */
void
nmg_offset_eu_vert(base, eu, face_normal, tip)
point_t			base;
CONST struct edgeuse	*eu;
CONST vect_t		face_normal;
int			tip;
{
	struct edgeuse	*prev_eu;
	CONST struct edgeuse	*this_eu;
	vect_t		prev_vec;	/* from cur_pt to prev_pt */
	vect_t		eu_vec;		/* from cur_pt to next_pt */
	vect_t		prev_left;
	vect_t		eu_left;
	vect_t		delta_vec;	/* offset vector from vertex */
	struct vertex_g	*this_vg, *mate_vg, *prev_vg;

	bzero(delta_vec, sizeof(vect_t)),
	prev_eu = RT_LIST_PLAST_CIRC( edgeuse, eu ); 
	this_eu = eu;

	NMG_CK_EDGEUSE(this_eu);
	NMG_CK_VERTEXUSE(this_eu->vu_p);
	NMG_CK_VERTEX(this_eu->vu_p->v_p);
    	this_vg = this_eu->vu_p->v_p->vg_p;
	NMG_CK_VERTEX_G(this_vg);

	NMG_CK_EDGEUSE(this_eu->eumate_p);
	NMG_CK_VERTEXUSE(this_eu->eumate_p->vu_p);
	NMG_CK_VERTEX(this_eu->eumate_p->vu_p->v_p);
    	mate_vg = this_eu->eumate_p->vu_p->v_p->vg_p;
	NMG_CK_VERTEX_G(mate_vg);

	NMG_CK_EDGEUSE(prev_eu);
	NMG_CK_VERTEXUSE(prev_eu->vu_p);
	NMG_CK_VERTEX(prev_eu->vu_p->v_p);
    	prev_vg = prev_eu->vu_p->v_p->vg_p;
	NMG_CK_VERTEX_G(prev_vg);

	/* get "left" vector for edgeuse */
	VSUB2(eu_vec, mate_vg->coord, this_vg->coord); 
	VUNITIZE(eu_vec);
	VCROSS(eu_left, face_normal, eu_vec);


	/* get "left" vector for previous edgeuse */
	VSUB2(prev_vec, this_vg->coord, prev_vg->coord);
	VUNITIZE(prev_vec);
	VCROSS(prev_left, face_normal, prev_vec);

	/* get "delta" vector to apply to vertex */
	VADD2(delta_vec, prev_left, eu_left);

	if (MAGSQ(delta_vec) > VDIVIDE_TOL) {
		VUNITIZE(delta_vec);
		VJOIN2(base, this_vg->coord,
			(nmg_eue_dist*1.3),delta_vec,
			(nmg_eue_dist*0.8),face_normal);
		
	} else if (tip) {
		VJOIN2(base, this_vg->coord,
			(nmg_eue_dist*1.3),prev_left,
			(nmg_eue_dist*0.8),face_normal);
	} else {
		VJOIN2(base, this_vg->coord,
			(nmg_eue_dist*1.3),eu_left,
			(nmg_eue_dist*0.8),face_normal);
	}
}



/*			N M G _ E U _ C O O R D S
 *
 *  Get the two (offset and shrunken) endpoints that represent
 *  an edgeuse.
 *  Return the base point, and a point 60% along the way towards the
 *  other end.
 */
static void nmg_eu_coords(eu, base, tip60)
CONST struct edgeuse *eu;
point_t base, tip60;
{
	point_t	tip;

	NMG_CK_EDGEUSE(eu);

	if (*eu->up.magic_p == NMG_SHELL_MAGIC ||
	    (*eu->up.magic_p == NMG_LOOPUSE_MAGIC &&
	    *eu->up.lu_p->up.magic_p == NMG_SHELL_MAGIC) ) {
	    	/* Wire edge, or edge in wire loop */
	    	VMOVE( base, eu->vu_p->v_p->vg_p->coord );
		NMG_CK_EDGEUSE(eu->eumate_p);
		VMOVE( tip, eu->eumate_p->vu_p->v_p->vg_p->coord );
	}
	else if (*eu->up.magic_p == NMG_LOOPUSE_MAGIC &&
	    *eu->up.lu_p->up.magic_p == NMG_FACEUSE_MAGIC) {
	    	/* Loop in face */
	    	struct faceuse	*fu;
	    	vect_t		face_normal;

	    	fu = eu->up.lu_p->up.fu_p;
	    	NMG_GET_FU_NORMAL( face_normal, fu );


	    	nmg_offset_eu_vert(base, eu, face_normal, 0);
		nmg_offset_eu_vert(tip, RT_LIST_PNEXT_CIRC(edgeuse, eu),
			face_normal, 1);

	} else
		rt_bomb("nmg_eu_coords: bad edgeuse up. What's going on?\n");

	VBLEND2( tip60, 0.4, base, 0.6, tip );
}

/*
 *			N M G _ E U _ R A D I A L
 *
 *  Find location for 80% tip on edgeuse's radial edgeuse.
 */
static void nmg_eu_radial(eu, tip)
CONST struct edgeuse *eu;
point_t tip;
{
	point_t	b2, t2;

	NMG_CK_EDGEUSE(eu->radial_p);
	NMG_CK_VERTEXUSE(eu->radial_p->vu_p);
	NMG_CK_VERTEX(eu->radial_p->vu_p->v_p);
	NMG_CK_VERTEX_G(eu->radial_p->vu_p->v_p->vg_p);

	nmg_eu_coords(eu->radial_p, b2, t2);

	/* find point 80% along other eu where radial pointer should touch */
	VCOMB2( tip, 0.8, t2, 0.2, b2 );
}

/*
 *			N M G _ E U _ L A S T
 *
 *  Find the tip of the last (previous) edgeuse from 'eu'.
 */
static void nmg_eu_last( eu, tip_out )
CONST struct edgeuse	*eu;
point_t		tip_out;
{
	point_t		radial_base;
	point_t		radial_tip;
	point_t		last_base;
	point_t		last_tip;
	point_t		p;
	struct edgeuse	*eulast;

	NMG_CK_EDGEUSE(eu);
	eulast = RT_LIST_PLAST_CIRC( edgeuse, eu );
	NMG_CK_EDGEUSE(eulast);
	NMG_CK_VERTEXUSE(eulast->vu_p);
	NMG_CK_VERTEX(eulast->vu_p->v_p);
	NMG_CK_VERTEX_G(eulast->vu_p->v_p->vg_p);

	nmg_eu_coords(eulast->radial_p, radial_base, radial_tip);

	/* find pt 80% along LAST eu's radial eu where radial ptr touches */
	VCOMB2( p, 0.8, radial_tip, 0.2, radial_base );

	/* get coordinates of last edgeuse */
	nmg_eu_coords(eulast, last_base, last_tip);

	/* Find pt 80% along other eu where last pointer should touch */
	VCOMB2( tip_out, 0.8, last_tip, 0.2, p );
}

/*
 *			N M G _ E U _ N E X T
 *
 *  Return the base of the next edgeuse
 */
static void nmg_eu_next_base( eu, next_base)
CONST struct edgeuse	*eu;
point_t		next_base;
{
	point_t	t2;
	register struct edgeuse	*nexteu;

	NMG_CK_EDGEUSE(eu);
	nexteu = RT_LIST_PNEXT_CIRC( edgeuse, eu );
	NMG_CK_EDGEUSE(nexteu);
	NMG_CK_VERTEXUSE(nexteu->vu_p);
	NMG_CK_VERTEX(nexteu->vu_p->v_p);
	NMG_CK_VERTEX_G(nexteu->vu_p->v_p->vg_p);

	nmg_eu_coords(nexteu, next_base, t2);
}

/************************************************************************
 *									*
 *		NMG to UNIX-Plot routines, for visualization		*
 *  XXX These should get replaced with calls to the vlblock routines	*
 *									*
 ************************************************************************/

/*
 *			N M G _ P L _ V
 */
void
nmg_pl_v(fp, v, b)
FILE			*fp;
CONST struct vertex	*v;
long			*b;
{
	pointp_t p;
	static char label[128];

	NMG_INDEX_RETURN_IF_SET_ELSE_SET(b, v->index);

	NMG_CK_VERTEX(v);
	NMG_CK_VERTEX_G(v->vg_p);
	p = v->vg_p->coord;

	pl_color(fp, 255, 255, 255);
	if (rt_g.NMG_debug & DEBUG_LABEL_PTS) {
		(void)sprintf(label, "%g %g %g", p[0], p[1], p[2]);
		pdv_3move( fp, p );
		pl_label(fp, label);
	}
	pdv_3point(fp, p);
}

/*
 *			N M G _ P L _ E
 */
void
nmg_pl_e(fp, e, b, red, green, blue)
FILE			*fp;
CONST struct edge	*e;
long			*b;
int			red, green, blue;
{
	pointp_t	p0, p1;
	point_t		end0, end1;
	vect_t		v;

	NMG_INDEX_RETURN_IF_SET_ELSE_SET(b, e->index);
	
	NMG_CK_EDGEUSE(e->eu_p);
	NMG_CK_VERTEXUSE(e->eu_p->vu_p);
	NMG_CK_VERTEX(e->eu_p->vu_p->v_p);
	NMG_CK_VERTEX_G(e->eu_p->vu_p->v_p->vg_p);
	p0 = e->eu_p->vu_p->v_p->vg_p->coord;

	NMG_CK_VERTEXUSE(e->eu_p->eumate_p->vu_p);
	NMG_CK_VERTEX(e->eu_p->eumate_p->vu_p->v_p);
	NMG_CK_VERTEX_G(e->eu_p->eumate_p->vu_p->v_p->vg_p);
	p1 = e->eu_p->eumate_p->vu_p->v_p->vg_p->coord;

	/* leave a little room between the edge endpoints and the vertex
	 * compute endpoints by forming a vector between verets, scale vector
	 * and modify points
	 */
	VSUB2SCALE(v, p1, p0, 0.95);
	VADD2(end0, p0, v);
	VSUB2(end1, p1, v);

	pl_color(fp, red, green, blue);
	pdv_3line( fp, end0, end1 );

	nmg_pl_v(fp, e->eu_p->vu_p->v_p, b);
	nmg_pl_v(fp, e->eu_p->eumate_p->vu_p->v_p, b);
}

/*
 *			M N G _ P L _ E U
 */
void
nmg_pl_eu(fp, eu, b, red, green, blue)
FILE			*fp;
CONST struct edgeuse	*eu;
long			*b;
int			red, green, blue;
{
	point_t base, tip;
	point_t	radial_tip;
	point_t	next_base;
	point_t	last_tip;

	NMG_CK_EDGEUSE(eu);
	NMG_CK_EDGE(eu->e_p);
	NMG_CK_VERTEXUSE(eu->vu_p);
	NMG_CK_VERTEX(eu->vu_p->v_p);
	NMG_CK_VERTEX_G(eu->vu_p->v_p->vg_p);

	NMG_CK_VERTEXUSE(eu->eumate_p->vu_p);
	NMG_CK_VERTEX(eu->eumate_p->vu_p->v_p);
	NMG_CK_VERTEX_G(eu->eumate_p->vu_p->v_p->vg_p);

	NMG_INDEX_RETURN_IF_SET_ELSE_SET(b, eu->index);

	nmg_pl_e(fp, eu->e_p, b, red, green, blue);

	if (*eu->up.magic_p == NMG_LOOPUSE_MAGIC &&
	    *eu->up.lu_p->up.magic_p == NMG_FACEUSE_MAGIC) {

	    	nmg_eu_coords(eu, base, tip);
	    	if (eu->up.lu_p->up.fu_p->orientation == OT_SAME)
	    		red += 50;
		else if (eu->up.lu_p->up.fu_p->orientation == OT_OPPOSITE)
			red -= 50;
	    	else
	    		red = green = blue = 255;

		pl_color(fp, red, green, blue);
	    	pdv_3line( fp, base, tip );

	    	nmg_eu_radial( eu, radial_tip );
		pl_color(fp, red, green-20, blue);
	    	pdv_3line( fp, tip, radial_tip );

		pl_color(fp, 0, 100, 0);
	    	nmg_eu_next_base( eu, next_base );
	    	pdv_3line( fp, tip, next_base );

/*** presently unused ***
	    	nmg_eu_last( eu, last_tip );
		pl_color(fp, 0, 200, 0);
	    	pdv_3line( fp, base, last_tip );
****/
	    }
}

/*
 *			N M G _ P L _ L U
 */
void
nmg_pl_lu(fp, lu, b, red, green, blue)
FILE			*fp;
CONST struct loopuse	*lu;
long			*b;
int			red, green, blue;
{
#if 0
	struct edgeuse	*eu;
	long		magic1;

	NMG_CK_LOOPUSE(lu);
	NMG_INDEX_RETURN_IF_SET_ELSE_SET(b, lu->index);

	magic1 = RT_LIST_FIRST_MAGIC( &lu->down_hd );
	if (magic1 == NMG_VERTEXUSE_MAGIC &&
	    lu->orientation != OT_BOOLPLACE) {
	    	nmg_pl_v(fp, RT_LIST_PNEXT(vertexuse, &lu->down_hd)->v_p, b);
	} else if (magic1 == NMG_EDGEUSE_MAGIC) {
		for( RT_LIST_FOR( eu, edgeuse, &lu->down_hd ) )  {
			nmg_pl_eu(fp, eu, b, red, green, blue);
		}
	}
#else
	struct rt_vlblock	*vbp;

	vbp = rt_vlblock_init();
	nmg_vlblock_lu(vbp, lu, b, red, green, blue, 0, 0);
	rt_plot_vlblock(fp, vbp);
	rt_vlblock_free(vbp);
#endif
}

/*
 *			M N G _ P L _ F U
 */
void
nmg_pl_fu(fp, fu, b, red, green, blue)
FILE			*fp;
CONST struct faceuse	*fu;
long			*b;
int			red, green, blue;
{
#if 0
	struct loopuse *lu;

	NMG_CK_FACEUSE(fu);

	NMG_INDEX_RETURN_IF_SET_ELSE_SET(b, fu->index);

	for( RT_LIST_FOR( lu, loopuse, &fu->lu_hd ) )  {

		nmg_pl_lu(fp, lu, b, red, green, blue);
	}
#else
	struct loopuse		*lu;
	struct rt_vlblock	*vbp;
	int 		loopnum = 0;

	NMG_CK_FACEUSE(fu);
	NMG_INDEX_RETURN_IF_SET_ELSE_SET(b, fu->index);

	vbp = rt_vlblock_init();

	for( RT_LIST_FOR( lu, loopuse, &fu->lu_hd ) )  {
		nmg_vlblock_lu(vbp, lu, b, red, green, blue, 1, loopnum++);
	}

	rt_plot_vlblock(fp, vbp);
	rt_vlblock_free(vbp);
#endif
}

/*
 *			N M G _ P L _ S
 *
 *  Note that "b" should probably be defined a level higher,
 *  to reduce malloc/free calls when plotting multiple shells.
 */
void
nmg_pl_s(fp, s)
FILE			*fp;
CONST struct shell	*s;
{
#if 0
	struct faceuse *fu;
	struct loopuse *lu;
	struct edgeuse *eu;
	long		*b;

	NMG_CK_SHELL(s);
	if( s->sa_p )  {
		NMG_CK_SHELL_A( s->sa_p );
		pdv_3space( fp, s->sa_p->min_pt, s->sa_p->max_pt );
	}

	/* get space for flag array, to ensure each item is done once */
	b = (long *)rt_calloc( s->r_p->m_p->maxindex, sizeof(long),
		"nmg_pl_s flag[]" );

	for( RT_LIST_FOR( fu, faceuse, &s->fu_hd ) )  {
		NMG_CK_FACEUSE(fu);
		nmg_pl_fu(fp, fu, b, 80, 100, 170);
	}

	for( RT_LIST_FOR( lu, loopuse, &s->lu_hd ) )  {
		NMG_CK_LOOPUSE(lu);
		nmg_pl_lu(fp, lu, b, 255, 0, 0);
	}

	for( RT_LIST_FOR( eu, edgeuse, &s->eu_hd ) )  {
		NMG_CK_EDGEUSE(eu);
		NMG_CK_EDGE(eu->e_p);

		nmg_pl_eu(fp, eu, b, 200, 200, 0 );
	}
	if (s->vu_p) {
		nmg_pl_v(fp, s->vu_p->v_p, b );
	}

	if( RT_LIST_IS_EMPTY( &s->fu_hd ) &&
	    RT_LIST_IS_EMPTY( &s->lu_hd ) &&
	    RT_LIST_IS_EMPTY( &s->eu_hd ) && !s->vu_p) {
	    	rt_log("WARNING nmg_pl_s(): shell has no children\n");
	}

	rt_free( (char *)b, "nmg_pl_s flag[]" );
#else
	struct rt_vlblock	*vbp;

	vbp = rt_vlblock_init();
	nmg_vlblock_s(vbp, s, 0);
	rt_plot_vlblock(fp, vbp);
	rt_vlblock_free(vbp);
#endif
}

/*
 *			N M G _ P L _ R
 */
void
nmg_pl_r(fp, r)
FILE			*fp;
CONST struct nmgregion	*r;
{
	struct rt_vlblock	*vbp;

	vbp = rt_vlblock_init();
	nmg_vlblock_r(vbp, r, 0);
	rt_plot_vlblock(fp, vbp);
	rt_vlblock_free(vbp);
}

/*
 *			N M G _ P L _ M
 */
void
nmg_pl_m(fp, m)
FILE			*fp;
CONST struct model	*m;
{
	struct rt_vlblock	*vbp;

	vbp = rt_vlblock_init();
	nmg_vlblock_m(vbp, m, 0);
	rt_plot_vlblock(fp, vbp);
	rt_vlblock_free(vbp);
}

/************************************************************************
 *									*
 *		Visualization of fancy edgeuses into VLBLOCKs		*
 *									*
 *  This is the preferred method of obtaining fancy NMG displays.	*
 *									*
 ************************************************************************/

/*
 *			N M G _ V L B L O C K _ V
 */
void
nmg_vlblock_v(vbp, v, tab)
struct rt_vlblock		*vbp;
CONST struct vertex		*v;
long				*tab;
{
	pointp_t p;
	static char label[128];
	struct rt_list	*vh;

	RT_CK_VLBLOCK(vbp);
	NMG_CK_VERTEX(v);
	NMG_INDEX_RETURN_IF_SET_ELSE_SET( tab, v->index );

	NMG_CK_VERTEX_G(v->vg_p);
	p = v->vg_p->coord;

	vh = rt_vlblock_find( vbp, 255, 255, 255 );
#if 0
	if (rt_g.NMG_debug & DEBUG_LABEL_PTS) {
		mat_t	mat;
		mat_idn(mat);
		(void)sprintf(label, "%g %g %g", p[0], p[1], p[2]);
		/* XXX What size characters to use? */
		rt_vlist_3string( vh, label, p, mat, scale );
	}
#endif
	RT_ADD_VLIST( vh, p, RT_VLIST_LINE_MOVE );
	RT_ADD_VLIST( vh, p, RT_VLIST_LINE_DRAW );
}

/*
 *			N M G _ V L B L O C K _ E
 */
void
nmg_vlblock_e(vbp, e, tab, red, green, blue, fancy)
struct rt_vlblock	*vbp;
CONST struct edge	*e;
long			*tab;
int			red, green, blue;
int			fancy;
{
	pointp_t p0, p1;
	point_t end0, end1;
	vect_t v;
	struct rt_list	*vh;

	RT_CK_VLBLOCK(vbp);
	NMG_CK_EDGE(e);
	NMG_INDEX_RETURN_IF_SET_ELSE_SET( tab, e->index );
	
	NMG_CK_EDGEUSE(e->eu_p);
	NMG_CK_VERTEXUSE(e->eu_p->vu_p);
	NMG_CK_VERTEX(e->eu_p->vu_p->v_p);
	NMG_CK_VERTEX_G(e->eu_p->vu_p->v_p->vg_p);
	p0 = e->eu_p->vu_p->v_p->vg_p->coord;

	NMG_CK_VERTEXUSE(e->eu_p->eumate_p->vu_p);
	NMG_CK_VERTEX(e->eu_p->eumate_p->vu_p->v_p);
	NMG_CK_VERTEX_G(e->eu_p->eumate_p->vu_p->v_p->vg_p);
	p1 = e->eu_p->eumate_p->vu_p->v_p->vg_p->coord;

	/* leave a little room between the edge endpoints and the vertex
	 * compute endpoints by forming a vector between verets, scale vector
	 * and modify points
	 */
	VSUB2SCALE(v, p1, p0, 0.90);
	VADD2(end0, p0, v);
	VSUB2(end1, p1, v);

	vh = rt_vlblock_find( vbp, red, green, blue );
	RT_ADD_VLIST( vh, end0, RT_VLIST_LINE_MOVE );
	RT_ADD_VLIST( vh, end1, RT_VLIST_LINE_DRAW );

	nmg_vlblock_v(vbp, e->eu_p->vu_p->v_p, tab);
	nmg_vlblock_v(vbp, e->eu_p->eumate_p->vu_p->v_p, tab);
}

/*
 *			M N G _ V L B L O C K _ E U
 */
void
nmg_vlblock_eu(vbp, eu, tab, red, green, blue, fancy, loopnum)
struct rt_vlblock		*vbp;
CONST struct edgeuse		*eu;
long				*tab;
int				red, green, blue;
int				fancy;
int				loopnum;
{
	point_t base, tip;
	point_t	radial_tip;
	point_t	next_base;
	point_t	last_tip;
	struct rt_list	*vh;

	RT_CK_VLBLOCK(vbp);
	NMG_CK_EDGEUSE(eu);
	NMG_INDEX_RETURN_IF_SET_ELSE_SET( tab, eu->index );

	NMG_CK_EDGE(eu->e_p);
	NMG_CK_VERTEXUSE(eu->vu_p);
	NMG_CK_VERTEX(eu->vu_p->v_p);
	NMG_CK_VERTEX_G(eu->vu_p->v_p->vg_p);

	NMG_CK_VERTEXUSE(eu->eumate_p->vu_p);
	NMG_CK_VERTEX(eu->eumate_p->vu_p->v_p);
	NMG_CK_VERTEX_G(eu->eumate_p->vu_p->v_p->vg_p);

	nmg_vlblock_e(vbp, eu->e_p, tab, red, green, blue, fancy);

	if( !fancy )  return;

	if (*eu->up.magic_p == NMG_LOOPUSE_MAGIC &&
	    *eu->up.lu_p->up.magic_p == NMG_FACEUSE_MAGIC) {

	    	/* if "fancy" doesn't specify plotting edgeuses of this
	    	 * particular face orientation, return
	    	 */
	    	if ( (eu->up.lu_p->up.fu_p->orientation == OT_SAME &&
	    	     (fancy & 1) == 0) ||
		     (eu->up.lu_p->up.fu_p->orientation == OT_OPPOSITE &&
		     (fancy & 2) == 0) )
	    		return;

	    	nmg_eu_coords(eu, base, tip);
	    	/* draw edgeuses of an OT_SAME faceuse in bright green,
	    	 * and  edgeuses of an OT_OPPOSITE faceuse in cyan.
	    	 * WIRE/UNSPEC edgeuses are drawn white.
	    	 */
	    	if (eu->up.lu_p->up.fu_p->orientation == OT_SAME) {
	    		if (eu->up.lu_p->orientation == OT_SAME) {
	    			/* green */
	    			red = 75;
	    			green = 250;
	    			blue = 75;
	    		} else if (eu->up.lu_p->orientation == OT_OPPOSITE) {
	    			/* yellow */
	    			red = 250;
	    			green = 250;
	    			blue = 75;
	    		} else {
	    			red = 250;
	    			green = 50;
	    			blue = 250;
	    		}
		} else if (eu->up.lu_p->up.fu_p->orientation == OT_OPPOSITE) {
	    		if (eu->up.lu_p->orientation == OT_SAME) {
	    			/* blue */
	    			red = 100;
	    			green = 100;
	    			blue = 250;
	    		} else if (eu->up.lu_p->orientation == OT_OPPOSITE) {
	    			/* cyan */
	    			red = 200;
	    			green = 100;
	    			blue = 250;
	    		} else {
	    			/* dark magenta */
	    			red = 125;
	    			green = 0;
	    			blue = 125;
	    		}
	    	} else
	    		red = green = blue = 255;

	    	/* draw the portion from the vertexuse to just beyond the
	    	 * midway point to represent the edgeuse
	    	 */
		vh = rt_vlblock_find( vbp, red, green, blue );
		RT_ADD_VLIST( vh, base, RT_VLIST_LINE_MOVE );
		RT_ADD_VLIST( vh, tip, RT_VLIST_LINE_DRAW );

	    	/* draw a line from the tip of the edgeuse part to a point
	    	 * behind the tip of the radial edgeuse.  This provides 2
	    	 * visual cues.  First it allows us to identify the radial
	    	 * edgeuse, and second, it makes a "half arrowhead" on the
	    	 * edgeuse, making it easier to recognize the direction
	    	 * of the edgeuse
	    	 */
	    	nmg_eu_radial( eu, radial_tip );
		vh = rt_vlblock_find( vbp, red, green-20, blue );
		RT_ADD_VLIST( vh, tip, RT_VLIST_LINE_MOVE );
		RT_ADD_VLIST( vh, radial_tip, RT_VLIST_LINE_DRAW );

	    	/* we draw a line from the tip of the edgeuse line
	    	 * to the vertexuse/start of the next edgeuse in the loop.
	    	 * This helps us to visually trace the loop from edgeuse to
	    	 * edgeuse.  The color of this part encodes the loopuse
	    	 * orientation.
	    	 */
	    	nmg_eu_next_base( eu, next_base );
	    	red *= 0.5;
	    	green *= 0.5;
	    	blue *= 0.5;
		vh = rt_vlblock_find( vbp, red, green, blue );
		RT_ADD_VLIST( vh, tip, RT_VLIST_LINE_MOVE );
		RT_ADD_VLIST( vh, next_base, RT_VLIST_LINE_DRAW );
	}
}

/*
 *			N M G _ V L B L O C K _ E U L E F T
 *
 *  Draw the left vector for this edgeuse.
 *  At the tip, write the angle around the edgeuse, in degrees.
 *
 *  Color is determined by caller.
 */
void
nmg_vlblock_euleft( vh, eu, center, mat, xvec, yvec, len, tol )
struct rt_list			*vh;
CONST struct edgeuse		*eu;
CONST point_t			center;
CONST mat_t			mat;
CONST vect_t			xvec;
CONST vect_t			yvec;
double				len;
CONST struct rt_tol		*tol;
{
	vect_t		left;
	point_t		tip;
	fastf_t		fan_len;
	fastf_t		char_scale;
	double		ang;
	char		str[128];

	NMG_CK_EDGEUSE(eu);
	RT_CK_TOL(tol);

	if( nmg_find_eu_leftvec( left, eu ) < 0 )  return;

	/* fan_len is baed on length of eu */
	fan_len = len * 0.2;
	VJOIN1( tip, center, fan_len, left );

	RT_ADD_VLIST( vh, center, RT_VLIST_LINE_MOVE );
	RT_ADD_VLIST( vh, tip, RT_VLIST_LINE_DRAW );

	ang = rt_angle_measure( left, xvec, yvec ) * rt_radtodeg;
	sprintf( str, "%g", ang );

	/* char_scale is based on length of eu */
	char_scale = len * 0.05;
	rt_vlist_3string( vh, str, tip, mat, char_scale );
}

/*
 *			N M G _ V L B L O C K _ A R O U N D _ E U
 *
 *  Given an edgeuse, plot all the edgeuses around the common edge.
 *  A graphical parallel to nmg_pr_fu_around_eu_vecs().
 *
 *  If the "fancy" flag is set, draw an angle fan around the edge midpoint,
 *  using the same angular reference as nmg_pr_fu_around_eu_vecs(), so
 *  that the printed output can be cross-referenced to this display.
 */
void
nmg_vlblock_around_eu(vbp, arg_eu, tab, fancy, tol )
struct rt_vlblock		*vbp;
CONST struct edgeuse		*arg_eu;
long				*tab;
int				fancy;
CONST struct rt_tol		*tol;
{
	CONST struct edgeuse		*orig_eu;
	register CONST struct edgeuse	*eu;
	vect_t			xvec, yvec, zvec;
	point_t			center;
	mat_t			mat;
	struct rt_list		*vh;
	fastf_t			len;

	RT_CK_VLBLOCK(vbp);
	NMG_CK_EDGEUSE(arg_eu);
	RT_CK_TOL(tol);

	if( fancy )  {
		VSUB2( xvec, arg_eu->eumate_p->vu_p->v_p->vg_p->coord,
			arg_eu->vu_p->v_p->vg_p->coord );
		len = MAGNITUDE( xvec );

		/* Erect coordinate system around eu */
		nmg_eu_2vecs_perp( xvec, yvec, zvec, arg_eu, tol );

		/*  Construct matrix to rotate characters from 2D drawing space
		 *  into model coordinates, oriented in plane perpendicular to eu.
		 */
		mat_zero( mat );
		mat[0] = xvec[X];
		mat[4] = xvec[Y];
		mat[8] = xvec[Z];

		mat[1] = yvec[X];
		mat[5] = yvec[Y];
		mat[9] = yvec[Z];

		mat[2] = zvec[X];
		mat[6] = zvec[Y];
		mat[10] = zvec[Z];
		mat[15] = 1;

		VADD2SCALE( center, arg_eu->vu_p->v_p->vg_p->coord,
			arg_eu->eumate_p->vu_p->v_p->vg_p->coord, 0.5 );

		/* Yellow, for now */
		vh = rt_vlblock_find( vbp, 255, 200, 0 );
	}

	orig_eu = arg_eu->eumate_p;

	eu = orig_eu;
	do {
		if(fancy) nmg_vlblock_euleft( vh, eu, center, mat, xvec, yvec, len, tol );

		nmg_vlblock_eu(vbp, eu, tab, 80, 100, 170, 3, 0);
		eu = eu->eumate_p;

		nmg_vlblock_eu(vbp, eu, tab, 80, 100, 170, 3, 0);
		eu = eu->radial_p;
	} while( eu != orig_eu );
}

/*
 *			N M G _ V L B L O C K _ L U
 */
void
nmg_vlblock_lu(vbp, lu, tab, red, green, blue, fancy, loopnum)
struct rt_vlblock	*vbp;
CONST struct loopuse	*lu;
long			*tab;
int			red, green, blue;
int			fancy;
int	 		loopnum;
{
	struct edgeuse	*eu;
	long		magic1;
	struct vertexuse *vu;

	RT_CK_VLBLOCK(vbp);
	NMG_CK_LOOPUSE(lu);
	NMG_INDEX_RETURN_IF_SET_ELSE_SET( tab, lu->index );

	magic1 = RT_LIST_FIRST_MAGIC( &lu->down_hd );
	if (magic1 == NMG_VERTEXUSE_MAGIC &&
	    lu->orientation != OT_BOOLPLACE) {
	    	vu = RT_LIST_PNEXT(vertexuse, &lu->down_hd);
	    	NMG_CK_VERTEXUSE(vu);
	    	nmg_vlblock_v(vbp, vu->v_p, tab);
	} else if (magic1 == NMG_EDGEUSE_MAGIC) {
		for( RT_LIST_FOR( eu, edgeuse, &lu->down_hd ) )  {
			nmg_vlblock_eu(vbp, eu, tab, red, green, blue,
				fancy, loopnum);
		}
	}
}

/*
 *			M N G _ V L B L O C K _ F U
 */
void
nmg_vlblock_fu(vbp, fu, tab, fancy)
struct rt_vlblock	*vbp;
CONST struct faceuse	*fu;
long			*tab;
int			fancy;
{
	struct loopuse *lu;
	int 		loopnum = 0;

	RT_CK_VLBLOCK(vbp);
	NMG_CK_FACEUSE(fu);
	NMG_INDEX_RETURN_IF_SET_ELSE_SET( tab, fu->index );

	for( RT_LIST_FOR( lu, loopuse, &fu->lu_hd ) )  {
		/* Draw in pale blue / purple */
		if( fancy )  {
			nmg_vlblock_lu(vbp, lu, tab, 80, 100, 170, fancy, loopnum++ );
		} else {
			/* Non-fancy */
			nmg_vlblock_lu(vbp, lu, tab, 80, 100, 170, 0, loopnum++ );
		}
	}
}

/*
 *			N M G _ V L B L O C K _ S
 */
void
nmg_vlblock_s(vbp, s, fancy)
struct rt_vlblock	*vbp;
CONST struct shell	*s;
int			fancy;
{
	struct faceuse *fu;
	struct loopuse *lu;
	struct edgeuse *eu;
	struct model	*m;
	long		*tab;

	RT_CK_VLBLOCK(vbp);
	NMG_CK_SHELL(s);
	NMG_CK_REGION(s->r_p);
	m = s->r_p->m_p;
	NMG_CK_MODEL(m);

	/* get space for list of items processed */
	tab = (long *)rt_calloc( m->maxindex+1, sizeof(long),
		"nmg_vlblock_s tab[]");

	for( RT_LIST_FOR( fu, faceuse, &s->fu_hd ) )  {
		NMG_CK_FACEUSE(fu);
		nmg_vlblock_fu(vbp, fu, tab, fancy );
	}

	for( RT_LIST_FOR( lu, loopuse, &s->lu_hd ) )  {
		NMG_CK_LOOPUSE(lu);
		if( fancy ) {
			nmg_vlblock_lu(vbp, lu, tab, 255, 0, 0, fancy, 0);
		} else {
			/* non-fancy, wire loops in red */
			nmg_vlblock_lu(vbp, lu, tab, 200, 0, 0, 0, 0);
		}
	}

	for( RT_LIST_FOR( eu, edgeuse, &s->eu_hd ) )  {
		NMG_CK_EDGEUSE(eu);
		NMG_CK_EDGE(eu->e_p);

		if( fancy )  {
			nmg_vlblock_eu(vbp, eu, tab, 200, 200, 0, fancy, 0 );
		} else {
			/* non-fancy, wire edges in yellow */
			nmg_vlblock_eu(vbp, eu, tab, 200, 200, 0, 0, 0); 
		}
	}
	if (s->vu_p) {
		nmg_vlblock_v(vbp, s->vu_p->v_p, tab );
	}

	rt_free( (char *)tab, "nmg_vlblock_s tab[]" );
}

/*
 *			N M G _ V L B L O C K _ R
 */
void
nmg_vlblock_r(vbp, r, fancy)
struct rt_vlblock	*vbp;
CONST struct nmgregion	*r;
int			fancy;
{
	struct shell *s;

	RT_CK_VLBLOCK(vbp);
	NMG_CK_REGION(r);

	for( RT_LIST_FOR( s, shell, &r->s_hd ) )  {
		nmg_vlblock_s(vbp, s, fancy);
	}
}

/*
 *			N M G _ V L B L O C K _ M
 */
void
nmg_vlblock_m(vbp, m, fancy)
struct rt_vlblock	*vbp;
CONST struct model	*m;
int			fancy;
{
	struct nmgregion *r;

	RT_CK_VLBLOCK(vbp);
	NMG_CK_MODEL(m);

	for( RT_LIST_FOR( r, nmgregion, &m->r_hd ) )  {
		nmg_vlblock_r(vbp, r, fancy);
	}
}

/************************************************************************
 *									*
 *		Visualization helper routines				*
 *									*
 ************************************************************************/

/*
 *  If another use of this edge is in another shell, plot all the
 *  uses around this edge.
 */
void
nmg_pl_edges_in_2_shells(vbp, b, eu, fancy, tol)
struct rt_vlblock	*vbp;
long			*b;
CONST struct edgeuse	*eu;
int			fancy;
CONST struct rt_tol	*tol;
{
	CONST struct edgeuse	*eur;
	CONST struct shell	*s;

	RT_CK_TOL(tol);
	eur = eu;
	NMG_CK_EDGEUSE(eu);
	NMG_CK_LOOPUSE(eu->up.lu_p);
	NMG_CK_FACEUSE(eu->up.lu_p->up.fu_p);
	s = eu->up.lu_p->up.fu_p->s_p;
	NMG_CK_SHELL(s);

	do {
		NMG_CK_EDGEUSE(eur);

		if (*eur->up.magic_p == NMG_LOOPUSE_MAGIC &&
		    *eur->up.lu_p->up.magic_p == NMG_FACEUSE_MAGIC &&
		    eur->up.lu_p->up.fu_p->s_p != s) {
		    	nmg_vlblock_around_eu(vbp, eu, b, fancy, tol);
		    	break;
		    }

		eur = eur->radial_p->eumate_p;
	} while (eur != eu);
}

/*
 *			N M G _ P L _ I S E C T
 *
 *  Called by nmg_bool.c
 */
void
nmg_pl_isect(filename, s, tol)
CONST char		*filename;
CONST struct shell	*s;
CONST struct rt_tol	*tol;
{
	struct faceuse		*fu;
	struct loopuse		*lu;
	struct edgeuse		*eu;
	long			*b;
	FILE			*fp;
	long			magic1;
	struct rt_vlblock	*vbp;

	NMG_CK_SHELL(s);
	RT_CK_TOL(tol);

	if ((fp=fopen(filename, "w")) == (FILE *)NULL) {
		(void)perror(filename);
		exit(-1);
	}

	b = (long *)rt_calloc( s->r_p->m_p->maxindex+1, sizeof(long),
		"nmg_pl_isect flags[]" );

	vbp = rt_vlblock_init();

	rt_log("Plotting to \"%s\"\n", filename);
	if( s->sa_p )  {
		NMG_CK_SHELL_A( s->sa_p );
#if 0
		pdv_3space( fp, s->sa_p->min_pt, s->sa_p->max_pt );
#endif
	}

	for( RT_LIST_FOR( fu, faceuse, &s->fu_hd ) )  {
		NMG_CK_FACEUSE(fu);
		for( RT_LIST_FOR( lu, loopuse, &fu->lu_hd ) )  {
			NMG_CK_LOOPUSE(lu);
			magic1 = RT_LIST_FIRST_MAGIC( &lu->down_hd );
			if (magic1 == NMG_EDGEUSE_MAGIC) {
				for( RT_LIST_FOR( eu, edgeuse, &lu->down_hd ) )  {
					NMG_CK_EDGEUSE(eu);
					nmg_pl_edges_in_2_shells(vbp, b, eu, 0, tol);
				}
			} else if (magic1 == NMG_VERTEXUSE_MAGIC) {
				;
			} else {
				rt_bomb("nmg_pl_isect() bad loopuse down\n");
			}
		}
	}

	rt_plot_vlblock(fp, vbp);
	rt_vlblock_free(vbp);

	rt_free( (char *)b, "nmg_pl_isect flags[]" );

	(void)fclose(fp);
}

/*
 *			N M G _ P L _ C O M B _ F U
 *
 *  Called from nmg_bool.c/nmg_face_combine()
 */
void
nmg_pl_comb_fu( num1, num2, fu1 )
int	num1;
int	num2;
CONST struct faceuse	*fu1;
{
	FILE			*fp;
	char			name[64];
	int			do_plot = 0;
	int			do_anim = 0;
	struct model		*m;
	long			*tab;
	struct rt_vlblock	*vbp;

	if(rt_g.NMG_debug & DEBUG_PLOTEM &&
	   rt_g.NMG_debug & DEBUG_FCUT ) do_plot = 1;
	if( rt_g.NMG_debug & DEBUG_PL_ANIM )  do_anim = 1;

	if( !do_plot && !do_anim )  return;

	m = nmg_find_model( &fu1->l.magic );
	NMG_CK_MODEL(m);
	/* get space for list of items processed */
	tab = (long *)rt_calloc( m->maxindex+1, sizeof(long),
		"nmg_pl_comb_fu tab[]");

	vbp = rt_vlblock_init();

	nmg_vlblock_fu(vbp, fu1, tab, 3);

	if( do_plot )  {
	    	(void)sprintf(name, "comb%d.%d.pl", num1, num2);
		if ((fp=fopen(name, "w")) == (FILE *)NULL) {
			(void)perror(name);
			return;
		}
		rt_log("Plotting %s\n", name);

		rt_plot_vlblock(fp, vbp);

		(void)fclose(fp);
	}

	if( do_anim )  {
		if( nmg_vlblock_anim_upcall )  {
			(*nmg_vlblock_anim_upcall)( vbp,
				(rt_g.NMG_debug&DEBUG_PL_SLOW) ? US_DELAY : 0,
				0 );
		} else {
			rt_log("null nmg_vlblock_anim_upcall, no animation\n");
		}
	}
	rt_vlblock_free(vbp);
	rt_free( (char *)tab, "nmg_pl_comb_fu tab[]" );
}

/*
 *			N M G _ P L _ 2 F U
 *
 *  Note that 'str' is expected to contain a %d to place the frame number.
 *
 *  Called from nmg_isect_2faces and other places.
 */
void
nmg_pl_2fu( str, unused, fu1, fu2, show_mates )
CONST char		*str;
int			unused;
CONST struct faceuse	*fu1;
CONST struct faceuse	*fu2;
int			show_mates;
{
	FILE		*fp;
	char		name[32];
	struct model	*m;
	long		*tab;
	static int	num = 1;
	struct rt_vlblock	*vbp;

	if( rt_g.NMG_debug & (DEBUG_PLOTEM|DEBUG_PL_ANIM) == 0 )  return;

	m = nmg_find_model( &fu1->l.magic );
	NMG_CK_MODEL(m);
	/* get space for list of items processed */
	tab = (long *)rt_calloc( m->maxindex+1, sizeof(long),
		"nmg_pl_comb_fu tab[]");

	/* Create the vlblock */
	vbp = rt_vlblock_init();

	nmg_vlblock_fu( vbp, fu1, tab, 3);
	if( show_mates )
		nmg_vlblock_fu( vbp, fu1->fumate_p, tab, 3);

	nmg_vlblock_fu( vbp, fu2, tab, 3);
	if( show_mates )
		nmg_vlblock_fu( vbp, fu2->fumate_p, tab, 3);

	if( rt_g.NMG_debug & DEBUG_PLOTEM )  {
		(void)sprintf(name, str, num++);
		rt_log("plotting to %s\n", name);
		if ((fp=fopen(name, "w")) == (FILE *)NULL)  {
			perror(name);
			return;
		}
		rt_plot_vlblock( fp, vbp );
		(void)fclose(fp);
	}

	if( rt_g.NMG_debug & DEBUG_PL_ANIM )  {
		/* Cause animation of boolean operation as it proceeds! */
		if( nmg_vlblock_anim_upcall )  {
			(*nmg_vlblock_anim_upcall)( vbp,
				(rt_g.NMG_debug&DEBUG_PL_SLOW) ? US_DELAY : 0,
				0 );
		}
	}

	rt_vlblock_free(vbp);
	rt_free( (char *)tab, "nmg_pl_2fu tab[]" );
}

/************************************************************************
 *									*
 *			Graphical display of classifier results		*
 *									*
 ************************************************************************/

int		nmg_class_nothing_broken=1;
static long	**global_classlist;
static struct rt_vlblock *vbp_old;
static long	*broken_tab;
static int 	broken_color;
static unsigned char broken_colors[][3] = {
	{ 100, 100, 255 },	/* NMG_CLASS_AinB (bright blue) */
	{ 255,  50,  50 },	/* NMG_CLASS_AonBshared (red) */
	{ 255,  50, 255 }, 	/* NMG_CLASS_AonBanti (magenta) */
	{  50, 255,  50 },	/* NMG_CLASS_AoutB (bright green) */
	{ 255, 255, 255 },	/* UNKNOWN (white) */
	{ 255, 255, 125 }	/* no classification list (cyan) */
};
#define PICK_BROKEN_COLOR(p) { \
	if (global_classlist == (long **)NULL) { \
		broken_color = 5; \
	} else if( NMG_INDEX_TEST(global_classlist[NMG_CLASS_AinB], (p)) ) \
		broken_color = NMG_CLASS_AinB; \
	else if( NMG_INDEX_TEST(global_classlist[NMG_CLASS_AonBshared], (p)) ) \
		broken_color = NMG_CLASS_AonBshared; \
	else if( NMG_INDEX_TEST(global_classlist[NMG_CLASS_AonBanti], (p)) ) \
		broken_color = NMG_CLASS_AonBanti; \
	else if ( NMG_INDEX_TEST(global_classlist[NMG_CLASS_AoutB], (p)) ) \
		broken_color = NMG_CLASS_AoutB; \
	else \
		broken_color = 4;}

/*
 *			S H O W _ B R O K E N _ V U
 */
static void
show_broken_vu(vbp, vu, fancy)
struct rt_vlblock *vbp;
int fancy;
CONST struct vertexuse *vu;
{
	pointp_t p;
	static char label[128];
	struct rt_list	*vh;
	struct vertex *v;
	point_t pt;

	NMG_CK_VERTEXUSE(vu);
	v = vu->v_p;
	NMG_CK_VERTEX(v);
	NMG_CK_VERTEX_G(v->vg_p);

	NMG_INDEX_RETURN_IF_SET_ELSE_SET( broken_tab, v->index );

	NMG_CK_VERTEX_G(v->vg_p);
	p = v->vg_p->coord;

	PICK_BROKEN_COLOR(vu->v_p);
	if (broken_color == 4) {
/*		fprintf(stderr, "vertex broken_color %d...", broken_color); */
		PICK_BROKEN_COLOR(vu);
/*		fprintf(stderr, "vertexuse broken_color %d\n", broken_color); */
	}
	vh = rt_vlblock_find( vbp, 
		broken_colors[broken_color][0], broken_colors[broken_color][1], broken_colors[broken_color][2]);

	RT_ADD_VLIST( vh, p, RT_VLIST_LINE_MOVE );
	RT_ADD_VLIST( vh, p, RT_VLIST_LINE_DRAW );


	VMOVE(pt, p);
	pt[0] += 0.05;
	RT_ADD_VLIST( vh, pt, RT_VLIST_LINE_MOVE );
	VMOVE(pt, p);
	pt[0] -= 0.05;
	RT_ADD_VLIST( vh, pt, RT_VLIST_LINE_DRAW );

	VMOVE(pt, p);
	pt[1] += 0.05;
	RT_ADD_VLIST( vh, pt, RT_VLIST_LINE_MOVE );
	VMOVE(pt, p);
	pt[1] -= 0.05;
	RT_ADD_VLIST( vh, pt, RT_VLIST_LINE_DRAW );

	VMOVE(pt, p);
	pt[2] += 0.05;
	RT_ADD_VLIST( vh, pt, RT_VLIST_LINE_MOVE );
	VMOVE(pt, p);
	pt[2] -= 0.05;
	RT_ADD_VLIST( vh, pt, RT_VLIST_LINE_DRAW );

	RT_ADD_VLIST( vh, p, RT_VLIST_LINE_MOVE );
}

static void
show_broken_e(vbp, eu, fancy)
struct rt_vlblock *vbp;
int fancy;
CONST struct edgeuse *eu;
{
	pointp_t p0, p1;
	point_t end0, end1;
	vect_t v;
	struct rt_list	*vh;

	NMG_CK_VERTEXUSE(eu->vu_p);
	NMG_CK_VERTEX(eu->vu_p->v_p);
	NMG_CK_VERTEX_G(eu->vu_p->v_p->vg_p);
	NMG_CK_VERTEXUSE(eu->eumate_p->vu_p);
	NMG_CK_VERTEX(eu->eumate_p->vu_p->v_p);
	NMG_CK_VERTEX_G(eu->eumate_p->vu_p->v_p->vg_p);

	NMG_INDEX_RETURN_IF_SET_ELSE_SET( broken_tab, eu->e_p->index );

	p0 = eu->vu_p->v_p->vg_p->coord;
	p1 = eu->eumate_p->vu_p->v_p->vg_p->coord;

	/* leave a little room between the edge endpoints and the vertex
	 * compute endpoints by forming a vector between verts, scale vector,
	 * and modify points
	 */
	VSUB2SCALE(v, p1, p0, 0.90);
	VADD2(end0, p0, v);
	VSUB2(end1, p1, v);


	PICK_BROKEN_COLOR(eu->e_p);
	if (broken_color == 4) {
/*		fprintf(stderr, "edge broken_color %d... ", broken_color); */
		PICK_BROKEN_COLOR(eu);
/*		fprintf(stderr, "edgeuse broken_color %d\n", broken_color); */
	}

	vh = rt_vlblock_find( vbp, 
		broken_colors[broken_color][0], broken_colors[broken_color][1], broken_colors[broken_color][2]);

	RT_ADD_VLIST( vh, end0, RT_VLIST_LINE_MOVE );
	RT_ADD_VLIST( vh, end1, RT_VLIST_LINE_DRAW );

	show_broken_vu(vbp, eu->vu_p, fancy);
	show_broken_vu(vbp, eu->eumate_p->vu_p, fancy);

}


static void
show_broken_eu(vbp, eu, fancy)
struct rt_vlblock *vbp;
int fancy;
CONST struct edgeuse *eu;
{
	vect_t v;
	struct rt_list	*vh;
    	int red, green, blue;
	point_t base, tip;
	point_t	radial_tip;
	point_t	next_base;
	point_t	last_tip;

	NMG_CK_EDGEUSE(eu);
	NMG_CK_EDGE(eu->e_p);

	show_broken_e(vbp, eu, fancy);

	if (!fancy) return;

	/* paint the edgeuse lines */
	if (*eu->up.magic_p == NMG_LOOPUSE_MAGIC &&
	    *eu->up.lu_p->up.magic_p == NMG_FACEUSE_MAGIC) {

	    	red = broken_colors[broken_color][0];
	    	green = broken_colors[broken_color][1];
	    	blue = broken_colors[broken_color][2];

	    	nmg_eu_coords(eu, base, tip);
	    	if (eu->up.lu_p->up.fu_p->orientation == OT_SAME)
	    		red += 50;
		else if (eu->up.lu_p->up.fu_p->orientation == OT_OPPOSITE)
			red -= 50;
	    	else
	    		red = green = blue = 255;

		vh = rt_vlblock_find( vbp, red, green, blue );
		RT_ADD_VLIST( vh, base, RT_VLIST_LINE_MOVE );
		RT_ADD_VLIST( vh, tip, RT_VLIST_LINE_DRAW );

	    	nmg_eu_radial( eu, radial_tip );
		vh = rt_vlblock_find( vbp, red, green-20, blue );
		RT_ADD_VLIST( vh, tip, RT_VLIST_LINE_MOVE );
		RT_ADD_VLIST( vh, radial_tip, RT_VLIST_LINE_DRAW );

	    	nmg_eu_next_base( eu, next_base );
		vh = rt_vlblock_find( vbp, 0, 100, 0 );
#if 0
		RT_ADD_VLIST( vh, tip, RT_VLIST_LINE_MOVE );
		RT_ADD_VLIST( vh, next_base, RT_VLIST_LINE_DRAW );
#else
	    	{
	    		register struct rt_vlist *_vp = RT_LIST_LAST( rt_vlist, (vh) );
	    		if( RT_LIST_IS_HEAD( _vp, (vh) ) || _vp->nused >= RT_VLIST_CHUNK ) {
#if 0
	    			RT_GET_VLIST(_vp);
#else
	    			_vp = RT_LIST_FIRST( rt_vlist, &rt_g.rtg_vlfree );
	    			if( RT_LIST_IS_HEAD( _vp, &rt_g.rtg_vlfree ) ) {
					_vp = (struct rt_vlist *)rt_malloc(sizeof(struct rt_vlist), "rt_vlist");
					_vp->l.magic = RT_VLIST_MAGIC;
				} else {
					RT_LIST_DEQUEUE( &(_vp->l) );
				}
				_vp->nused = 0;
#endif
	    			RT_LIST_INSERT( (vh), &(_vp->l) );
	    		}
	    		VMOVE( _vp->pt[_vp->nused], (tip) );
	    		_vp->cmd[_vp->nused++] = (RT_VLIST_LINE_MOVE);
	    	}
	    	{
	    		register struct rt_vlist *_vp = RT_LIST_LAST( rt_vlist, (vh) );
	    		if( RT_LIST_IS_HEAD( _vp, (vh) ) || _vp->nused >= RT_VLIST_CHUNK ) {
	    			RT_GET_VLIST(_vp);
	    			RT_LIST_INSERT( (vh), &(_vp->l) );
	    		}
	    		VMOVE( _vp->pt[_vp->nused], (next_base) );
	    		_vp->cmd[_vp->nused++] = (RT_VLIST_LINE_DRAW);
	    	}
#endif
	}

}

static void
show_broken_lu(vbp, lu, fancy)
struct rt_vlblock *vbp;
int fancy;
CONST struct loopuse *lu;
{
	register struct edgeuse *eu;
	struct rt_list	*vh;
	vect_t		n;

	NMG_CK_LOOPUSE(lu);

	if( RT_LIST_FIRST_MAGIC(&lu->down_hd)==NMG_VERTEXUSE_MAGIC )  {
		register struct vertexuse *vu;
		vu = RT_LIST_FIRST(vertexuse, &lu->down_hd);
		show_broken_vu(vbp, vu, fancy);
		return;
	}

	if (rt_g.NMG_debug & DEBUG_GRAPHCL)  {
		for (RT_LIST_FOR(eu, edgeuse, &lu->down_hd))
			show_broken_eu(vbp, eu, fancy);
	}

	/* Draw colored polygons for the actual face loops */
	/* Faces are not classified, only loops */
	/* This can obscure the edge/vertex info */
	PICK_BROKEN_COLOR(lu->l_p);
	vh = rt_vlblock_find( vbp, 
		broken_colors[broken_color][0], broken_colors[broken_color][1], broken_colors[broken_color][2]);

	if( *lu->up.magic_p == NMG_FACEUSE_MAGIC )  {
		NMG_GET_FU_NORMAL( n, lu->up.fu_p );
	} else {
		/* For wire loops, use a constant normal */
		VSET( n, 0, 0, 1 );
	}

	if ((rt_g.NMG_debug & (DEBUG_GRAPHCL|DEBUG_PL_LOOP)) == (DEBUG_PL_LOOP) ) {
		/* If only DEBUG_PL_LOOP set, just draw lu as wires */
		nmg_lu_to_vlist( vh, lu, 0, n );
	} else if ((rt_g.NMG_debug & (DEBUG_GRAPHCL|DEBUG_PL_LOOP)) == (DEBUG_GRAPHCL|DEBUG_PL_LOOP) ) {
		/* Draw as polygons if both set */
		nmg_lu_to_vlist( vh, lu, 1, n );
	} else {
		/* If only DEBUG_GRAPHCL set, don't draw lu's at all */
	}
}



static void
show_broken_fu(vbp, fu, fancy)
struct rt_vlblock *vbp;
int fancy;
CONST struct faceuse *fu;
{
	register struct loopuse *lu;

	NMG_CK_FACEUSE(fu);
	for (RT_LIST_FOR(lu, loopuse, &fu->lu_hd)) {
		show_broken_lu(vbp, lu, fancy);
	}
}

static void
show_broken_s(vbp, s, fancy)
struct rt_vlblock *vbp;
CONST struct shell *s;
int fancy;
{
	struct faceuse *fu;
	struct loopuse *lu;
	struct edgeuse *eu;
	struct vertexuse *vu;

	NMG_CK_SHELL(s);
	for ( RT_LIST_FOR(fu, faceuse, &s->fu_hd ))
		show_broken_fu(vbp, fu, fancy);
	for ( RT_LIST_FOR(lu, loopuse, &s->lu_hd ))
		show_broken_lu(vbp, lu, fancy);
	for ( RT_LIST_FOR(eu, edgeuse, &s->eu_hd ))
		show_broken_eu(vbp, eu, fancy);
	if ( s->vu_p )
		show_broken_vu(vbp, s->vu_p, fancy);
}
static void
show_broken_r(vbp, r, fancy)
struct rt_vlblock *vbp;
CONST struct nmgregion *r;
int fancy;
{
	register struct shell *s;

	NMG_CK_REGION(r);
	for ( RT_LIST_FOR(s, shell, & r->s_hd))
		show_broken_s(vbp, s, fancy);
}

static void
show_broken_m(vbp, m, fancy)
struct rt_vlblock *vbp;
CONST struct model *m;
int fancy;
{
	register struct nmgregion *r;

	NMG_CK_MODEL(m);
	for (RT_LIST_FOR(r, nmgregion, &m->r_hd))
		show_broken_r(vbp, r, fancy);
}

static struct rt_vlblock *vbp = (struct rt_vlblock *)NULL;
static int stepalong = 0;

void
sigstepalong(i)
int i;
{
	stepalong=1;
}

/*
 *			S H O W _ B R O K E N _ S T U F F
 *
 * XXX Needs new name, with nmg_ prefix, and a stronger indication
 * that this is a graphical display of classifier operation.
 */
void
show_broken_stuff(p, classlist, all_new, fancy, a_string)
long	*classlist[4];
long	*p;
int	all_new;
CONST char	*a_string;
{
	struct model *m;

/*	printf("showing broken stuff\n"); */

	global_classlist = classlist;

	nmg_class_nothing_broken = 0;

	if (!vbp)
		vbp = rt_vlblock_init();
	else if (all_new) {
		rt_vlblock_free(vbp);
		vbp = (struct rt_vlblock *)NULL;
		vbp = rt_vlblock_init();
	}

	m = nmg_find_model(p);
	/* get space for list of items processed */
	if (!broken_tab) {
		broken_tab = (long *)rt_calloc( m->maxindex+1, sizeof(long),
			"nmg_vlblock_s tab[]");
	} else if (all_new) {
		bzero(broken_tab,  (m->maxindex+1) * sizeof(long));
	}


	switch (*p) {
	case NMG_MODEL_MAGIC:
		show_broken_m( vbp, (struct model *)p, fancy);
		break;
	case NMG_REGION_MAGIC:
		show_broken_r( vbp, (struct nmgregion *)p, fancy);
		break;
	case NMG_SHELL_MAGIC:
		show_broken_s( vbp, (struct shell *)p, fancy);
		break;
	case NMG_FACE_MAGIC:
		show_broken_fu( vbp, ((struct face *)p)->fu_p, fancy);
		break;
	case NMG_FACEUSE_MAGIC:
		show_broken_fu( vbp, (struct faceuse *)p, fancy);
#if 0
		{
			struct rt_vlblock *vbp2 = vbp;
			register struct loopuse *lu;
			struct faceuse *fu = (struct faceuse *)p;
			int i;
			void            (*cur_sigint)();

			cur_sigint = signal(SIGINT, sigstepalong);
			for (stepalong=0;!stepalong;) {
				for (RT_LIST_FOR(lu, loopuse, &fu->lu_hd)) {
					show_broken_stuff(lu, classlist, 1, fancy);
					for (i=0 ; ++i ; );
				}
			}
			signal(SIGINT, cur_sigint);

			show_broken_fu( vbp, (struct faceuse *)p, fancy);
		}
#endif
		break;
	case NMG_LOOPUSE_MAGIC:
		show_broken_lu( vbp, (struct loopuse *)p, fancy);
		break;
	case NMG_EDGE_MAGIC:
		 show_broken_eu( vbp, ((struct edge *)p)->eu_p, fancy);
		 break;
	case NMG_EDGEUSE_MAGIC:
		show_broken_eu( vbp, (struct edgeuse *)p, fancy);
		break;
	case NMG_VERTEXUSE_MAGIC:
		show_broken_vu( vbp, (struct vertexuse *)p, fancy);
		break;
	default: fprintf(stderr, "Unknown magic number %ld %0x %ld %0x\n", *p, *p, p, p);
				break;
	}

	/* Cause animation of boolean operation as it proceeds! */
	/* The "copy" flag on nmg_vlblock_anim_upcall() means that
	 * the vlist will remain, undisturbed, for further use. */
	if( nmg_vlblock_anim_upcall )  {
		struct rt_vlblock *vbp2 = vbp;
		register struct loopuse *lu;
		struct faceuse *fu = (struct faceuse *)p;
		int i;
		void            (*cur_sigint)();

		if (!a_string) {
			(*nmg_vlblock_anim_upcall)( vbp,
				(rt_g.NMG_debug&DEBUG_PL_SLOW) ? US_DELAY : 0,
				1 );
		} else {

			fprintf(stderr, "NMG Intermediate display Ctrl-C to continue (%s)\n", a_string);
			cur_sigint = signal(SIGINT, sigstepalong);
			(*nmg_vlblock_anim_upcall)( vbp,
				(rt_g.NMG_debug&DEBUG_PL_SLOW) ? US_DELAY : 0,
				1 );
			for (stepalong = 0; !stepalong ; ) {
				(*nmg_mged_debug_display_hack)();
			}
			signal(SIGINT, cur_sigint);
			fprintf(stderr, "Continuing\n");
		}
	} else {
		rt_vlblock_free(vbp);
		vbp = (struct rt_vlblock *)NULL;
		rt_free((char *)broken_tab, "broken_tab");
		broken_tab = (long *)NULL;
	}
}

/*
 *			N M G _ F A C E _ P L O T
 */
void
nmg_face_plot( fu )
CONST struct faceuse	*fu;
{
	FILE		*fp;
	char		name[32];
	extern void (*nmg_vlblock_anim_upcall)();
	struct model		*m;
	struct rt_vlblock	*vbp;
	struct face_g_plane	*fg;
	long		*tab;
	int		fancy;
	static int	num = 1;

	if( rt_g.NMG_debug & (DEBUG_PLOTEM|DEBUG_PL_ANIM) == 0 )  return;

	NMG_CK_FACEUSE(fu);

	m = nmg_find_model( (long *)fu );
	NMG_CK_MODEL(m);

	/* get space for list of items processed */
	tab = (long *)rt_calloc( m->maxindex+1, sizeof(long),
		"nmg_face_plot tab[]");

	vbp = rt_vlblock_init();

	fancy = 3;	/* show both types of edgeuses */
	nmg_vlblock_fu(vbp, fu, tab, fancy );

	if( rt_g.NMG_debug & DEBUG_PLOTEM )  {
		(void)sprintf(name, "face%d.pl", num++);
		rt_log("plotting to %s\n", name);
		if ((fp=fopen(name, "w")) == (FILE *)NULL)  {
			perror(name);
			return;
		}
		rt_plot_vlblock( fp, vbp );
		(void)fclose(fp);
	}

	if( rt_g.NMG_debug & DEBUG_PL_ANIM )  {
		/* Cause animation of boolean operation as it proceeds! */
		if( nmg_vlblock_anim_upcall )  {
			/* if requested, delay 3/4 second */
			(*nmg_vlblock_anim_upcall)( vbp,
				(rt_g.NMG_debug&DEBUG_PL_SLOW) ? 750000 : 0,
				0 );
		} else {
			rt_log("null nmg_vlblock_anim_upcall, no animation\n");
		}
	}
	rt_vlblock_free(vbp);
	rt_free( (char *)tab, "nmg_face_plot tab[]" );

}

/*
 *			N M G _ 2 F A C E _ P L O T
 *
 *  Just like nmg_face_plot, except it draws two faces each iteration.
 */
void
nmg_2face_plot( fu1, fu2 )
CONST struct faceuse	*fu1, *fu2;
{
	extern void (*nmg_vlblock_anim_upcall)();
	struct model		*m;
	struct rt_vlblock	*vbp;
	struct face_g_plane	*fg;
	long		*tab;
	int		fancy;

	if( ! (rt_g.NMG_debug & DEBUG_PL_ANIM) )  return;

	NMG_CK_FACEUSE(fu1);
	NMG_CK_FACEUSE(fu2);

	m = nmg_find_model( (long *)fu1 );
	NMG_CK_MODEL(m);

	/* get space for list of items processed */
	tab = (long *)rt_calloc( m->maxindex+1, sizeof(long),
		"nmg_2face_plot tab[]");

	vbp = rt_vlblock_init();

	fancy = 3;	/* show both types of edgeuses */
	nmg_vlblock_fu(vbp, fu1, tab, fancy );
	nmg_vlblock_fu(vbp, fu2, tab, fancy );

	/* Cause animation of boolean operation as it proceeds! */
	if( nmg_vlblock_anim_upcall )  {
		/* if requested, delay 3/4 second */
		(*nmg_vlblock_anim_upcall)( vbp,
			(rt_g.NMG_debug&DEBUG_PL_SLOW) ? 750000 : 0,
			0 );
	} else {
		rt_log("null nmg_vlblock_anim_upcall, no animation\n");
	}
	rt_vlblock_free(vbp);
	rt_free( (char *)tab, "nmg_2face_plot tab[]" );

}

/*
 *			N M G _ F A C E _ L U _ P L O T
 *
 *  Plot the loop, and a ray from vu1 to vu2.
 */
void
nmg_face_lu_plot( lu, vu1, vu2 )
CONST struct loopuse		*lu;
CONST struct vertexuse		*vu1, *vu2;
{
	FILE	*fp;
	struct model	*m;
	long		*b;
	char		buf[128];
	static int	num = 0;
	vect_t		dir;
	point_t		p1, p2;

	if(!(rt_g.NMG_debug&DEBUG_PLOTEM)) return;

	NMG_CK_LOOPUSE(lu);
	NMG_CK_VERTEXUSE(vu1);
	NMG_CK_VERTEXUSE(vu2);

	m = nmg_find_model((long *)lu);
	sprintf(buf, "loop%d.pl", num++ );

	if( (fp = fopen(buf, "w")) == NULL )  {
		perror(buf);
		return;
	}
	b = (long *)rt_calloc( m->maxindex, sizeof(long), "nmg_face_lu_plot flag[]" );
	nmg_pl_lu(fp, lu, b, 255, 0, 0);

	/*
	 *  Two yellow lines for the ray.
	 *  Overshoot edge by +/-10%, for visibility.
	 *  Don't draw over top of the actual edge, it might hide verts.
	 */
	pl_color(fp, 255, 255, 0);
	VSUB2( dir, vu2->v_p->vg_p->coord, vu1->v_p->vg_p->coord );
	VJOIN1( p1, vu1->v_p->vg_p->coord, -0.1, dir );
	pdv_3line(fp, p1, vu1->v_p->vg_p->coord );
	VJOIN1( p2, vu1->v_p->vg_p->coord,  1.1, dir );
	pdv_3line(fp, vu2->v_p->vg_p->coord, p2 );

	fclose(fp);
	rt_log("wrote %s\n", buf);
	rt_free( (char *)b, "nmg_face_lu_plot flag[]" );
}

/*
 *			N M G _ P L O T _ L U _ R A Y
 *
 *  Plot the loop, a ray from vu1 to vu2, and the left vector.
 */
void
nmg_plot_lu_ray( lu, vu1, vu2, left )
CONST struct loopuse		*lu;
CONST struct vertexuse		*vu1, *vu2;
CONST vect_t			left;
{
	FILE	*fp;
	struct model	*m;
	long		*b;
	char		buf[128];
	static int	num = 0;
	vect_t		dir;
	point_t		p1, p2;
	fastf_t		left_mag;

	if(!(rt_g.NMG_debug&DEBUG_PLOTEM)) return;

	NMG_CK_LOOPUSE(lu);
	NMG_CK_VERTEXUSE(vu1);
	NMG_CK_VERTEXUSE(vu2);

	m = nmg_find_model((long *)lu);
	sprintf(buf, "loop%d.pl", num++ );

	if( (fp = fopen(buf, "w")) == NULL )  {
		perror(buf);
		return;
	}
	b = (long *)rt_calloc( m->maxindex, sizeof(long), "nmg_plot_lu_ray flag[]" );
	nmg_pl_lu(fp, lu, b, 255, 0, 0);

	/*
	 *  Two yellow lines for the ray, and a third for the left vector.
	 *  Overshoot edge by +/-10%, for visibility.
	 *  Don't draw over top of the actual edge, it might hide verts.
	 */
	pl_color(fp, 255, 255, 0);
	VSUB2( dir, vu2->v_p->vg_p->coord, vu1->v_p->vg_p->coord );
	VJOIN1( p1, vu1->v_p->vg_p->coord, -0.1, dir );
	pdv_3line(fp, p1, vu1->v_p->vg_p->coord );
	VJOIN1( p2, vu1->v_p->vg_p->coord,  1.1, dir );
	pdv_3line(fp, vu2->v_p->vg_p->coord, p2 );

	/* The left vector */
	left_mag = 0.1 * MAGNITUDE(dir);
	VJOIN1( p2, p1, left_mag, left );
	pdv_3line(fp, p1, p2);

	fclose(fp);
	rt_log("wrote %s\n", buf);
	rt_free( (char *)b, "nmg_plot_lu_ray flag[]" );
}

/*
 *			N M G _ P L O T _ R A Y _ F A C E
 */
void
nmg_plot_ray_face(fname, pt, dir, fu)
CONST char *fname;
point_t pt;
CONST vect_t dir;
CONST struct faceuse *fu;
{
	FILE *fd;
	long *b;
	point_t pp;
	static int i=0;
	char name[1024];

	if ( ! (rt_g.NMG_debug & DEBUG_NMGRT) )
		return;

	sprintf(name, "%s%0d.pl", fname, i++);
	if ((fd = fopen(name, "w")) == (FILE *)NULL) {
		perror(name);
		rt_log("plot_ray_face cannot open %s", name);
		rt_bomb("aborting");
	}

	b = (long *)rt_calloc( fu->s_p->r_p->m_p->maxindex, sizeof(long), "bit vec");

	nmg_pl_fu(fd, fu, b, 200, 200, 200);

	rt_free((char *)b, "bit vec");

	VSCALE(pp, dir, 1000.0);
	VADD2(pp, pt, pp);
	pdv_3line( fd, pt, pp );
	(void)fclose(fd);
}

/*
 *			N M G _ P L O T _ L U _ A R O U N D _ E U
 *
 *  Draw and label all the loopuses gathered around this edgeuse.
 *
 *  Called by nmg_radial_join_eu().
 */
void
nmg_plot_lu_around_eu( prefix, eu, tol )
CONST char		*prefix;
CONST struct edgeuse	*eu;
CONST struct rt_tol	*tol;
{
	char			file[256];
	static int		num=0;
	struct model		*m;
	struct rt_vlblock	*vbp;
	long			*tab;
	CONST struct edgeuse	*eur;
	FILE			*fp;

	NMG_CK_EDGEUSE(eu);
	RT_CK_TOL(tol);

	sprintf(file, "%s%0d.pl", prefix, num++);
	rt_log("plotting to %s\n", file);
	if ((fp = fopen(file, "w")) == (FILE *)NULL) {
		rt_log("plot_lu_around_eu() cannot open %s", file);
		return;
	}

	m = nmg_find_model( (long *)eu );
	NMG_CK_MODEL(m);
	tab = (long *)rt_calloc( m->maxindex, sizeof(long), "bit vec");

	vbp = rt_vlblock_init();

	/* Draw all the left vectors, and a fancy edgeuse plot */
	nmg_vlblock_around_eu(vbp, eu, tab, 3, tol );

	eur = eu;
	do {
		NMG_CK_EDGEUSE(eur);

		if (*eur->up.magic_p == NMG_LOOPUSE_MAGIC )  {
			/* Draw this loop in non-fancy format, for context */
			nmg_vlblock_lu(vbp, eur->up.lu_p, tab, 80, 100, 170, 0, 0 );
		}
		eur = eur->radial_p->eumate_p;
	} while (eur != eu);

	rt_plot_vlblock( fp, vbp );
	(void)fclose(fp);
	rt_vlblock_free(vbp);
	rt_free((char *)tab, "bit vec");
}

/*
 *			N M G _ H A C K _ S N U R B
 *
 *  Convert a new NMG format snurb to the older LIBNURB format,
 *  by copying data and pointers.
 *  Under no circumstances should the output of this routine be freed,
 *  or it will corrupt the NMG original!  Just discard it.
 *
 *  XXX Temporary hack until LIBNURB is updated to new data structures.
 */
void
nmg_hack_snurb( old, fg )
struct snurb	*old;
CONST struct face_g_snurb	*fg;
{
	bzero( (char *)old, sizeof(struct snurb) );

	RT_LIST_INIT( &old->l );
	old->l.magic = RT_SNURB_MAGIC;

	old->order[0] = fg->order[0];
	old->order[1] = fg->order[1];
	old->u_knots = fg->u;		/* struct copy, including pointers! */
	old->v_knots = fg->v;
	old->s_size[0] = fg->s_size[0];
	old->s_size[1] = fg->s_size[1];
	old->pt_type = fg->pt_type;
	old->ctl_points = fg->ctl_points;
}

/*
 *			N M G _ S N U R B _ T O _ V L I S T
 *
 *  A routine to draw the entire surface of a face_g_snurb.
 *  No handling of trimming curves is done.
 */
int
nmg_snurb_to_vlist( vhead, fg, n_interior )
struct rt_list			*vhead;
CONST struct face_g_snurb	*fg;
int				n_interior;	/* typ. 10 */
{
	register int		i;
	register int		j;
	register fastf_t	* vp;
	int			s;
	struct knot_vector 	tkv1,
				tkv2,
				tau1,
				tau2;
	struct snurb		n;	/* XXX hack, don't free! */
	struct snurb 	*r, *c;
	int 		coords;

	RT_CK_LIST( vhead );
	NMG_CK_FACE_G_SNURB(fg);

	rt_nurb_kvgen( &tkv1,
		fg->u.knots[0],
		fg->u.knots[fg->u.k_size-1], n_interior);

	rt_nurb_kvgen( &tkv2,
		fg->v.knots[0],
		fg->v.knots[fg->v.k_size-1], n_interior);
		
	rt_nurb_kvmerge(&tau1, &tkv1, &fg->u);
	rt_nurb_kvmerge(&tau2, &tkv2, &fg->v);

	nmg_hack_snurb( &n, fg );	/* XXX */
	NMG_CK_SNURB(&n);
	r = (struct snurb *) rt_nurb_s_refine( &n, RT_NURB_SPLIT_COL, &tau2);
	NMG_CK_SNURB(r);
	c = (struct snurb *) rt_nurb_s_refine( r, RT_NURB_SPLIT_ROW, &tau1);
	NMG_CK_SNURB(c);

	coords = RT_NURB_EXTRACT_COORDS(c->pt_type);
	
	if( RT_NURB_IS_PT_RATIONAL(c->pt_type))
	{
		vp = c->ctl_points;
		for(i= 0; 
			i < c->s_size[0] * c->s_size[1]; 
			i++)
		{
			FAST fastf_t	div;
			vp[0] *= (div = 1/vp[3]);
			vp[1] *= div;
			vp[2] *= div;
			vp[3] *= div;
			vp += coords;
		}
	}

	vp = c->ctl_points;
	for( i = 0; i < c->s_size[0]; i++)
	{
		RT_ADD_VLIST( vhead, vp, RT_VLIST_LINE_MOVE );
		vp += coords;
		for( j = 1; j < c->s_size[1]; j++)
		{
			RT_ADD_VLIST( vhead, vp, RT_VLIST_LINE_DRAW );
			vp += coords;
		}
	}

	for( j = 0; j < c->s_size[1]; j++)
	{
		int stride;
			
		stride = c->s_size[1] * coords;
		vp = &c->ctl_points[j * coords];
		RT_ADD_VLIST( vhead, vp, RT_VLIST_LINE_MOVE );
		for( i = 0; i < c->s_size[0]; i++)
		{
			RT_ADD_VLIST( vhead, vp, RT_VLIST_LINE_DRAW );
			vp += stride;
		}
	}
	rt_nurb_free_snurb(c);
	rt_nurb_free_snurb(r);

	rt_free( (char *) tau1.knots, "rt_nurb_plot:tau1.knots");
	rt_free( (char *) tau2.knots, "rt_nurb_plot:tau2.knots");
	rt_free( (char *) tkv1.knots, "rt_nurb_plot:tkv1>knots");
	rt_free( (char *) tkv2.knots, "rt_nurb_plot:tkv2.knots");

	return(0);
}
