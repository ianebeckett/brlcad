/*
 *			P R E D I C T O R . C
 *
 *  Put a predictor frame into view, as an aid to velocity-based
 *  navigation through an MGED model.
 *
 *  Inspired by the paper "Manipulating the Future:  Predictor Based
 *  Feedback for Velocity Control in Virtual Environment Navigation"
 *  by Dale Chapman and Colin Ware, <cware@unb.ca>, in
 *  ACM SIGGRAPH Computer Graphics Special Issue on 1992 Symposium
 *  on Interactive 3D Graphics.
 *
 *  Author -
 *	Michael John Muuss
 *  
 *  Source -
 *	SECAD/VLD Computing Consortium
 *	The U. S. Army Ballistic Research Laboratory
 *	Aberdeen Proving Ground, Maryland  21005-5066
 *  
 *  Copyright Notice -
 *	This software is Copyright (C) 1992 by the United States Army.
 *	All rights reserved.
 */
#ifndef lint
static char RCSid[] = "@(#)$Header$ (BRL)";
#endif

#include <stdio.h>
#ifdef BSD
#include <strings.h>
#else
#include <string.h>
#endif

#include "machine.h"
#include "externs.h"
#include "vmath.h"
#include "rtstring.h"
#include "raytrace.h"
#include "./ged.h"

extern double	frametime;		/* time needed to draw last frame */
extern mat_t	ModelDelta;		/* Changed to Viewrot this frame */

#define MAX_TRAIL	32
struct trail {
	int	cur_index;	/* index of first free entry */
	int	nused;		/* max index in use */
	point_t	pt[MAX_TRAIL];
};

/*
 *			I N I T _ T R A I L
 */
static void
init_trail(tp)
struct trail	*tp;
{
	tp->cur_index = 0;
	tp->nused = 0;
}

/*
 *			P U S H _ T R A I L
 *
 *  Add a new point to the end of the trail.
 */
static void
push_trail(tp, pt)
struct trail	*tp;
point_t		pt;
{
	VMOVE( tp->pt[tp->cur_index], pt );
	if( tp->cur_index >= tp->nused )  tp->nused++;
	tp->cur_index++;
	if( tp->cur_index >= MAX_TRAIL )  tp->cur_index = 0;
}

/*
 *			D R A W _ T R A I L
 *
 *  Draw from the most recently added point, backwards.
 */
static void
draw_trail(vhead, tp)
struct rt_list	*vhead;
struct trail	*tp;
{
	int	i;
	int	todo = tp->nused;

	RT_LIST_INIT( vhead );
	if( tp->nused <= 0 )  return;
	if( (i = tp->cur_index-1) < 0 )  i = tp->nused-1;
	for( ; todo > 0; todo-- )  {
		if( todo == tp->nused )  {
			RT_ADD_VLIST( vhead, tp->pt[i], RT_VLIST_LINE_MOVE );
		}  else  {
			RT_ADD_VLIST( vhead, tp->pt[i], RT_VLIST_LINE_DRAW );
		}
		if( (--i) < 0 )  i = tp->nused-1;
	}
}

static struct trail	tA, tB, tC, tD;
static struct trail	tE, tF, tG, tH;

#define PREDICTOR_NAME	"_PREDIC_FRAME_"

/*
 *			P R E D I C T O R _ K I L L
 */
void
predictor_kill()
{
	struct rt_vls	str;

	rt_vls_init( &str );
	rt_vls_printf( &str, "kill %s\n", PREDICTOR_NAME );
	(void)cmdline( &str );
	rt_vls_trunc( &str, 0 );

	rt_vls_strcat( &str, "kill _PREDIC_TRAIL_*\n" );
	(void)cmdline( &str );
	rt_vls_free( &str );

	init_trail( &tA );
	init_trail( &tB );
	init_trail( &tC );
	init_trail( &tD );

	init_trail( &tE );
	init_trail( &tF );
	init_trail( &tG );
	init_trail( &tH );
}

#define TF_BORD	0.01
#define TF_X	0.14
#define TF_Y	0.07
#define TF_Z	(1.0-0.15)	/* To prevent Z clipping of TF_X */

#define TF_VL( _m, _v ) \
	{ vect_t edgevect_m; \
	MAT4X3VEC( edgevect_m, predictorXv2m, _v ); \
	VADD2( _m, framecenter_m, edgevect_m ); }

/*
 *			P R E D I C T O R _ F R A M E
 *
 *  Draw the frame itself as four polygons:
 *	ABFE, HGCD, EILH, and JFGK.
 *  The streamers will attach at edges AE, BF, GC, and HD.
 *	
 *		D --------------- C
 *		|                 |
 *		H -L-----------K- G
 *		|  |           |  |
 *		|  |           |  |
 *		|  |           |  |
 *		E -I-----------J- F
 *		|                 |
 *		A --------------- B
 */
void
predictor_frame()
{
	int	i;
	int	nframes;
	mat_t	predictor;
	mat_t	predictorXv2m;
	point_t	v;		/* view coords */
	point_t	m;		/* model coords */
	point_t	mA,mB,mC,mD,mE,mF,mG,mH,mI,mJ,mK,mL;
	struct rt_list	vhead;
	struct rt_list	trail;
	point_t	framecenter_m;
	point_t	framecenter_v;
	point_t	center_m;
	vect_t	delta_v;

	if( rateflag_rotate == 0 && rateflag_slew == 0 && rateflag_zoom == 0 )  {
		/* If no motion, and predictor is drawn, get rid of it */
		if( db_lookup( dbip, PREDICTOR_NAME, LOOKUP_QUIET ) != DIR_NULL )  {
			predictor_kill();
			dmaflag = 1;
		}
		return;
	}

	/* Advance into the future */
	nframes = (int)(mged_variables.predictor_advance / frametime);
	if( nframes < 1 )  nframes = 1;

	/* Build view2model matrix for the future time */
	mat_idn( predictor );
	for( i=0; i < nframes; i++ )  {
		mat_mul2( ModelDelta, predictor );
	}
	mat_mul( predictorXv2m, predictor, view2model );

	MAT_DELTAS_GET_NEG( center_m, toViewcenter );
	MAT4X3PNT( framecenter_m, predictor, center_m );
	MAT4X3PNT( framecenter_v, model2view, framecenter_m );

	/*
	 * Draw the frame around the point framecenter_v.
	 */
	RT_LIST_INIT( &vhead );

	/* Centering dot */
	VSETALL( delta_v, 0 );
	TF_VL( m, delta_v );
	RT_ADD_VLIST( &vhead, m, RT_VLIST_LINE_MOVE );
	RT_ADD_VLIST( &vhead, m, RT_VLIST_LINE_DRAW );

	/* The exterior rectangle */
	VSET( delta_v, -TF_X, -TF_Y, 0 );
	TF_VL( mA, delta_v );

	VSET( delta_v,  TF_X, -TF_Y, 0 );
	TF_VL( mB, delta_v );

	VSET( delta_v,  TF_X,  TF_Y, 0 );
	TF_VL( mC, delta_v );

	VSET( delta_v, -TF_X,  TF_Y, 0 );
	TF_VL( mD, delta_v );

	/* The EFGH rectangle */
	VSET( delta_v, -TF_X, -TF_Y+TF_BORD, 0 );
	TF_VL( mE, delta_v );

	VSET( delta_v,  TF_X, -TF_Y+TF_BORD, 0 );
	TF_VL( mF, delta_v );

	VSET( delta_v,  TF_X,  TF_Y-TF_BORD, 0 );
	TF_VL( mG, delta_v );

	VSET( delta_v, -TF_X,  TF_Y-TF_BORD, 0 );
	TF_VL( mH, delta_v );

	/* The IJKL rectangle */
	VSET( delta_v, -TF_X+TF_BORD, -TF_Y+TF_BORD, 0 );
	TF_VL( mI, delta_v );

	VSET( delta_v,  TF_X-TF_BORD, -TF_Y+TF_BORD, 0 );
	TF_VL( mJ, delta_v );

	VSET( delta_v,  TF_X-TF_BORD,  TF_Y-TF_BORD, 0 );
	TF_VL( mK, delta_v );

	VSET( delta_v, -TF_X+TF_BORD,  TF_Y-TF_BORD, 0 );
	TF_VL( mL, delta_v );

	RT_ADD_VLIST( &vhead, mA, RT_VLIST_LINE_MOVE );
	RT_ADD_VLIST( &vhead, mB, RT_VLIST_LINE_DRAW );
	RT_ADD_VLIST( &vhead, mF, RT_VLIST_LINE_DRAW );
	RT_ADD_VLIST( &vhead, mE, RT_VLIST_LINE_DRAW );
	RT_ADD_VLIST( &vhead, mA, RT_VLIST_LINE_DRAW );

	RT_ADD_VLIST( &vhead, mE, RT_VLIST_LINE_MOVE );
	RT_ADD_VLIST( &vhead, mI, RT_VLIST_LINE_DRAW );
	RT_ADD_VLIST( &vhead, mL, RT_VLIST_LINE_DRAW );
	RT_ADD_VLIST( &vhead, mH, RT_VLIST_LINE_DRAW );
	RT_ADD_VLIST( &vhead, mE, RT_VLIST_LINE_DRAW );

	RT_ADD_VLIST( &vhead, mH, RT_VLIST_LINE_MOVE );
	RT_ADD_VLIST( &vhead, mG, RT_VLIST_LINE_DRAW );
	RT_ADD_VLIST( &vhead, mC, RT_VLIST_LINE_DRAW );
	RT_ADD_VLIST( &vhead, mD, RT_VLIST_LINE_DRAW );
	RT_ADD_VLIST( &vhead, mH, RT_VLIST_LINE_DRAW );

	RT_ADD_VLIST( &vhead, mJ, RT_VLIST_LINE_MOVE );
	RT_ADD_VLIST( &vhead, mF, RT_VLIST_LINE_DRAW );
	RT_ADD_VLIST( &vhead, mG, RT_VLIST_LINE_DRAW );
	RT_ADD_VLIST( &vhead, mK, RT_VLIST_LINE_DRAW );
	RT_ADD_VLIST( &vhead, mJ, RT_VLIST_LINE_DRAW );

	invent_solid( PREDICTOR_NAME, &vhead, 0x00FFFFFFL );

	push_trail( &tA, mA );
	push_trail( &tB, mB );
	push_trail( &tC, mC );
	push_trail( &tD, mD );

	push_trail( &tE, mE );
	push_trail( &tF, mF );
	push_trail( &tG, mG );
	push_trail( &tH, mH );

	/* Draw the trails */

	draw_trail( &trail, &tA );
	invent_solid( "_PREDIC_TRAIL_LL_", &trail, 0x00FF00FFL );

	draw_trail( &trail, &tB );
	invent_solid( "_PREDIC_TRAIL_LR_", &trail, 0x0000FFFFL );

	draw_trail( &trail, &tC );
	invent_solid( "_PREDIC_TRAIL_UR_", &trail, 0x00FF00FFL );

	draw_trail( &trail, &tD );
	invent_solid( "_PREDIC_TRAIL_UL_", &trail, 0x0000FFFFL );

	/* Done */
	mat_idn( ModelDelta );
}

/*
 *			P R E D I C T O R _ H O O K
 *
 *  Called from set.c when the predictor variables are modified.
 */
void
predictor_hook()
{
	if( mged_variables.predictor > 0 )  {
		/* Allocate storage? */
	} else {
		/* Release storage? */
		predictor_kill();
	}
	dmaflag = 1;
}
