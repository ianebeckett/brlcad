/*                           N M G . H
 * BRL-CAD
 *
 * Copyright (c) 2004-2022 United States Government as represented by
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
/** @addtogroup libnmg
 *
 * Definition of data structures for "Non-Manifold Geometry
 * Modelling."  Developed from "Non-Manifold Geometric Boundary
 * Modeling" by Kevin Weiler, 5/7/87 (SIGGraph 1989 Course #20 Notes)
 *
 * See also "Topological Structures for Geometric Modeling"
 * by Kevin J. Weiler - RPI Phd thesis from 1986.
 *
 */
/** @{ */
/** @file include/nmg.h */

/*
 * TODO: since this original work, there have been a number of subsequent
 * papers published on ways of representing non-manifold geometry.  In
 * particular, the "partial entity" structure is purported to have essentially
 * the same representation power as Weiler's radial edge approach but with
 * more compact storage:
 *
 *  Sang Hun Lee and Kunwoo Lee, Partial Entity Structure: A Compact Boundary
 *  Representation for Non-Manifold Geometric Modeling, J. Comput. Inf. Sci.
 *  Eng 1(4), 356-365 (Nov 01, 2001)
 *
 *  Sang Hun Lee and Kunwoo Lee, Partial entity structure: a compact
 *  non-manifold boundary representation based on partial topological entities,
 *  SMA '01 Proceedings of the sixth ACM symposium on Solid modeling and
 *  applications Pages 159-170
 *
 * Also, from a design perspective, it would be nice if we could set this up so
 * that the more general non-manifold data structures include less general
 * containers and just extend them with the necessary extra information...
 * Don't know if that's possible, but it would be really nice from a data
 * conversion standpoint...
 *
 * TODO:  This paper may be worth a look from an API design perspective:
 * Topological Operators for Non-Manifold Modeling (1995) 
 * http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.50.1961
 *
 * also potentially useful:
 * https://www.cs.purdue.edu/homes/cmh/distribution/books/geo.html
 * https://cs.nyu.edu/faculty/yap/book/egc/
 *
 */

#ifndef NMG_H
#define NMG_H

#include "common.h"

/* system headers */
#include <stdio.h> /* for FILE */

/* interface headers */
#include "vmath.h"
#include "bu/list.h"
#include "bu/log.h"
#include "bu/magic.h"
#include "bu/ptbl.h"
#include "bg/plane.h"
#include "bn/tol.h"
#include "bv/vlist.h"
#include "vmath.h"

__BEGIN_DECLS

#include "nmg/defines.h"

#include "nmg/debug.h"
#include "nmg/globals.h"
#include "nmg/vertex.h"
#include "nmg/edge.h"
#include "nmg/loop.h"
#include "nmg/face.h"
#include "nmg/shell.h"
#include "nmg/region.h"
#include "nmg/model.h"
#include "nmg/nurb.h"
#include "nmg/ray.h"
#include "nmg/plot.h"
#include "nmg/print.h"
#include "nmg/index.h"
#include "nmg/radial.h"
#include "nmg/visit.h"


/*
 * NOTE: We rely on the fact that the first 32 bits in a struct is the
 * magic number (which is used to identify the struct type).  This may
 * be either a magic value, or an rt_list structure, which starts with
 * a magic number.
 *
 * To these ends, there is a standard ordering for fields in
 * "object-use" structures.  That ordering is:
 *
 *   1) magic number, or rt_list structure
 *   2) pointer to parent
 *   5) pointer to mate
 *   6) pointer to geometry
 *   7) pointer to attributes
 *   8) pointer to child(ren)
 */


struct nmg_inter_struct {
    uint32_t            magic;
    struct bu_ptbl      *l1;            /**< @brief  vertexuses on the line of */
    struct bu_ptbl      *l2;            /**< @brief  intersection between planes */
    fastf_t             *mag1;          /**< @brief  Distances along intersection line */
    fastf_t             *mag2;          /**< @brief  for each vertexuse in l1 and l2. */
    size_t              mag_len;        /**< @brief  Array size of mag1 and mag2 */
    struct shell        *s1;
    struct shell        *s2;
    struct faceuse      *fu1;           /**< @brief  null if l1 comes from a wire */
    struct faceuse      *fu2;           /**< @brief  null if l2 comes from a wire */
    struct bn_tol       tol;
    int                 coplanar;       /**< @brief  a flag */
    struct edge_g_lseg  *on_eg;         /**< @brief  edge_g for line of intersection */
    point_t             pt;             /**< @brief  3D line of intersection */
    vect_t              dir;
    point_t             pt2d;           /**< @brief  2D projection of isect line */
    vect_t              dir2d;
    fastf_t             *vert2d;        /**< @brief  Array of 2d vertex projections [index] */
    size_t              maxindex;       /**< @brief  size of vert2d[] */
    mat_t               proj;           /**< @brief  Matrix to project onto XY plane */
    const uint32_t      *twod;          /**< @brief  ptr to face/edge of 2d projection */
};
#define NMG_CK_INTER_STRUCT(_p) NMG_CKMAG(_p, NMG_INTER_STRUCT_MAGIC, "nmg_inter_struct")


/* From file nmg_mk.c */
/*      MAKE routines */
NMG_EXPORT extern struct model *nmg_mm(void);
NMG_EXPORT extern struct model *nmg_mmr(void);
NMG_EXPORT extern struct nmgregion *nmg_mrsv(struct model *m);
NMG_EXPORT extern struct shell *nmg_msv(struct nmgregion *r_p);
NMG_EXPORT extern struct faceuse *nmg_mf(struct loopuse *lu1);
NMG_EXPORT extern struct loopuse *nmg_mlv(uint32_t *magic,
					  struct vertex *v,
					  int orientation);
NMG_EXPORT extern struct edgeuse *nmg_me(struct vertex *v1,
					 struct vertex *v2,
					 struct shell *s);
NMG_EXPORT extern struct edgeuse *nmg_meonvu(struct vertexuse *vu);
NMG_EXPORT extern struct loopuse *nmg_ml(struct shell *s);
/*      KILL routines */
NMG_EXPORT extern int nmg_keg(struct edgeuse *eu);
NMG_EXPORT extern int nmg_kvu(struct vertexuse *vu);
NMG_EXPORT extern int nmg_kfu(struct faceuse *fu1);
NMG_EXPORT extern int nmg_klu(struct loopuse *lu1);
NMG_EXPORT extern int nmg_keu(struct edgeuse *eu);
NMG_EXPORT extern int nmg_keu_zl(struct shell *s,
				 const struct bn_tol *tol);
NMG_EXPORT extern int nmg_ks(struct shell *s);
NMG_EXPORT extern int nmg_kr(struct nmgregion *r);
NMG_EXPORT extern void nmg_km(struct model *m);
/*      Geometry and Attribute routines */
NMG_EXPORT extern void nmg_vertex_gv(struct vertex *v,
				     const point_t pt);
NMG_EXPORT extern void nmg_vertex_g(struct vertex *v,
				    fastf_t x,
				    fastf_t y,
				    fastf_t z);
NMG_EXPORT extern void nmg_vertexuse_nv(struct vertexuse *vu,
					const vect_t norm);
NMG_EXPORT extern void nmg_vertexuse_a_cnurb(struct vertexuse *vu,
					     const fastf_t *uvw);
NMG_EXPORT extern void nmg_edge_g(struct edgeuse *eu);
NMG_EXPORT extern void nmg_edge_g_cnurb(struct edgeuse *eu,
					int order,
					int n_knots,
					fastf_t *kv,
					int n_pts,
					int pt_type,
					fastf_t *points);
NMG_EXPORT extern void nmg_edge_g_cnurb_plinear(struct edgeuse *eu);
NMG_EXPORT extern int nmg_use_edge_g(struct edgeuse *eu,
				     uint32_t *eg);
NMG_EXPORT extern void nmg_loop_g(struct loop *l,
				  const struct bn_tol *tol);
NMG_EXPORT extern void nmg_face_g(struct faceuse *fu,
				  const plane_t p);
NMG_EXPORT extern void nmg_face_new_g(struct faceuse *fu,
				      const plane_t pl);
NMG_EXPORT extern void nmg_face_g_snurb(struct faceuse *fu,
					int u_order,
					int v_order,
					int n_u_knots,
					int n_v_knots,
					fastf_t *ukv,
					fastf_t *vkv,
					int n_rows,
					int n_cols,
					int pt_type,
					fastf_t *mesh);
NMG_EXPORT extern void nmg_face_bb(struct face *f,
				   const struct bn_tol *tol);
NMG_EXPORT extern void nmg_shell_a(struct shell *s,
				   const struct bn_tol *tol);
NMG_EXPORT extern void nmg_region_a(struct nmgregion *r,
				    const struct bn_tol *tol);
/*      DEMOTE routines */
NMG_EXPORT extern int nmg_demote_lu(struct loopuse *lu);
NMG_EXPORT extern int nmg_demote_eu(struct edgeuse *eu);
/*      MODIFY routines */
NMG_EXPORT extern void nmg_movevu(struct vertexuse *vu,
				  struct vertex *v);
NMG_EXPORT extern void nmg_je(struct edgeuse *eudst,
			      struct edgeuse *eusrc);
NMG_EXPORT extern void nmg_unglueedge(struct edgeuse *eu);
NMG_EXPORT extern void nmg_jv(struct vertex *v1,
			      struct vertex *v2);
NMG_EXPORT extern void nmg_jfg(struct face *f1,
			       struct face *f2);
NMG_EXPORT extern void nmg_jeg(struct edge_g_lseg *dest_eg,
			       struct edge_g_lseg *src_eg);



NMG_EXPORT extern void nmg_count_shell_kids(const struct model *m,
                                            size_t *total_wires,
                                            size_t *total_faces,
                                            size_t *total_points);

/* From nmg_mod.c */

/*      SHELL Routines */
NMG_EXPORT extern void nmg_shell_coplanar_face_merge(struct shell *s,
						     const struct bn_tol *tol,
						     const int simplify,
						     struct bu_list *vlfree);
NMG_EXPORT extern int nmg_simplify_shell(struct shell *s, struct bu_list *vlfree);
NMG_EXPORT extern void nmg_rm_redundancies(struct shell *s,
					   struct bu_list *vlfree,
					   const struct bn_tol *tol);
NMG_EXPORT extern void nmg_sanitize_s_lv(struct shell *s,
					 int orient);
NMG_EXPORT extern void nmg_s_split_touchingloops(struct shell *s,
						 const struct bn_tol *tol);
NMG_EXPORT extern void nmg_s_join_touchingloops(struct shell             *s,
						const struct bn_tol      *tol);
NMG_EXPORT extern void nmg_js(struct shell       *s1,
			      struct shell       *s2,
			      struct bu_list *vlfree,
			      const struct bn_tol        *tol);
NMG_EXPORT extern void nmg_invert_shell(struct shell             *s);

/*      FACE Routines */
NMG_EXPORT extern struct faceuse *nmg_cmface(struct shell *s,
					     struct vertex **vt[],
					     int n);
NMG_EXPORT extern struct faceuse *nmg_cface(struct shell *s,
					    struct vertex **vt,
					    int n);
NMG_EXPORT extern struct faceuse *nmg_add_loop_to_face(struct shell *s,
						       struct faceuse *fu,
						       struct vertex **verts,
						       int n,
						       int dir);
NMG_EXPORT extern int nmg_fu_planeeqn(struct faceuse *fu,
				      const struct bn_tol *tol);
NMG_EXPORT extern void nmg_gluefaces(struct faceuse *fulist[],
				     int n,
				     struct bu_list *vlfree,
				     const struct bn_tol *tol);
NMG_EXPORT extern int nmg_simplify_face(struct faceuse *fu, struct bu_list *vlfree);
NMG_EXPORT extern void nmg_reverse_face(struct faceuse *fu);
NMG_EXPORT extern void nmg_mv_fu_between_shells(struct shell *dest,
						struct shell *src,
						struct faceuse *fu);
NMG_EXPORT extern void nmg_jf(struct faceuse *dest_fu,
			      struct faceuse *src_fu);
NMG_EXPORT extern struct faceuse *nmg_dup_face(struct faceuse *fu,
					       struct shell *s);

/*      VERTEX Routines */
NMG_EXPORT extern void nmg_mv_vu_between_shells(struct shell *dest,
						struct shell *src,
						struct vertexuse *vu);

/* From nmg_info.c */
/* Model routines */
NMG_EXPORT extern struct model *nmg_find_model(const uint32_t *magic_p);
NMG_EXPORT extern struct shell *nmg_find_shell(const uint32_t *magic_p);
NMG_EXPORT extern void nmg_model_bb(point_t min_pt,
				    point_t max_pt,
				    const struct model *m);


/* Edge routines */
NMG_EXPORT extern struct edgeuse *nmg_find_matching_eu_in_s(const struct edgeuse *eu1,
							    const struct shell   *s2);
NMG_EXPORT extern struct edgeuse *nmg_findeu(const struct vertex *v1,
					     const struct vertex *v2,
					     const struct shell *s,
					     const struct edgeuse *eup,
					     int dangling_only);
NMG_EXPORT extern struct edgeuse *nmg_find_eu_in_face(const struct vertex *v1,
						      const struct vertex *v2,
						      const struct faceuse *fu,
						      const struct edgeuse *eup,
						      int dangling_only);
NMG_EXPORT extern struct edgeuse *nmg_find_e(const struct vertex *v1,
					     const struct vertex *v2,
					     const struct shell *s,
					     const struct edge *ep);
NMG_EXPORT extern struct edgeuse *nmg_find_eu_of_vu(const struct vertexuse *vu);
NMG_EXPORT extern struct edgeuse *nmg_find_eu_with_vu_in_lu(const struct loopuse *lu,
							    const struct vertexuse *vu);
NMG_EXPORT extern const struct edgeuse *nmg_faceradial(const struct edgeuse *eu);
NMG_EXPORT extern const struct edgeuse *nmg_radial_face_edge_in_shell(const struct edgeuse *eu);
NMG_EXPORT extern const struct edgeuse *nmg_find_edge_between_2fu(const struct faceuse *fu1,
								  const struct faceuse *fu2,
								  struct bu_list *vlfree,
								  const struct bn_tol *tol);
NMG_EXPORT extern struct edge *nmg_find_e_nearest_pt2(uint32_t *magic_p,
						      const point_t pt2,
						      const mat_t mat,
						      struct bu_list *vlfree,
						      const struct bn_tol *tol);
NMG_EXPORT extern struct edgeuse *nmg_find_matching_eu_in_s(const struct edgeuse *eu1,
							    const struct shell *s2);
NMG_EXPORT extern void nmg_eu_2vecs_perp(vect_t xvec,
					 vect_t yvec,
					 vect_t zvec,
					 const struct edgeuse *eu,
					 const struct bn_tol *tol);
NMG_EXPORT extern int nmg_find_eu_leftvec(vect_t left,
					  const struct edgeuse *eu);
NMG_EXPORT extern int nmg_find_eu_left_non_unit(vect_t left,
						const struct edgeuse     *eu);
NMG_EXPORT extern struct edgeuse *nmg_find_ot_same_eu_of_e(const struct edge *e);

/* Vertex routines */
NMG_EXPORT extern struct vertexuse *nmg_find_v_in_face(const struct vertex *,
						       const struct faceuse *);
NMG_EXPORT extern struct vertexuse *nmg_find_v_in_shell(const struct vertex *v,
							const struct shell *s,
							int edges_only);
NMG_EXPORT extern struct vertexuse *nmg_find_pnt_in_lu(const struct loopuse *lu,
						      const point_t pt,
						      const struct bn_tol *tol);
NMG_EXPORT extern struct vertexuse *nmg_find_pnt_in_face(const struct faceuse *fu,
							const point_t pt,
							const struct bn_tol *tol);
NMG_EXPORT extern struct vertex *nmg_find_pnt_in_shell(const struct shell *s,
						      const point_t pt,
						      const struct bn_tol *tol);
NMG_EXPORT extern struct vertex *nmg_find_pnt_in_model(const struct model *m,
						      const point_t pt,
						      const struct bn_tol *tol);
NMG_EXPORT extern int nmg_is_vertex_in_edgelist(const struct vertex *v,
						const struct bu_list *hd);
NMG_EXPORT extern int nmg_is_vertex_in_looplist(const struct vertex *v,
						const struct bu_list *hd,
						int singletons);
NMG_EXPORT extern struct vertexuse *nmg_is_vertex_in_face(const struct vertex *v,
							  const struct face *f);
NMG_EXPORT extern int nmg_is_vertex_a_selfloop_in_shell(const struct vertex *v,
							const struct shell *s);
NMG_EXPORT extern int nmg_is_vertex_in_facelist(const struct vertex *v,
						const struct bu_list *hd);
NMG_EXPORT extern int nmg_is_edge_in_edgelist(const struct edge *e,
					      const struct bu_list *hd);
NMG_EXPORT extern int nmg_is_edge_in_looplist(const struct edge *e,
					      const struct bu_list *hd);
NMG_EXPORT extern int nmg_is_edge_in_facelist(const struct edge *e,
					      const struct bu_list *hd);
NMG_EXPORT extern int nmg_is_loop_in_facelist(const struct loop *l,
					      const struct bu_list *fu_hd);

/* Tabulation routines */
NMG_EXPORT extern void nmg_vertex_tabulate(struct bu_ptbl *tab,
					   const uint32_t *magic_p,
					   struct bu_list *vlfree);
NMG_EXPORT extern void nmg_vertexuse_normal_tabulate(struct bu_ptbl *tab,
						     const uint32_t *magic_p,
						     struct bu_list *vlfree);
NMG_EXPORT extern void nmg_edgeuse_tabulate(struct bu_ptbl *tab,
					    const uint32_t *magic_p,
					    struct bu_list *vlfree);
NMG_EXPORT extern void nmg_edge_tabulate(struct bu_ptbl *tab,
					 const uint32_t *magic_p,
					 struct bu_list *vlfree);
NMG_EXPORT extern void nmg_edge_g_tabulate(struct bu_ptbl *tab,
					   const uint32_t *magic_p,
					   struct bu_list *vlfree);
NMG_EXPORT extern void nmg_face_tabulate(struct bu_ptbl *tab,
					 const uint32_t *magic_p,
					 struct bu_list *vlfree);
NMG_EXPORT extern void nmg_edgeuse_with_eg_tabulate(struct bu_ptbl *tab,
						    const struct edge_g_lseg *eg);
NMG_EXPORT extern void nmg_edgeuse_on_line_tabulate(struct bu_ptbl *tab,
						    const uint32_t *magic_p,
						    const point_t pt,
						    const vect_t dir,
						    struct bu_list *vlfree,
						    const struct bn_tol *tol);
NMG_EXPORT extern void nmg_e_and_v_tabulate(struct bu_ptbl *eutab,
					    struct bu_ptbl *vtab,
					    const uint32_t *magic_p,
					    struct bu_list *vlfree);
NMG_EXPORT extern int nmg_2edgeuse_g_coincident(const struct edgeuse     *eu1,
						const struct edgeuse     *eu2,
						const struct bn_tol      *tol);

/* From nmg_extrude.c */
NMG_EXPORT extern void nmg_translate_face(struct faceuse *fu, const vect_t Vec, struct bu_list *vlfree);
NMG_EXPORT extern int nmg_extrude_face(struct faceuse *fu, const vect_t Vec, struct bu_list *vlfree, const struct bn_tol *tol);
NMG_EXPORT extern struct vertexuse *nmg_find_vertex_in_lu(const struct vertex *v, const struct loopuse *lu);
NMG_EXPORT extern void nmg_fix_overlapping_loops(struct shell *s, struct bu_list *vlfree, const struct bn_tol *tol);
NMG_EXPORT extern void nmg_break_crossed_loops(struct shell *is, const struct bn_tol *tol);
NMG_EXPORT extern struct shell *nmg_extrude_cleanup(struct shell *is, const int is_void, struct bu_list *vlfree, const struct bn_tol *tol);
NMG_EXPORT extern void nmg_hollow_shell(struct shell *s, const fastf_t thick, const int approximate, struct bu_list *vlfree, const struct bn_tol *tol);
NMG_EXPORT extern struct shell *nmg_extrude_shell(struct shell *s, const fastf_t dist, const int normal_ward, const int approximate, struct bu_list *vlfree, const struct bn_tol *tol);

/* From nmg_misc.c */
NMG_EXPORT extern int nmg_snurb_calc_lu_uv_orient(const struct loopuse *lu);
NMG_EXPORT extern void nmg_snurb_fu_eval(const struct faceuse *fu,
					 const fastf_t u,
					 const fastf_t v,
					 point_t pt_on_srf);
NMG_EXPORT extern void nmg_snurb_fu_get_norm(const struct faceuse *fu,
					     const fastf_t u,
					     const fastf_t v,
					     vect_t norm);
NMG_EXPORT extern void nmg_snurb_fu_get_norm_at_vu(const struct faceuse *fu,
						   const struct vertexuse *vu,
						   vect_t norm);
NMG_EXPORT extern void nmg_find_zero_length_edges(const struct model *m, struct bu_list *vlfree);
NMG_EXPORT extern struct face *nmg_find_top_face_in_dir(const struct shell *s,
							int dir, long *flags);
NMG_EXPORT extern struct face *nmg_find_top_face(const struct shell *s,
						 int *dir,
						 long *flags);
NMG_EXPORT extern int nmg_find_outer_and_void_shells(struct nmgregion *r,
						     struct bu_ptbl ***shells,
						     struct bu_list *vlfree,
						     const struct bn_tol *tol);
NMG_EXPORT extern int nmg_mark_edges_real(const uint32_t *magic_p, struct bu_list *vlfree);
NMG_EXPORT extern void nmg_tabulate_face_g_verts(struct bu_ptbl *tab,
						 const struct face_g_plane *fg);
NMG_EXPORT extern void nmg_isect_shell_self(struct shell *s,
					    struct bu_list *vlfree,
					    const struct bn_tol *tol);
NMG_EXPORT extern struct edgeuse *nmg_next_radial_eu(const struct edgeuse *eu,
						     const struct shell *s,
						     const int wires);
NMG_EXPORT extern struct edgeuse *nmg_prev_radial_eu(const struct edgeuse *eu,
						     const struct shell *s,
						     const int wires);
NMG_EXPORT extern int nmg_radial_face_count(const struct edgeuse *eu,
					    const struct shell *s);
NMG_EXPORT extern int nmg_check_closed_shell(const struct shell *s,
					     const struct bn_tol *tol);
NMG_EXPORT extern int nmg_move_lu_between_fus(struct faceuse *dest,
					      struct faceuse *src,
					      struct loopuse *lu);
NMG_EXPORT extern void nmg_loop_plane_newell(const struct loopuse *lu,
					     plane_t pl);
NMG_EXPORT extern fastf_t nmg_loop_plane_area(const struct loopuse *lu,
					      plane_t pl);
NMG_EXPORT extern fastf_t nmg_loop_plane_area2(const struct loopuse *lu,
					       plane_t pl,
					       const struct bn_tol *tol);
NMG_EXPORT extern int nmg_calc_face_plane(struct faceuse *fu_in,
					  plane_t pl, struct bu_list *vlfree);
NMG_EXPORT extern int nmg_calc_face_g(struct faceuse *fu, struct bu_list *vlfree);
NMG_EXPORT extern fastf_t nmg_faceuse_area(const struct faceuse *fu);
NMG_EXPORT extern fastf_t nmg_shell_area(const struct shell *s);
NMG_EXPORT extern fastf_t nmg_region_area(const struct nmgregion *r);
NMG_EXPORT extern fastf_t nmg_model_area(const struct model *m);
/* Some stray rt_ plane functions here */
NMG_EXPORT extern void nmg_purge_unwanted_intersection_points(struct bu_ptbl *vert_list,
							      fastf_t *mag,
							      const struct faceuse *fu,
							      const struct bn_tol *tol);
NMG_EXPORT extern int nmg_in_or_ref(struct vertexuse *vu,
				    struct bu_ptbl *b);
NMG_EXPORT extern void nmg_rebound(struct model *m,
				   const struct bn_tol *tol);
NMG_EXPORT extern struct edgeuse *nmg_pop_eu(struct bu_ptbl *stack);
NMG_EXPORT extern void nmg_reverse_radials(struct faceuse *fu,
					   const struct bn_tol *tol);
NMG_EXPORT extern void nmg_reverse_face_and_radials(struct faceuse *fu,
						    const struct bn_tol *tol);
NMG_EXPORT extern int nmg_shell_is_void(const struct shell *s);
NMG_EXPORT extern void nmg_propagate_normals(struct faceuse *fu_in,
					     long *flags,
					     const struct bn_tol *tol);
NMG_EXPORT extern void nmg_connect_same_fu_orients(struct shell *s);
NMG_EXPORT extern void nmg_fix_decomposed_shell_normals(struct shell *s,
							const struct bn_tol *tol);
NMG_EXPORT extern struct model *nmg_mk_model_from_region(struct nmgregion *r,
							 int reindex, struct bu_list *vlfree);
NMG_EXPORT extern void nmg_fix_normals(struct shell *s_orig,
				       struct bu_list *vlfree,
				       const struct bn_tol *tol);
NMG_EXPORT extern int nmg_break_long_edges(struct shell *s,
					   const struct bn_tol *tol);
NMG_EXPORT extern struct faceuse *nmg_mk_new_face_from_loop(struct loopuse *lu);
NMG_EXPORT extern int nmg_split_loops_into_faces(uint32_t *magic_p, struct bu_list *vlfree,
						 const struct bn_tol *tol);
NMG_EXPORT extern int nmg_decompose_shell(struct shell *s, struct bu_list *vlfree,
					  const struct bn_tol *tol);
NMG_EXPORT extern int nmg_unbreak_region_edges(uint32_t *magic_p, struct bu_list *vlfree);
NMG_EXPORT extern void nmg_vlist_to_eu(struct bu_list *vlist,
				       struct shell *s);
NMG_EXPORT extern int nmg_mv_shell_to_region(struct shell *s,
					     struct nmgregion *r);
NMG_EXPORT extern int nmg_find_isect_faces(const struct vertex *new_v,
					   struct bu_ptbl *faces,
					   int *free_edges,
					   const struct bn_tol *tol);
NMG_EXPORT extern int nmg_simple_vertex_solve(struct vertex *new_v,
					      const struct bu_ptbl *faces,
					      const struct bn_tol *tol);
NMG_EXPORT extern int nmg_ck_vert_on_fus(const struct vertex *v,
					 const struct bn_tol *tol);
NMG_EXPORT extern void nmg_make_faces_at_vert(struct vertex *new_v,
					      struct bu_ptbl *int_faces,
					      struct bu_list *vlfree,
					      const struct bn_tol *tol);
NMG_EXPORT extern void nmg_kill_cracks_at_vertex(const struct vertex *vp);
NMG_EXPORT extern int nmg_complex_vertex_solve(struct vertex *new_v,
					       const struct bu_ptbl *faces,
					       const int free_edges,
					       const int approximate,
					       struct bu_list *vlfree,
					       const struct bn_tol *tol);
NMG_EXPORT extern int nmg_bad_face_normals(const struct shell *s,
					   const struct bn_tol *tol);
NMG_EXPORT extern int nmg_faces_are_radial(const struct faceuse *fu1,
					   const struct faceuse *fu2);
NMG_EXPORT extern int nmg_move_edge_thru_pnt(struct edgeuse *mv_eu,
					     const point_t pt,
					     const struct bn_tol *tol);
NMG_EXPORT extern void nmg_vlist_to_wire_edges(struct shell *s,
					       const struct bu_list *vhead);
NMG_EXPORT extern void nmg_follow_free_edges_to_vertex(const struct vertex *vpa,
						       const struct vertex *vpb,
						       struct bu_ptbl *bad_verts,
						       const struct shell *s,
						       const struct edgeuse *eu,
						       struct bu_ptbl *verts,
						       int *found);
NMG_EXPORT extern int nmg_in_vert(struct vertex *new_v,
				  const int approximate,
				  struct bu_list *vlfree,
				  const struct bn_tol *tol);
NMG_EXPORT extern void nmg_mirror_model(struct model *m, struct bu_list *vlfree);
NMG_EXPORT extern int nmg_kill_cracks(struct shell *s);
NMG_EXPORT extern int nmg_kill_zero_length_edgeuses(struct model *m);
NMG_EXPORT extern void nmg_make_faces_within_tol(struct shell *s,
						 struct bu_list *vlfree,
						 const struct bn_tol *tol);
NMG_EXPORT extern void nmg_intersect_loops_self(struct shell *s,
						const struct bn_tol *tol);
NMG_EXPORT extern struct edge_g_cnurb *nmg_join_cnurbs(struct bu_list *crv_head);
NMG_EXPORT extern struct edge_g_cnurb *nmg_arc2d_to_cnurb(point_t i_center,
							  point_t i_start,
							  point_t i_end,
							  int point_type,
							  const struct bn_tol *tol);
NMG_EXPORT extern int nmg_break_edge_at_verts(struct edge *e,
					      struct bu_ptbl *verts,
					      const struct bn_tol *tol);
NMG_EXPORT extern fastf_t nmg_loop_plane_area(const struct loopuse *lu,
					      plane_t pl);
NMG_EXPORT extern int nmg_break_edges(uint32_t *magic_p, struct bu_list *vlfree,
				      const struct bn_tol *tol);
NMG_EXPORT extern int nmg_lu_is_convex(struct loopuse *lu,
				       struct bu_list *vlfree,
				       const struct bn_tol *tol);
NMG_EXPORT extern int nmg_edge_collapse(struct model *m,
					const struct bn_tol *tol,
					const fastf_t tol_coll,
					const fastf_t min_angle,
					struct bu_list *vlfree);

/* From nmg_copy.c */
NMG_EXPORT extern struct model *nmg_clone_model(const struct model *original);

/* From nmg_tri.c */
NMG_EXPORT extern void nmg_triangulate_shell(struct shell *s,
					     struct bu_list *vlfree,
					     const struct bn_tol  *tol);


NMG_EXPORT extern void nmg_triangulate_model(struct model *m,
					     struct bu_list *vlfree,
					     const struct bn_tol *tol);
NMG_EXPORT extern int nmg_triangulate_fu(struct faceuse *fu,
					 struct bu_list *vlfree,
					 const struct bn_tol *tol);
NMG_EXPORT extern void nmg_dump_model(struct model *m);

/* nmg_manif.c */
NMG_EXPORT extern int nmg_dangling_face(const struct faceuse *fu,
					const char *manifolds);
/* static paint_face */
/* static set_edge_sub_manifold */
/* static set_loop_sub_manifold */
/* static set_face_sub_manifold */
NMG_EXPORT extern char *nmg_manifolds(struct model *m);

/* From nmg_fuse.c */
NMG_EXPORT extern int nmg_is_common_bigloop(const struct face *f1,
					    const struct face *f2);
NMG_EXPORT extern void nmg_region_v_unique(struct nmgregion *r1, struct bu_list *vlfree,
					   const struct bn_tol *tol);
NMG_EXPORT extern int nmg_ptbl_vfuse(struct bu_ptbl *t,
				     const struct bn_tol *tol);
NMG_EXPORT extern int	nmg_region_both_vfuse(struct bu_ptbl *t1,
					      struct bu_ptbl *t2,
					      const struct bn_tol *tol);
/* nmg_two_region_vertex_fuse replaced with nmg_vertex_fuse */
NMG_EXPORT extern int nmg_vertex_fuse(const uint32_t *magic_p,struct bu_list *vlfree,
				      const struct bn_tol *tol);
NMG_EXPORT extern int nmg_cnurb_is_linear(const struct edge_g_cnurb *cnrb);
NMG_EXPORT extern int nmg_snurb_is_planar(const struct face_g_snurb *srf,
					  const struct bn_tol *tol);
NMG_EXPORT extern void nmg_eval_linear_trim_curve(const struct face_g_snurb *snrb,
						  const fastf_t uvw[3],
						  point_t xyz);
NMG_EXPORT extern void nmg_eval_trim_curve(const struct edge_g_cnurb *cnrb,
					   const struct face_g_snurb *snrb,
					   const fastf_t t, point_t xyz);
/* nmg_split_trim */
NMG_EXPORT extern void nmg_eval_trim_to_tol(const struct edge_g_cnurb *cnrb,
					    const struct face_g_snurb *snrb,
					    const fastf_t t0,
					    const fastf_t t1,
					    struct bu_list *head,
					    const struct bn_tol *tol);
/* nmg_split_linear_trim */
NMG_EXPORT extern void nmg_eval_linear_trim_to_tol(const struct edge_g_cnurb *cnrb,
						   const struct face_g_snurb *snrb,
						   const fastf_t uvw1[3],
						   const fastf_t uvw2[3],
						   struct bu_list *head,
						   const struct bn_tol *tol);
NMG_EXPORT extern int nmg_cnurb_lseg_coincident(const struct edgeuse *eu1,
						const struct edge_g_cnurb *cnrb,
						const struct face_g_snurb *snrb,
						const point_t pt1,
						const point_t pt2,
						const struct bn_tol *tol);
NMG_EXPORT extern int	nmg_cnurb_is_on_crv(const struct edgeuse *eu,
					    const struct edge_g_cnurb *cnrb,
					    const struct face_g_snurb *snrb,
					    const struct bu_list *head,
					    const struct bn_tol *tol);
NMG_EXPORT extern int nmg_edge_fuse(const uint32_t *magic_p,struct bu_list *vlfree,
				    const struct bn_tol *tol);
NMG_EXPORT extern int nmg_edge_g_fuse(const uint32_t *magic_p,struct bu_list *vlfree,
				      const struct bn_tol	*tol);
NMG_EXPORT extern int nmg_ck_fu_verts(struct faceuse *fu1,
				      struct face *f2,
				      const struct bn_tol *tol);
NMG_EXPORT extern int nmg_ck_fg_verts(struct faceuse *fu1,
				      struct face *f2,
				      const struct bn_tol *tol);
NMG_EXPORT extern int	nmg_two_face_fuse(struct face	*f1,
					  struct face *f2,
					  const struct bn_tol *tol);
NMG_EXPORT extern int nmg_model_face_fuse(struct model *m,struct bu_list *vlfree,
					  const struct bn_tol *tol);
NMG_EXPORT extern int nmg_break_all_es_on_v(uint32_t *magic_p,
					    struct vertex *v,struct bu_list *vlfree,
					    const struct bn_tol *tol);
NMG_EXPORT extern int nmg_break_e_on_v(const uint32_t *magic_p,struct bu_list *vlfree,
				       const struct bn_tol *tol);
/* DEPRECATED: use nmg_break_e_on_v */
NMG_EXPORT extern int nmg_model_break_e_on_v(const uint32_t *magic_p,
					     struct bu_list *vlfree,
					     const struct bn_tol *tol);
NMG_EXPORT extern int nmg_model_fuse(struct model *m,
				     struct bu_list *vlfree,
				     const struct bn_tol *tol);


NMG_EXPORT extern struct edge_g_lseg	*nmg_pick_best_edge_g(struct edgeuse *eu1,
							      struct edgeuse *eu2,
							      const struct bn_tol *tol);

/* nmg_class.c */
NMG_EXPORT extern int nmg_classify_pnt_loop(const point_t pt,
					   const struct loopuse *lu,
					   struct bu_list *vlfree,
					   const struct bn_tol *tol);

NMG_EXPORT extern int nmg_classify_s_vs_s(struct shell *s,
					  struct shell *s2,
					  struct bu_list *vlfree,
					  const struct bn_tol *tol);

NMG_EXPORT extern int nmg_classify_lu_lu(const struct loopuse *lu1,
					 const struct loopuse *lu2,
					 struct bu_list *vlfree,
					 const struct bn_tol *tol);

NMG_EXPORT extern int nmg_class_pnt_f(const point_t pt,
				     const struct faceuse *fu,
				     const struct bn_tol *tol);
NMG_EXPORT extern int nmg_class_pnt_s(const point_t pt,
				     const struct shell *s,
				     const int in_or_out_only,
				     struct bu_list *vlfree,
				     const struct bn_tol *tol);

/* From nmg_pt_fu.c */
NMG_EXPORT extern int nmg_eu_is_part_of_crack(const struct edgeuse *eu);

NMG_EXPORT extern int nmg_class_pnt_lu_except(point_t             pt,
					     const struct loopuse        *lu,
					     const struct edge           *e_p,
					     struct bu_list *vlfree,
					     const struct bn_tol *tol);

NMG_EXPORT extern int nmg_class_pnt_fu_except(const point_t pt,
					     const struct faceuse *fu,
					     const struct loopuse *ignore_lu,
					     void (*eu_func)(struct edgeuse *, point_t, const char *, struct bu_list *),
					     void (*vu_func)(struct vertexuse *, point_t, const char *),
					     const char *priv,
					     const int call_on_hits,
					     const int in_or_out_only,
					     struct bu_list *vlfree,
					     const struct bn_tol *tol);

/* from nmg_mesh.c */
NMG_EXPORT extern int nmg_mesh_two_faces(struct faceuse *fu1,
					 struct faceuse *fu2,
					 const struct bn_tol     *tol);
NMG_EXPORT extern void nmg_radial_join_eu(struct edgeuse *eu1,
					  struct edgeuse *eu2,
					  const struct bn_tol *tol);
NMG_EXPORT extern void nmg_mesh_faces(struct faceuse *fu1,
				      struct faceuse *fu2,
				      struct bu_list *vlfree,
				      const struct bn_tol *tol);
NMG_EXPORT extern int nmg_mesh_face_shell(struct faceuse *fu1,
					  struct shell *s,
					  const struct bn_tol *tol);
NMG_EXPORT extern int nmg_mesh_shell_shell(struct shell *s1,
					   struct shell *s2,
					   struct bu_list *vlfree,
					   const struct bn_tol *tol);
NMG_EXPORT extern double nmg_measure_fu_angle(const struct edgeuse *eu,
					      const vect_t xvec,
					      const vect_t yvec,
					      const vect_t zvec);

/* from nmg_bool.c */
NMG_EXPORT extern struct nmgregion *nmg_do_bool(struct nmgregion *s1,
						struct nmgregion *s2,
						const int oper, struct bu_list *vlfree, const struct bn_tol *tol);
NMG_EXPORT extern int nmg_two_region_vertex_fuse(struct nmgregion *r1,
						 struct nmgregion *r2,
						 const struct bn_tol *tol);


/* from nmg_class.c */
NMG_EXPORT extern void nmg_class_shells(struct shell *sA,
					struct shell *sB,
					char **classlist,
					struct bu_list *vlfree,
					const struct bn_tol *tol);

/* from nmg_fcut.c */
/* static void ptbl_vsort */
NMG_EXPORT extern int nmg_ck_vu_ptbl(struct bu_ptbl      *p,
				     struct faceuse      *fu);
NMG_EXPORT extern double nmg_vu_angle_measure(struct vertexuse   *vu,
					      vect_t x_dir,
					      vect_t y_dir,
					      int assessment,
					      int in);
NMG_EXPORT extern int nmg_wedge_class(int        ass,    /* assessment of two edges forming wedge */
				      double     a,
				      double     b);
NMG_EXPORT extern void nmg_sanitize_fu(struct faceuse    *fu);
NMG_EXPORT extern void nmg_unlist_v(struct bu_ptbl       *b,
				    fastf_t *mag,
				    struct vertex        *v);
NMG_EXPORT extern struct edge_g_lseg *nmg_face_cutjoin(struct bu_ptbl *b1,
						       struct bu_ptbl *b2,
						       fastf_t *mag1,
						       fastf_t *mag2,
						       struct faceuse *fu1,
						       struct faceuse *fu2,
						       point_t pt,
						       vect_t dir,
						       struct edge_g_lseg *eg,
						       struct bu_list *vlfree,
						       const struct bn_tol *tol);
NMG_EXPORT extern void nmg_fcut_face_2d(struct bu_ptbl *vu_list,
					fastf_t *mag,
					struct faceuse *fu1,
					struct faceuse *fu2,
					struct bu_list *vlfree,
					struct bn_tol *tol);
NMG_EXPORT extern int nmg_insert_vu_if_on_edge(struct vertexuse *vu1,
					       struct vertexuse *vu2,
					       struct edgeuse *new_eu,
					       struct bn_tol *tol);
/* nmg_face_state_transition */

#define nmg_mev(_v, _u) nmg_me((_v), (struct vertex *)NULL, (_u))

/* From nmg_eval.c */
NMG_EXPORT extern void nmg_ck_lu_orientation(struct loopuse *lu,
					     const struct bn_tol *tolp);
NMG_EXPORT extern const char *nmg_class_name(int class_no);
NMG_EXPORT extern void nmg_evaluate_boolean(struct shell *sA,
					    struct shell *sB,
					    int          op,
					    char         **classlist,
					    struct bu_list *vlfree,
					    const struct bn_tol  *tol);

/* The following functions cannot be publicly declared because struct
 * nmg_bool_state is private to nmg_eval.c
 */
/* nmg_eval_shell */
/* nmg_eval_action */
/* nmg_eval_plot */

/* From nmg_ck.c */
NMG_EXPORT extern void nmg_vvg(const struct vertex_g *vg);
NMG_EXPORT extern void nmg_vvertex(const struct vertex *v,
				   const struct vertexuse *vup);
NMG_EXPORT extern void nmg_vvua(const uint32_t *vua);
NMG_EXPORT extern void nmg_vvu(const struct vertexuse *vu,
			       const uint32_t *up_magic_p);
NMG_EXPORT extern void nmg_veg(const uint32_t *eg);
NMG_EXPORT extern void nmg_vedge(const struct edge *e,
				 const struct edgeuse *eup);
NMG_EXPORT extern void nmg_veu(const struct bu_list      *hp,
			       const uint32_t *up_magic_p);
NMG_EXPORT extern void nmg_vlg(const struct loop_g *lg);
NMG_EXPORT extern void nmg_vloop(const struct loop *l,
				 const struct loopuse *lup);
NMG_EXPORT extern void nmg_vlu(const struct bu_list      *hp,
			       const uint32_t *up);
NMG_EXPORT extern void nmg_vfg(const struct face_g_plane *fg);
NMG_EXPORT extern void nmg_vface(const struct face *f,
				 const struct faceuse *fup);
NMG_EXPORT extern void nmg_vfu(const struct bu_list      *hp,
			       const struct shell *s);
NMG_EXPORT extern void nmg_vsshell(const struct shell *s,
				   const struct nmgregion *r);
NMG_EXPORT extern void nmg_vshell(const struct bu_list *hp,
				  const struct nmgregion *r);
NMG_EXPORT extern void nmg_vregion(const struct bu_list *hp,
				   const struct model *m);
NMG_EXPORT extern void nmg_vmodel(const struct model *m);

/* checking routines */
NMG_EXPORT extern void nmg_ck_e(const struct edgeuse *eu,
				const struct edge *e,
				const char *str);
NMG_EXPORT extern void nmg_ck_vu(const uint32_t *parent,
				 const struct vertexuse *vu,
				 const char *str);
NMG_EXPORT extern void nmg_ck_eu(const uint32_t *parent,
				 const struct edgeuse *eu,
				 const char *str);
NMG_EXPORT extern void nmg_ck_lg(const struct loop *l,
				 const struct loop_g *lg,
				 const char *str);
NMG_EXPORT extern void nmg_ck_l(const struct loopuse *lu,
				const struct loop *l,
				const char *str);
NMG_EXPORT extern void nmg_ck_lu(const uint32_t *parent,
				 const struct loopuse *lu,
				 const char *str);
NMG_EXPORT extern void nmg_ck_fg(const struct face *f,
				 const struct face_g_plane *fg,
				 const char *str);
NMG_EXPORT extern void nmg_ck_f(const struct faceuse *fu,
				const struct face *f,
				const char *str);
NMG_EXPORT extern void nmg_ck_fu(const struct shell *s,
				 const struct faceuse *fu,
				 const char *str);
NMG_EXPORT extern int nmg_ck_eg_verts(const struct edge_g_lseg *eg,
				      const struct bn_tol *tol);
NMG_EXPORT extern size_t nmg_ck_geometry(const struct model *m,
					 struct bu_list *vlfree,
					 const struct bn_tol *tol);
NMG_EXPORT extern int nmg_ck_face_worthless_edges(const struct faceuse *fu);
NMG_EXPORT extern void nmg_ck_lueu(const struct loopuse *lu, const char *s);
NMG_EXPORT extern int nmg_check_radial(const struct edgeuse *eu, const struct bn_tol *tol);
NMG_EXPORT extern int nmg_eu_2s_orient_bad(const struct edgeuse  *eu,
					   const struct shell    *s1,
					   const struct shell    *s2,
					   const struct bn_tol   *tol);
NMG_EXPORT extern int nmg_ck_closed_surf(const struct shell *s,
					 const struct bn_tol *tol);
NMG_EXPORT extern int nmg_ck_closed_region(const struct nmgregion *r,
					   const struct bn_tol *tol);
NMG_EXPORT extern void nmg_ck_v_in_2fus(const struct vertex *vp,
					const struct faceuse *fu1,
					const struct faceuse *fu2,
					const struct bn_tol *tol);
NMG_EXPORT extern void nmg_ck_vs_in_region(const struct nmgregion *r,
					   struct bu_list *vlfree,
					   const struct bn_tol *tol);


/* From nmg_inter.c */
NMG_EXPORT extern struct vertexuse *nmg_make_dualvu(struct vertex *v,
						    struct faceuse *fu,
						    const struct bn_tol *tol);
NMG_EXPORT extern struct vertexuse *nmg_enlist_vu(struct nmg_inter_struct        *is,
						  const struct vertexuse *vu,
						  struct vertexuse *dualvu,
						  fastf_t dist);
NMG_EXPORT extern void nmg_isect2d_prep(struct nmg_inter_struct *is,
					const uint32_t *assoc_use);
NMG_EXPORT extern void nmg_isect2d_cleanup(struct nmg_inter_struct *is);
NMG_EXPORT extern void nmg_isect2d_final_cleanup(void);
NMG_EXPORT extern int nmg_isect_2faceuse(point_t pt,
					 vect_t dir,
					 struct faceuse *fu1,
					 struct faceuse *fu2,
					 const struct bn_tol *tol);
NMG_EXPORT extern void nmg_isect_vert2p_face2p(struct nmg_inter_struct *is,
					       struct vertexuse *vu1,
					       struct faceuse *fu2);
NMG_EXPORT extern struct edgeuse *nmg_break_eu_on_v(struct edgeuse *eu1,
						    struct vertex *v2,
						    struct faceuse *fu,
						    struct nmg_inter_struct *is);
NMG_EXPORT extern void nmg_break_eg_on_v(const struct edge_g_lseg        *eg,
					 struct vertex           *v,
					 const struct bn_tol     *tol);
NMG_EXPORT extern int nmg_isect_2colinear_edge2p(struct edgeuse  *eu1,
						 struct edgeuse  *eu2,
						 struct faceuse          *fu,
						 struct nmg_inter_struct *is,
						 struct bu_ptbl          *l1,
						 struct bu_ptbl          *l2);
NMG_EXPORT extern int nmg_isect_edge2p_edge2p(struct nmg_inter_struct    *is,
					      struct edgeuse             *eu1,
					      struct edgeuse             *eu2,
					      struct faceuse             *fu1,
					      struct faceuse             *fu2);
NMG_EXPORT extern int nmg_isect_construct_nice_ray(struct nmg_inter_struct       *is,
						   struct faceuse                *fu2);
NMG_EXPORT extern void nmg_enlist_one_vu(struct nmg_inter_struct *is,
					 const struct vertexuse  *vu,
					 fastf_t                 dist);
NMG_EXPORT extern int    nmg_isect_line2_edge2p(struct nmg_inter_struct  *is,
						struct bu_ptbl           *list,
						struct edgeuse           *eu1,
						struct faceuse           *fu1,
						struct faceuse           *fu2);
NMG_EXPORT extern void nmg_isect_line2_vertex2(struct nmg_inter_struct   *is,
					       struct vertexuse  *vu1,
					       struct faceuse            *fu1);
NMG_EXPORT extern int nmg_isect_two_ptbls(struct nmg_inter_struct        *is,
					  const struct bu_ptbl   *t1,
					  const struct bu_ptbl   *t2);
NMG_EXPORT extern struct edge_g_lseg     *nmg_find_eg_on_line(const uint32_t *magic_p,
							      const point_t              pt,
							      const vect_t               dir,
							      struct bu_list *vlfree,
							      const struct bn_tol        *tol);
NMG_EXPORT extern int nmg_k0eu(struct vertex     *v);
NMG_EXPORT extern struct vertex *nmg_repair_v_near_v(struct vertex               *hit_v,
						     struct vertex               *v,
						     const struct edge_g_lseg    *eg1,
						     const struct edge_g_lseg    *eg2,
						     int                 bomb,
						     const struct bn_tol *tol);
NMG_EXPORT extern struct vertex *nmg_search_v_eg(const struct edgeuse    *eu,
						 int                     second,
						 const struct edge_g_lseg        *eg1,
						 const struct edge_g_lseg        *eg2,
						 struct vertex           *hit_v,
						 const struct bn_tol     *tol);
NMG_EXPORT extern struct vertex *nmg_common_v_2eg(struct edge_g_lseg     *eg1,
						  struct edge_g_lseg     *eg2,
						  const struct bn_tol    *tol);
NMG_EXPORT extern int nmg_is_vertex_on_inter(struct vertex *v,
					     struct faceuse *fu1,
					     struct faceuse *fu2,
					     struct nmg_inter_struct *is,
					     struct bu_list *vlfree);
NMG_EXPORT extern void nmg_isect_eu_verts(struct edgeuse *eu,
					  struct vertex_g *vg1,
					  struct vertex_g *vg2,
					  struct bu_ptbl *verts,
					  struct bu_ptbl *inters,
					  const struct bn_tol *tol);
NMG_EXPORT extern void nmg_isect_eu_eu(struct edgeuse *eu1,
				       struct vertex_g *vg1a,
				       struct vertex_g *vg1b,
				       vect_t dir1,
				       struct edgeuse *eu2,
				       struct bu_ptbl *verts,
				       struct bu_ptbl *inters,
				       const struct bn_tol *tol);
NMG_EXPORT extern void nmg_isect_eu_fu(struct nmg_inter_struct *is,
				       struct bu_ptbl            *verts,
				       struct edgeuse            *eu,
				       struct faceuse          *fu,
				       struct bu_list *vlfree);
NMG_EXPORT extern void nmg_isect_fu_jra(struct nmg_inter_struct  *is,
					struct faceuse           *fu1,
					struct faceuse           *fu2,
					struct bu_ptbl           *eu1_list,
					struct bu_ptbl           *eu2_list,
					struct bu_list *vlfree);
NMG_EXPORT extern void nmg_isect_line2_face2pNEW(struct nmg_inter_struct *is,
						 struct faceuse *fu1, struct faceuse *fu2,
						 struct bu_ptbl *eu1_list,
						 struct bu_ptbl *eu2_list,
						 struct bu_list *vlfree);
NMG_EXPORT extern int    nmg_is_eu_on_line3(const struct edgeuse *eu,
					    const point_t                pt,
					    const vect_t         dir,
					    const struct bn_tol  *tol);
NMG_EXPORT extern struct edge_g_lseg     *nmg_find_eg_between_2fg(const struct faceuse   *ofu1,
								  const struct faceuse   *fu2,
								  struct bu_list *vlfree,
								  const struct bn_tol    *tol);
NMG_EXPORT extern struct edgeuse *nmg_does_fu_use_eg(const struct faceuse        *fu1,
						     const uint32_t      *eg);
NMG_EXPORT extern int rt_line_on_plane(const point_t     pt,
				       const vect_t      dir,
				       const plane_t     plane,
				       const struct bn_tol       *tol);
NMG_EXPORT extern void nmg_cut_lu_into_coplanar_and_non(struct loopuse *lu,
							plane_t pl,
							struct nmg_inter_struct *is,
							struct bu_list *vlfree);
NMG_EXPORT extern void nmg_check_radial_angles(char *str,
					       struct shell *s,
					       struct bu_list *vlfree,
					       const struct bn_tol *tol);
NMG_EXPORT extern int nmg_faces_can_be_intersected(struct nmg_inter_struct *bs,
						   const struct faceuse *fu1,
						   const struct faceuse *fu2,
						   struct bu_list *vlfree,
						   const struct bn_tol *tol);
NMG_EXPORT extern void nmg_isect_two_generic_faces(struct faceuse                *fu1,
						   struct faceuse                *fu2,
						   struct bu_list *vlfree,
						   const struct bn_tol   *tol);
NMG_EXPORT extern void nmg_crackshells(struct shell *s1,
				       struct shell *s2,
				       struct bu_list *vlfree,
				       const struct bn_tol *tol);
NMG_EXPORT extern int nmg_fu_touchingloops(const struct faceuse *fu);


/* From nmg_index.c */
NMG_EXPORT extern void nmg_m_set_high_bit(struct model *m);
NMG_EXPORT extern void nmg_m_reindex(struct model *m, long newindex);
NMG_EXPORT extern void nmg_vls_struct_counts(struct bu_vls *str,
					     const struct nmg_struct_counts *ctr);
NMG_EXPORT extern void nmg_pr_struct_counts(const struct nmg_struct_counts *ctr,
					    const char *str);
NMG_EXPORT extern uint32_t **nmg_m_struct_count(struct nmg_struct_counts *ctr,
						const struct model *m);
NMG_EXPORT extern void nmg_pr_m_struct_counts(const struct model *m,
					      const char         *str);
NMG_EXPORT extern void nmg_merge_models(struct model *m1,
					struct model *m2);
NMG_EXPORT extern long nmg_find_max_index(const struct model *m);

/* From nmg_rt_isect.c */
NMG_EXPORT extern int nmg_class_ray_vs_shell(struct nmg_ray *rp,
					     const struct shell *s,
					     const int in_or_out_only,
					     struct bu_list *vlfree,
					     const struct bn_tol *tol);

NMG_EXPORT extern void nmg_isect_ray_model(struct nmg_ray_data *rd, struct bu_list *vlfree);

__END_DECLS

#endif /* NMG_H */

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
