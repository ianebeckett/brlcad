#include "bucket-pr.h"
#include "bucket-pr_impl.hpp"

#include "rt/primitives/bot.h"
#include "rt/rt_instance.h"
#include "rt/tie.h"
#include "bio.h"
#include "rt/geom.h"

template<typename Float>
using tree_t = jk::tree::KDTree<tie_tri_s, 3, 32UL, jk::tree::SquaredL2, Float>;

template<typename Float>
void bucketpr_push(struct tie_s *tie, TIE_3 **tlist, unsigned int tnum, void *plist, unsigned int pstride, tree_t<Float> tree) {
    unsigned int i;

    /* expand the tri buffer if needed */
    if (tnum + tie->tri_num > tie->tri_num_alloc) {
	tie->tri_list = (struct tie_tri_s *)bu_realloc(tie->tri_list, sizeof(struct tie_tri_s) * (tie->tri_num + tnum), "tri_list during tie_push");
	tie->tri_num_alloc += tnum;
    }

    for (i = 0; i < tnum; i++) {

	if (tie_check_degenerate) {
	    TIE_3 u, v, w;

	    VSUB2(u.v,  (*tlist[i*3+1]).v,  (*tlist[i*3+0]).v);
	    VSUB2(v.v,  (*tlist[i*3+2]).v,  (*tlist[i*3+0]).v);
	    VCROSS(w.v,  u.v,  v.v);

	    if (MAGNITUDE(w.v) < 0.0001 * 0.0001) {
		bu_log("WARNING: degenerate triangle found: %f %f %f | %f %f %f | %f %f %f\n",
		       V3ARGS((*tlist[i*3+0]).v),  V3ARGS((*tlist[i*3+1]).v), V3ARGS((*tlist[i*3+2]).v));
		continue;
	    }
	}

	/* pack pack pack */
	tie->tri_list[tie->tri_num].data[0] = *tlist[i*3+0];
	tie->tri_list[tie->tri_num].data[1] = *tlist[i*3+1];
	tie->tri_list[tie->tri_num].data[2] = *tlist[i*3+2];

	/* set the association pointer */
	tie->tri_list[tie->tri_num].ptr = plist;
	if (plist)
	    plist = (void *)((intptr_t)plist + pstride);

	V2SETALL(tie->tri_list[tie->tri_num].v, 0.0);
	tie->tri_list[tie->tri_num].b = 0;
	tie->tri_num++;
    }
    return;
}

template<typename Float>
int bucketpr_build( struct soltab *stp, struct rt_bot_internal *bot_ip, struct rt_i *rtip ) {
    

    struct tie_s *tie;
    struct bot_specific *bot;
    size_t tri_index, i;
    TIE_3 *tribuf = NULL, **tribufp = NULL;

    tree_t kd_tree;

    RT_BOT_CK_MAGIC(bot_ip);

    BU_GET(bot, struct bot_specific);
    stp->st_specific = (void *)bot;
    bot->bot_mode = bot_ip->mode;
    bot->bot_orientation = bot_ip->orientation;
    bot->bot_flags = bot_ip->bot_flags;
    if (bot_ip->thickness) {
	bot->bot_thickness = (fastf_t *)bu_calloc(bot_ip->num_faces, sizeof(fastf_t), "bot_thickness");
	for (tri_index = 0; tri_index < bot_ip->num_faces; tri_index++)
	    bot->bot_thickness[tri_index] = bot_ip->thickness[tri_index];
    } else {
	bot->bot_thickness = NULL;
    }

    if (bot_ip->face_mode) {
	bot->bot_facemode = bu_bitv_dup(bot_ip->face_mode);
    } else {
	bot->bot_facemode = BU_BITV_NULL;
    }
    bot->bot_facelist = NULL;

    tie = (struct tie_s *)bottie_allocn_double(bot_ip->num_faces);
    if (tie != NULL) {
	bot_ip->tie = bot->tie = tie;
    } else {
	return -1;
    }

    if ((tribuf = (TIE_3 *)bu_malloc(sizeof(TIE_3) * 3 * bot_ip->num_faces, "triangle tribuffer")) == NULL) {
	tie_free_double(tie);
	return -1;
    }
    if ((tribufp = (TIE_3 **)bu_malloc(sizeof(TIE_3*) * 3 * bot_ip->num_faces, "triangle tribuffer pointer")) == NULL) {
	tie_free_double(tie);
	bu_free(tribuf, "tribuf");
	return -1;
    }

    for (i = 0; i < bot_ip->num_faces*3; i++) {
	tribufp[i] = &tribuf[i];
	VMOVE(tribuf[i].v, (bot_ip->vertices+3*bot_ip->faces[i]));
    }

    /* tie_pushX sig: (struct tie_s *,
     *                 TIE_3 **,
     *                 unsigned int,
     *                 void *,
     *                 unsigned int);
     */
    //prbucket_push_double((struct tie_s *)bot_ip->tie, tribufp, bot_ip->num_faces, bot, 0);
    tie_push_double((struct tie_s *)bot_ip->tie, tribufp, bot_ip->num_faces, bot, 0);

    bu_free(tribuf, "tribuffer");
    bu_free(tribufp, "tribufp");

    tie_prep_double((struct tie_s *)bot->tie);

    VMOVE(stp->st_min, tie->amin);
    VMOVE(stp->st_max, tie->amax);

    /* zero thickness will get missed by the raytracer */
    BBOX_NONDEGEN(stp->st_min, stp->st_max, rtip->rti_tol.dist);

    VMOVE(stp->st_center, tie->mid);
    stp->st_aradius = tie->radius;
    stp->st_bradius = tie->radius;

    return 0;
}