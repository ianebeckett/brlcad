#include "bvh_impl.h"

#include <bits/chrono.h>
#include <chrono>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "rt/primitives/bot.h"
#include "rt/rt_instance.h"
#include "rt/tie.h"
#include "bio.h"
#include "rt/geom.h"
#include "rt/seg.h"
#include "common.h"
#include "vmath.h"

#include "src/bvh.h"
#include "../../../librt_private.h"

template< typename Float, typename ... Floats >
Float ilist_func_apply( Float(*func)(std::initializer_list<Float>), Float f, Floats ... floats ) {
  static_assert( sizeof...( Floats ) >= 1, "ilist_func_apply Can only be called with at least 2 arguments...");
  std::initializer_list<Float> ilist = { std::decay_t<Float>(f), std::forward< Float >( floats ) ... };
  return func( ilist );
}

/**
 * @brief Max taking a variable number of floating-point arguments
 */
template< typename F, typename ... Float >
F max( F f, Float ... floats ) {
  return ilist_func_apply( std::max, f, floats ... );
}

template< typename F, typename ... Float >
F min( F f, Float ... floats ) {
  return ilist_func_apply( std::min, f, floats ... );
}

/**
 * @brief Implementation of NanoRT BVH Building
 *
 * Based on (read: copied from) bottie_prep and tie_prep
 */
template< typename Float >
int bvh_build( struct soltab *stp, struct rt_bot_internal *bot_ip, struct rt_i *rtip ) {
  struct tie_s *tie;
  struct bot_specific *bot;
  size_t tri_index, i;
  TIE_3 *tribuf = NULL, **tribufp = NULL;

  // nanort::BVHBuildOptions<Float> build_options;
  // build_options.min_primitives_for_parallel_build = 0;
  // nanort::BVHAccel<Float> * accel = new nanort::BVHAccel<Float>();

  RT_BOT_CK_MAGIC(bot_ip);

  BU_GET(bot, struct bot_specific);
  stp->st_specific = (void *)bot;
  bot->bot_mode = bot_ip->mode;
  bot->bot_orientation = bot_ip->orientation;
  bot->bot_flags = bot_ip->bot_flags;
  bot_ip->nanort = bot->nanort = nullptr;
  bot_ip->tie = NULL;
  printf("Mode: %d\nOrientation: %d\n", bot_ip->mode, bot_ip->orientation );
  if( stp->st_meth->ft_shot != rt_bot_shot ) {
    throw std::runtime_error("NanoRT ft_shot is not rt_bot_shot!");
  }
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


  bot->bot_facelist = bot_ip->faces;
  bot->bot_facearray = (void**)bot_ip->vertices;
  bot->bot_ntri = bot_ip->num_faces * 3;


  // std::cerr << "Building triangle mesh and pred..." << std::endl;
  // nanort::TriangleMesh<Float> triangle_mesh( bot_ip->vertices, (const unsigned *) bot_ip->faces, sizeof(fastf_t)*3 );
  // nanort::TriangleSAHPred<Float> triangle_pred( bot_ip->vertices, (const unsigned *) bot_ip->faces, sizeof(fastf_t)*3 );


  // std::cerr << "Building accelerator..." << std::endl;
  // std::cerr << bot_ip->num_faces * 3 << " triangles" << std::endl;
  // auto start_time = std::chrono::system_clock::now().time_since_epoch();
  // auto ret = accel->Build( bot_ip->num_faces, triangle_mesh, triangle_pred, build_options );
  // auto end_time = std::chrono::system_clock::now().time_since_epoch() - start_time;
  // std::cerr << "Accelerator built! Took " << std::chrono::duration_cast< std::chrono::milliseconds >( end_time ).count() << " ms " << std::endl;
  // nanort::BVHBuildStatistics stats = accel->GetStatistics();

  // if( not ret ) throw std::runtime_error("Failed to build NanoRT BVH!");

  // printf("  BVH statistics:\n");
  // printf("    # of leaf   nodes: %d\n", stats.num_leaf_nodes);
  // printf("    # of branch nodes: %d\n", stats.num_branch_nodes);
  // printf("  Max tree depth     : %d\n", stats.max_tree_depth);

  // // Set the min and max bounding boxes.
  // accel->BoundingBox( stp->st_min, stp->st_max );
  // printf("  BVH Bounding Box Min: %.3f %.3f %.3f\n", stp->st_min[0], stp->st_min[1], stp->st_min[2]);
  // printf("  BVH Bounding Box Max: %.3f %.3f %.3f\n", stp->st_max[0], stp->st_max[1], stp->st_max[2]);

  // VMOVE(stp->st_min, tie->amin);
  // VMOVE(stp->st_max, tie->amax);

  /* zero thickness will get missed by the raytracer */
  BBOX_NONDEGEN(stp->st_min, stp->st_max, rtip->rti_tol.dist);

  // VMOVE(stp->st_center, tie->mid);
  // Calculate center of bounding box
  for( int i = 0; i < 3; i++ ) {
    stp->st_center[i] =  (stp->st_max[i] - stp->st_min[i]) / (Float)2.f;
  }
  // TODO: How to calculate this?
  // FIXME: Actually calculate. Currently just the max of the distances from the center...
  auto radius = max( std::abs( stp->st_max[0] - stp->st_center[0] ),
                                std::abs( stp->st_max[1] - stp->st_center[1] ),
                                std::abs( stp->st_max[2] - stp->st_center[2] ),
                                std::abs( stp->st_center[0] - stp->st_min[0] ),
                                std::abs( stp->st_center[1] - stp->st_min[1] ),
                                std::abs( stp->st_center[2] - stp->st_min[2] ) );
  stp->st_aradius = radius;
  stp->st_bradius = radius;

  return 0;
}

#define MAXHITS 128

struct hitdata_s {
    int nhits;
    struct hit hits[MAXHITS];
    struct tri_specific ts[MAXHITS];
    struct xray *rp;
};

static void *
hitfunc(struct tie_ray_s *ray, struct tie_id_s *id, struct tie_tri_s *UNUSED(tri), void *ptr)
{
  struct hitdata_s *h = (struct hitdata_s *)ptr;
  struct tri_specific *tsp;
  struct hit *hp;

  if (h->nhits > (MAXHITS-1)) {
    bu_log("Too many hits!\n");
    return (void *)1;
  }

  hp = &h->hits[h->nhits];
  hp->hit_private = &h->ts[h->nhits];
  tsp = (struct tri_specific *)hp->hit_private;
  h->nhits++;


  hp->hit_magic = RT_HIT_MAGIC;
  hp->hit_dist = id->dist;

  /* hit_vpriv is used later to clean up odd hits, exit before entrance, or
   * dangling entrance in bot_makesegs_(). When TIE was leaving this
   * unset, BOT hits were disappearing from the segment depending on the
   * random hit_vpriv[X] uninitialized values - set it using id->norm and the
   * ray direction. */
  VMOVE(tsp->tri_N, id->norm);
  hp->hit_vpriv[X] = VDOT(tsp->tri_N, ray->dir);

  /* Of the hit_vpriv assignments added in commit 50164, only hit_vpriv[X]
   * was based on initialized calculations.  Rather than leaving the other
   * assignments (which were based on calculations using uninitialized values
   * in the tsp struct) we simply assign 0 values.
   *
   * NOTE: bot_norm_ does use Y and Z in the normal calculations, so those
   * still won't work.  Likewise, bot_plate_segs_ uses hit_surfno which is
   * also not set correctly.  Perhaps the currently unused tri would have
   * the needed info? */
  hp->hit_vpriv[Y] = 0.0;
  hp->hit_vpriv[Z] = 0.0;
  hp->hit_surfno = 0;

  /* add hitdist into array and add one to nhits */
  return NULL;	/* continue firing */
}

template< typename Float >
int  bvh_shot(struct soltab *stp, struct xray *rp, struct application *ap, struct seg *seghead) {
  struct bot_specific * bot = (bot_specific*)stp->st_specific;
  // nanort::BVHAccel<Float> * accel = (nanort::BVHAccel<Float>*) bot->nanort;

  // printf("Performing NanoRT shot!\n");

  struct tie_s *tie;
  // struct hitdata_s hitdata;
  struct tie_id_s id;
  struct tie_ray_s ray;
  int i;
  fastf_t dirlen;


  tie = (struct tie_s *)bot->tie;

  // hitdata.nhits = 0;
  // hitdata.rp = &ap->a_ray;
  /* do not need to init 'hits' and 'ts', tracked by 'nhits' */

  /* small backout applied to ray origin */
  dirlen = MAGSQ(rp->r_dir);
  VSUB2(ray.pos, rp->r_pt, rp->r_dir);	/* step back one dirlen */
  VMOVE(ray.dir, rp->r_dir);
  ray.depth = ray.kdtree_depth = 0;


  // nanort::TriangleIntersector< Float, nanort::TriangleIntersection< Float > > intersector( (Float*)bot->bot_facearray, (const unsigned*)bot->bot_facelist, sizeof(Float)*3 );

  // // Single-point intersection
  // nanort::TriangleIntersection<Float> isect;
  // nanort::Ray<Float> nrt_ray( ray.pos, ray.dir, rp->r_min, rp->r_max );
  // // nanort::StackVector<nanort::NodeHit<Float>, 128> hits;
  // // nanort::BVHTraceOptions options;
  // // auto hit = accel->ListNodeIntersections( nrt_ray, 128, intersector, &hits );
  // bool hit = accel->Traverse( nrt_ray, intersector, &isect );

  // if( hit ) {
  //   struct hitdata_s hitdata;
  //   hitdata.rp = rp;
  //   struct tie_ray_s tie_ray;
  //   struct tie_id_s tie_id;

  //   Float *A = &((Float*)bot->bot_facearray)[ isect.prim_id + 0 ];
  //   Float *B = &((Float*)bot->bot_facearray)[ isect.prim_id + 1 ];
  //   Float *C = &((Float*)bot->bot_facearray)[ isect.prim_id + 2 ];

  //   Float AC[3], AB[3], NORM;
  //   VSUB2(AC, C, A);
  //   VSUB2(AB, B, A);
  //   VCROSS(tie_id.norm, AC, AB);
  //   VUNITIZE( tie_id.norm );
  //   tie_id.dist = isect.t;

  //   VMOVE(tie_ray.dir, rp->r_dir);
  //   VMOVE(tie_ray.pos, rp->r_pt );

  //   hitfunc( &tie_ray, &tie_id, nullptr, &hitdata );
  //   return rt_bot_makesegs( hitdata.hits, hitdata.nhits, stp, rp, ap, seghead, NULL );

  //   // printf("Ray hit some nodes!\n");
  //   // printf("Ray hit tri%d @ t = %f\n", isect.prim_id, isect.t);
  // }

  return -1;
}
