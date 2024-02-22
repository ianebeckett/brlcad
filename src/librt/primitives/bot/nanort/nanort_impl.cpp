#include "nanort_impl.h"

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

#include "nanort.h"
#include "../../../librt_private.h"

extern "C" {

  int nanort_build_double( struct soltab *stp, struct rt_bot_internal *bot_ip, struct rt_i *rtip ) {
    return nanort_build<double>(stp, bot_ip, rtip);
  }



  void nanort_push_double(void *vtie, TIE_3 **tri, unsigned int ntri, void *usr, unsigned int pstride);
  int  nanort_prep_double(struct soltab *stp, struct rt_bot_internal *bot, struct rt_i *rtip) {
    std::cerr << "NANORT IMPLEMENTATION ENABLED!" << std::endl;

    return -1;
  }
  int  nanort_shot_double(struct soltab *stp, struct xray *rp, struct application *ap, struct seg *seghead) {
    return nanort_shot< double >( stp, rp, ap, seghead );
  }
  void nanort_free_double(void *vtie);
}

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
int nanort_build( struct soltab *stp, struct rt_bot_internal *bot_ip, struct rt_i *rtip ) {
  struct tie_s *tie;
  struct bot_specific *bot;
  size_t tri_index, i;
  TIE_3 *tribuf = NULL, **tribufp = NULL;

  nanort::BVHBuildOptions<Float> build_options;
  nanort::BVHAccel<Float> * accel = new nanort::BVHAccel<Float>();

  RT_BOT_CK_MAGIC(bot_ip);

  BU_GET(bot, struct bot_specific);
  stp->st_specific = (void *)bot;
  bot->bot_mode = bot_ip->mode;
  bot->bot_orientation = bot_ip->orientation;
  bot->bot_flags = bot_ip->bot_flags;
  bot_ip->nanort = bot->nanort = accel;
  bot_ip->tie = NULL;
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

  std::cerr << "Building triangle mesh and pred..." << std::endl;
  static_assert( std::is_same_v< Float, std::decay_t<decltype(*bot_ip->vertices)>>, "Invalid floating point type!" );

  // 3 triangles per face; 3 floats per triangle
  Float * tri_vert = (Float*)bu_malloc( sizeof(Float) * 3 * 3 * bot_ip->num_faces, "NanoRT vertex vector" );
  unsigned int * tri_faces = (unsigned int*)bu_malloc( sizeof(unsigned int) * 3 * bot_ip->num_faces, "NanoRT faces vector" );

  // Copy data around
  for (i = 0; i < bot_ip->num_faces*3; i++) {
    // tribufp[i] = &tribuf[i];
    // VMOVE(tribuf[i].v, (bot_ip->vertices+3*bot_ip->faces[i]));
    VMOVE( &tri_vert[i*3], (bot_ip->vertices + 3*bot_ip->faces[i] ) );

    tri_faces[i] = bot_ip->faces[i];
  }

  std::cerr << "Memcpy done..." << std::endl;

  // nanort::TriangleMesh<Float> triangle_mesh( tri_vert, tri_faces, sizeof(Float)*3 );
  // nanort::TriangleSAHPred<Float> triangle_pred( tri_vert, tri_faces, sizeof(Float)*3 );
  nanort::TriangleMesh<Float> triangle_mesh( bot_ip->vertices, (const unsigned *) bot_ip->faces, sizeof(Float)*3 );
  nanort::TriangleSAHPred<Float> triangle_pred( bot_ip->vertices, (const unsigned *) bot_ip->faces, sizeof(Float)*3 );


  std::cerr << "Building accelerator..." << std::endl;
  std::cerr << bot_ip->num_faces * 3 << " triangles" << std::endl;
  auto start_time = std::chrono::system_clock::now().time_since_epoch();
  auto ret = accel->Build( bot_ip->num_faces, triangle_mesh, triangle_pred, build_options );
  auto end_time = std::chrono::system_clock::now().time_since_epoch() - start_time;
  std::cerr << "Accelerator built! Took " << std::chrono::duration_cast< std::chrono::milliseconds >( end_time ).count() << " ms " << std::endl;
  nanort::BVHBuildStatistics stats = accel->GetStatistics();

  printf("  BVH statistics:\n");
  printf("    # of leaf   nodes: %d\n", stats.num_leaf_nodes);
  printf("    # of branch nodes: %d\n", stats.num_branch_nodes);
  printf("  Max tree depth     : %d\n", stats.max_tree_depth);

  // Set the min and max bounding boxes.
  accel->BoundingBox( stp->st_min, stp->st_max );
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

static void *
hitfunc(struct tie_ray_s *ray, struct tie_id_s *id, struct tie_tri_s *UNUSED(tri), void *ptr) {

  return nullptr;
}

template< typename Float >
int  nanort_shot(struct soltab *stp, struct xray *rp, struct application *ap, struct seg *seghead) {
  struct bot_specific * bot = (bot_specific*)stp->st_specific;
  nanort::BVHAccel<Float> * accel = (nanort::BVHAccel<Float>*) bot->nanort;

  // printf("Performing NanoRT shot!\n");

  struct tie_s *tie;
  // struct hitdata_s hitdata;
  struct tie_id_s id;
  struct tie_ray_s ray;
  int i;
  fastf_t dirlen;

  bot = (struct bot_specific *)stp->st_specific;
  tie = (struct tie_s *)bot->tie;

  // hitdata.nhits = 0;
  // hitdata.rp = &ap->a_ray;
  /* do not need to init 'hits' and 'ts', tracked by 'nhits' */

  /* small backout applied to ray origin */
  dirlen = MAGSQ(rp->r_dir);
  VSUB2(ray.pos, rp->r_pt, rp->r_dir);	/* step back one dirlen */
  VMOVE(ray.dir, rp->r_dir);
  ray.depth = ray.kdtree_depth = 0;


  nanort::TriangleIntersector< Float, nanort::TriangleIntersection< Float > > intersector( (Float*)bot->bot_facearray, (const unsigned*)bot->bot_facelist, sizeof(Float)*3 );

  // Single-point intersection
  nanort::TriangleIntersection<Float> isect;
  nanort::Ray<Float> nrt_ray( ray.pos, ray.dir, 0, std::numeric_limits<Float>::max() );
  bool hit = accel->Traverse( nrt_ray, intersector, &isect );

  if( hit ) {
    printf("Ray hit @ t = %f\n", isect.t);
  }

  // bool hit = accel->Traverse(

  // tie_work_double(tie, &ray, &id, hitfunc, &hitdata);

  // /* use hitfunc to build the hit list */
  // if (hitdata.nhits == 0)
  //   return 0;

  // /* adjust hit distances to initial ray origin */
  // for (i = 0; i < hitdata.nhits; i++)
  //   hitdata.hits[i].hit_dist = hitdata.hits[i].hit_dist - dirlen;

  // /* FIXME: we don't have the hit_surfno but at least initialize it */
  // for (i = 0; i < hitdata.nhits; i++)
  //   hitdata.hits[i].hit_surfno = 0;

  // return rt_bot_makesegs(hitdata.hits, hitdata.nhits, stp, rp, ap, seghead, NULL);
  return 0;
  return rt_bot_makesegs(0, 0, stp, rp, ap, seghead, NULL);


  return -1;
}
