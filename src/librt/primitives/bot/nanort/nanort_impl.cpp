#include "nanort_impl.hpp"

#include <bits/chrono.h>
#include <chrono>
#include <iostream>
#include <stdexcept>
#include <type_traits>

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
  bot->bot_facelist = NULL;

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

  nanort::TriangleMesh<Float> triangle_mesh( tri_vert, tri_faces, sizeof(Float)*3 );
  nanort::TriangleSAHPred<Float> triangle_pred( tri_vert, tri_faces, sizeof(Float)*3 );


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
  // Allocate triangles and buffer space...

  // tie = (struct tie_s *)bottie_allocn_double(bot_ip->num_faces);
  // if (tie != NULL) {
  //   bot_ip->tie = bot->tie = tie;
  // } else {
  //   return -1;
  // }

  // if ((tribuf = (TIE_3 *)bu_malloc(sizeof(TIE_3) * 3 * bot_ip->num_faces, "triangle tribuffer")) == NULL) {
  //   tie_free_double(tie);
  //   return -1;
  // }
  // if ((tribufp = (TIE_3 **)bu_malloc(sizeof(TIE_3*) * 3 * bot_ip->num_faces, "triangle tribuffer pointer")) == NULL) {
  //   tie_free_double(tie);
  //   bu_free(tribuf, "tribuf");
  //   return -1;
  // }

  // Copy data around
  // for (i = 0; i < bot_ip->num_faces*3; i++) {
  //   tribufp[i] = &tribuf[i];
  //   VMOVE(tribuf[i].v, (bot_ip->vertices+3*bot_ip->faces[i]));
  // }


  /* tie_pushX sig: (struct tie_s *,
   *                 TIE_3 **,
   *                 unsigned int,
   *                 void *,
   *                 unsigned int);
   */
  // Add to KD tree
  // tie_push_double((struct tie_s *)bot_ip->tie, tribufp, bot_ip->num_faces, bot, 0);

  // bu_free(tribuf, "tribuffer");
  // bu_free(tribufp, "tribufp");

  // // Perform KD tree build
  // tie_prep_double((struct tie_s *)bot->tie);

  // Set the min and max bounding boxes.
  accel->BoundingBox( stp->st_min, stp->st_max );
  // VMOVE(stp->st_min, tie->amin);
  // VMOVE(stp->st_max, tie->amax);

  /* zero thickness will get missed by the raytracer */
  BBOX_NONDEGEN(stp->st_min, stp->st_max, rtip->rti_tol.dist);

  // VMOVE(stp->st_center, tie->mid);
  for( int i = 0; i < 3; i++ ) {
    stp->st_center[i] =  (stp->st_max[i] - stp->st_min[i]) / (Float)2.f;
  }
  // TODO: How to calculate this?
  stp->st_aradius = 2;
  stp->st_bradius = 2;

  return 0;
}

template< typename Float >
int  nanort_shot(struct soltab *stp, struct xray *rp, struct application *ap, struct seg *seghead) {
  struct bot_specific * bot = (bot_specific*)stp->st_specific;
  nanort::BVHAccel<Float> * accel = (nanort::BVHAccel<Float>*) bot->nanort;



  return -1;
}
