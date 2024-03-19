#include "bvh_impl.h"

#include <bits/chrono.h>
#include <chrono>
#include <initializer_list>
#include <iomanip>
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
#include "src/bbox.h"
#include "src/executor.h"
#include "src/thread_pool.h"
#include "src/vec.h"
#include "vmath.h"

#include "src/bvh.h"
#include "src/default_builder.h"
#include "src/tri.h"
#include "src/stack.h"
#include "../../../librt_private.h"

#include "bu/parallel.h"

#define PERMUTE_PRIMS true

namespace BVH = bvh::v2;

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

template< typename Float >
struct Accel {
  using Scalar  = Float;
  using Vec3    = bvh::v2::Vec<Scalar, 3>;
  using BBox    = bvh::v2::BBox<Scalar, 3>;
  using Tri     = bvh::v2::Tri<Scalar, 3>;
  using Node    = bvh::v2::Node<Scalar, 3>;
  using Bvh     = bvh::v2::Bvh<Node>;
  using Ray     = bvh::v2::Ray<Scalar, 3>;
  using PtrTri  = bvh::v2::PtrTri<Scalar, 3>;
  using PrecomputedTri = bvh::v2::PrecomputedTri<Scalar>;

  Bvh bvh;
  std::vector<PrecomputedTri> tris;
  std::vector<Tri> ptrTris;
};

/**
 * @brief Implementation of NanoRT BVH Building
 *
 * Based on (read: copied from) bottie_prep and tie_prep
 */
template< typename Float >
int bvh_build( struct soltab *stp, struct rt_bot_internal *bot_ip, struct rt_i *rtip ) {
  using Scalar  = Float;
  using Vec3    = bvh::v2::Vec<Scalar, 3>;
  using BBox    = bvh::v2::BBox<Scalar, 3>;
  using Tri     = bvh::v2::Tri<Scalar, 3>;
  using Node    = bvh::v2::Node<Scalar, 3>;
  using Bvh     = bvh::v2::Bvh<Node>;
  using Ray     = bvh::v2::Ray<Scalar, 3>;
  using PtrTri  = BVH::PtrTri<Scalar, 3>;

  struct tie_s *tie;
  struct bot_specific *bot;
  size_t tri_index, i;
  TIE_3 *tribuf = NULL, **tribufp = NULL;

  RT_BOT_CK_MAGIC(bot_ip);

  BU_GET(bot, struct bot_specific);
  stp->st_specific = (void *)bot;
  bot->bot_mode = bot_ip->mode;
  bot->bot_orientation = bot_ip->orientation;
  bot->bot_flags = bot_ip->bot_flags;
  bot_ip->nanort = bot->nanort = nullptr;
  bot_ip->tie = NULL;
  printf("Mode: %d\nOrientation: %d\n", bot_ip->mode, bot_ip->orientation );
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

  bot->nanort = bot_ip->nanort = new Accel<Float>();
  Accel<Float> & accel = *(Accel<Float>*)bot_ip->nanort;

  // nanort::BVHBuildOptions<Float> build_options;
  // build_options.min_primitives_for_parallel_build = 0;
  // nanort::BVHAccel<Float> * accel = new nanort::BVHAccel<Float>();
  typename BVH::DefaultBuilder<Node>::Config config;
  config.quality = BVH::DefaultBuilder<Node>::Quality::Low;


  BVH::ThreadPool threadpool( 1 );
  // BVH::ThreadPool threadpool( bu_avail_cpus() );
  BVH::ParallelExecutor executor( threadpool );

  BBox model_bbox = BBox::make_empty();
  std::vector<BBox> bboxes( bot_ip->num_faces );
  std::vector<Vec3> centers( bot_ip->num_faces );
  std::vector<PtrTri> tris( bot_ip->num_faces );

  std::cout << "Made threadpool and executor!" << std::endl;

  auto vecString = [](Float * vec) -> std::string {
    std::stringstream ss;
    ss << std::setprecision( 3 );
    ss << vec[0] << " " << vec[1] << " " << vec[2];
    return ss.str();
  };

  accel.ptrTris.resize( bot_ip->num_faces );
  accel.tris.resize( bot_ip->num_faces );

  for( unsigned i = 0; i < bot_ip->num_faces; ++i ) {
    auto face0 = bot_ip->faces[i + 0];
    auto face1 = bot_ip->faces[i + 1];
    auto face2 = bot_ip->faces[i + 2];
    Float * v0 = bot_ip->vertices + face0;
    Float * v1 = bot_ip->vertices + face1;
    Float * v2 = bot_ip->vertices + face2;
    PtrTri tri( v0, v1, v2 );

    std::cerr << "Triangle ( face IDX " << i << " / " << bot_ip->num_faces << " )" << std::endl;
    std::cerr << "\tV0: " << vecString( v0 ) << std::endl;
    std::cerr << "\tV1: " << vecString( v1 ) << std::endl;
    std::cerr << "\tV2: " << vecString( v2 ) << std::endl;

    centers[i] = tri.get_center();
    bboxes[i] = tri.get_bbox();
    model_bbox.extend( bboxes[i] );
    accel.ptrTris[i] = std::move( tri );
    accel.tris[i] = accel.ptrTris[i];
    accel.tris[i].prim_id = i;
  }
  //executor.for_each( 0, bot_ip->num_faces, [&bot_ip, &centers, &bboxes, &tris, &model_bbox, &accel, vecString](size_t begin, size_t end) -> void {
  //  for( unsigned i = begin; i < end; ++i ) {
  //    auto face = bot_ip->faces[i];
  //    auto tri_idx = face * 3 * 3;
  //    Float * v0 = bot_ip->vertices + tri_idx;
  //    Float * v1 = bot_ip->vertices + tri_idx + 3;
  //    Float * v2 = bot_ip->vertices + tri_idx + 6;
  //    PtrTri tri( v0, v1, v2 );

  //    std::cerr << "Triangle " << face << ": ( face IDX " << i << " / " << bot_ip->num_faces << " )" << std::endl;
  //    std::cerr << "\tV0: " << vecString( v0 ) << std::endl;
  //    std::cerr << "\tV1: " << vecString( v1 ) << std::endl;
  //    std::cerr << "\tV2: " << vecString( v2 ) << std::endl;

  //    centers[i] = tri.get_center();
  //    bboxes[i] = tri.get_bbox();
  //    model_bbox.extend( bboxes[i] );
  //    accel.ptrTris[i] = std::move( tri );
  //  }
  //});


  std::cout << "Bounding boxes and centers calculated" << std::endl;

  accel.bvh = BVH::DefaultBuilder<Node>::build( threadpool, bboxes, centers, config );

  std::cerr << "BVH Built! Whoop!" << std::endl;
  std::cerr << "Num faces: " << bot_ip->num_faces << std::endl;

  // Precompute Triangles...
  // accel.tris.resize( bot_ip->num_faces );
  // executor.for_each( 0, bot_ip->num_faces, [&accel, &bot_ip, &tris]( size_t begin, size_t end ) {
  //   for( int i = begin; i < end; i++ ) {
  //     int j;
  //     if constexpr ( PERMUTE_PRIMS ) {
  //       j = bot_ip->faces[i];
  //     }
  //     else {
  //       j = i;
  //     }

  //     accel.tris[i] = accel.ptrTris[i];
  //     accel.tris[i].prim_id = j;
  //   }
  // });

  std::cerr << "Triangles Computed!" << std::endl;


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
  VMOVE( stp->st_min, model_bbox.min.values );
  VMOVE( stp->st_max, model_bbox.max.values );

  // VMOVE(stp->st_min, tie->amin);
  // VMOVE(stp->st_max, tie->amax);

  /* zero thickness will get missed by the raytracer */
  BBOX_NONDEGEN(stp->st_min, stp->st_max, rtip->rti_tol.dist);

  // VMOVE(stp->st_center, tie->mid);
  // Calculate center of bounding box
  VMOVE( stp->st_center, model_bbox.get_center().values );

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

  // ADDED BY CAPSTONE TEAM
  // bvh::v2::Vec<fastf_t, 3, false> r_dir( ray->dir ), r_org( ray->pos );
  // auto pos = r_dir * id->dist + r_org;
  // VMOVE( hp->hit_point, pos.values );
  VJOIN1( hp->hit_point, ray->pos, id->dist, ray->dir );
  VMOVE( hp->hit_normal, id->norm );

  /* hit_vpriv is used later to clean up odd hits, exit before entrance, or
   * dangling entrance in bot_makesegs_(). When TIE was leaving this
   * unset, BOT hits were disappearing from the segment depending on the
   * random hit_vpriv[X] uninitialized values - set it using id->norm and the
   * ray direction. */
  VMOVE(tsp->tri_N, id->norm);
  hp->hit_vpriv[X] = VDOT(tsp->tri_N, ray->dir);
  tsp->tri_normals = tsp->tri_N;

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
  static int shot_number = 0;
  static int hit_number = 0;
  shot_number++;

  using Scalar  = Float;
  using Vec3    = bvh::v2::Vec<Scalar, 3>;
  using BBox    = bvh::v2::BBox<Scalar, 3>;
  using Tri     = bvh::v2::Tri<Scalar, 3>;
  using Node    = bvh::v2::Node<Scalar, 3>;
  using Bvh     = bvh::v2::Bvh<Node>;
  using Ray     = bvh::v2::Ray<Scalar, 3, false>;
  using PtrTri  = BVH::PtrTri<Scalar, 3>;
  constexpr size_t invalid_prim_id = std::numeric_limits<size_t>::max();

  struct bot_specific * bot = (bot_specific*)stp->st_specific;
  // nanort::BVHAccel<Float> * accel = (nanort::BVHAccel<Float>*) bot->nanort;
  Accel<Float> & accel = *(Accel<Float>*)bot->nanort;

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

  struct hitdata_s hitdata;
  hitdata.nhits = 0;
  hitdata.rp = rp;

  bvh::v2::SmallStack<typename Bvh::Index, 128 > stack;
  Ray bvhray( ray.pos, ray.dir, rp->r_min, rp->r_max );
  Float hit_isect = std::numeric_limits<Float>::max();
  size_t hit_prim = -1;

  auto vecString = [](Float * vec) -> std::string {
    std::stringstream ss;
    ss << std::setprecision( 3 );
    ss << vec[0] << " " << vec[1] << " " << vec[2];
    return ss.str();
  };

  accel.bvh.template intersect<false, false>( bvhray, accel.bvh.get_root().index, stack, [&](size_t begin, size_t end) {
    bool hit = false;
    for( int i = begin; i < end; ++i ) {
      int j = i;
      if constexpr( PERMUTE_PRIMS ) {
        j = accel.bvh.prim_ids[i];
      }

      // HIT
      Float isect = 0;
      auto uv = accel.tris[j].intersect( bvhray, isect, -std::numeric_limits<Float>::epsilon() );
      if( uv ) {
        if( isect < hit_isect ) {
          hit_prim = accel.tris[j].prim_id;
          hit_isect = isect;
        }
        hit = true;
        // printf("Hit(%5d): %4d %.4f\n", shot_number, accel.tris[j].prim_id, isect );
        // auto & [u, v] = *uv;
        tie_id_s tie_id;
        tie_ray_s tie_ray;

        auto prim_id = accel.tris[j].prim_id;
        auto tri = accel.ptrTris[j];


        // auto face0 = ((unsigned int*)bot->bot_facelist)[i + 0];
        // auto face1 = ((unsigned int*)bot->bot_facelist)[i + 1];
        // auto face2 = ((unsigned int*)bot->bot_facelist)[i + 2];

        Float *A = tri.p0.values; // &((Float*)bot->bot_facearray)[ face0 ];
        Float *B = tri.p1.values; // &((Float*)bot->bot_facearray)[ face1 ];
        Float *C = tri.p2.values; // &((Float*)bot->bot_facearray)[ face2 ];
        // std::cerr << "Triangle " << j << ":" << std::endl;
        // std::cerr << "\tA: " << vecString( A ) << std::endl;
        // std::cerr << "\tB: " << vecString( B ) << std::endl;
        // std::cerr << "\tC: " << vecString( C ) << std::endl;

        Float AC[3], AB[3];
        VSUB2(AC, C, A);
        VSUB2(AB, B, A);
        VCROSS(tie_id.norm, AC, AB);
        VUNITIZE( tie_id.norm );
        tie_id.dist = isect;

        VMOVE(tie_ray.dir, bvhray.dir.values );
        VMOVE(tie_ray.pos, bvhray.org.values );

        struct tri_specific * ts = hitdata.ts + hitdata.nhits;

        ts->tri_normals = ts->tri_N;

        auto hfret = hitfunc( &tie_ray, &tie_id, nullptr, &hitdata );
        if( hfret != nullptr ) {
          throw std::runtime_error("Too many hits!!!");
        }
        break;
      }

    }
    return hit;
  });


  // if( hit_prim != -1 ) {
  //   printf("Found hit(%4d). T = %.4f  ID=%5d\n", shot_number, hit_isect, hit_prim );
  //   printf("Hit number %d\n", ++hit_number );
  // }

  if( hitdata.nhits > 0 ) {
    for (i = 0; i < hitdata.nhits; i++)
      hitdata.hits[i].hit_dist = hitdata.hits[i].hit_dist - dirlen;
    for (i = 1; i < hitdata.nhits; i++) {
        hitdata.hits[i].hit_surfno = 0;
    }

    std::sort( hitdata.hits, hitdata.hits + hitdata.nhits, []( struct hit const & h1, struct hit const & h2 ) {
        return h1.hit_dist < h2.hit_dist;
    });

    // printf("Closest Hit: %.3f : %.3f %.3f %.3f\n", hitdata.hits[0].hit_dist, hitdata.hits[0].hit_point[0]
    //                                                                        , hitdata.hits[0].hit_point[1]
    //                                                                        , hitdata.hits[0].hit_point[2] );

    return rt_bot_makesegs( hitdata.hits, hitdata.nhits, stp, rp, ap, seghead, nullptr );
  }

  return 0;

}
