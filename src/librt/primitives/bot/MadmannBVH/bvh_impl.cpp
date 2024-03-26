#include "bvh_impl.h"

#include <bits/chrono.h>
#include <initializer_list>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "rt/primitives/bot.h"
#include "rt/rt_instance.h"
#include "rt/tie.h"
#include "rt/geom.h"
#include "rt/seg.h"

#include "bu/parallel.h"
#include "common.h"
#include "vmath.h"
#include "../../../librt_private.h"

#include "src/bvh.h"
#include "src/default_builder.h"
#include "src/tri.h"
#include "src/stack.h"
#include "src/bbox.h"
#include "src/executor.h"
#include "src/thread_pool.h"
#include "src/vec.h"


#define PERMUTE_PRIMS false

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

  bvh::v2::SmallStack<typename Bvh::Index, 128 > * stacks;

  Accel()
    : stacks( new std::decay_t<decltype(*stacks)>[ bu_avail_cpus() ] ) {}

  ~Accel() {
    delete[] stacks;
  }
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

  // BVH::ThreadPool threadpool( 1 );
  BVH::ThreadPool threadpool( bu_avail_cpus() );
  BVH::ParallelExecutor executor( threadpool );

  BBox model_bbox = BBox::make_empty();
  std::vector<BBox> bboxes( bot_ip->num_faces );
  std::vector<Vec3> centers( bot_ip->num_faces );

  auto vecString = [](Float * vec) -> std::string {
    std::stringstream ss;
    ss << std::setprecision( 3 );
    ss << vec[0] << " " << vec[1] << " " << vec[2];
    return ss.str();
  };

  accel.ptrTris.resize( bot_ip->num_faces );
  accel.tris.resize( bot_ip->num_faces );

  executor.for_each( 0, bot_ip->num_faces, [&](size_t start, size_t end) {
    for( unsigned i = start; i < end; ++i ) {
      auto face0 = 3 * bot_ip->faces[3*i + 0];
      auto face1 = 3 * bot_ip->faces[3*i + 1];
      auto face2 = 3 * bot_ip->faces[3*i + 2];
      Float * v0 = bot_ip->vertices + face0;
      Float * v1 = bot_ip->vertices + face1;
      Float * v2 = bot_ip->vertices + face2;
      PtrTri tri( v0, v1, v2 );

      centers[i] = tri.get_center();
      bboxes[i] = tri.get_bbox();
      model_bbox.extend( bboxes[i] );
      accel.ptrTris[i] = std::move( tri );
    }
  });

  typename BVH::DefaultBuilder<Node>::Config config;
  config.quality = BVH::DefaultBuilder<Node>::Quality::High;
  accel.bvh = BVH::DefaultBuilder<Node>::build( threadpool, bboxes, centers, config );

  // Precompute Triangles...
  executor.for_each( 0, bot_ip->num_faces, [&](size_t begin, size_t end) {
    for( int i = 0; i < bot_ip->num_faces; ++i ) {
      int j = PERMUTE_PRIMS ? accel.bvh.prim_ids[i] : i;
      accel.tris[i] = accel.ptrTris[j];
      accel.tris[i].prim_id = j;
    }
  });

  bu_log( "Built BVH with %d triangles\n", bot_ip->num_faces );


  // // Set the min and max bounding boxes.
  VMOVE( stp->st_min, model_bbox.min.values );
  VMOVE( stp->st_max, model_bbox.max.values );

  /* zero thickness will get missed by the raytracer */
  BBOX_NONDEGEN(stp->st_min, stp->st_max, rtip->rti_tol.dist);

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

#ifndef MAXHITS
#define MAXHITS 128
#endif

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

  auto ID = bu_parallel_id() - 1;
  constexpr size_t invalid_prim_id = std::numeric_limits<size_t>::max();

  using Scalar  = Float;
  using Vec3    = bvh::v2::Vec<Scalar, 3>;
  using BBox    = bvh::v2::BBox<Scalar, 3>;
  using Tri     = bvh::v2::Tri<Scalar, 3>;
  using Node    = bvh::v2::Node<Scalar, 3>;
  using Bvh     = bvh::v2::Bvh<Node>;
  using Ray     = bvh::v2::Ray<Scalar, 3, false>;
  using PtrTri  = BVH::PtrTri<Scalar, 3>;

  struct bot_specific * bot = (bot_specific*)stp->st_specific;
  // nanort::BVHAccel<Float> * accel = (nanort::BVHAccel<Float>*) bot->nanort;
  Accel<Float> & accel = *(Accel<Float>*)bot->nanort;

  struct tie_s *tie;
  // struct hitdata_s hitdata;
  struct tie_id_s id;
  struct tie_ray_s ray;
  int i;
  fastf_t dirlen;
  struct hitdata_s hitdata;

  tie = (struct tie_s *)bot->tie;

  hitdata.nhits = 0;
  hitdata.rp = &ap->a_ray;
  /* do not need to init 'hits' and 'ts', tracked by 'nhits' */

  /* small backout applied to ray origin */
  dirlen = MAGSQ(rp->r_dir);
  VSUB2(ray.pos, rp->r_pt, rp->r_dir);	/* step back one dirlen */
  VMOVE(ray.dir, rp->r_dir);

  Ray bvhray( ray.pos, ray.dir, rp->r_min, rp->r_max );
  Float hit_isect = std::numeric_limits<Float>::max();
  size_t hit_prim = -1;

  auto vecString = [](Float * vec) -> std::string {
    std::stringstream ss;
    ss << std::setprecision( 3 );
    ss << vec[0] << " " << vec[1] << " " << vec[2];
    return ss.str();
  };

  accel.stacks[ ID ].size = 0;
  accel.bvh.template intersect<false, false>( bvhray, accel.bvh.get_root().index, accel.stacks[ ID ], [&](size_t begin, size_t end) {
    bool hit = false;
    for( int i = begin; i < end; ++i ) {
      int j = PERMUTE_PRIMS ? i : accel.bvh.prim_ids[i];

      // HIT
      Float isect = 0;
      auto uv = accel.tris[j].intersect( bvhray, bvhray.tmax, -std::numeric_limits<Float>::epsilon() );
      if( uv ) {
        hit = true;
        tie_id_s tie_id;
        tie_ray_s tie_ray;

        std::tie( tie_id.alpha, tie_id.beta ) = *uv;

        if( isect < hit_isect ) {
          hit_prim = j;
          hit_isect = bvhray.tmax;
        }

        auto const & tri = accel.ptrTris[j];

        Float const* A = tri.p0.values; // &((Float*)bot->bot_facearray)[ face0 ];
        Float const* B = tri.p1.values; // &((Float*)bot->bot_facearray)[ face1 ];
        Float const* C = tri.p2.values; // &((Float*)bot->bot_facearray)[ face2 ];

        Float AC[3], AB[3];
        VSUB2(AC, C, A);
        VSUB2(AB, B, A);
        VCROSS(tie_id.norm, AC, AB);
        VUNITIZE( tie_id.norm );
        tie_id.dist = bvhray.tmax;

        VMOVE(tie_ray.dir, bvhray.dir.values );
        VMOVE(tie_ray.pos, bvhray.org.values );

        struct tri_specific * ts = hitdata.ts + hitdata.nhits;

        ts->tri_normals = nullptr; // ts->tri_N;

        auto hfret = hitfunc( &tie_ray, &tie_id, nullptr, &hitdata );
        if( hfret != nullptr ) {
          throw std::runtime_error("Too many hits!!!");
        }
      }

    }
    return hit;
  });

  if( hitdata.nhits > 0 ) {
    for (i = 1; i < hitdata.nhits; i++) {
      hitdata.hits[i].hit_dist = hitdata.hits[i].hit_dist - dirlen;
      hitdata.hits[i].hit_surfno = 0;
    }

    std::sort( hitdata.hits, hitdata.hits + hitdata.nhits, []( struct hit const & h1, struct hit const & h2 ) {
        return h1.hit_dist < h2.hit_dist;
    });


    return rt_bot_makesegs( hitdata.hits, hitdata.nhits, stp, rp, ap, seghead, nullptr );
  }

  return 0;

}

template<typename Float>
void bvh_free( struct bot_specific * bot ) {
  auto * accel = (Accel<Float>*)bot->nanort;
  bot->nanort = nullptr;

  delete accel;
}
