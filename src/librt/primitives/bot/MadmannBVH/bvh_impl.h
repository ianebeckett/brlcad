#ifndef LIBRT_PRIMITIVE_BOT_BVH_IMPL_H
#define LIBRT_PRIMITIVE_BOT_BVH_IMPL_H

/**
 * @file librt/primitives/bot/nanort/nanort_impl.h
 *
 * Interface header between NanoRT and LibRT. Meant to be included by both consumers of NanoRT
 * (librt - c sources and headers) and by the C++ implementation
 */

#include <rt/tie.h>

#ifdef __cplusplus
// Forward declarations of rt_bot_* functions for implementation
extern "C" {
  int rt_bot_shot(struct soltab *, struct xray *, struct application *, struct seg *);
  int rt_bot_makesegs(struct hit *hits, size_t nhits, struct soltab *stp, struct xray *rp, struct application *ap, struct seg *seghead, struct rt_piecestate *psp);
}

template< typename Float >
int bvh_build( struct soltab *stp, struct rt_bot_internal *bot_ip, struct rt_i *rtip );

template< typename Float >
int bvh_shot(struct soltab *stp, struct xray *rp, struct application *ap, struct seg *seghead);

template< typename Float >
void bvh_free( struct bot_specific * bot );
#endif


// EXTERN 'C' DECLARATIONS FOR C INCLUDERS
//
#ifdef __cplusplus
extern "C" {
#endif
  int bvh_build_double( struct soltab *stp, struct rt_bot_internal *bot, struct rt_i *rtip );
  int bvh_shot_double(struct soltab *stp, struct xray *rp, struct application *ap, struct seg *seghead);

  void bvh_free_double( struct bot_specific * );
  void bvh_free_float( struct bot_specific * );

  int bvh_build_float( struct soltab *stp, struct rt_bot_internal *bot, struct rt_i *rtip );
  int bvh_shot_float(struct soltab *stp, struct xray *rp, struct application *ap, struct seg *seghead);


#ifdef __cplusplus
}
#endif

//
// EXTERN "C" DEFINITIONS FOR C++ IMPLEMENTATION
//
#ifdef __cplusplus
extern "C" {
   int bvh_build_double( struct soltab *stp, struct rt_bot_internal *bot, struct rt_i *rtip ) {
    return bvh_build<fastf_t>(stp,bot,rtip);
  }
   int bvh_shot_double(struct soltab *stp, struct xray *rp, struct application *ap, struct seg *seghead) {
    return bvh_shot<fastf_t>(stp,rp,ap,seghead);
  }
   int bvh_build_float( struct soltab *stp, struct rt_bot_internal *bot, struct rt_i *rtip ) {
    return bvh_build<fastf_t>(stp,bot,rtip);
  }
   int bvh_shot_float(struct soltab *stp, struct xray *rp, struct application *ap, struct seg *seghead) {
    return bvh_shot<fastf_t>(stp,rp,ap,seghead);
  }

  void bvh_free_double( struct bot_specific * bot ) {
    bvh_free<fastf_t>( bot );
  }
  void bvh_free_float( struct bot_specific * bot ) {
    bvh_free<fastf_t>( bot );
  }

}
#endif

#endif
