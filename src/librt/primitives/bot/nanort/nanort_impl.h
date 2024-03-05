#ifndef LIBRT_PRIMITIVE_BOT_NANORT_IMPL_H
#define LIBRT_PRIMITIVE_BOT_NANORT_IMPL_H

/**
 * @file librt/primitives/bot/nanort/nanort_impl.h
 *
 * Interface header between NanoRT and LibRT. Meant to be included by both consumers of NanoRT
 * (librt - c sources and headers) and by the C++ implementation
 */


#include "rt/tie.h"

#ifdef __cplusplus
extern "C" {
#endif
  int nanort_build_double( struct soltab *stp, struct rt_bot_internal *bot, struct rt_i *rtip );
  void nanort_push_double(void *vtie, TIE_3 **tri, unsigned int ntri, void *usr, unsigned int pstride);
  int  nanort_prep_double(struct soltab *stp, struct rt_bot_internal *bot, struct rt_i *rtip);
  int  nanort_shot_double(struct soltab *stp, struct xray *rp, struct application *ap, struct seg *seghead);
  void nanort_free_double(void *vtie);

  void nanort_push_float(void *vtie, float **tri, unsigned int ntri, void *usr, unsigned int pstride);
  int  nanort_prep_float(struct soltab *stp, struct rt_bot_internal *bot, struct rt_i *rtip);
  int  nanort_shot_float(struct soltab *stp, struct xray *rp, struct application *ap, struct seg *seghead);
  void nanort_free_float(void *vtie);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
// Forward declarations of rt_bot_* functions for implementation
extern "C" {
  int rt_bot_shot(struct soltab *, struct xray *, struct application *, struct seg *);
  int rt_bot_makesegs(struct hit *hits, size_t nhits, struct soltab *stp, struct xray *rp, struct application *ap, struct seg *seghead, struct rt_piecestate *psp);
}

template< typename Float >
int nanort_build( struct soltab *stp, struct rt_bot_internal *bot_ip, struct rt_i *rtip );

template< typename Float >
int  nanort_shot(struct soltab *stp, struct xray *rp, struct application *ap, struct seg *seghead);
#endif

#endif
