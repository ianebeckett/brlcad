#pragma once

#include "rt/tie.h"

extern "C" {
  #include "nanort_interface.h"
  int rt_bot_shot(struct soltab *, struct xray *, struct application *, struct seg *);
}

template< typename Float >
int nanort_build( struct soltab *stp, struct rt_bot_internal *bot_ip, struct rt_i *rtip );

template< typename Float >
int  nanort_shot(struct soltab *stp, struct xray *rp, struct application *ap, struct seg *seghead);
