#pragma once

#include "rt/tie.h"

extern "C" {
    #include "nanort_interface.h"
    int rt_bot_shot(struct soltab *, struct xray *, struct application *, struct seg *);
}

template<typename Float>
int bucketpr_build(struct soltab *, struct rt_bot_internal *, struct rt_i *);