#pragma once

#include "rt/tie.h"

extern "C" {
  #include "nanort_interface.h"
}

template< typename Float >
int nanort_build( struct soltab *stp, struct rt_bot_internal *bot_ip, struct rt_i *rtip );
