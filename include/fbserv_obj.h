/*                    F B S E R V _ O B J . H
 * BRL-CAD
 *
 * Copyright (c) 2004-2007 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @addtogroup libfb */
/*@{*/
/** @file fbserv_obj.h
 *
 */
#ifndef __FBSERV_OBJ_H__
#define __FBSERV_OBJ_H__

#include "fb.h"
#include "pkg.h"

#define NET_LONG_LEN	4	/**< @brief #  bytes to network long */
#define MAX_CLIENTS 32
#define MAX_PORT_TRIES 100
#define FBS_CALLBACK_NULL	(void (*)())NULL

struct fbserv_listener {
  int			fbsl_fd;			/**< @brief socket to listen for connections */
  int			fbsl_port;			/**< @brief port number to listen on */
  int			fbsl_listen;			/**< @brief !0 means listen for connections */
  struct fbserv_obj	*fbsl_fbsp;			/**< @brief points to its fbserv object */
};

struct fbserv_client {
  int			fbsc_fd;
  struct pkg_conn	*fbsc_pkg;
  struct fbserv_obj	*fbsc_fbsp;			/**< @brief points to its fbserv object */
};

struct fbserv_obj {
  FBIO				*fbs_fbp;			/**< @brief framebuffer pointer */
  struct fbserv_listener	fbs_listener;			/**< @brief data for listening */
  struct fbserv_client		fbs_clients[MAX_CLIENTS];	/**< @brief connected clients */
  void				(*fbs_callback)();		/**< @brief callback function */
  genptr_t			fbs_clientData;
};

FB_EXPORT extern int fbs_open();
FB_EXPORT extern int fbs_close();

#endif  /* __FBSERV_OBJ_H__ */
/*@}*/
/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * c-basic-offset: 4
 * indent-tabs-mode: t
 * End:
 * ex: shiftwidth=4 tabstop=8
 */

