/*                      Q T C A D V I E W . C P P
 * BRL-CAD
 *
 * Copyright (c) 2021-2023 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file QtCADView.cpp
 *
 * Wrapper widget that handles the various widget types which may
 * constitute a Qt based geometry view.
 *
 */

#include "common.h"

#include "bg/polygon.h"
#include "bv.h"
#include "qtcad/QtCADView.h"
#include "qtcad/SignalFlags.h"

extern "C" {
#include "bu/malloc.h"
}


QtCADView::QtCADView(QWidget *parent, int type, struct fb *fbp)
    : QWidget(parent)
{
    this->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    l = new QBoxLayout(QBoxLayout::LeftToRight, this);
    l->setSpacing(0);
    l->setContentsMargins(0, 0, 0, 0);

    switch (type) {
#ifdef BRLCAD_OPENGL
	case QtCADView_GL:
	    canvas_gl = new QtGL(this, fbp);
	    canvas_gl->setMinimumSize(50,50);
	    canvas_gl->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
	    l->addWidget(canvas_gl);
	    QObject::connect(canvas_gl, &QtGL::changed, this, &QtCADView::do_view_changed);
	    QObject::connect(canvas_gl, &QtGL::init_done, this, &QtCADView::do_init_done);
	    break;
#endif
	case QtCADView_SW:
	    canvas_sw = new QtSW(this, fbp);
	    canvas_sw->setMinimumSize(50,50);
	    canvas_sw->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
	    l->addWidget(canvas_sw);
	    QObject::connect(canvas_sw, &QtSW::changed, this, &QtCADView::do_view_changed);
	    QObject::connect(canvas_sw, &QtSW::init_done, this, &QtCADView::do_init_done);
	    break;
	default:
#ifdef BRLCAD_OPENGL
	    canvas_gl = new QtGL(this, fbp);
	    canvas_gl->setMinimumSize(50,50);
	    canvas_gl->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
	    l->addWidget(canvas_gl);
	    QObject::connect(canvas_gl, &QtGL::changed, this, &QtCADView::do_view_changed);
	    QObject::connect(canvas_gl, &QtGL::init_done, this, &QtCADView::do_init_done);
#else
	    canvas_sw = new QtSW(this, fbp);
	    canvas_sw->setMinimumSize(50,50);
	    canvas_sw->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
	    l->addWidget(canvas_sw);
	    QObject::connect(canvas_sw, &QtSW::changed, this, &QtCADView::do_view_changed);
	    QObject::connect(canvas_sw, &QtSW::init_done, this, &QtCADView::do_init_done);
#endif
	    return;
    }
}

QtCADView::~QtCADView()
{
#ifdef BRLCAD_OPENGL
    if (canvas_gl)
	delete canvas_gl;
#endif
    if (canvas_sw)
	delete canvas_sw;
}

bool
QtCADView::isValid()
{
    if (canvas_sw)
	return true;

#ifdef BRLCAD_OPENGL
    if (canvas_gl)
	return canvas_gl->isValid();
#endif

    return false;
}

int
QtCADView::view_type()
{
#ifdef BRLCAD_OPENGL
    if (canvas_gl)
	return QtCADView_GL;
#endif
    if (canvas_sw)
	return QtCADView_SW;

    return -1;
}


void
QtCADView::save_image(int UNUSED(quad))
{
}

void
QtCADView::do_view_changed()
{
    QTCAD_SLOT("QtCADView::do_view_changed", 1);
    emit changed();
}

void
QtCADView::need_update(unsigned long long)
{
    bv_log(4, "QtCADView::need_update");
    QTCAD_SLOT("QtCADView::need_update", 1);
#ifdef BRLCAD_OPENGL
    if (canvas_gl) {
	canvas_gl->need_update();
	return;
    }
#endif
    if (canvas_sw) {
	canvas_sw->need_update();
	return;
    }
}

struct bview *
QtCADView::view()
{
#ifdef BRLCAD_OPENGL
    if (canvas_gl)
	return canvas_gl->v;
#endif
    if (canvas_sw)
	return canvas_sw->v;

    return NULL;
}

struct dm *
QtCADView::dmp()
{
#ifdef BRLCAD_OPENGL
    if (canvas_gl)
	return canvas_gl->dmp;
#endif
    if (canvas_sw)
	return canvas_sw->dmp;

    return NULL;
}

struct fb *
QtCADView::ifp()
{
#ifdef BRLCAD_OPENGL
    if (canvas_gl)
	return canvas_gl->ifp;
#endif
    if (canvas_sw)
	return canvas_sw->ifp;

    return NULL;
}

void
QtCADView::set_view(struct bview *nv)
{
#ifdef BRLCAD_OPENGL
    if (canvas_gl) {
	canvas_gl->v = nv;
	if (canvas_gl->dmp && canvas_gl->v) {
	    canvas_gl->v->dmp = canvas_gl->dmp;
	}
    }
#endif
    if (canvas_sw) {
	canvas_sw->v = nv;
    	if (canvas_sw->dmp && canvas_sw->v) {
	    canvas_sw->v->dmp = canvas_sw->dmp;
	}
    }
}

void
QtCADView::stash_hashes()
{
#ifdef BRLCAD_OPENGL
    if (canvas_gl) {
	canvas_gl->stash_hashes();
    }
#endif
    if (canvas_sw) {
	canvas_sw->stash_hashes();
    }
}

bool
QtCADView::diff_hashes()
{
#ifdef BRLCAD_OPENGL
    if (canvas_gl) {
	return canvas_gl->diff_hashes();
    }
#endif
    if (canvas_sw) {
	return canvas_sw->diff_hashes();
    }

    return false;
}

void
QtCADView::aet(double a, double e, double t)
{
#ifdef BRLCAD_OPENGL
    if (canvas_gl) {
	canvas_gl->aet(a, e, t);
    }
#endif
    if (canvas_sw) {
	canvas_sw->aet(a, e, t);
    }
}

void
QtCADView::set_current(int i)
{
#ifdef BRLCAD_OPENGL
    if (canvas_gl) {
	canvas_gl->current = i;
    }
#endif
    if (canvas_sw) {
	canvas_sw->current = i;
    }
}

int
QtCADView::current()
{
#ifdef BRLCAD_OPENGL
    if (canvas_gl) {
	return canvas_gl->current;
    }
#endif
    if (canvas_sw) {
	return canvas_sw->current;
    }

    return 0;
}

void
QtCADView::add_event_filter(QObject *o)
{
    curr_event_filter = o;
#ifdef BRLCAD_OPENGL
    if (canvas_gl) {
	canvas_gl->installEventFilter(o);
	return;
    }
#endif
    if (canvas_sw) {
	canvas_sw->installEventFilter(o);
	return;
    }
}

void
QtCADView::clear_event_filter(QObject *o)
{
    if (!o)
	return;
#ifdef BRLCAD_OPENGL
    if (canvas_gl) {
	canvas_gl->removeEventFilter(o);
    }
#endif
    if (canvas_sw) {
	canvas_sw->removeEventFilter(o);
    }
    curr_event_filter = NULL;
}

void
QtCADView::set_draw_custom(void (*draw_custom)(struct bview *, void *), void *draw_udata)
{

#ifdef BRLCAD_OPENGL
    if (canvas_gl) {
	canvas_gl->draw_custom = draw_custom;
	canvas_gl->draw_udata = draw_udata;
	return;
    }
#endif
    if (canvas_sw) {
	canvas_sw->draw_custom = draw_custom;
	canvas_sw->draw_udata = draw_udata;
	return;
    }
}

void
QtCADView::enableDefaultKeyBindings()
{
#ifdef BRLCAD_OPENGL
    if (canvas_gl) {
	canvas_gl->enableDefaultKeyBindings();
	return;
    }
#endif
    if (canvas_sw) {
	canvas_sw->enableDefaultKeyBindings();
	return;
    }
}

void
QtCADView::disableDefaultKeyBindings()
{
#ifdef BRLCAD_OPENGL
    if (canvas_gl) {
	canvas_gl->disableDefaultKeyBindings();
	return;
    }
#endif
    if (canvas_sw) {
	canvas_sw->disableDefaultKeyBindings();
	return;
    }
}

void
QtCADView::enableDefaultMouseBindings()
{
#ifdef BRLCAD_OPENGL
    if (canvas_gl) {
	canvas_gl->enableDefaultMouseBindings();
	return;
    }
#endif
    if (canvas_sw) {
	canvas_sw->enableDefaultMouseBindings();
	return;
    }


}

void
QtCADView::disableDefaultMouseBindings()
{
#ifdef BRLCAD_OPENGL
    if (canvas_gl) {
	canvas_gl->disableDefaultMouseBindings();
	return;
    }
#endif
    if (canvas_sw) {
	canvas_sw->disableDefaultMouseBindings();
	return;
    }
}


void
QtCADView::set_lmouse_move_default(int mm)
{
    QTCAD_SLOT("QtCADView::set_lmouse_move_default", 1);

#ifdef BRLCAD_OPENGL
    if (canvas_gl) {
	canvas_gl->set_lmouse_move_default(mm);
	return;
    }
#endif
    if (canvas_sw) {
	canvas_sw->set_lmouse_move_default(mm);
	return;
    }
}


void
QtCADView::do_init_done()
{
    QTCAD_SLOT("QtCADView::do_init_done", 1);
    emit init_done();
}

bool
QPolyFilter::eventFilter(QObject *, QEvent *e)
{
    if (!cv)
	return false;

    struct bview *v = cv->view();

    QMouseEvent *m_e = NULL;

    if (e->type() == QEvent::MouseButtonPress || e->type() == QEvent::MouseButtonRelease || e->type() == QEvent::MouseButtonDblClick || e->type() == QEvent::MouseMove) {

	m_e = (QMouseEvent *)e;

	v->gv_prevMouseX = v->gv_mouse_x;
	v->gv_prevMouseY = v->gv_mouse_y;
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
	v->gv_mouse_x = m_e->x();
	v->gv_mouse_y = m_e->y();
#else
	v->gv_mouse_x = m_e->position().x();
	v->gv_mouse_y = m_e->position().y();
#endif
    }

    if (!m_e)
	return false;

    // If we have modifiers, we're most likely doing shift grips
    if (m_e->modifiers() != Qt::NoModifier)
	return false;

    if (m_e->type() == QEvent::MouseButtonPress && m_e->buttons().testFlag(Qt::LeftButton)) {
	if (!p) {
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
	    p = bv_create_polygon(v, BV_VIEW_OBJS, ptype, m_e->x(), m_e->y());
#else
	    p = bv_create_polygon(v, BV_VIEW_OBJS, ptype, m_e->position().x(), m_e->position().y());
#endif
	    p->s_v = v;
	    struct bv_polygon *ip = (struct bv_polygon *)p->s_i_data;

	    // Get edge color
	    bu_color_to_rgb_chars(&edge_color, p->s_color);

	    // fill color
	    BU_COLOR_CPY(&ip->fill_color, &fill_color);

	    // fill settings
	    vect2d_t vdir = V2INIT_ZERO;
	    vdir[0] = fill_slope_x;
	    vdir[1] = fill_slope_y;
	    V2MOVE(ip->fill_dir, vdir);
	    ip->fill_delta = fill_density;

	    // Set fill
	    if (fill_poly) {
		ip->fill_flag = 1;
		bv_update_polygon(p, p->s_v, BV_POLYGON_UPDATE_PROPS_ONLY);
	    }

	    // Name appropriately
	    bu_vls_init(&p->s_uuid);

	    // It doesn't get a "proper" name until its finalized
	    bu_vls_printf(&p->s_uuid, "_tmp_view_polygon");

	    emit view_updated(QTCAD_VIEW_REFRESH);
	    return true;
	}

	// If we're creating a general polygon, we're appending points after
	// the initial creation
	struct bv_polygon *ip = (struct bv_polygon *)p->s_i_data;
	if (ip->type == BV_POLYGON_GENERAL) {
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
	    p->s_v->gv_mouse_x = m_e->x();
	    p->s_v->gv_mouse_y = m_e->y();
#else
	    p->s_v->gv_mouse_x = m_e->position().x();
	    p->s_v->gv_mouse_y = m_e->position().y();
#endif
	    bv_update_polygon(p, p->s_v, BV_POLYGON_UPDATE_PT_APPEND);

	    emit view_updated(QTCAD_VIEW_REFRESH);
	    return true;
	}

	// When we're dealing with polygons stray left clicks shouldn't zoom - just
	// consume them if we're not using them above.
	return true;
    }

    if (m_e->type() == QEvent::MouseButtonPress && m_e->buttons().testFlag(Qt::RightButton)) {
	// No-op if no current polygon is defined
	if (!p)
	    return true;

	// Non-general polygon creation doesn't use right click.
	struct bv_polygon *ip = (struct bv_polygon *)p->s_i_data;
	if (ip->type != BV_POLYGON_GENERAL) {
	    return true;
	}

	// General polygon, have right click - finish up.
	finalize(true);
	return true;
    }

    if (m_e->type() == QEvent::MouseButtonPress) {
	// We also don't want other stray mouse clicks to do something surprising
	return true;
    }

    // During initial add/creation of non-general polygons, mouse movement
    // adjusts the shape
    if (m_e->type() == QEvent::MouseMove) {
	// No-op if no current polygon is defined
	if (!p)
	    return true;

	// General polygon creation doesn't use mouse movement.
	struct bv_polygon *ip = (struct bv_polygon *)p->s_i_data;
	if (ip->type == BV_POLYGON_GENERAL) {
	    return true;
	}

	// For every other polygon type, call the libbv update routine
	// with the view's x,y coordinates
	if (m_e->buttons().testFlag(Qt::LeftButton) && m_e->modifiers() == Qt::NoModifier) {
	    bv_update_polygon(p, p->s_v, BV_POLYGON_UPDATE_DEFAULT);
	    emit view_updated(QTCAD_VIEW_REFRESH);
	    return true;
	}
    }

    if (m_e->type() == QEvent::MouseButtonRelease) {

	// No-op if no current polygon is defined
	if (!p)
	    return true;

	struct bv_polygon *ip = (struct bv_polygon *)p->s_i_data;
	if (ip->type == BV_POLYGON_GENERAL) {
	    // General polygons are finalized by an explicit close
	    // (either right mouse click or the close checkbox)
	    return true;
	}

	// For all non-general polygons, mouse release is the signal
	// to finish up.
	finalize(true);

	return true;
    }

    return false;
}

void
QPolyFilter::finalize(bool)
{
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
