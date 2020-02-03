/*                        C D T . H
 * BRL-CAD
 *
 * Copyright (c) 2007-2019 United States Government as represented by
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
/** @addtogroup libbrep */
/** @{ */
/** @file cdt.h
 *
 * Constrained Delaunay Triangulation of NURBS B-Rep objects.
 *
 */


/*
Design thoughts:

When performing triangle insertions from the CDT, can use the ordered
edges and unordered edges to valid the triangles as they go in, rather
than trying to do it after the fact:

1.  Any triangle with two vertices on a singularity is degenerate - skip.

2.  Before inserting a triangle into the mesh, make sure all three triangle
edges would reference unique unordered edges in 3D.  We've seen at least one
instance of a 2D triangle that ended up using the same interior edge for 2 of
its edges - that's a degenerate triangle, just filter it out.

3.  For a proposed insertion triangle, check it's three edges against the
unordered edge tree.   RTree uedge lookup will be used for this - in principle
it should be quicker than a std::map lookup by taking advantage of the spatial
grouping of the mesh.  There should be at least one pre existing edge (from the
mesh boundary if not from an already inserted triangle) that is already in the
containers and will be linked to the new triangle.  Check that unordered edge's
associated ordered edge that doesn't have an associated triangle against the
new triangle edges -  if there's not a match, flip the triangle ordering and
try again.

What this will do is force all triangles in 3D to form a valid topology,
regardless of whether or not the original CDT triangles would have formed such
a topology if naively translated.  Then we can check the triangle normals
against the brep normals and locally address those problem areas, but without
unmated edges being one of the potential issues.  I.e. the mesh will stay
as "valid" as possible in BoT terms at all times, even during the initial build.


During the initial edge breakdown, split all curved edges to a 45 degree normal
tolerance at an absolute minimum.  Do the same for the surface.  Drive curve
refinement using the size of the surface patches overlapping the edges, and
vice versa - needs to be handled so we don't get infinite looping, but the idea
is that curves, surfaces, and otherwise unsplit curves and surfaces in the
vicinity of curves and surfaces needing splitting will be (locally) refined
consistently.  Should terminate with those curve/surface areas that need it
(and only those) refined per tolerance settings, or to minimal standard.
Hopefully that will avoid the need for the various special case rules currently
defined based on linear/non-linear pairings and such, which don't always do
what we want/need locally in any event.

*/


#include "common.h"

#include <map>
#include <queue>
#include <string>

#include "../cdt/RTree.h"
#include "vmath.h"
#include "bn/tol.h"
#include "brep/defines.h"
#include "brep/cdt2.h"

#define BREP_PLANAR_TOL 0.05

class mesh_edge_t;
class mesh_uedge_t;
class mesh_tri_t;
class mesh_t;
class poly_edge_t;
class poly_uedge_t;
class polygon_t;

typedef enum {
    B_VERT,
    B_EDGE,
    B_SURF
} mesh_point_type;

class mesh_point_t {
    public:
	mesh_point_type type;

	ON_3dPoint p;
	ON_3dVector n;

	// Minimum associated edge length.  If < 0, no associated edges
	double min_len();

	// Moves the point to the new value, sets the new normal, and updates
	// all impacted data (vert bboxes, edge bboxes and lengths, triangle
	// bboxes, RTrees, etc.) using edge connectivity to find everything
	// that needs updating.
	void update(ON_3dPoint &np, ON_3dVector &nn);

	// Topology information: associated object and edges
	struct brep_cdt *cdt;
	std::set<size_t> uedges;

	void bbox_update();
	ON_BoundingBox bb;

	size_t vect_ind;

	// Origin information - specifics about what vertex, edge or surface on the
	// mesh is connected to this particular point.  Which set of information is
	// valid for this particular point depends on its type.
	// B_VERT
	int vert_index;
	bool singular;
	// B_EDGE
	int trim_index;
	int edge_index;
	double t;
	// B_SURF
	int face_index;
	double u;
	double v;


	// Variables used to store the associated index when assembling BoT
	// meshes
	int bot_vert_ind = -1;
	int bot_norm_ind = -1;

};

typedef enum {
    B_UNSET,
    B_SINGULAR,
    B_BOUNDARY,
    B_INTERIOR
} edge_type_t;

class mesh_edge_t {
    public:
	mesh_edge_t() {
	    v[0] = v[1] = -1;
	    m = NULL;
	    uedge = NULL;
	    type = B_UNSET;
	}

	void clear();
	bool set(long i, long j);

	bool operator<(mesh_edge_t other) const
	{
	    bool c1 = (v[0] < other.v[0]);
	    bool c1e = (v[0] == other.v[0]);
	    bool c2 = (v[1] < other.v[1]);
	    return (c1 || (c1e && c2));
	}

	bool operator==(mesh_edge_t other) const
	{
	    bool c0 = (type == other.type);
	    if (!c0) return false;
	    bool c1 = (v[0] == other.v[0]);
	    bool c2 = (v[1] == other.v[1]);
	    return (c1 && c2);
	}
	bool operator!=(mesh_edge_t other) const
	{
	    bool c0 = (type == other.type);
	    if (!c0) return true;
	    bool c1 = (v[0] != other.v[0]);
	    bool c2 = (v[1] != other.v[1]);
	    return (c1 || c2);
	}

	ON_BoundingBox &bbox();

	int vect_ind = -1;

	long v[2];
	edge_type_t type;
	mesh_t *m;
	mesh_uedge_t *uedge;

	mesh_tri_t *tri;

	double len = 0.0;
	bool unused = false;
	bool current = false;

    private:
	ON_BoundingBox bb;
};

class mesh_uedge_t {
    public:

	mesh_uedge_t() {
	    v[0] = v[1] = -1;
	    e[0] = e[1] = NULL;
	    pe[0] = pe[1] = NULL;
	    tri[0] = tri[1] = NULL;
	    type = B_UNSET;
	}

	void clear();
	bool set(mesh_edge_t *e1, mesh_edge_t *e2);

	bool operator<(mesh_uedge_t other) const
	{
	    bool c1 = (v[0] < other.v[0]);
	    bool c1e = (v[0] == other.v[0]);
	    bool c2 = (v[1] < other.v[1]);
	    return (c1 || (c1e && c2));
	}

	bool operator==(mesh_uedge_t other) const
	{
	    bool c0 = (type == other.type);
	    if (!c0) return false;
	    bool c1 = (v[0] == other.v[0]);
	    bool c2 = (v[1] == other.v[1]);
	    return (c1 && c2);
	}
	bool operator!=(mesh_uedge_t other) const
	{
	    bool c0 = (type == other.type);
	    if (!c0) return true;
	    bool c1 = (v[0] != other.v[0]);
	    bool c2 = (v[1] != other.v[1]);
	    return (c1 || c2);
	}

	bool split(ON_3dPoint &p);

	ON_BoundingBox *bbox();


	int vect_ind = -1;

	struct brep_cdt *cdt;
	edge_type_t type;
	mesh_edge_t *e[2];
	poly_edge_t *pe[2];
	mesh_tri_t *tri[2];
	long v[2];

	// If mesh_uedge_t is also a boundary edge, more information is needed
	int edge_ind = -1;
	double cp_len = DBL_MAX;
	ON_NurbsCurve *nc = NULL;
	double t_start = DBL_MAX;
	double t_end = DBL_MAX;

	int edge_start_pnt = -1;
	int edge_end_pnt = -1;
	ON_3dVector edge_tan_start = ON_3dVector::UnsetVector;
	ON_3dVector edge_tan_end = ON_3dVector::UnsetVector;

	void bbox_update();
	ON_BoundingBox bb;

	double len = 0.0;
	bool linear = false;
	bool unused = false;
	bool current = false;
};


class mesh_tri_t {
    public:

	mesh_tri_t()
	{
	    v[0] = v[1] = v[2] = -1;
	    ind = LONG_MAX;
	    m = NULL;
	    current = false;
	}

	mesh_tri_t(mesh_t *mesh, long i, long j, long k)
	{
	    v[0] = i;
	    v[1] = j;
	    v[2] = k;
	    ind = LONG_MAX;
	    m = mesh;
	    current = false;
	}

	mesh_tri_t(const mesh_tri_t &other)
	{
	    v[0] = other.v[0];
	    v[1] = other.v[1];
	    v[2] = other.v[2];
	    ind = other.ind;
	    m = other.m;
	    current = other.current;
	    bb = other.bb;
	}

	mesh_tri_t& operator=(const mesh_tri_t &other) = default;

	bool operator<(mesh_tri_t other) const
	{
	    std::vector<long> vca, voa;
	    for (int i = 0; i < 3; i++) {
		vca.push_back(v[i]);
		voa.push_back(other.v[i]);
	    }
	    std::sort(vca.begin(), vca.end());
	    std::sort(voa.begin(), voa.end());
	    bool c1 = (vca[0] < voa[0]);
	    bool c1e = (vca[0] == voa[0]);
	    bool c2 = (vca[1] < voa[1]);
	    bool c2e = (vca[1] == voa[1]);
	    bool c3 = (vca[2] < voa[2]);
	    return (c1 || (c1e && c2) || (c1e && c2e && c3));
	}
	bool operator==(mesh_tri_t other) const
	{
	    if (m != other.m) {
		return false;
	    }
	    bool c1 = ((v[0] == other.v[0]) || (v[0] == other.v[1]) || (v[0] == other.v[2]));
	    bool c2 = ((v[1] == other.v[0]) || (v[1] == other.v[1]) || (v[1] == other.v[2]));
	    bool c3 = ((v[2] == other.v[0]) || (v[2] == other.v[1]) || (v[2] == other.v[2]));
	    return (c1 && c2 && c3);
	}
	bool operator!=(mesh_tri_t other) const
	{
	    if (m != other.m) {
		return true;
	    }
	    bool c1 = ((v[0] != other.v[0]) && (v[0] != other.v[1]) && (v[0] != other.v[2]));
	    bool c2 = ((v[1] != other.v[0]) && (v[1] != other.v[1]) && (v[1] != other.v[2]));
	    bool c3 = ((v[2] != other.v[0]) && (v[2] != other.v[1]) && (v[2] != other.v[2]));
	    return (c1 || c2 || c3);
	}

	ON_3dPoint &center();
	ON_3dVector &bnorm();
	ON_3dVector &tnorm();

	ON_BoundingBox &bbox();

	long v[3];
	mesh_edge_t e[3];
	size_t ind;
	mesh_t *m;

	bool unused = false;
    private:
	bool current = false;
	ON_BoundingBox bb;
	ON_Plane tplane;
	ON_Plane bplane;
};

class poly_point_t {
    public:
	mesh_point_t *mp = NULL;
	double u;
	double v;
};

class polygon_t;

class poly_edge_t {
    public:

	poly_edge_t(polygon_t *p) {
	    v[0] = v[1] = -1;
	    polygon = p;
	}

	poly_edge_t(polygon_t *p, long i, long j) {
	    v[0] = i;
	    v[1] = j;
	    polygon = p;
	}

	void set(long i, long j) {
	    v[0] = i;
	    v[1] = j;
	}

	bool operator<(poly_edge_t other) const
	{
	    bool c1 = (v[0] < other.v[0]);
	    bool c1e = (v[0] == other.v[0]);
	    bool c2 = (v[1] < other.v[1]);
	    return (c1 || (c1e && c2));
	}

	bool operator==(poly_edge_t other) const
	{
	    bool c0 = (type == other.type);
	    if (!c0) return false;
	    bool c1 = (v[0] == other.v[0]);
	    bool c2 = (v[1] == other.v[1]);
	    return (c1 && c2);
	}
	bool operator!=(poly_edge_t other) const
	{
	    bool c0 = (type == other.type);
	    if (!c0) return true;
	    bool c1 = (v[0] != other.v[0]);
	    bool c2 = (v[1] != other.v[1]);
	    return (c1 || c2);
	}

	ON_BoundingBox *bbox();

	int vect_ind = -1;

	edge_type_t type;
	polygon_t *polygon;
	long v[2];
	mesh_edge_t *ue;
	double len = 0.0;
	bool current = false;
    private:
	ON_BoundingBox bb;
};

class poly_uedge_t {
    public:

	poly_uedge_t(polygon_t *p) {
	    v[0] = v[1] = -1;
	    polygon = p;
	    len = 0.0;
	    current = false;
	}

	poly_uedge_t(polygon_t *p, long i, long j) {
	    v[0] = i;
	    v[1] = j;
	    polygon = p;
	    len = 0.0;
	    current = false;
	}

	poly_uedge_t(poly_edge_t e) {
	    v[0] = (e.v[0] <= e.v[1]) ? e.v[0] : e.v[1];
	    v[1] = (e.v[0] > e.v[1]) ? e.v[0] : e.v[1];
	    len = e.len;
	    bb = *e.bbox();
	    polygon = e.polygon;
	    len = e.len;
	    current = e.current;
	}

	void set(long i, long j) {
	    v[0] = i;
	    v[1] = j;
	    current = false;
	}

	bool operator<(poly_uedge_t other) const
	{
	    bool c1 = (v[0] < other.v[0]);
	    bool c1e = (v[0] == other.v[0]);
	    bool c2 = (v[1] < other.v[1]);
	    return (c1 || (c1e && c2));
	}

	bool operator==(poly_uedge_t other) const
	{
	    bool c0 = (type == other.type);
	    if (!c0) return false;
	    bool c1 = (v[0] == other.v[0]);
	    bool c2 = (v[1] == other.v[1]);
	    return (c1 && c2);
	}
	bool operator!=(poly_uedge_t other) const
	{
	    bool c0 = (type == other.type);
	    if (!c0) return true;
	    bool c1 = (v[0] != other.v[0]);
	    bool c2 = (v[1] != other.v[1]);
	    return (c1 || c2);
	}

	ON_BoundingBox *bbox();

	int vect_ind = -1;

	double len;
	edge_type_t type;
	polygon_t *polygon;
	long v[2];
    private:
	bool current;
	ON_BoundingBox bb;
};

class polygon_t {
    std::vector<poly_point_t> p_pnts;
    std::map<long, long> o2p;
    RTree<size_t, double, 2> p_edges_tree; // 2D spatial lookup for polygon edges
};

class mesh_t
{
    public:
	// Primary container for points unique to this face
	std::vector<mesh_point_t> m_pnts;

	// Primary containers for face edges
	std::vector<mesh_edge_t> m_edges_vect;
	std::vector<mesh_uedge_t> m_uedges_vect;
	std::vector<poly_edge_t> m_pedges_vect;

	// Primary triangle container
	std::vector<mesh_tri_t> tris_vect;
	RTree<size_t, double, 3> tris_tree;

	// Polygonal approximation of face trimming loops
	std::vector<polygon_t> loops;

	// Find mesh_uedge_t in m_uedges_vect associated with v1,v2.  Returns
	// -1 if no such edge exists.
	long find_uedge(long v1, long v2);

	// Find mesh_edge_t associated with v1,v2 (looks up uedge,
	// then checks uedge's edges for the mesh_edge_t with the
	// correct ordering.  Returns -1 if no such edge exists
	long find_edge(long v1, long v2);

	// Create a new 3D edge
	mesh_edge_t &new_edge();
	// Create a new 3D triangle
	mesh_tri_t &new_tri();

	// Identify and return the unordered edge associated with
	// the ordered edge, or create such an unordered edge if
	// one does not already exist.
	mesh_uedge_t &uedge(mesh_edge_t &e);

	void delete_edge(mesh_edge_t &e);
	void delete_uedge(mesh_uedge_t &ue);
	void delete_tri(mesh_tri_t &t);

	struct brep_cdt *cdt;
	int f_id;
    private:

	std::queue<size_t> m_equeue; // Available (unused) entries in m_edges_vect
	std::queue<size_t> m_pequeue; // Available (unused) entries in m_pdges_vect
	std::queue<size_t> m_tqueue; // Available (unused) entries in tris_vect
};

class brep_cdt_state {
    public:
	// Primary container for points
	std::vector<mesh_point_t> b_pnts;
	// Spatial lookup tree for points
	RTree<size_t, double, 3> b_pnts_tree;

	// All brep unordered edges
	std::vector<mesh_uedge_t> b_uedges_vect;
	// Spatial lookup for unordered edges
	RTree<size_t, double, 3> b_uedges_tree;

	// Spatial lookup for face bboxes
	RTree<int, double, 3> b_faces_tree;

	mesh_uedge_t &new_uedge();
	void delete_uedge(mesh_uedge_t &ue);

	ON_Brep *orig_brep = NULL;
	ON_Brep *brep = NULL;

	std::string name;
	struct brep_cdt *cdt;

	void brep_cpy_init();
	void verts_init();
	void uedges_init();
	void edges_init();

    private:
	std::queue<size_t> b_equeue; // Available (unused) entries in edges_vect
	std::queue<size_t> b_uequeue; // Available (unused) entries in uedges_vect
};

struct brep_cdt_impl {
    brep_cdt_state s;
};


// Vertex processing routines
ON_3dVector singular_vert_norm(ON_Brep *brep, int index);

/** @} */

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8

