
#include "bg/plane.h"
#include "common.h"

#include "rt/hit.h"
#include "rt/tie.h"

#include "bu/assert.h"
#include "bu/malloc.h"

#include <stdbool.h>
#include <inttypes.h>

#define MAX(val1, val2) ((val1) < (val2) ? (val2) : (val1))
#define MIN(val1, val2) ((val1) > (val2) ? (val2) : (val1))

typedef struct _triangle {
    TIE_3 verts[3];
} triangle;

typedef struct _bounding_box {
    TIE_3 lows;
    TIE_3 highs;
} bounding_box;

TIE_3 div(TIE_3 a, TFLOAT b) {
    return (TIE_3){a.v[0] / b, a.v[1] / b, a.v[2] / b};
}

const int refs_that_fit_in_bounding_box = sizeof(bounding_box)/sizeof(triangle*);
typedef struct _bb_refs {
    triangle* refs[sizeof(bounding_box)/sizeof(triangle*)];
} bb_refs;

typedef struct _bvh_node {
    bool is_bounding_box;
    union {
	bounding_box bbox;
	bb_refs refs;
    } data;
} bvh_node;

void BuildBounds(triangle* trip, int trip_size, bvh_node* bounds_arr, int bounds_size, int offset, int cur_dim) {
    
    int offset_c1 = offset * 2 + 1;
    int offset_c2 = offset * 2 + 2;
    
    //if this assert fires, we have not allocated enough space
    //in the bounds array/tree
    BU_ASSERT(offset_c2 < bounds_size);
    
    bounds_arr[offset].is_bounding_box = true;
    bounding_box* bounds = &bounds_arr[offset].data.bbox;
    
    VSETALL(bounds->lows.v ,  INFINITY);
    VSETALL(bounds->highs.v, -INFINITY);
    
    //check each vertex in each triangle
    TIE_3 average = {0.0f};
    for (int t = 0; t < trip_size; t++) {
	for (int i = 0; i < 3; i++) {
	    TIE_3 pos = trip[t].verts[i];
	    VMIN(bounds->lows.v, pos.v);
	    VMAX(bounds->highs.v, pos.v);
	    VADD2(average.v, average.v, div(pos, (float)trip_size).v);
	}
    }
    average = div(average, 3.0f);
    
    if (trip_size <= refs_that_fit_in_bounding_box * 2) {
	//build leaf nodes and return
	bounds_arr[offset_c1].is_bounding_box = false;
	bounds_arr[offset_c2].is_bounding_box = false;
	
	for (int i = 0; i < refs_that_fit_in_bounding_box; i++) {
	    bounds_arr[offset_c1].data.refs.refs[i] = &trip[i];
	}
	for (int i = 0; i < trip_size -refs_that_fit_in_bounding_box; i++) {
	    bounds_arr[offset_c2].data.refs.refs[i] = &trip[i + refs_that_fit_in_bounding_box];
	}
	return;
    }
    // This is currently an insertion sort, because that was simple
    // However we want to convert this to a qsort pivot around the 
    // average, which we can calculate while we're checking the 
    // minimum and maximum bounds.
    
    //insertion sort
    for (int i = 1; i < trip_size; i++) {
	triangle* i_tri = &trip[i];
	float i_avg;
	i_avg  = i_tri->verts[0].v[cur_dim];
	i_avg += i_tri->verts[1].v[cur_dim];
	i_avg += i_tri->verts[2].v[cur_dim];
	i_avg /= 3.0;
	triangle swap = trip[i];
	int j;
	for (j = i-1; j >= 0; j--) {
	    triangle* j_tri = &trip[j];
	    float j_avg;
	    j_avg  = j_tri->verts[0].v[cur_dim];
	    j_avg += j_tri->verts[1].v[cur_dim];
	    j_avg += j_tri->verts[2].v[cur_dim];
	    j_avg /= 3.0;
	    if (j_avg < i_avg) {
		break;
	    }
	    trip[j+1] = trip[j];
	}
	trip[j+1] = swap;
    }
    
    int next_dim = (cur_dim + 1) %3;
    
    BuildBounds(trip                , trip_size /2            , bounds_arr, bounds_size, offset_c1, next_dim);
    BuildBounds(trip + trip_size /2, trip_size - trip_size/2, bounds_arr, bounds_size, offset_c2, next_dim);
    //we do the subtraction so we don't have to care if trip_size is a multiple of 2.
}

int bottie_prep_bvh(TIE_3* tri_arr, int num_tris, void* persistant_storage, int persistant_storage_size) {
    // We want to eventually convert this to a index into the triangle array
    // because that (u32) will take less space than a 64bit pointer
    BU_ASSERT(num_tris < (uint64_t)1 << 32);
    BU_ASSERT(num_tris > refs_that_fit_in_bounding_box *2);
    int num_leaf_nodes = num_tris / refs_that_fit_in_bounding_box;
    if (num_tris % refs_that_fit_in_bounding_box) {
	num_leaf_nodes += 1;
    }
    int log_count = 1;
    for (int num_nodes = num_leaf_nodes; num_nodes > 1; log_count++) {
	int add = num_nodes & 0x1;
	num_nodes = num_nodes >> 1;
	num_nodes += add;
    }
    int bounds_size = 1 << log_count;
    //m->bounds_size = bounds_size;
    //m->bounds = new BBnode[bounds_size];
    bu_log("Allocate %d nodes.\n", bounds_size);
    bvh_node* bounds = bu_calloc(bounds_size, sizeof(bvh_node), "BVH node preallocated arena, heap ordering");
    triangle* tris = bu_malloc(num_tris * sizeof(triangle), "Triangle Buffer for BVH");
    for (int i = 0; i < num_tris; i++) {
	VMOVE(tris[i].verts[0].v, tri_arr[i*3+0].v);
	VMOVE(tris[i].verts[1].v, tri_arr[i*3+1].v);
	VMOVE(tris[i].verts[2].v, tri_arr[i*3+2].v);
    }
    
    BuildBounds(tris, num_tris, bounds, bounds_size, 0, 0);
    
    struct tie_s * pass_up = (struct tie_s *)persistant_storage;
    pass_up->kdtree = (struct tie_kdtree_s*)bounds;
    pass_up->tri_list = (struct tie_tri_s*)tris;
    pass_up->rays_fired = 0;
    pass_up->max_depth = log_count;
    pass_up->tri_num = num_tris;
    pass_up->tri_num_alloc = num_tris;
    pass_up->stat = 0;
    pass_up->kdmethod = 0;
    VMOVE(pass_up->min, bounds[0].data.bbox.lows.v);
    VMOVE(pass_up->max, bounds[0].data.bbox.highs.v);
    VADD2(pass_up->mid, pass_up->min, pass_up->max);
    VSCALE(pass_up->mid, pass_up->mid, 0.5f);
    VSUB2(pass_up->amin, pass_up->min, pass_up->mid);
    VSUB2(pass_up->amax, pass_up->max, pass_up->mid);
    pass_up->radius = sqrt(VDOT(pass_up->amax, pass_up->amax));
}

void intersect_triangle(triangle* tri, 
			struct tie_ray_s* ray, 
			void* hitdata, 
			void *(*hitfunc)(struct tie_ray_s*, struct tie_id_s*, struct tie_tri_s*, void *ptr)) {
    if (!tri) return;
    TIE_3 E1, E2, S, S1, S2, norm;
    VSUB2(E1.v, tri->verts[1].v, tri->verts[0].v);
    VSUB2(E2.v, tri->verts[2].v, tri->verts[0].v);
    VSUB2( S.v, ray->pos, tri->verts[0].v);
    VCROSS(S1.v,ray->dir, E2.v);
    TFLOAT denom = VDOT(S1.v, E1.v);
    TFLOAT beta = VDOT(S1.v, S.v) / denom;
    //early out beta
    if ((beta < 0.0f) | (beta > 1.0f)) { return;}
    VCROSS(S2.v, S.v, E1.v);
    TFLOAT gamma = VDOT(S2.v, ray->dir) / denom;
    //early out gamma
    if ((gamma < 0.0f) | (gamma + beta > 1.0f)) { return;}
    TFLOAT alpha = 1 - beta - gamma;
    TFLOAT dist = VDOT(S2.v, E2.v) / denom;
    VCROSS(norm.v, E1.v, E2.v); 
    struct tie_id_s id;
    id.dist = dist;
    VMOVE(id.norm, norm.v);
    hitfunc(ray, &id, (void*)tri, hitdata);
}

TFLOAT minimum(TIE_3 point) {
    return MIN(point.v[0], MIN(point.v[1], point.v[2]));
}

TFLOAT maximum(TIE_3 point) {
    return MAX(point.v[0], MAX(point.v[1], point.v[2]));
}

bool intersect_bbox(bounding_box bbox, struct tie_ray_s* ray) {
    
    TIE_3  lows_t;
    TIE_3 highs_t;

    for (int i = 0; i < 3; i++) {
	 lows_t.v[i] = (bbox.lows.v[i]  - ray->pos[i]) / ray->dir[i];
	highs_t.v[i] = (bbox.highs.v[i] - ray->pos[i]) / ray->dir[i];
    }
    
    TIE_3 low_ts;
    TIE_3 high_ts;
    
    // VSETALL( low_ts.v,  INFINITY);
    // VSETALL(high_ts.v, -INFINITY);
    
    // VMINMAX(low_ts.v, high_ts.v,  lows_t.v);
    // VMINMAX(low_ts.v, high_ts.v, highs_t.v);
     low_ts.v[0] = MIN(lows_t.v[0], highs_t.v[0]);
     low_ts.v[1] = MIN(lows_t.v[1], highs_t.v[1]);
     low_ts.v[2] = MIN(lows_t.v[2], highs_t.v[2]);
    
    high_ts.v[0] = MAX(lows_t.v[0], highs_t.v[0]);
    high_ts.v[1] = MAX(lows_t.v[1], highs_t.v[1]);
    high_ts.v[2] = MAX(lows_t.v[2], highs_t.v[2]);
    
    TFLOAT high_t = minimum(high_ts);
    TFLOAT low_t  = maximum(low_ts);
    
    return (!ZERO(high_t)) & (low_t < high_t);
}
	
void recursive_bvh_intersect(bvh_node* bounds, 
			     int offset, 
			     struct tie_ray_s* ray, 
			     void* hitdata, 
			     void *(*hitfunc)(struct tie_ray_s*, struct tie_id_s*, struct tie_tri_s*, void*)) {

    if (bounds[offset].is_bounding_box) {
	if (intersect_bbox(bounds[offset].data.bbox, ray)) {
	    recursive_bvh_intersect(bounds, offset*2 +1, ray, hitdata, hitfunc);
	    recursive_bvh_intersect(bounds, offset*2 +2, ray, hitdata, hitfunc);
	}
    } else {
	for (int i = 0; i < refs_that_fit_in_bounding_box; i++) {
	    intersect_triangle(bounds[offset].data.refs.refs[i], ray, hitdata, hitfunc);
	}
    }
}

void bottie_shot_bvh(void* persistant_storage, 
		     struct tie_ray_s* ray, 
		     void *(*hitfunc)(struct tie_ray_s*, struct tie_id_s*, struct tie_tri_s*, void *ptr),
		     void* hitdata) {
    struct tie_s * pass_down = (struct tie_s *)persistant_storage;
    recursive_bvh_intersect((bvh_node*)pass_down->kdtree, 0, ray, hitdata, hitfunc);
}

void bottie_free_bvh(struct tie_s* tie) {
    bu_free(tie->kdtree, "BVH nodes");
    bu_free(tie->tri_list, "Triangle list");
}










