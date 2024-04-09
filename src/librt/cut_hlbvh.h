

struct bvh_build_node {
    fastf_t bounds[6];
    struct bvh_build_node *children[2];
    long first_prim_offset, n_primitives;
    uint8_t split_axis;
};

struct bvh_flat_node {
    fastf_t bounds[6];
    long n_primitives;
    union {
	long first_prim_offset;
	struct bvh_flat_node *other_child;
    } data;
};

#ifndef HLBVH_IMPLEMENTATION

RT_EXPORT extern struct bu_pool *
hlbvh_init_pool(size_t n_primatives);

RT_EXPORT extern struct bvh_build_node *
hlbvh_create(long max_prims_in_node, struct bu_pool *pool, const fastf_t *centroids_prims,
	     const fastf_t *bounds_prims, long *total_nodes,
	     const long n_primitives, long **ordered_prims);

RT_EXPORT void
hlbvh_shot(struct bvh_build_node* root, struct xray* rp, long** check_tris, size_t* num_check_tris);

#endif // HLBVH_IMPLEMENTATION