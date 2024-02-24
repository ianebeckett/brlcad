
#ifndef LIBRT_PRIMATIVES_BOT_BVH_H
#define LIBRT_PRIMATIVES_BOT_BVH_H
extern int bottie_prep_bvh(TIE_3*, int, void*, int);
extern void bottie_shot_bvh(void*, struct tie_ray_s*, void *(*func)(struct tie_ray_s*, struct tie_id_s*, struct tie_tri_s*, void*), void*);
extern void bottie_free_bvh(struct tie_s*);
#endif /* LIBRT_PRIMATIVES_BOT_BVH_H */