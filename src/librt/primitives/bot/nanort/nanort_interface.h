#ifndef LIBRT_PRIMITIVES_BOT_NANORT_IFACE_H
#define LIBRT_PRIMITIVES_BOT_NANORT_IFACE_H

void nanort_push_double(void *vtie, TIE_3 **tri, unsigned int ntri, void *usr, unsigned int pstride);
int  nanort_prep_double(struct soltab *stp, struct rt_bot_internal *bot, struct rt_i *rtip);
int  nanort_shot_double(struct soltab *stp, register struct xray *rp, struct application *ap, struct seg *seghead);
void nanort_free_double(void *vtie);

void nanort_push_float(void *vtie, float **tri, unsigned int ntri, void *usr, unsigned int pstride);
int  nanort_prep_float(struct soltab *stp, struct rt_bot_internal *bot, struct rt_i *rtip);
int  nanort_shot_float(struct soltab *stp, register struct xray *rp, struct application *ap, struct seg *seghead);
void nanort_free_float(void *vtie);

#endif
