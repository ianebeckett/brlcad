#ifndef LIBRT_PRIMITIVES_BOT_PRBUCKET_IFACE_H
#define LIBRT_PRIMITIVES_BOT_PRBUCKET_IFACE_H



int bucketpr_build_double(struct soltab *stp, struct rt_bot_internal *bot_ip, struct rt_i *rtip);

void bucketpr_push_double(void *vtie, TIE_3 **tri, unsigned int ntri, void *usr, unsigned int pstride);
#endif