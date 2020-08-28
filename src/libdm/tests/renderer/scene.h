#ifndef SCENE_H
#define SCENE_H

#include "graphics.h"
#include "maths.h"
#include "mesh.h"
//#include "texture.h"

typedef struct {
    float frame_time;
    float delta_time;
    vec3_t light_dir;
    vec3_t camera_pos;
    mat4_t light_view_matrix;
    mat4_t light_proj_matrix;
    mat4_t camera_view_matrix;
    mat4_t camera_proj_matrix;
    float ambient_intensity;
    float punctual_intensity;
    int layer_view;
} perframe_t;

typedef struct model {
    mesh_t *mesh;
    program_t *program;
    mat4_t transform;
    /* for sorting */
    int opaque;
    float distance;
    /* polymorphism */
    void (*update)(struct model *model, perframe_t *perframe);
    void (*draw)(struct model *model, framebuffer_t *framebuffer,
                 int shadow_pass);
    void (*release)(struct model *model);
} model_t;

typedef struct {
    vec4_t background;
    model_t **models;
    /* light intensity */
    float ambient_intensity;
    float punctual_intensity;
} scene_t;

scene_t *scene_create(vec3_t background, model_t **models,
                      float ambient_intensity, float punctual_intensity
                      );
void scene_release(scene_t *scene);

#endif
