#ifndef GLOBALS
#define GLOBALS

#define MARCH

// Camera stuff
#define MIN_R (0.1f)
#define MAX_R (5.0f)

// Profiling
#define PROFILE 500

// Shared between all files
#define WORK_GROUP_SIZE 64
#define WORK_GROUP_SIZE_X 8
#define WORK_GROUP_SIZE_Y 8

// Shared between renderers
#define NCLIP 0
#define FCLIP 4
#define DIFFUSE_BSDF (float3)(1, 1, 1) //includes color
#define AREA_LIGHT_POS (float3)(0, 1, 0)
#define AREA_LIGHT_DIR (float3)(0, -1, 0)
#define AREA_LIGHT_DIMX (float3)(0.3, 0, 0)
#define AREA_LIGHT_DIMY (float3)(0, 0, 0.3)
#define AREA_LIGHT_AREA 1 //make sure this matches dimx * dimy
#define AREA_LIGHT_RADIANCE (float3)(1, 1, 1)
#define EPS_F (0.00001f)

#ifdef MARCH

// Shared
#define WIND_WIDTH 800
#define WIND_HEIGHT 600
#define SPHERE_RADIUS (0.015f)
#define SPHERE_COUNT 1024
#define NUM_RAYS 1

// Ray marching
#define MIN_THRESHOLD (0.001f)
#define MAX_THRESHOLD 10
#define NUM_ITERATIONS 40
#define K_SMOOTH (0.06f)
#define EPS_GRAD (0.001f)
#define EPS_BBOX (0.01f)

#define LIGHT_SAMPLES_M 1
#define GLOBAL_ILLUMINATION_M (float3)(0.2, 0.2, 0.2)

// Ray tracing with light blocking
#define LIGHT_SAMPLES 16
#define GLOBAL_ILLUMINATION (float3)(0.05, 0.05, 0.05)

#else
#define WIND_WIDTH 1280
#define WIND_HEIGHT 960
#define SPHERE_RADIUS (0.0075f)
#define SPHERE_COUNT 4096

// Ray tracing without light blocking
#define NUM_RAYS 1
#define LIGHT_SAMPLES 8
#define GLOBAL_ILLUMINATION (float3)(0.05, 0.05, 0.05)

#endif

#endif
