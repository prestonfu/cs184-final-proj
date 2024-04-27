#include <colormaps.h>

#define H_FOV 50
#define V_FOV 35

#define SPHERE_RADIUS 0.1
#define SPHERE_COUNT 4096
#define MIN_THRESHOLD 0.01
#define MAX_THRESHOLD 9999
#define NUM_ITERATIONS 50

#define DIFFUSE_BSDF (float3)(1, 1, 1) //includes color
#define AREA_LIGHT_POS (float3)(0, 0, 1)
#define AREA_LIGHT_DIR (float3)(0, 0, -1)
#define AREA_LIGHT_DIMX (float3)(1, 0, 0)
#define AREA_LIGHT_DIMY (float3)(0, 1, 0)
#define AREA_LIGHT_AREA 1 //make sure this matches dimx * dimy
#define AREA_LIGHT_RADIANCE (float3)(1, 1, 1)

#define NUM_RAYS 8
#define LIGHT_SAMPLES 8

typedef struct {
    float3 o;
    float3 d;
    float mint;
    float maxt;
} Ray;

typedef struct {
    float t; 
    float3 n;
} Intersection;


int rand(int* seed) // 1 <= *seed < m
{
    int const a = 16807; //ie 7**5
    int const m = 2147483647; //ie 2**31-1

    *seed = ((long)(*seed) * a) % m;
    return(*seed);
}

float randf(int* seed) // rand from [0, 1]
{
    float const reciprocal = 1.0 / 2147483646;
    return rand(seed) * reciprocal;
}

float signed_dist(Ray r, float x, float y, float z) 
{
    //distance of ray point to center of circle
    float distance = sqrt((r.o.x - x) * (r.o.x - x) + (r.o.y - y) * (r.o.y - y) + (r.o.z - z) * (r.o.z - z));
    return distance - SPHERE_RADIUS;
}


//this function is for when two spheres intersect with each other, other options available
// https://iquilezles.org/articles/smin/
float smoothMin(float dist_first, float dist_second, float k) 
{
    float h = max(k - fabs(dist_first - dist_second), 0.0)/k;
    float minimum = min(dist_first, dist_second);
    return minimum - h*h*h*k*1/6.0;
}




bool intersect(Ray r, Intersection *intersection, global const float* sphere_list)
{
    // 100 is the max num steps we "grow" our ray out until giving up
    // for (int i = 0; i < NUM_ITERATIONS; i++) {
    //     float dist = INFINITY;
    //     for (uint j = 0; j < SPHERE_COUNT; j++) {
    //         // Define a sphere at point x_coord, y_coord, z_coord 
    //         float x_coord = sphere_list[3 * j];
    //         float y_coord = sphere_list[3 * j + 1];
    //         float z_coord = sphere_list[3 * j + 2];
            
    //         float new_dist = signed_dist(r, x_coord, y_coord, z_coord);
    //         if (new_dist < dist) {
    //             dist = new_dist;
    //         }
            
    //     }
    //     if (dist < MIN_THRESHOLD && length(r.o + r.d * dist) > r.mint && length(r.o + r.d * dist) < r.maxt) {
    //         return true;
    //     }  else if (dist > MAX_THRESHOLD) {
    //         return false;
    //     }
    //     r.o = r.o + r.d * dist;
    // }
    return false;
}

float3 mat_mul(const float* mat, float3 v)
{
    float3 result;
    result.x = mat[0] * v.x + mat[3] * v.y + mat[6] * v.z;
    result.y = mat[1] * v.x + mat[4] * v.y + mat[7] * v.z;
    result.z = mat[2] * v.x + mat[5] * v.y + mat[8] * v.z;
    return result;
}

Ray generate_ray(float x, float y, const float* c2w, float3 camPos) 
{

  // TODO (Part 1.1):
  // compute position of the input sensor sample coordinate on the
  // canonical sensor plane one unit away from the pinhole.
  // Note: hFov and vFov are in degrees.
  //
  float new_x = x - 0.5;
  new_x = new_x * tan(radians(H_FOV/2.0));
  new_x = new_x * 2.0;

  double new_y = y - 0.5;
  new_y = new_y * tan(radians(V_FOV/2.0));
  new_y = new_y * 2.0;

  float3 vect = (float3)(new_x, new_y, -1);

  vect = mat_mul(c2w, vect);

  Ray r = { camPos, normalize(vect) };

  return r;
}

kernel void raytrace
(
    write_only image2d_t out,
    const uint width,
    const uint height,
    const float3 camPos,
    global const float* c2w,
    global const float* spheres,
    global int* seedMemory
)
{
    const int xi = get_global_id(0);
    const int yi = get_global_id(1);
    if (xi >= width || yi >= height)
        return;

    uint id = get_global_id(1) * get_global_size(0) + get_global_id(0);
    int seed = seedMemory[id];

    float private_c2w[9];
    for (int i = 0; i < 9; i++)
        private_c2w[i] = c2w[i];


    Ray ray = generate_ray(xi + 0.5, yi + 0.5, private_c2w, camPos);
    ray.mint = FLT_EPSILON;
    ray.maxt = MAX_THRESHOLD;
    Intersection isect;
    
    float3 radiance = (float3)(0, 0, 0);
    for (int i = 0; i < NUM_RAYS; i++)
    {
        if (intersect(ray, &isect, spheres))
        {
            float3 hit_p = ray.o + ray.d * isect.t;
            for (int j = 0; j < LIGHT_SAMPLES; j++)
            {
                float2 sample = (float2)(randf(&seed) - 0.5, randf(&seed) - 0.5);
                float3 d = AREA_LIGHT_POS + sample.x * AREA_LIGHT_DIMX + sample.y * AREA_LIGHT_DIMY - hit_p;
                float cosTheta = dot(d, AREA_LIGHT_DIR);
                float sqDist = d.x * d.x + d.y * d.y + d.z * d.z;
                float dist = sqrt(sqDist);

                if (cosTheta < 0)
                    continue;

                float3 wi = d / dist;
                float distToLight = dist;
                float pdf = sqDist / (AREA_LIGHT_AREA * fabs(cosTheta));
                float3 rad = AREA_LIGHT_RADIANCE;

                Ray ray;
                ray.o = hit_p;
                ray.d = wi;
                ray.mint = FLT_EPSILON;
                ray.maxt = distToLight - FLT_EPSILON;
                Intersection temp;
                if (intersect(ray, &temp, spheres))
                    continue;
                radiance += DIFFUSE_BSDF * rad * dot(wi, isect.n) / pdf / LIGHT_SAMPLES;
            }
        }
    }
    radiance /= NUM_RAYS;
    // radiance.x /= 10;
    // radiance.y /= 10;
    // radiance.z /= 10;
    write_imagef(out, (int2)(xi, yi), (float4)(radiance.x, radiance.y, radiance.z, 1.0));

    seedMemory[id] = seed;
}