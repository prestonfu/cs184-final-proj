#define SPHERE_RADIUS 0.01
#define SPHERE_COUNT 400
#define MIN_THRESHOLD (0.00001f)
#define MAX_THRESHOLD 100
#define NUM_ITERATIONS 100

#define NCLIP 0
#define FCLIP 5

#define DIFFUSE_BSDF (float3)(1, 1, 1) //includes color
#define AREA_LIGHT_POS (float3)(0, 2, 0)
#define AREA_LIGHT_DIR (float3)(0, -1, 0)
#define AREA_LIGHT_DIMX (float3)(2, 0, 0)
#define AREA_LIGHT_DIMY (float3)(0, 0, 2)
#define AREA_LIGHT_AREA 4 //make sure this matches dimx * dimy
#define AREA_LIGHT_RADIANCE (float3)(1, 1, 1)

#define GLOBAL_ILLUMINATION (float3)(0.1, 0.1, 0.1)

#define NUM_RAYS 1
#define LIGHT_SAMPLES 1

#define EPS_F (0.00001f)

typedef struct {
    float3 o;
    float3 d; //should be normalized
    float mint;
    float maxt;
} Ray;

typedef struct {
    float t; 
    float3 n;
    bool hit;
} Intersection;

#define randf(seed) ((float)(seed = (((long)(seed) * 16807) % 2147483647)) / 2147483647)

kernel void raytrace
(
    write_only image2d_t out,
    const uint width,
    const uint height,
    const float hFov,
    const float vFov,
    const float3 camPos,
    global const float* c2w,
    global const float* spheres,
    global int* seedMemory
)
{
    const uint xi = get_global_id(0);
    const uint yi = get_global_id(1);
    if (xi >= width || yi >= height)
        return; //if we use local memory we may need these threads still to copy data

    uint id = yi * width + xi;
    int seed = id;
    float3 c2w_0 = (float3)(c2w[0], c2w[3], c2w[6]);
    float3 c2w_1 = (float3)(c2w[1], c2w[4], c2w[7]);
    float3 c2w_2 = (float3)(c2w[2], c2w[5], c2w[8]);
    
    float3 radiance = (float3)(0, 0, 0);
    for (int i = 0; i < NUM_RAYS; i++)
    {
        float x = (xi + randf(seed)) / width;
        float y = 1 - (yi + randf(seed)) / height;
        float3 vec = (float3)(2 * (x - 0.5) * tan(0.5 * hFov * M_PI / 180), 2 * (y - 0.5) * tan(0.5 * vFov * M_PI / 180), -1);
        float3 res = vec.x * c2w_0 + vec.y * c2w_1 + vec.z * c2w_2;

        Ray r;
        r.o = camPos;
        r.d = normalize(res);
        r.mint = NCLIP;
        r.maxt = FCLIP;

        Intersection isect;
        isect.t = r.mint;
        isect.hit = false;
    
        float3 hit_p = r.o + isect.t * r.d;
        for (int j = 0; j < NUM_ITERATIONS; j++)
        {
            float dist = MAX_THRESHOLD;
            float3 c;
            for (uint k = 0; k < SPHERE_COUNT; k++)
            {
                float3 center = (float3)(spheres[3 * k], spheres[3 * k + 1], spheres[3 * k + 2]);
                float new_dist = length(hit_p - center) - SPHERE_RADIUS;

                if (new_dist < dist)
                {
                    dist = new_dist;
                    c = center;
                }
            }
            isect.t += dist;
            
            hit_p = r.o + isect.t * r.d;
            if (dist < MIN_THRESHOLD)
            {
                isect.n = normalize(hit_p - c);
                isect.hit = true;
            }
        
            if (isect.t > r.maxt || dist > MAX_THRESHOLD )
            {
                isect.hit = false;
            }
        }
       
        if (isect.hit)
        {
            float3 hit_p = r.o + r.d * isect.t;
            for (int j = 0; j < LIGHT_SAMPLES; j++)
            {
                float2 sample = (float2)(randf(seed) - 0.5, randf(seed) - 0.5);
                float3 d = AREA_LIGHT_POS + sample.x * AREA_LIGHT_DIMX + sample.y * AREA_LIGHT_DIMY - hit_p;
                float cosTheta = dot(d, AREA_LIGHT_DIR);
                float sqDist = d.x * d.x + d.y * d.y + d.z * d.z;
                float dist = sqrt(sqDist);

                if (cosTheta > 0)
                    continue;

                float3 wi = d / dist;
                float pdf = sqDist / (AREA_LIGHT_AREA * -cosTheta);
                float3 rad = AREA_LIGHT_RADIANCE;

                Ray ray;
                ray.o = hit_p;
                ray.d = wi;
                ray.mint = EPS_F;
                ray.maxt = dist - EPS_F;
                
                // Intersection temp = intersect(ray, spheres);
                // if (temp.hit)
                //     continue;
                if (dot(wi, isect.n) < 0)
                    continue;
                radiance += DIFFUSE_BSDF * rad * dot(wi, isect.n) / pdf / LIGHT_SAMPLES;
            }
            radiance += GLOBAL_ILLUMINATION;
        }
    }
    radiance /= NUM_RAYS;
    write_imagef(out, (int2)(xi, yi), (float4)(radiance, 1.0));
    seedMemory[id] = seed;
}