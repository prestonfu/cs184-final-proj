#define SPHERE_RADIUS 0.015
#define SPHERE_COUNT 512
#define WORK_GROUP_SIZE 64

#define NCLIP 0
#define FCLIP 4

#define DIFFUSE_BSDF (float3)(0.6, 0.6, 0.6) //includes color
#define AREA_LIGHT_POS (float3)(0, 1.5, 0)
#define AREA_LIGHT_DIR (float3)(0, -1, 0)
#define AREA_LIGHT_DIMX (float3)(0.5, 0, 0)
#define AREA_LIGHT_DIMY (float3)(0, 0, 0.5)
#define AREA_LIGHT_AREA 4 //make sure this matches dimx * dimy
#define AREA_LIGHT_RADIANCE (float3)(1, 1, 1)

#define GLOBAL_ILLUMINATION (float3)(0.05, 0.05, 0.05)
#define GLOBAL_ILLUMINATION_M (float3)(0.2, 0.2, 0.2)

#define NUM_RAYS 1
#define LIGHT_SAMPLES 16
#define LIGHT_SAMPLES_M 1
#define NUM_BOUNCES 0

#define EPS_F (0.00001f)

#define randf(seed) ((float)(seed = (((uint)(seed) * 16807) % 2147483647)) / 2147483647)
//#define randf(seed) ((float)(seed = abs(seed + 2147483647 / 4)) / 2147483647)

#define RAY_CAMERA 0
#define RAY_LIGHT 1
#define RAY_BOUNCE 2

#define MIN_THRESHOLD (0.001f)
#define MAX_THRESHOLD 10
#define NUM_ITERATIONS 20

#define K_SMOOTH (0.06f)
#define EPS_GRAD (0.001f)
#define EPS_BBOX (0.01f)

typedef struct {
    float3 o;
    float3 d; //should be normalized
    float mint;
    float maxt;
    float3 rad;
    int type;
} Ray;

typedef struct {
    float t;
    float3 n;
    bool hit;
    float3 color;
    bool light;
} Intersection;

bool intersect_bbox(Ray r, float3 min_vals, float3 max_vals)
{
    float t0 = r.mint;
    float t1 = r.maxt;
    for (int i = 0; i < 3; i++)
    {
        const float recip = 1 / r.d[i];
        t0 = max(t0, min((min_vals[i] - r.o[i]) * recip, (max_vals[i] - r.o[i]) * recip));
        t1 = min(t1, max((min_vals[i] - r.o[i]) * recip, (max_vals[i] - r.o[i]) * recip));
    }
    return t1 >= t0;
}

inline float smooth_min(float distA, float distB)
{
    float h = max(K_SMOOTH - fabs(distA - distB), 0.0f) * 1 / K_SMOOTH;
    return min(distA, distB) - h * h * h * K_SMOOTH * 1 / 6.0f;
}

kernel void raytrace
(
    write_only image2d_t out,
    const uint width,
    const uint height,
    const float hFov_expr,
    const float vFov_expr,
    const float3 camPos,
    const int selectedIndex,
    global const float* c2w,
    global const float* spheres,
    global const int* permutation, //for sphere positions
    global const float* bboxes,
    global const float* colors,
    global int* seedMemory
)
{
    local float3 localBuffer[WORK_GROUP_SIZE];
    local float3 localColor[WORK_GROUP_SIZE];
    local ushort localFlags[WORK_GROUP_SIZE];

    const uint xi = get_global_id(0);
    const uint yi = get_global_id(1);
    const uint lid = get_local_id(1) * get_local_size(0) + get_local_id(0);

    uint id = yi * width + xi;
    int seed = id;
    float3 c2w_0 = (float3)(c2w[0], c2w[3], c2w[6]);
    float3 c2w_1 = (float3)(c2w[1], c2w[4], c2w[7]);
    float3 c2w_2 = (float3)(c2w[2], c2w[5], c2w[8]);
    
    float3 radiance = (float3)(0, 0, 0);
    Ray r;
    float3 hit_p;
    float3 normal;
    float3 color;
    bool active = false;

    int raysLeft = NUM_RAYS;
    int bouncesLeft = 0;
    int lightSamplesLeft = 0;

    while(raysLeft > 0 || bouncesLeft > 0 || lightSamplesLeft > 0)
    {
        if (bouncesLeft == 0 && lightSamplesLeft == 0)
        {
            raysLeft--;
            bouncesLeft = NUM_BOUNCES;
            lightSamplesLeft = LIGHT_SAMPLES;

            float x = (xi + randf(seed)) / width;
            float y = 1 - (yi + randf(seed)) / height;
            float3 vec = (float3)((x - 0.5) * hFov_expr, (y - 0.5) * vFov_expr, -1);
            float3 res = vec.x * c2w_0 + vec.y * c2w_1 + vec.z * c2w_2;

            r.o = camPos;
            r.d = normalize(res);
            r.mint = NCLIP;
            r.maxt = FCLIP;
            r.type = RAY_CAMERA;
            
            active = true;
        }
        else if (lightSamplesLeft > 0)
        {
            lightSamplesLeft--;
            
            if (active)
            {
                float2 sample = (float2)(randf(seed) - 0.5, randf(seed) - 0.5);
                float3 d = AREA_LIGHT_POS + sample.x * AREA_LIGHT_DIMX + sample.y * AREA_LIGHT_DIMY - hit_p;
                float cosTheta = dot(d, AREA_LIGHT_DIR);
                float sqDist = dot(d, d);
                float dist = sqrt(sqDist);

                if (cosTheta > 0)
                    active = false;

                float3 wi = d / dist;
                float pdf = sqDist / (AREA_LIGHT_AREA * -cosTheta);

                r.o = hit_p;
                r.d = wi;
                r.mint = EPS_F;
                r.maxt = dist - EPS_F;
                r.rad = color * AREA_LIGHT_RADIANCE * dot(wi, normal) / pdf / LIGHT_SAMPLES;
                r.type = RAY_LIGHT;
            }
        }

        // Start intersection
        Intersection isect;
        isect.t = r.mint;
        isect.hit = false;

        localFlags[lid] = 0;
        if (active)
        {
            for (int k = 0; k < SPHERE_COUNT / WORK_GROUP_SIZE; k++)
            {
                if (intersect_bbox(r, (float3)(bboxes[6 * k], bboxes[6 * k + 1], bboxes[6 * k + 2]), (float3)(bboxes[6 * k + 3], bboxes[6 * k + 4], bboxes[6 * k + 5])))
                    localFlags[lid] |= (1 << k);
            }
        }

        barrier(CLK_LOCAL_MEM_FENCE);

        ushort flags = 0;
        for (int k = 0; k < WORK_GROUP_SIZE; k++)
            flags |= localFlags[k];
            
        for (int k = 0; k < SPHERE_COUNT / WORK_GROUP_SIZE; k++)
        {
            if ((flags & (1 << k)) == 0)// || (selectedIndex != -1 && k != selectedIndex))
                continue;
            int index = permutation[k * WORK_GROUP_SIZE + lid];
            localBuffer[lid] = (float3)(spheres[3 * index], spheres[3 * index + 1], spheres[3 * index + 2]);
            localColor[lid] = (float3)(colors[3 * index], colors[3 * index + 1], colors[3 * index + 2]);

            barrier(CLK_LOCAL_MEM_FENCE);

            if (active)
            {
                for (uint l = 0; l < WORK_GROUP_SIZE; l++)
                {
                    // test for intersections
                    float a = dot(r.d, r.d);
                    float b = 2.0 * dot(r.o - localBuffer[l], r.d);
                    float c = dot(r.o - localBuffer[l], r.o - localBuffer[l]) - SPHERE_RADIUS * SPHERE_RADIUS;

                    if (b*b - 4.0 * a * c >= 0.0)
                    {
                        float t_min = (-1.0 * b - sqrt(b*b - 4.0 * a * c)) / (2.0 * a);
                        float t_max = (-1.0 * b + sqrt(b*b - 4.0 * a * c)) / (2.0 * a);
                        if (t_max >= r.mint) 
                        {
                            if (r.mint <= t_min && t_min <= r.maxt)
                                isect.t = t_min;
                            else if (t_max <= r.maxt) 
                                isect.t = t_max;
                            else 
                                continue;

                            r.maxt = isect.t;
                            isect.n = normalize(r.o + (isect.t) * r.d - localBuffer[l]);
                            isect.color = localColor[l];
                            isect.hit = true;
                        }
                    }
                }
            }

            barrier(CLK_LOCAL_MEM_FENCE);
        }

        if (!isect.hit)
        {
            float t = (AREA_LIGHT_POS.y - r.o.y) / r.d.y;
            float3 p = r.o + r.d * t;
            if (t > 0 && fabs(p.x) < AREA_LIGHT_DIMX.x && fabs(p.z) < AREA_LIGHT_DIMY.y)
            {
                //isect.hit = true;
                active = false;
                radiance += AREA_LIGHT_RADIANCE;
            }
        }

        // Finish intersection

        if (r.type == RAY_CAMERA)
        {
            if (isect.hit)
            {
                hit_p = r.o + r.d * isect.t;
                normal = isect.n;
                color = isect.color;
                radiance += GLOBAL_ILLUMINATION;
            }
            else
                active = false;
        }
        else if (r.type == RAY_LIGHT)
        {
            if (active && !isect.hit)
            {
                radiance += r.rad;
                //radiance += isect.color * AREA_LIGHT_RADIANCE * dot(r.d, isect.n) / pdf / LIGHT_SAMPLES;
            }
        }
        // if (isect.hit)
        // {
        //     float3 hit_p = r.o + r.d * isect.t;

        //     float2 sample = (float2)(randf(seed) - 0.5, randf(seed) - 0.5);
        //     float3 d = AREA_LIGHT_POS + sample.x * AREA_LIGHT_DIMX + sample.y * AREA_LIGHT_DIMY - hit_p;
        //     float cosTheta = dot(d, AREA_LIGHT_DIR);
        //     float sqDist = dot(d, d);
        //     float dist = sqrt(sqDist);

        //     // if (cosTheta > 0)
        //     //     continue;

        //     float3 wi = d / dist;
        //     float pdf = sqDist / (AREA_LIGHT_AREA * -cosTheta);

        //     r.o = hit_p;
        //     r.d = wi;
        //     r.mint = EPS_F;
        //     r.maxt = dist - EPS_F;
        //     r.rad = isect.color * AREA_LIGHT_RADIANCE * dot(wi, isect.n) / pdf / LIGHT_SAMPLES;
                
        //         // Intersection temp = intersect(ray, spheres);
        //         // if (temp.hit)
        //         //     continue;
        //         // if (dot(wi, isect.n) < 0)
        //         // {
        //         //     continue;
        //         // }
        //         // radiance += isect.rad * rad * dot(wi, isect.n) / pdf / LIGHT_SAMPLES;
        //     //radiance += GLOBAL_ILLUMINATION;
        // }
    }
    radiance /= NUM_RAYS;
    
    write_imagef(out, (int2)(xi, yi), (float4)(radiance, 1.0));
    seedMemory[id] = seed;
}

kernel void raymarch
(
    write_only image2d_t out,
    const uint width,
    const uint height,
    const float hFov_expr,
    const float vFov_expr,
    const float3 camPos,
    const int selectedIndex,
    global const float* c2w,
    global const float* spheres,
    global const int* permutation, //for sphere positions
    global const float* bboxes,
    global const float* color,
    global int* seedMemory
)
{
    local float3 localBuffer[WORK_GROUP_SIZE];
    local ushort localFlags[WORK_GROUP_SIZE];

    const uint xi = get_global_id(0);
    const uint yi = get_global_id(1);
    const uint lid = get_local_id(1) * get_local_size(0) + get_local_id(0);

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
        float3 vec = (float3)((x - 0.5) * hFov_expr, (y - 0.5) * vFov_expr, -1);
        float3 res = vec.x * c2w_0 + vec.y * c2w_1 + vec.z * c2w_2;

        Ray r;
        r.o = camPos;
        r.d = normalize(res);
        r.mint = NCLIP;
        r.maxt = FCLIP;

        localFlags[lid] = 0;
        for (int k = 0; k < SPHERE_COUNT / WORK_GROUP_SIZE; k++)
        {
            if (intersect_bbox(r, (float3)(bboxes[6 * k] - EPS_BBOX, bboxes[6 * k + 1] - EPS_BBOX, bboxes[6 * k + 2] - EPS_BBOX), 
                                  (float3)(bboxes[6 * k + 3] + EPS_BBOX, bboxes[6 * k + 4] + EPS_BBOX, bboxes[6 * k + 5] + EPS_BBOX)))
                localFlags[lid] |= (1 << k);
        }

        barrier(CLK_LOCAL_MEM_FENCE);

        ushort flags = 0;
        for (int k = 0; k < WORK_GROUP_SIZE; k++)
            flags |= localFlags[k];

        Intersection isect;
        isect.t = r.mint;
        isect.hit = false;
    
        float3 hit_p = r.o + isect.t * r.d;
        bool done = false;
        for (int j = 0; j < NUM_ITERATIONS; j++)
        {
            float dist = MAX_THRESHOLD;
            
            for (int k = 0; k < SPHERE_COUNT / WORK_GROUP_SIZE; k++)
            {
                if ((flags & (1 << k)) == 0 || (selectedIndex != -1 && k != selectedIndex))
                    continue;
                int index = permutation[k * WORK_GROUP_SIZE + lid];

                localBuffer[lid] = (float3)(spheres[3 * index], spheres[3 * index + 1], spheres[3 * index + 2]);

                barrier(CLK_LOCAL_MEM_FENCE);

                if (!done)
                {
                    for (uint l = 0; l < WORK_GROUP_SIZE; l++)
                    {
                            dist = smooth_min(dist, fast_length(hit_p - localBuffer[l]) - SPHERE_RADIUS);
                    }
                }

                barrier(CLK_LOCAL_MEM_FENCE);
            }
            if (done)
                continue;
            isect.t += dist;
            hit_p = r.o + isect.t * r.d;

            if (dist < MIN_THRESHOLD)
            {
                isect.hit = true;
                done = true;
            }
        
            if (isect.t > r.maxt)
            {
                isect.hit = false;
                done = true;
            }
        }
        // Calculate normal
        float3 hits[6] = {hit_p, hit_p, hit_p, hit_p, hit_p, hit_p};
        hits[0].x += EPS_GRAD;
        hits[1].x -= EPS_GRAD;
        hits[2].y += EPS_GRAD;
        hits[3].y -= EPS_GRAD;
        hits[4].z += EPS_GRAD;
        hits[5].z -= EPS_GRAD;
        float dists[6] = {MAX_THRESHOLD, MAX_THRESHOLD, MAX_THRESHOLD, MAX_THRESHOLD, MAX_THRESHOLD, MAX_THRESHOLD};
        for (int k = 0; k < SPHERE_COUNT / WORK_GROUP_SIZE; k++)
        {
            int index = k * WORK_GROUP_SIZE + lid;
            localBuffer[lid] = (float3)(spheres[3 * index], spheres[3 * index + 1], spheres[3 * index + 2]);

            barrier(CLK_LOCAL_MEM_FENCE);

            if (isect.hit)
            {
                for (uint l = 0; l < WORK_GROUP_SIZE; l++)
                {
                    for (int j = 0; j < 6; j++)
                    {
                        float3 disp = hits[j] - localBuffer[l];
                        float new_dist = fast_length(disp) - SPHERE_RADIUS;

                        dists[j] = smooth_min(dists[j], new_dist);
                    }
                }
            }

            barrier(CLK_LOCAL_MEM_FENCE);
        }

        isect.n = normalize((float3)(dists[0] - dists[1], dists[2] - dists[3], dists[4] - dists[5]));

        if (isect.hit)
        {
            float3 hit_p = r.o + r.d * isect.t;
            for (int j = 0; j < LIGHT_SAMPLES_M; j++)
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

                if (dot(wi, isect.n) < 0)
                    continue;
                radiance += DIFFUSE_BSDF * rad * dot(wi, isect.n) / pdf / LIGHT_SAMPLES_M;
            }
            radiance += GLOBAL_ILLUMINATION;
        }
    }
    radiance /= NUM_RAYS;
    write_imagef(out, (int2)(xi, yi), (float4)(radiance, 1.0));
    seedMemory[id] = seed;
}
