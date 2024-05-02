#define randf(seed) ((float)(seed = (((long)(seed) * 16807) % 2147483647)) / 2147483647)

constant float3 debugColors[16] =
{
    (float3)(0, 0, 1),
    (float3)(0, 1, 0),
    (float3)(0, 1, 1),
    (float3)(1, 0, 0),
    (float3)(1, 0, 1),
    (float3)(1, 1, 0),
    (float3)(1, 1, 1),
    (float3)(0, 0.5, 0.5),
    (float3)(0.5, 0, 0.5),
    (float3)(0.5, 0.5, 0),
    (float3)(1, 0.5, 0),
    (float3)(1, 0, 0.5),
    (float3)(0.5, 1, 0),
    (float3)(0.5, 0, 1),
    (float3)(0, 0.5, 1),
    (float3)(0, 1, 0.5)
};

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
    float3 rad;
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
    local ulong localFlags[WORK_GROUP_SIZE];

    const uint xi = get_global_id(0);
    const uint yi = get_global_id(1);
    const uint lid = get_local_id(1) * get_local_size(0) + get_local_id(0);

    uint id = yi * width + xi;
    int seed = seedMemory[id];
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
            if (intersect_bbox(r, (float3)(bboxes[6 * k], bboxes[6 * k + 1], bboxes[6 * k + 2]), (float3)(bboxes[6 * k + 3], bboxes[6 * k + 4], bboxes[6 * k + 5])))
                localFlags[lid] |= (1 << k);
        }

        barrier(CLK_LOCAL_MEM_FENCE);

        ulong flags = 0;
        for (int k = 0; k < WORK_GROUP_SIZE; k++)
            flags |= localFlags[k];

        Intersection isect;
        isect.t = r.mint;
        isect.hit = false;
        
        for (int k = 0; k < SPHERE_COUNT / WORK_GROUP_SIZE; k++)
        {
            if ((flags & (1 << k)) == 0)// || (selectedIndex != -1 && k != selectedIndex))
               continue;
            int index = permutation[k * WORK_GROUP_SIZE + lid];
            localBuffer[lid] = (float3)(spheres[3 * index], spheres[3 * index + 1], spheres[3 * index + 2]);
            localColor[lid] = (float3)(colors[3 * index], colors[3 * index + 1], colors[3 * index + 2]);

            barrier(CLK_LOCAL_MEM_FENCE);

            for (uint l = 0; l < WORK_GROUP_SIZE; l++)
            {
                float t1, t2;

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
                        t1 = t_min;
                        t2 = t_max;
                        if (r.mint <= t1 && t1 <= r.maxt)
                        {
                            isect.t = t1;
                        } 
                        else if (t2 <= r.maxt) 
                        {
                            isect.t = t2;
                        } 
                        else 
                        {
                            continue;
                        }

                        isect.n = normalize(r.o + (isect.t) * r.d - localBuffer[l]);
                        isect.rad = localColor[l];
                        if (selectedIndex == -2)
                            isect.rad = debugColors[k % 16];
                        r.maxt = isect.t;
                        isect.hit = true;
                    }
                }
            }

            barrier(CLK_LOCAL_MEM_FENCE);
        }
/* global ver
        for (int k = 0; k < SPHERE_COUNT; k++)
        {
            if ((flags & (1 << (k / WORK_GROUP_SIZE))) == 0)// || (selectedIndex != -1 && k != selectedIndex))
            {
                k += WORK_GROUP_SIZE;
                continue;
            }
            int index = permutation[k];
            float3 pos = (float3)(spheres[3 * index], spheres[3 * index + 1], spheres[3 * index + 2]);
            float3 color = (float3)(colors[3 * index], colors[3 * index + 1], colors[3 * index + 2]);

            float t1, t2;

            // test for intersections
            float a = dot(r.d, r.d);
            float b = 2.0 * dot(r.o - pos, r.d);
            float c = dot(r.o - pos, r.o - pos) - SPHERE_RADIUS * SPHERE_RADIUS;

            if (b*b - 4.0 * a * c >= 0.0)
            {
                float t_min = (-1.0 * b - sqrt(b*b - 4.0 * a * c)) / (2.0 * a);
                float t_max = (-1.0 * b + sqrt(b*b - 4.0 * a * c)) / (2.0 * a);
                if (t_max >= r.mint) 
                {
                    t1 = t_min;
                    t2 = t_max;
                    if (r.mint <= t1 && t1 <= r.maxt)
                    {
                        isect.t = t1;
                    } 
                    else if (t2 <= r.maxt) 
                    {
                        isect.t = t2;
                    } 
                    else 
                    {
                        continue;
                    }

                    isect.n = normalize(r.o + (isect.t) * r.d - pos);
                    isect.rad = color;
                    r.maxt = isect.t;
                    isect.hit = true;
                }
            }
        }
*/
        float t = (AREA_LIGHT_POS.y - r.o.y) / r.d.y;
        float3 p = r.o + r.d * t;
        if (t > 0 && fabs(p.x) < AREA_LIGHT_DIMX.x && fabs(p.z) < AREA_LIGHT_DIMY.z)
        {
            if (!isect.hit || t < isect.t)
            {
                isect.hit = false;
                radiance += r.o.y < AREA_LIGHT_POS.y ? AREA_LIGHT_RADIANCE : (float3)(0.17, 0.17, 0.17);
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
                float sqDist = dot(d, d);
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
                radiance += isect.rad * rad * dot(wi, isect.n) / pdf / LIGHT_SAMPLES;
            }
            radiance += GLOBAL_ILLUMINATION;
        }
    }
    radiance /= NUM_RAYS;
    write_imagef(out, (int2)(xi, yi), (float4)(radiance, 1.0));
    seedMemory[id] = seed;
}
