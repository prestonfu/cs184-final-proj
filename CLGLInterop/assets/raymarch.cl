#define SPHERE_RADIUS 0.01
#define SPHERE_COUNT 4096
#define MIN_THRESHOLD (0.001f)
#define MAX_THRESHOLD 10
#define NUM_ITERATIONS 40
#define WORK_GROUP_SIZE 256

#define NCLIP 0
#define FCLIP 4

#define DIFFUSE_BSDF (float3)(1, 1, 1) //includes color
#define AREA_LIGHT_POS (float3)(0, 2, 0)
#define AREA_LIGHT_DIR (float3)(0, -1, 0)
#define AREA_LIGHT_DIMX (float3)(2, 0, 0)
#define AREA_LIGHT_DIMY (float3)(0, 0, 2)
#define AREA_LIGHT_AREA 4 //make sure this matches dimx * dimy
#define AREA_LIGHT_RADIANCE (float3)(1, 1, 1)

#define GLOBAL_ILLUMINATION (float3)(0.2, 0.2, 0.2)

#define NUM_RAYS 1
#define LIGHT_SAMPLES 4

#define EPS_F (0.00001f)
#define EPS_GRAD (0.001f)

//#define K_SMOOTH (0.0005f)
#define K_SMOOTH (0.02f)

#define randf(seed) ((float)(seed = (((long)(seed) * 16807) % 2147483647)) / 2147483647)

#define SMOOTHING true

constant float3 colors[16] =
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

//bbox_1[0] = min_x, bbox_1[1] = min_y, bbox_1[2] = min_z
//bbox_2[0] = max_x, bbox_2[1] = max_y, bbox_2[2] = max_z
bool intersect_bbox(Ray r, float3 min_vals, float3 max_vals)
{
    // float3 t_min_per_dim = (float3)(INFINITY, INFINITY, INFINITY);
    // float3 t_max_per_dim = (float3)(-INFINITY, -INFINITY, -INFINITY);;
    // for (int i = 0; i < 3; i++)
    // {
    //     t_min_per_dim[i] = min(t_min_per_dim[i], (ray.d[i] == 0) ? INFINITY : (min_vals[i] - ray.o[i]) / ray.d[i]);
    //     t_max_per_dim[i] = max(t_max_per_dim[i], (ray.d[i] == 0) ? -INFINITY : (max_vals[i] - ray.o[i]) / ray.d[i]);
    // }
    // float t_min = max(t_min_per_dim.x, min(t_min_per_dim.y, t_min_per_dim.z));
    // float t_max = min(t_max_per_dim.x, min(t_max_per_dim.y, t_max_per_dim.z));
    
    // return t_min != INFINITY && t_max != -INFINITY;

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
    global const float* targetColor,
    global int* seedMemory
)
{
    local float3 localBuffer[WORK_GROUP_SIZE];
    local ushort localFlags[WORK_GROUP_SIZE];
    //local int cnt;

    const uint xi = get_global_id(0);
    const uint yi = get_global_id(1);
    const uint lid = get_local_id(1) * get_local_size(0) + get_local_id(0);
    //const uint localSize = get_local_size(0) * get_local_size(1);
    // if (xi >= width || yi >= height)
    //     return; //if we use local memory we may need these threads still to copy data

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
            if (intersect_bbox(r, (float3)(bboxes[6 * k], bboxes[6 * k + 1], bboxes[6 * k + 2]), (float3)(bboxes[6 * k + 3], bboxes[6 * k + 4], bboxes[6 * k + 5])))
                localFlags[lid] |= (1 << k);
        }
        // if (lid == 0)
        //     cnt = WORK_GROUP_SIZE;

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
                //int index = k * WORK_GROUP_SIZE + lid;
                //if (index < SPHERE_COUNT)
                localBuffer[lid] = (float3)(spheres[3 * index], spheres[3 * index + 1], spheres[3 * index + 2]);
                //localFlags[lid] = done;
                //else
                //    localBuffer[lid] = (float3)(INFINITY, INFINITY, INFINITY);

                barrier(CLK_LOCAL_MEM_FENCE);

                if (!done)
                {
                    for (uint l = 0; l < WORK_GROUP_SIZE; l++)
                    {
                        if (SMOOTHING)
                        {
                            dist = smooth_min(dist, fast_length(hit_p - localBuffer[l]) - SPHERE_RADIUS);
                        }
                        else
                        { 
                            float3 disp = hit_p - localBuffer[l];
                            float new_dist = disp.x * disp.x + disp.y * disp.y + disp.z * disp.z;
                            if (new_dist < dist)
                            {
                                dist = new_dist;
                                isect.n = localBuffer[l];
                                isect.rad = colors[k];
                            }
                        }
                    }
                }

                barrier(CLK_LOCAL_MEM_FENCE);
            }
            if (done)
                continue;            
            if (!SMOOTHING)
                dist = sqrt(dist) - SPHERE_RADIUS;
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
        if (SMOOTHING)
        {
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
                            //float new_dist = disp.x * disp.x + disp.y * disp.y + disp.z * disp.z;
                            float new_dist = fast_length(disp) - SPHERE_RADIUS;

                            dists[j] = smooth_min(dists[j], new_dist);
                        }
                    }
                }

                barrier(CLK_LOCAL_MEM_FENCE);
            }
            //for (int j = 0; j < 6; j++)
            //    dists[j] = sqrt(dists[j]) - SPHERE_RADIUS;

            isect.n = normalize((float3)(dists[0] - dists[1], dists[2] - dists[3], dists[4] - dists[5]));
        }
//#endif
        else
        {
            isect.n = normalize(hit_p - isect.n);
        }
        // if (xi >= width || yi >= height)
        //     continue;

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
                //radiance += isect.rad * rad * dot(wi, isect.n) / pdf / LIGHT_SAMPLES;
                radiance += DIFFUSE_BSDF * rad * dot(wi, isect.n) / pdf / LIGHT_SAMPLES;
            }
            radiance += GLOBAL_ILLUMINATION;
        }
    }
    radiance /= NUM_RAYS;
    write_imagef(out, (int2)(xi, yi), (float4)(radiance, 1.0));
    seedMemory[id] = seed;
}