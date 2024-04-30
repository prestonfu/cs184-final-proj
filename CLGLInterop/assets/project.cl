#define NUM_SPHERES 4096
#define SPHERE_RADIUS 0.01

kernel void project
(    
    const uint width,
    const uint height,
    const float hFov_expr,
    const float vFov_expr,
    const float3 camPos,
    global const float* c2w, 
    global const float* spheres,
    global float* projection
)
{
    float3 w2c_0 = (float3)(c2w[0], c2w[1], c2w[2]);
    float3 w2c_1 = (float3)(c2w[3], c2w[4], c2w[5]);
    float3 w2c_2 = (float3)(c2w[6], c2w[7], c2w[8]);
    
    int i = get_global_id(0);
    float sphere_offset_x = spheres[3 * i] - camPos.x;
    float sphere_offset_y = spheres[3 * i + 1] - camPos.y;
    float sphere_offset_z = spheres[3 * i + 2] - camPos.z;
    float3 vec = w2c_0 * sphere_offset_x + w2c_1 * sphere_offset_y + w2c_2 * sphere_offset_z;

    float min_x = INFINITY;
    float min_y = INFINITY;
    float max_x = -INFINITY;
    float max_y = -INFINITY;
    for (int j = 0; j < 8; j++) {
        //float3 projection = vec + (2*(j & 4) - 1) * w2c_0 + (2*(j & 2) - 1) * w2c_1 + (2*(j & 1) - 1) * w2c_2;
        float3 projection = vec + SPHERE_RADIUS * (float3)((j & 4) / 2 - 1, (j & 2) - 1, (j & 1) * 2 - 1);
        projection = - projection / projection.z;
        float x = 0.5 + projection[0] / hFov_expr;
        float y = 0.5 + projection[1] / vFov_expr;
        min_x = min(min_x, x);
        min_y = min(min_y, y);
        max_x = max(max_x, x);
        max_y = max(max_y, y);
    }
    
    projection[4 * i + 0] = min_x;
    projection[4 * i + 1] = min_y;
    projection[4 * i + 2] = max_x;
    projection[4 * i + 3] = max_y;
}
