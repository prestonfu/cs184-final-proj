__kernel void ragen
(
    // Input
    uint width,
    uint height,
    Camera camera,
    __global uint* sample_counter,
    // Output
    __global Ray*    rays,
    __global uint*   ray_counter,
    __global uint*   pixel_indices,
    __global float3* throughputs,
    __global float3* diffuse_albedo,
    __global float*  depth_buffer,
    __global float3* normal_buffer,
    __global float2* velocity_buffer
)


{
    uint ray_idx = get_global_id(0);

    if (ray_idx >= width * height)
    {
        return;
    }

    uint pixel_idx = ray_idx;
    uint pixel_x = pixel_idx % width;
    uint pixel_y = pixel_idx / width;

    float inv_width = 1.0f / (float)(width);
    float inv_height = 1.0f / (float)(height);

    uint sample_idx = sample_counter[0];
    unsigned int seed = pixel_idx + HashUInt32(sample_idx);

#if 1
    float x = (pixel_x + GetRandomFloat(&seed)) * inv_width;
    float y = (pixel_y + GetRandomFloat(&seed)) * inv_height;
#else
    float x = (pixel_x + 0.5f) * inv_width;
    float y = (pixel_y + 0.5f) * inv_height;
#endif

    float angle = tan(0.5f * camera.fov);
    x = (x * 2.0f - 1.0f) * angle * camera.aspect_ratio;
    y = (y * 2.0f - 1.0f) * angle;

    float3 dir = normalize(x * cross(camera.front, camera.up) + y * camera.up + camera.front);

    // Depth of Field, as a primer for camera position
    float3 point_aimed = camera.position + camera.focus_distance * dir;
    float2 dof_dir = PointInHexagon(&seed);
    float r = camera.aperture;
    float3 new_pos = camera.position + dof_dir.x * r * cross(camera.front, camera.up) + dof_dir.y * r * camera.up;

    Ray ray;
    ray.origin.xyz = new_pos;
    ray.origin.w = 0.0;
    ray.direction.xyz = normalize(point_aimed - new_pos);
    // not too sure what distance should be added here
    ray.direction.w = MAX_RENDER_DIST; 

    rays[ray_idx] = ray;
    pixel_indices[ray_idx] = pixel_idx;
    throughputs[pixel_idx] = (float3)(1.0f, 1.0f, 1.0f);
    diffuse_albedo[pixel_idx] = (float3)(0.0f, 0.0f, 0.0f);
    // same case of not sure what distance to add here, will leave as is for now
    depth_buffer[pixel_idx] = MAX_RENDER_DIST;
    normal_buffer[pixel_idx] = (float3)(0.0f, 0.0f, 0.0f);
    velocity_buffer[pixel_idx] = (float2)(0.0f, 0.0f);

    // We are writing to global ray counter
    if (ray_idx == 0)
    {
        ray_counter[0] = width * height;
    }
}

