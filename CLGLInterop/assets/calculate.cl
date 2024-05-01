#define SPHERE_COUNT 256
#define WORK_GROUP_SIZE 256

#define RADIUS_SEPARATION 0.1
#define RADIUS_COHESION 0.15
#define RADIUS_ALIGNMENT 0.08
#define RADIUS_CONTAINMENT 1
#define RADIUS_ARRIVAL 0.05
#define RADIUS_CENTER 0.5

#define K_SEPARATION 7
#define K_COHESION 7.5
#define K_ALIGNMENT 4
#define K_CONTAINMENT 20
#define K_ARRIVAL 70
#define K_CENTER 2

#define T_CONTAINMENT 1

#define MAX_VELOCITY 0.5

kernel void calculate
(
    global float* positions, 
    global const float* velocities, 
    global float* accelerations,
    global float* targetPos,
    const int selectedIndex
)
{
    local float3 otherPos[WORK_GROUP_SIZE];
    local float3 otherVel[WORK_GROUP_SIZE];

    int id = get_global_id(0);
    int lid = get_local_id(0);
    float3 pos = (float3)(positions[3 * id], positions[3 * id + 1], positions[3 * id + 2]);
    float3 vel = (float3)(velocities[3 * id], velocities[3 * id + 1], velocities[3 * id + 2]);

    float3 accel = (float3)(0.0f, 0.0f, 0.0f);
    float3 separation, cohesion, alignment;
    separation = cohesion = alignment = accel;
    int nCohesion = 0, nAlignment = 0;

    for (int i = 0; i < SPHERE_COUNT; i++)
    {
        if (i == id)
            continue;
            
        
            
        float3 otherPos = (float3)(positions[3 * i], positions[3 * i + 1], positions[3 * i + 2]);
        float distance = length(pos - otherPos);
        if (dot(pos - otherPos, normalize(vel)) < -0.5 * distance)
            continue;
        float3 otherVel = (float3)(velocities[3 * i], velocities[3 * i + 1], velocities[3 * i + 2]);
        if (distance < RADIUS_SEPARATION)
        {
            float3 displacement = pos - otherPos;
            separation += displacement / distance / distance;
        }
        if (distance < RADIUS_COHESION)
        {
            cohesion += otherPos;
            nCohesion++;
        }
        if (distance < RADIUS_ALIGNMENT)
        {
            alignment += otherVel;
            nAlignment++;
        }
    }
    // for (int i = 0; i < SPHERE_COUNT; i += WORK_GROUP_SIZE)
    // {
    //     int index = i + lid;
    //     otherPos[lid] = (float3)(positions[3 * index], positions[3 * index + 1], positions[3 * index + 2]);
    //     otherVel[lid] = (float3)(velocities[3 * index], velocities[3 * index + 1], velocities[3 * index + 2]);

    //     barrier(CLK_LOCAL_MEM_FENCE);

    //     if (index != id)
    //     {
    //         for (int j = 0; j < WORK_GROUP_SIZE; j++)
    //         {
    //             float3 otherPos = otherPos[j];//(float3)(positions[3 * i], positions[3 * i + 1], positions[3 * i + 2]);
    //             float distance = length(pos - otherPos);
    //             if (dot(pos - otherPos, normalize(vel)) < -0.5 * distance)
    //                 continue;
    //             float3 otherVel = otherVel[j];//(float3)(velocities[3 * i], velocities[3 * i + 1], velocities[3 * i + 2]);
    //             if (distance < RADIUS_SEPARATION)
    //             {
    //                 float3 displacement = pos - otherPos;
    //                 separation += displacement / distance / distance;
    //             }
    //             if (distance < RADIUS_COHESION)
    //             {
    //                 cohesion += otherPos;
    //                 nCohesion++;
    //             }
    //             if (distance < RADIUS_ALIGNMENT)
    //             {
    //                 alignment += otherVel;
    //                 nAlignment++;
    //             }
    //         }
    //     }

    //     barrier(CLK_LOCAL_MEM_FENCE);
    // }

    accel += K_SEPARATION * normalize(separation);
    if (nCohesion > 0)
    {
        float3 desiredPosition = cohesion / nCohesion;
        float3 desiredVelocity = normalize(desiredPosition - pos) * MAX_VELOCITY;
        accel += K_COHESION * normalize(desiredVelocity - vel);
    }
    if (nAlignment > 0)
    {
        float3 desiredVelocity = normalize(alignment / nAlignment) * MAX_VELOCITY;
        accel += K_ALIGNMENT * normalize(desiredVelocity - vel);
    }

    float a = vel.x * vel.x + vel.y * vel.y + vel.z * vel.z;
    float b = 2 * dot(pos, vel);
    float c = pos.x * pos.x + pos.y * pos.y + pos.z * pos.z - RADIUS_CONTAINMENT * RADIUS_CONTAINMENT;

    float t = (sqrt(b * b - 4 * a * c) - b) / (2 * a);
    if (t > 0 && t < T_CONTAINMENT)
    {
        float3 normal = -normalize(pos + t * vel);
        float3 dir = normalize(vel);
        accel += K_CONTAINMENT * normalize(normal - dot(normal, dir) * dir);
    }

    if (length(pos) > RADIUS_CENTER)
        accel += -normalize(pos) * K_CENTER;

    if (selectedIndex != -1)
    {
        float3 targetOffset = (float3)(targetPos[3 * id], targetPos[3 * id + 1], targetPos[3 * id + 2]) - pos;
        float3 desiredVelocity = min(MAX_VELOCITY, length(targetOffset) / RADIUS_ARRIVAL) * normalize(targetOffset);
        accel += K_ARRIVAL * normalize(desiredVelocity - vel);
    }

    accelerations[3 * id] = accel.x;
    accelerations[3 * id + 1] = accel.y;
    accelerations[3 * id + 2] = accel.z;
}

/*
 Vector3[] forces = new Vector3[num];
        for (int i = 0; i < num; i++)
        {
            Vector3 separation = new Vector3();
            (Vector3 pos, int n) cohesion = (new Vector3(), 0);
            (Vector3 vel, int n) alignment = (new Vector3(), 0);
            for (int j = 0; j < num; j++)
            {
                if (i == j)
                    continue;
                if (Vector3.Dot((boids[j].position - boids[i].position).normalized, velocity[i].normalized) < -0.5)
                    continue;
                float distance = (boids[i].position - boids[j].position).magnitude;
                if (distance < radiusSeparation)
                {
                    Vector3 displacement = boids[i].position - boids[j].position;
                    separation += displacement / displacement.sqrMagnitude;
                    //separation += displacement.normalized;
                }
                if (distance < radiusCohesion)
                {
                    cohesion.pos += boids[j].position;
                    cohesion.n++;
                }
                if (distance < radiusAlignment)
                {
                    alignment.vel += velocity[j];
                    alignment.n++;
                }
            }
            
            forces[i] += kSeparation * separation.normalized;
            //if (separation != new Vector3())
            //    forces[i] += kSeparation * (separation.normalized * maxVelocity - velocity[i]).normalized;

            if (cohesion.n > 0)
            {
                Vector3 desiredPosition = cohesion.pos / cohesion.n;
                Vector3 desiredVelocity = (desiredPosition - boids[i].position).normalized * maxVelocity;
                forces[i] += kCohesion * (desiredVelocity - velocity[i]).normalized;
            }
            if (alignment.n > 0)
            {
                Vector3 desiredVelocity = (alignment.vel / alignment.n).normalized * maxVelocity;
                forces[i] += kAlignment * (desiredVelocity - velocity[i]).normalized;
            }

            float a = velocity[i].sqrMagnitude;
            float b = 2 * Vector3.Dot(boids[i].position, velocity[i]);
            float c = boids[i].position.sqrMagnitude - radiusContainment * radiusContainment;

            float t = (Mathf.Sqrt(b * b - 4 * a * c) - b) / (2 * a);
            if (t > 0 && t < tContainment)
            {
                Vector3 normal = -(boids[i].position + t * velocity[i]).normalized;
                Vector3 dir = velocity[i].normalized;
                forces[i] += kContainment * (normal - Vector3.Dot(normal, dir) * dir).normalized;
            }

            if (targets != null && targets.Length > 0)
            {
                Vector3 targetOffset = targets[i] - boids[i].position;
                Vector3 desiredVelocity = Mathf.Min(maxVelocity, targetOffset.magnitude / radiusArrival) * targetOffset.normalized;
                forces[i] += kArrival * (desiredVelocity - velocity[i]).normalized;
            }
        }

        for (int i = 0; i < num; i++)
        {
            if (forces[i].magnitude > maxForce)
                forces[i] *= maxForce / forces[i].magnitude;
            boids[i].position += Time.deltaTime * velocity[i] + Time.deltaTime * Time.deltaTime / 2 * forces[i];
            velocity[i] += Time.deltaTime * forces[i];
            if (velocity[i].magnitude > maxVelocity)
                velocity[i] *= maxVelocity / velocity[i].magnitude;

            boids[i].forward = velocity[i];
            //boids[i].forward = boids[i].right;
            //boids[i].LookAt(velocity[i], Vector3.right);
        }
    }
*/