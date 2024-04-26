#define RADIUS_SEPARATION 0.5
#define RADIUS_COHESION 0.5
#define RADIUS_ALIGNMENT 0.5
#define RADIUS_CONTAINMENT 0.5

#define K_SEPARATION 0.5

kernel void calculate(const int n, global float* positions, global float* velocities, global float* accelerations, const float deltaTime)
{
    int id = get_global_id(0);
    float2 pos = (float2)(positions[2 * id], positions[2 * id + 1]);
    float2 vel = (float2)(velocities[2 * id], velocities[2 * id + 1]);

    float2 accel = (float2)(0.0, 0.0);
    float2 separation = (float2)(0.0, 0.0);
    // float2 cohesion = (float2)(0.0, 0.0);
    // float2 alignment = (float2)(0.0, 0.0);

    for (int i = 0; i < n; i++)
    {
        if (i == id)
            continue;
        float2 otherPos = (float2)(positions[2 * i], positions[2 * i + 1]);
        float distance = length(pos - otherPos);
        if (distance < RADIUS_SEPARATION) {
            float2 displacement = pos - otherPos;
            separation += displacement / distance / distance;
        }
        // if (distance < RADIUS_COHESION) {
        //     cohesion += 
        // }
    }

    accel += K_SEPARATION * normalize(separation);
    accelerations[2 * id] = accel.x;
    accelerations[2 * id + 1] = accel.y;
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