#include "pathtracer.h"

#include "misc.h"
#include "scene/light.h"
#include "scene/sphere.h"
#include "scene/triangle.h"
#include "util/random_util.h"
#include "vector3D.h"


using namespace CGL::SceneObjects;

namespace CGL {

PathTracer::PathTracer() {
  gridSampler = new UniformGridSampler2D();
  hemisphereSampler = new UniformHemisphereSampler3D();

  tm_gamma = 2.2f;
  tm_level = 1.0f;
  tm_key = 0.18;
  tm_wht = 5.0f;
}

PathTracer::~PathTracer() {
  delete gridSampler;
  delete hemisphereSampler;
}

void PathTracer::set_frame_size(size_t width, size_t height) {
  sampleBuffer.resize(width, height);
  sampleCountBuffer.resize(width * height);
}

void PathTracer::clear() {
  bvh = NULL;
  scene = NULL;
  camera = NULL;
  sampleBuffer.clear();
  sampleCountBuffer.clear();
  sampleBuffer.resize(0, 0);
  sampleCountBuffer.resize(0, 0);
}

void PathTracer::write_to_framebuffer(ImageBuffer &framebuffer, size_t x0,
                                      size_t y0, size_t x1, size_t y1) {
  sampleBuffer.toColor(framebuffer, x0, y0, x1, y1);
}

Vector3D
PathTracer::estimate_direct_lighting_hemisphere(const Ray &r,
                                                const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // For this function, sample uniformly in a hemisphere.

  // Note: When comparing Cornel Box (CBxxx.dae) results to importance sampling, you may find the "glow" around the light source is gone.
  // This is totally fine: the area lights in importance sampling has directionality, however in hemisphere sampling we don't model this behaviour.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d);

  // This is the same number of total samples as
  // estimate_direct_lighting_importance (outside of delta lights). We keep the
  // same number of samples for clarity of comparison.
  //cout << "check" << endl;
  int num_samples = scene->lights.size() * ns_area_light;
  Vector3D L_out(0,0,0);

  double float_val;
  Vector3D w;

  // TODO (Part 3): Write your sampling loop here
  // TODO BEFORE YOU BEGIN
  // UPDATE `est_radiance_global_illumination` to return direct lighting instead of normal shading 
  

  for(int i = 0; i <= num_samples; i++) {

    Vector3D bsdf = isect.bsdf->sample_f(w_out, &w, &float_val);
    //we want the depth to be equal to 1, from a TA
    int depth = 1;

    Intersection i_2;

    Ray ray_2 = Ray(hit_p,o2w * w, depth);
    ray_2.min_t = EPS_F;
    ray_2.max_t = 100000;

    if(bvh->intersect(ray_2, &i_2) &&  w.z > 0) {
      L_out = L_out + (i_2.bsdf->get_emission() * bsdf * w.z) / float_val;
    }


  }
  return  L_out / num_samples / 4;
  //cout << L_out << endl;
  return L_out;

}

Vector3D
PathTracer::estimate_direct_lighting_importance(const Ray &r,
                                                const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // To implement importance sampling, sample only from lights, not uniformly in
  // a hemisphere.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d);
  Vector3D L_out(0,0,0);

  auto end = scene->lights.end();



  for(auto begin = scene->lights.begin(); begin != end; begin++) {
    bool is_light = (*begin)->is_delta_light();
    Vector3D summer(0,0,0);
    int area = ns_area_light;
    Vector3D w;
    double distance;
    double p;

    if(is_light) {
      //cout << "checker" << endl;
      area = 1;
    }
    //cout << area << endl;
    for (int i = 0; i < area; i++) {
      
      Vector3D spec = (*begin)->sample_L(hit_p, &w, &distance, &p);
      if (dot(w, isect.n) > 0) {
        Ray ray = Ray(hit_p, w, distance - EPS_F);
        ray.depth = max_ray_depth;
        Intersection intersection;
        ray.min_t = EPS_F;
        if (!bvh->intersect(ray, &intersection)) {
          summer = summer + (spec * isect.bsdf->f(w_out, w) * dot(w, isect.n)) / p;
          //cout << summer << endl;
        }
      }
    }

    L_out = L_out  + summer/area;


  }
  return L_out/4;


  return Vector3D(1.0);

}

Vector3D PathTracer::zero_bounce_radiance(const Ray &r,
                                          const Intersection &isect) {
  // TODO: Part 3, Task 2
  // Returns the light that results from no bounces of light

  Vector3D emission = isect.bsdf->get_emission();
  return emission;


}

Vector3D PathTracer::one_bounce_radiance(const Ray &r,
                                         const Intersection &isect) {
  // TODO: Part 3, Task 3
  // Returns either the direct illumination by hemisphere or importance sampling
  // depending on `direct_hemisphere_sample`

  if (direct_hemisphere_sample) {
    //cout << "check" << endl;
    return estimate_direct_lighting_hemisphere(r, isect);
  }
  //cout << "checkerrr" << endl;
  return estimate_direct_lighting_importance(r, isect);


  return Vector3D(1.0);


}

Vector3D PathTracer::at_least_one_bounce_radiance(const Ray &r,
                                                  const Intersection &isect) {
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  Vector3D hit_p = r.o + r.d * isect.t;
  Vector3D w_out = w2o * (-r.d);

  Vector3D L_out(0, 0, 0);

  // TODO: Part 4, Task 2
  // Returns the one bounce radiance + radiance from extra bounces at this point.
  // Should be called recursively to simulate extra bounces.
  double value;
  Vector3D in_w;
  Vector3D bsdf = isect.bsdf->sample_f(w_out, &in_w, &value);

  int the_depth = r.depth;
  bool coin_flips = coin_flip(0.35);

  if (the_depth == 0) {
    return L_out;
  }


  //cout << "asdfasdf: " << the_depth << endl;
  //isAccumBounces = 1;
  Ray new_ray = Ray(hit_p + EPS_F * o2w*in_w, o2w * in_w);
  new_ray.depth = r.depth - 1;
  new_ray.min_t = EPS_F;
  new_ray.max_t = 99999;
  double value_2 = dot(o2w * in_w, isect.n);
  Intersection intersect;

  if(value == 0) {
    return L_out;
  }
  

  bool intersect_bool = bvh->intersect(new_ray, &intersect);
  //cout << "intersect_bool: " << intersect_bool << endl;

  if (intersect_bool) {
    L_out += one_bounce_radiance(r, isect);
  } else {
    L_out = Vector3D(0,0,0);
    return L_out;
  }
  if(intersect_bool && value_2 != 0 && coin_flips && new_ray.depth > 0) {
    

    //cout << "the depth: " << the_depth << endl;
    //cout << "the max_ray_depth: " << max_ray_depth << endl;

    if(isAccumBounces) {
      //cout << "check" << endl;
      Vector3D new_one = at_least_one_bounce_radiance(new_ray, intersect) * bsdf * value_2 / value;
      L_out += (new_one);
      //L_out = L_out; 
      return L_out;
    } else {

      //cout << "check" << endl;
      L_out = L_out + 0;
      return L_out;
    }
  }


  return L_out;

}

Vector3D PathTracer::est_radiance_global_illumination(const Ray &r) {
  Intersection isect;
  Vector3D L_out;

  // You will extend this in assignment 3-2.
  // If no intersection occurs, we simply return black.
  // This changes if you implement hemispherical lighting for extra credit.

  // The following line of code returns a debug color depending
  // on whether ray intersection with triangles or spheres has
  // been implemented.
  //
  // REMOVE THIS LINE when you are ready to begin Part 3.
  
  if (!bvh->intersect(r, &isect))
    return L_out;


  //L_out = (isect.t == INF_D) ? debug_shading(r.d) : normal_shading(isect.n);

  // TODO (Part 3): Return the direct illumination.

  // TODO (Part 4): Accumulate the "direct" and "indirect"
  // parts of global illumination into L_out rather than just direct
  /*
  if(r.depth == 0) {
    return zero_bounce_radiance(r, isect);
  } else if (r.depth == 1) {
    return one_bounce_radiance(r, isect);
  }
  */
  //cout << "check" << endl;
  //cout << isAccumBounces << endl;
  if (!isAccumBounces && max_ray_depth != 0 && 1 == 0) {
    L_out = at_least_one_bounce_radiance(r, isect);
    return L_out;
  }
  L_out = L_out + zero_bounce_radiance(r, isect) + one_bounce_radiance(r, isect);

  return L_out;
}

void PathTracer::raytrace_pixel(size_t x, size_t y) {
  // TODO (Part 1.2):
  // Make a loop that generates num_samples camera rays and traces them
  // through the scene. Return the average Vector3D.
  // You should call est_radiance_global_illumination in this function.
  double value_1 = 0;
  double value_2 = 0;
  int actual_count_samples = 0;
  // TODO (Part 5):
  // Modify your implementation to include adaptive sampling.
  // Use the command line parameters "samplesPerBatch" and "maxTolerance"
  int num_samples = ns_aa;          // total samples to evaluate

  Vector2D origin = Vector2D(x, y); // bottom left corner of the pixel
  Vector3D avg(0,0,0);
  int num_batches = 0;

  for (double i = 0; i < num_samples; i++) {


    if (num_batches == samplesPerBatch) {
      double the_value_1 = value_1 / i;

      double firstly = 1.96 * sqrt(((1.0 / (i - 1.0)) * (value_2 - (value_1 * value_1) / i)) / i);
      num_batches = 0;
      if (firstly <= maxTolerance * the_value_1) {
        actual_count_samples = i;
        //cout << "check";
        break;
      }
    }

    Vector2D random_sample = gridSampler->get_sample() + origin;
    double value_x = random_sample.x / sampleBuffer.w;
    double value_y = random_sample.y /  sampleBuffer.h;
    Ray ray = camera->generate_ray(value_x, value_y);
    ray.depth = max_ray_depth;
    avg += est_radiance_global_illumination(ray);

    double ilum = est_radiance_global_illumination(ray).illum();

    value_1 += ilum;
    value_2 += ilum * ilum;

    actual_count_samples = actual_count_samples + 1;
    num_batches = num_batches + 1;

  }
  avg = avg / (double) actual_count_samples;
  //avg = avg * num_samples;


  sampleBuffer.update_pixel(avg, x, y);
  sampleCountBuffer[x + y * sampleBuffer.w] = actual_count_samples;


}

void PathTracer::autofocus(Vector2D loc) {
  Ray r = camera->generate_ray(loc.x / sampleBuffer.w, loc.y / sampleBuffer.h);
  Intersection isect;

  bvh->intersect(r, &isect);

  camera->focalDistance = isect.t;
}

} // namespace CGL
