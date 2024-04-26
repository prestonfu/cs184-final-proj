#include "sphere.h"

#include <algorithm>
#include <cmath>
#include <ostream>

#include "pathtracer/bsdf.h"
#include "util/sphere_drawing.h"
#include "vector3D.h"

namespace CGL {
namespace SceneObjects {

const int MAX_STEPS = 50;
const float MAX_DIST = 1000.0;
const float MIN_DIST = 0.05;


// distance between two vectors
float dist(Vector3D first, Vector3D second) {
  return sqrt((first.x - second.x)*(first.x - second.x) + (first.y - second.y)*(first.y - second.y) 
  + (first.z - second.z)*(first.z - second.z));
}

// signed distance between a ray endpoint and a sphere's radius (value is negative if ray endpoint is within sphere)
float Sphere::sdf_sphere(Vector3D vector, Vector3D center, float radius) const {
  return dist(vector, center) - radius;
}


// raymarching algorithm
float Sphere::raymarching(const Ray &r, Intersection *i) const {

  float distance_traveled = 0.0;

  for (int i = 0; i < MAX_STEPS; i++) {
    Vector3D ray = r.o + distance_traveled * r.d;

    float dist_to_object = sdf_sphere(ray, o, o.r);

    // if we get very close to our object
    if (dist_to_object < MIN_DIST) {
      return dist_to_object;
    }

    //if we don't git our object and therefore are very far away from our target
    if (dist_to_object > MAX_DIST) {
      return dist_to_object;
    }

  }

  //if we don't intersect
  return MAX_DIST;
}


// IGNORE
bool Sphere::test(const Ray &r, double &t1, double &t2) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.


  return true;

}


//IGNORE 
bool checker(const Ray &r, Vector3D vect, double val) {
  if (2.0 * dot((r.o - vect), r.d) * 2.0 * dot((r.o - vect), r.d) - 4.0 * r.d.norm2() * ((r.o - vect).norm2() - val) >= 0) {
    return true;
  }
  return false;
}


// IGNORE
bool Sphere::has_intersection(const Ray &r) const {
  double first_val = r.min_t - 10;
  double second_val = r.max_t + 10;
  if (checker(r, o, r2)) {
    second_val = sqrt(2.0 * dot((r.o - o), r.d) * 2.0 * dot((r.o - o), r.d) - 4.0 * r.d.norm2() * ((r.o - o).norm2() - r2))
     - 2.0 * dot((r.o - o), r.d) ;
    second_val = second_val / (2.0 * r.d.norm2());

    first_val =  -1.0 * (sqrt(2.0 * dot((r.o - o), r.d) * 2.0 * dot((r.o - o), r.d) - 4.0 * r.d.norm2() * ((r.o - o).norm2() - r2)))
     - 2.0 * dot((r.o - o), r.d) ;
    first_val = first_val / (2.0 * r.d.norm2());
  } 
  double real_val = min(first_val, second_val);
  cout << "check" << endl;
  if (((real_val >= r.min_t) && (real_val <= r.max_t))) {
    return true;
  }
  return false;
  

  
  
  // TODO (Part 1.4):
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.
}

bool Sphere::intersect(const Ray &r, Intersection *i) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.

  //cout << "check" << endl;

  float value = raymarching(r, i);
  if (value >= MAX_DIST) {
    return false;
  } else {
    r.max_t = value;
    i->t = value;
    Vector3D val = r.o + value * r.d - o;
    val.normalize();
    i->n = val;
    i->primitive = this;
    i->bsdf = this->get_bsdf();
    return true;

  }

  /*

  double first_val = r.min_t - 10;
  double second_val = r.max_t + 10;
  //cout << "hello" << endl;
  if (checker(r, o, r2)) {
    //cout << "hello" << endl;
    second_val = sqrt(2.0 * dot((r.o - o), r.d) * 2.0 * dot((r.o - o), r.d) - 4.0 * r.d.norm2() * ((r.o - o).norm2() - r2))
     - 2.0 * dot((r.o - o), r.d) ;
    second_val = second_val / (2.0 * r.d.norm2());

    first_val =  -1.0 * (sqrt(2 * dot((r.o - o), r.d) * 2.0 * dot((r.o - o), r.d) - 4.0 * r.d.norm2() * ((r.o - o).norm2() - r2)))
     - 2.0 * dot((r.o - o), r.d) ;
    first_val = first_val / (2.0 * r.d.norm2());

    double real_val = min(first_val, second_val);
    //cout << "check" << endl;
    if (real_val == second_val) {
      //cout << "checker" << endl;
    }

    if((real_val >= r.min_t) && (real_val <= r.max_t)) {
        r.max_t = real_val;
        i->t = real_val;
        Vector3D val = r.o + first_val * r.d - o;
        val.normalize();
        i->n = val;
        i->primitive = this;
        i->bsdf = this->get_bsdf();
        return true;
    } else if ((second_val >= r.min_t) && (second_val <= r.max_t)) {
        cout << "check" << endl;
        r.max_t = second_val;
        i->t = second_val;
        Vector3D val = r.o + second_val * r.d - o;
        val.normalize();
        i->n = val;
        i->primitive = this;
        i->bsdf = this->get_bsdf();
        return true;
    }
    return false;
  }
  return false;

  */

}

void Sphere::draw(const Color &c, float alpha) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color &c, float alpha) const {
  // Misc::draw_sphere_opengl(o, r, c);
}

} // namespace SceneObjects
} // namespace CGL
