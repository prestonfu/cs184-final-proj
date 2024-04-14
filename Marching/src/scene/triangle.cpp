#include "triangle.h"

#include "CGL/CGL.h"
#include "GL/glew.h"

namespace CGL {
namespace SceneObjects {

Triangle::Triangle(const Mesh *mesh, size_t v1, size_t v2, size_t v3) {
  p1 = mesh->positions[v1];
  p2 = mesh->positions[v2];
  p3 = mesh->positions[v3];
  n1 = mesh->normals[v1];
  n2 = mesh->normals[v2];
  n3 = mesh->normals[v3];
  bbox = BBox(p1);
  bbox.expand(p2);
  bbox.expand(p3);

  bsdf = mesh->get_bsdf();
}

BBox Triangle::get_bbox() const { return bbox; }

bool checker(double x, double y, double z, const Ray &r) {
    if(x >= r.min_t && x <= r.max_t) {
      if(y >= 0 && z >= 0 && (1 - y - z) >= 0) {
        if(y <= 1 && z <= 1 && (1 - y - z) <= 1) {
          return true;
        }
      }
    }
  return false;
}

bool Triangle::has_intersection(const Ray &r) const {
  // Part 1, Task 3: implement ray-triangle intersection
  // The difference between this function and the next function is that the next
  // function records the "intersection" while this function only tests whether
  // there is a intersection.
  double first_val = dot(cross(r.d, p3 - p1), p2 - p1);
  double second_val = dot(cross(r.o - p1, p2 - p1), p3 - p1);
  double third_val = dot(cross(r.d, p3 - p1), r.o - p1);
  double fourth_val = dot(cross(r.o-p1, p2 - p1), r.d);
  Vector3D vect(0,0,0);
  Vector3D adder(second_val,third_val,fourth_val);
  vect += adder;
  vect = vect / first_val;

    if(checker(vect.x, vect.y, vect.z, r)) {
      return true;
    }
    return false;

  }

bool Triangle::intersect(const Ray &r, Intersection *isect) const {
  // Part 1, Task 3:
  // implement ray-triangle intersection. When an intersection takes
  // place, the Intersection data should be updated accordingly

  double first_val = dot(cross(r.d, p3 - p1), p2 - p1);
  double second_val = dot(cross(r.o - p1, p2 - p1), p3 - p1);
  double third_val = dot(cross(r.d, p3 - p1), r.o - p1);
  double fourth_val = dot(cross(r.o-p1, p2 - p1), r.d);
  Vector3D vect(0,0,0);
  Vector3D adder(second_val,third_val,fourth_val);
  vect += adder;
  vect = vect / first_val;

  if(checker(vect.x, vect.y, vect.z, r)) {
    //cout << "check!" << endl;
      r.max_t = vect.x;
      isect->t = vect.x;
      isect->n = n1 * (1 - vect.y - vect.z) + n2 * vect.y + n3 * vect.z;
      isect->primitive = this;
      isect->bsdf = get_bsdf();
      return true;
  }
    
  
  return false;


}



void Triangle::draw(const Color &c, float alpha) const {
  glColor4f(c.r, c.g, c.b, alpha);
  glBegin(GL_TRIANGLES);
  glVertex3d(p1.x, p1.y, p1.z);
  glVertex3d(p2.x, p2.y, p2.z);
  glVertex3d(p3.x, p3.y, p3.z);
  glEnd();
}

void Triangle::drawOutline(const Color &c, float alpha) const {
  glColor4f(c.r, c.g, c.b, alpha);
  glBegin(GL_LINE_LOOP);
  glVertex3d(p1.x, p1.y, p1.z);
  glVertex3d(p2.x, p2.y, p2.z);
  glVertex3d(p3.x, p3.y, p3.z);
  glEnd();
}

} // namespace SceneObjects
} // namespace CGL
