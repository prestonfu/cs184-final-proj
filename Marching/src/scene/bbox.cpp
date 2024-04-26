#include "bbox.h"

#include "GL/glew.h"
#include "vector3D.h"

#include <algorithm>
#include <iostream>

namespace CGL {

bool BBox::intersect(const Ray& r, double& t0, double& t1) const {

  // TODO (Part 2.2):
  // Implement ray - bounding box intersection test
  // If the ray intersected the bouding box within the range given by
  // t0, t1, update t0 and t1 with the new intersection times.

  //The three vectors
  Vector3D xy(0,0,1);
  Vector3D yz(1,0,0);
  Vector3D xz(0,1,0);

  //values
  double xy0 = dot((min - r.o), xy)/dot(r.d, xy);
  double xy1 = dot((max - r.o), xy)/dot(r.d, xy);

  double yz0 = dot((min - r.o), yz)/dot(r.d, yz);
  double yz1 = dot((max - r.o), yz)/dot(r.d, yz);

  double xz0 = dot((min - r.o), xz)/dot(r.d, xz);
  double xz1 = dot((max - r.o), xz)/dot(r.d, xz);

  //the swaps
  if(xz0 > xz1) {
    double temp = xz0;
    xz0 = xz1;
    xz1 = temp;
  }

  if(yz0 > yz1) {
    double temp = yz0;
    yz0 = yz1;
    yz1 = temp;
  }

  if(xy0 > xy1) {
    double temp = xy0;
    xy0 = xy1;
    xy1 = temp;
  }


  //the intersections
  double max_t = std::min(std::min(xz1, yz1), xy1);
  double min_t = std::max(std::max(xz0, yz0), xy0);

  if (min_t <= max_t) {
    if ((min_t < t1) && (max_t >t0)) {
          return true;
        }
  }


  return false;

}

void BBox::draw(Color c, float alpha) const {

  glColor4f(c.r, c.g, c.b, alpha);

  // top
  glBegin(GL_LINE_STRIP);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(max.x, max.y, max.z);
  glEnd();

  // bottom
  glBegin(GL_LINE_STRIP);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, min.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glEnd();

  // side
  glBegin(GL_LINES);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(min.x, min.y, max.z);
  glEnd();

}

std::ostream& operator<<(std::ostream& os, const BBox& b) {
  return os << "BBOX(" << b.min << ", " << b.max << ")";
}

} // namespace CGL
