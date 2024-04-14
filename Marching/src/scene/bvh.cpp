#include "bvh.h"

#include "CGL/CGL.h"
#include "scene/primitive.h"
#include "triangle.h"

#include <cmath>
#include <iostream>
#include <stack>
#include <vector>

using namespace std;

namespace CGL {
namespace SceneObjects {

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {

  primitives = std::vector<Primitive *>(_primitives);
  root = construct_bvh(primitives.begin(), primitives.end(), max_leaf_size);
}

BVHAccel::~BVHAccel() {
  if (root)
    delete root;
  primitives.clear();
}

BBox BVHAccel::get_bbox() const { return root->bb; }

void BVHAccel::draw(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->draw(c, alpha);
    }
  } else {
    draw(node->l, c, alpha);
    draw(node->r, c, alpha);
  }
}

void BVHAccel::drawOutline(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->drawOutline(c, alpha);
    }
  } else {
    drawOutline(node->l, c, alpha);
    drawOutline(node->r, c, alpha);
  }
}

BVHNode *BVHAccel::construct_bvh(std::vector<Primitive *>::iterator start,
                                 std::vector<Primitive *>::iterator end,
                                 size_t max_leaf_size) {

  // TODO (Part 2.1):
  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.


  BBox bbox;

  for (auto p = start; p != end; p++) {
    BBox bb = (*p)->get_bbox();
    bbox.expand(bb);
  }

  BVHNode *node = new BVHNode(bbox);

  if(max_leaf_size >= (end - start)) {
    node->start = start;
    node->end = end;
    node->l = NULL;
    node->r = NULL;
    return node;
  }
  int best_axis = 0;

  Vector3D summer(0,0,0);
  auto it = start;
  float best_heuristic = INFINITY;
  for (auto p = start; p != end; p++) {
    summer += (*p)->get_bbox().centroid();
  }
  summer = summer / (end - start);

  int num_axes = 3;
  for (int i = 0; i < num_axes; i++) {
    std::list<Primitive*> left_values = {};
    std::list<Primitive*> right_values = {};

    for (auto p = start; p != end; p++) {
      double value = (*p)->get_bbox().centroid()[i];
      if (value <= summer[i]) {
        left_values.push_back(*p);
      } else {
        right_values.push_back(*p);
      }
    }

    vector<Primitive*> left_prim;
    vector<Primitive*> right_prim;


    BBox left_box;
    BBox right_box;

    for (Primitive* v_1: left_values) {
      left_prim.push_back(v_1);
      left_box.expand(v_1->get_bbox());
    }

    for (Primitive* v_2: right_values) {
      right_prim.push_back(v_2);
      right_box.expand(v_2->get_bbox());
    }

    float first_left_val = left_box.extent.x * left_box.extent.y + left_box.extent.y * left_box.extent.z
     + left_box.extent.z * left_box.extent.x;

    float second_right_val = right_box.extent.x * right_box.extent.y + right_box.extent.y * right_box.extent.z
     + right_box.extent.z * right_box.extent.x;

    float final_val = first_left_val * left_prim.size() + second_right_val * right_prim.size();

    if (final_val < best_heuristic) {
      best_heuristic = final_val;
      best_axis = i;
    }

  }

  vector<Primitive*> left_prim;
  vector<Primitive*> right_prim;

  std::list<Primitive*> left_values = {};
  std::list<Primitive*> right_values = {};

  for (auto p = start; p != end; p++) {
    if ((*p)->get_bbox().centroid()[best_axis] <= summer[best_axis]) {
      left_values.push_back(*p);
    } else {
      right_values.push_back(*p);
    }
  }

  for (Primitive* v_1: left_values) {
      left_prim.push_back(v_1);
    }

    for (Primitive* v_2: right_values) {
      right_prim.push_back(v_2);
    }

  //now we check where to "divide" up the start and end for the two recursive calls
  vector<Primitive *>::iterator left = start;
  vector<Primitive *>::iterator right = start;

  int size_of_left = left_prim.size();
  int size_of_right = right_prim.size();

  for (int i = 0; i < (end - start); i++) {
    if (i < size_of_left) {
      *left = left_prim[i];
      left = left + 1;
      right = right + 1;
    } else {
      *right = right_prim[i - size_of_left];
      right = right + 1;
    }
  }

  node->l = construct_bvh(start, left, max_leaf_size);
  node->r = construct_bvh(left, end, max_leaf_size);




  return node;


}

bool BVHAccel::has_intersection(const Ray &ray, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.
  // Take note that this function has a short-circuit that the
  // Intersection version cannot, since it returns as soon as it finds
  // a hit, it doesn't actually have to find the closest hit.


  if (!node->bb.intersect(ray, ray.min_t, ray.max_t)) {
    return false;
  }

  if(!node->isLeaf()) {
    bool left = has_intersection(ray, node->l);
    bool right = has_intersection(ray, node->r);
    return left || right;
  }

  for (auto p = node->start; p != node->end; p++) {
    total_isects++;
    if ((*p)->has_intersection(ray))
      return true;
  }
  return false;


}

bool BVHAccel::intersect(const Ray &ray, Intersection *i, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.

  if (!node->bb.intersect(ray, ray.min_t, ray.max_t)) {
    return false;
  }

  if(!node->isLeaf()) {
    bool left = intersect(ray, i, node->l);
    bool right = intersect(ray, i, node->r);
    return left || right;
  }



  bool hit = false;
  for (auto p = node->start; p != node->end; p++) {
    total_isects++;
    hit = (*p)->intersect(ray, i) || hit;
  }
  return hit;


}

} // namespace SceneObjects
} // namespace CGL
