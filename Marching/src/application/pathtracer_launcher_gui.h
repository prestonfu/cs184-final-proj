#pragma once
#include "GLFW/glfw3.h"
#include <cstddef>
#include <cstring>
#include <string>
#include <sys/stat.h>
namespace PathtracerLauncherGUI {
/**
 * @brief Container for all settings tweakable from the GUI. The GUI settings are then read by the 
 * main application and used to configure the pathtracer.
 * The settings can be serialized via serialize() and deserialize().
 */
struct GUISettings {
  GUISettings() {
    pathtracer_ns_aa = 1;
    pathtracer_max_ray_depth = 1;
    pathtracer_ns_area_light = 1;

    pathtracer_ns_diff = 1;
    pathtracer_ns_glsy = 1;
    pathtracer_ns_refr = 1;

    pathtracer_num_threads = 1;

    pathtracer_samples_per_patch = 32;
    pathtracer_max_tolerance = 0.05f;
    pathtracer_direct_hemisphere_sample = false;

    pathtracer_lensRadius = 0.0;
    pathtracer_focalDistance = 4.7;
  }

  void serialize(const std::string &a_file_path);

  void deserialize(const std::string &a_file_path);

  size_t pathtracer_ns_aa;
  size_t pathtracer_max_ray_depth;
  size_t pathtracer_ns_area_light;

  size_t pathtracer_ns_diff;
  size_t pathtracer_ns_glsy;
  size_t pathtracer_ns_refr;

  size_t pathtracer_num_threads;

  float pathtracer_max_tolerance;
  size_t pathtracer_samples_per_patch;


  bool pathtracer_direct_hemisphere_sample;
  bool render_custom_region = false;
  bool pathtracer_accumulate_bounces = true;
  double pathtracer_lensRadius;
  double pathtracer_focalDistance;

  size_t settings_window_width = 600;
  size_t settings_window_height = 800;

  // settings but not for the app.
  bool write_to_file = false;
  int w = 1024, h = 768, x = -1, y = 0, dx = 0, dy = 0;
  std::string output_file_name, cam_settings = "";
  std::string scene_file_path = "";
};
void render_loop(GLFWwindow *a_window, GUISettings &a_settings);
int draw(GUISettings &a_settings);

bool file_exists(const std::string &name);
bool dae_exists(const std::string &name);
};
