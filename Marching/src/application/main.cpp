#include "CGL/CGL.h"
#include "CGL/viewer.h"

#define TINYEXR_IMPLEMENTATION
#include "CGL/tinyexr.h"

#include "application.h"
typedef uint32_t gid_t;
#include "util/image.h"
typedef uint32_t gid_t;

#include <iostream>
#ifdef _WIN32
#include "util/win32/getopt.h"
#else
#include <unistd.h>
#endif
#include "pathtracer_launcher_gui.h"

using namespace std;
using namespace CGL;

#define msg(s) cerr << "[PathTracer] " << s << endl;

void usage(const char *binaryName) {
  printf("Usage: %s [options] <scenefile>\n", binaryName);
  printf("Program Options:\n");
  printf("  -s  <INT>        Number of camera rays per pixel\n");
  printf("  -l  <INT>        Number of samples per area light\n");
  printf("  -t  <INT>        Number of render threads\n");
  printf("  -m  <INT>        Maximum ray depth\n");
  printf("  -o  <INT>        Accumulate Bounces of Light \n");
  printf("  -e  <PATH>       Path to environment map\n");
  printf("  -b  <FLOAT>      The size of the aperture\n");
  printf("  -d  <FLOAT>      The focal distance\n");
  printf("  -f  <FILENAME>   Image (.png) file to save output to in windowless "
         "mode\n");
  printf(
      "  -r  <INT> <INT>  Width and height of output image (if windowless)\n");
  printf("  -h               Print this help message\n");
  printf("\n");
}

HDRImageBuffer *load_exr(const char *file_path) {

  const char *err;

  EXRImage exr;
  InitEXRImage(&exr);

  int ret = ParseMultiChannelEXRHeaderFromFile(&exr, file_path, &err);
  if (ret != 0) {
    msg("Error parsing OpenEXR file: " << err);
    return NULL;
  }

  for (int i = 0; i < exr.num_channels; i++) {
    if (exr.pixel_types[i] == TINYEXR_PIXELTYPE_HALF) {
      exr.requested_pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT;
    }
  }

  ret = LoadMultiChannelEXRFromFile(&exr, file_path, &err);
  if (ret != 0) {
    msg("Error loading OpenEXR file: " << err);
    exit(EXIT_FAILURE);
  }

  HDRImageBuffer *envmap = new HDRImageBuffer();
  envmap->resize(exr.width, exr.height);
  float *channel_r = (float *)exr.images[2];
  float *channel_g = (float *)exr.images[1];
  float *channel_b = (float *)exr.images[0];
  for (size_t i = 0; i < exr.width * exr.height; i++) {
    envmap->data[i] = Vector3D(channel_r[i], channel_g[i], channel_b[i]);
  }

  return envmap;
}

int main(int argc, char **argv) {

  // get the options
  AppConfig config;
  int opt;
  bool write_to_file = false;
  size_t w = 0, h = 0, x = -1, y = 0, dx = 0, dy = 0;
  string output_file_name, cam_settings = "";
  string sceneFilePath;
  if (argc == 1) { // no argument specifiers, launch GUI to get the settings
#define SETTINGSFILE_PATH "settings.txt"
    PathtracerLauncherGUI::GUISettings settings;
    settings.deserialize(SETTINGSFILE_PATH);
    PathtracerLauncherGUI::draw(settings);
    // done drawing
    { // extract settings that don't belong to config
      w = settings.w;
      h = settings.h;
      x = settings.x;
      y = settings.y;
      dx = settings.dx;
      dy = settings.dy;
      write_to_file = settings.write_to_file;
      sceneFilePath = settings.scene_file_path;
      output_file_name = settings.output_file_name;
      cam_settings = settings.cam_settings;
    }
    { // extract settings that belong to config
      config.pathtracer_direct_hemisphere_sample =
          settings.pathtracer_direct_hemisphere_sample;
      config.pathtracer_focalDistance = settings.pathtracer_focalDistance;
      config.pathtracer_lensRadius = settings.pathtracer_lensRadius;
      config.pathtracer_max_ray_depth = settings.pathtracer_max_ray_depth;
      config.pathtracer_max_tolerance = settings.pathtracer_max_tolerance;
      config.pathtracer_ns_aa = settings.pathtracer_ns_aa;
      config.pathtracer_ns_area_light = settings.pathtracer_ns_area_light;
      config.pathtracer_ns_diff = settings.pathtracer_ns_diff;
      config.pathtracer_ns_glsy = settings.pathtracer_ns_glsy;
      config.pathtracer_ns_refr = settings.pathtracer_ns_refr;
      config.pathtracer_num_threads = settings.pathtracer_num_threads;
      config.pathtracer_samples_per_patch =
          settings.pathtracer_samples_per_patch;
      config.pathtracer_accumulate_bounces = settings.pathtracer_accumulate_bounces;
    }
  } else {
    while ((opt = getopt(argc, argv, "s:l:t:m:o:e:h:H:f:r:c:b:d:a:p:")) !=
           -1) { // for each option...
      switch (opt) {
      case 'f':
        write_to_file = true;
        output_file_name = string(optarg);
        break;
      case 'r':
        w = atoi(argv[optind - 1]);
        h = atoi(argv[optind]);
        optind++;
        break;
      case 'p':
        x = atoi(argv[optind - 1]);
        y = atoi(argv[optind - 0]);
        dx = atoi(argv[optind + 1]);
        dy = atoi(argv[optind + 2]);
        optind += 3;
        break;
      case 's':
        config.pathtracer_ns_aa = atoi(optarg);
        break;
      case 'l':
        config.pathtracer_ns_area_light = atoi(optarg);
        break;
      case 't':
        config.pathtracer_num_threads = atoi(optarg);
        break;
      case 'm':
        config.pathtracer_max_ray_depth = atoi(optarg);
        break;
      case 'o':
        config.pathtracer_accumulate_bounces = atoi(optarg) > 0;
        break;
      case 'e':
        std::cout << "[PathTracer] Loading environment map " << optarg
                  << std::endl;
        config.pathtracer_envmap = load_exr(optarg);
        break;
      case 'c':
        cam_settings = string(optarg);
        break;
      case 'b':
        config.pathtracer_lensRadius = atof(optarg);
        break;
      case 'd':
        config.pathtracer_focalDistance = atof(optarg);
        break;
      case 'a':
        config.pathtracer_samples_per_patch = atoi(argv[optind - 1]);
        config.pathtracer_max_tolerance = atof(argv[optind]);
        optind++;
        break;
      case 'H':
        config.pathtracer_direct_hemisphere_sample = true;
        optind--;
        break;
      default:
        usage(argv[0]);
        return 1;
      }
    }

    // print usage if no argument given
    if (optind >= argc) {
      usage(argv[0]);
      return 1;
    }

    sceneFilePath = argv[optind];
  }
  msg("Input scene file: " << sceneFilePath);
  string sceneFile = sceneFilePath.substr(sceneFilePath.find_last_of('/') + 1);
  sceneFile = sceneFile.substr(0, sceneFile.find(".dae"));
  config.pathtracer_filename = sceneFile;

  // parse scene
  Collada::SceneInfo *sceneInfo = new Collada::SceneInfo();
  if (Collada::ColladaParser::load(sceneFilePath.c_str(), sceneInfo) < 0) {
    delete sceneInfo;
    exit(0);
  }

  // create application
  Application *app = new Application(config, !write_to_file);

  msg("Rendering using " << config.pathtracer_num_threads << " threads");

  // write straight to file without opening a window if -f option provided
  if (write_to_file) {
    app->init();
    app->load(sceneInfo);
    delete sceneInfo;

    if (w && h)
      app->resize(w, h);

    if (cam_settings != "")
      app->load_camera(cam_settings);

    app->render_to_file(output_file_name, x, y, dx, dy);
    return 0;
  }

  // create viewer
  Viewer viewer = Viewer();

  // set renderer
  viewer.set_renderer(app);

  // init viewer
  viewer.init();

  // load scene
  app->load(sceneInfo);

  delete sceneInfo;

  if (w && h)
    viewer.resize(w, h);

  if (cam_settings != "")
    app->load_camera(cam_settings);

  // start viewer
  viewer.start();

  return 0;
}
