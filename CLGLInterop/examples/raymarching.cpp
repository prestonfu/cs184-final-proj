#include <glad/glad.h>
#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#else
#include <GL/gl.h>
#endif

#include <common/OpenCLUtil.h>

#ifdef OS_WIN
#define GLFW_EXPOSE_NATIVE_WIN32
#define GLFW_EXPOSE_NATIVE_WGL
#endif

#ifdef OS_LNX
#define GLFW_EXPOSE_NATIVE_X11
#define GLFW_EXPOSE_NATIVE_GLX
#endif

#include <GLFW/glfw3.h>
#include <GLFW/glfw3native.h>

#include <common/OpenGLUtil.h>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <map>
#include <vector>
#include <random>
#include <algorithm>

#include <common/read_npz.h>

using namespace std;
using namespace cl;

#define EPS_F (0.00001f)
#define MIN_R (0.1f)
#define MAX_R (5.0f)
#define LOCAL_WORK_SIZE 256
#define LOCAL_WORK_SIZE_X 16
#define LOCAL_WORK_SIZE_Y 16
#define SPHERE_RADIUS (0.02f)
#define SPHERE_COUNT 256

//#define PROFILE

typedef unsigned int uint;

static const uint NUM_JSETS = 9;

static const float matrix[16] =
{
    1.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 1.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 1.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 1.0f
};

static const float vertices[12] =
{
    -1.0f,-1.0f, 0.0,
     1.0f,-1.0f, 0.0,
     1.0f, 1.0f, 0.0,
    -1.0f, 1.0f, 0.0
};

static const float texcords[8] =
{
    0.0, 1.0,
    1.0, 1.0,
    1.0, 0.0,
    0.0, 0.0
};

static const uint indices[6] = {0,1,2,0,2,3};

static const float CJULIA[] = {
    -0.700f, 0.270f,
    -0.618f, 0.000f,
    -0.400f, 0.600f,
     0.285f, 0.000f,
     0.285f, 0.010f,
     0.450f, 0.143f,
    -0.702f,-0.384f,
    -0.835f,-0.232f,
    -0.800f, 0.156f,
     0.279f, 0.000f
};

static map<string, pair<std::vector<float>, std::vector<float>>> pointClouds;

static int wind_width = 1280;//640;
static int wind_height= 960;//480;
static int const nparticles = SPHERE_COUNT;
static int selectedIndex = -1;
static bool paused = false;

static std::vector<int> permutation(nparticles);
static std::vector<uint> codes(nparticles);
static std::vector<float> spheres(3 * nparticles);
static std::vector<float> bboxes(6 * nparticles / LOCAL_WORK_SIZE);

typedef struct {
    Device d;
    CommandQueue q;
    Program p;
    Kernel raytrace;
    Kernel raymarch;
    ImageGL tex;
    cl_float3 cameraPos;
    Buffer c2w;
    Buffer spheres;
    Buffer permutation;
    Buffer bboxes;
    Buffer seed;
    Buffer colors;

    Buffer targetPos;
    Buffer targetColor;

    Kernel kc; // calculate.cl kernel (calculate acceleration given positions and velocities)
    Kernel ki; // integrate.cl kernel (calculate position and velocity given acceleration)
    Buffer i; // position buffer
    Buffer v; // velocity buffer
    Buffer a; // acceleration buffer
    //Buffer 
} process_params;

typedef struct {
    GLuint prg;
    GLuint vao;
    GLuint tex;
} render_params;

typedef struct {
    float targetPos[3];
    float hFov;
    float vFov;
    float phi;
    float theta;
    float r;
} camera_state;

typedef struct {
    float x;
    float y;
    bool pressed; //left mouse button only
    bool pan;
} mouse_state;

process_params params;
render_params rparams;
camera_state cam;
mouse_state mouse;

static std::vector<string> pointCloudNames;

void loadPointCloud(int index)
{
    if (index == -1)
    {
        std::vector<float> data(3 * nparticles, 1);
        params.q.enqueueWriteBuffer(params.targetColor, CL_TRUE, 0, sizeof(float) * 3 * nparticles, data.data());
        cout << "model unloaded" << endl;
        return;
    }
    string name = pointCloudNames[index];
    if (pointClouds.count(name) == 0)
    {
        cout << name << " not found" << endl;
        return;
    }
    params.q.enqueueWriteBuffer(params.targetPos, CL_TRUE, 0, sizeof(float) * 3 * nparticles, pointClouds[name].first.data());
    params.q.enqueueWriteBuffer(params.targetColor, CL_TRUE, 0, sizeof(float) * 3 * nparticles, pointClouds[name].second.data());
    cout << name << " loaded" << endl;
    // std::vector<float> cloud = pointClouds[name];
    // for (int i = cloud.size() / 3 - 1; i > 0; i--)
    // {
    //     int j = rand() % i;
    //     swap(cloud[3 * i], cloud[3 * j]);
    //     swap(cloud[3 * i + 1], cloud[3 * j + 1]);
    //     swap(cloud[3 * i + 2], cloud[3 * j + 2]);
    // }
    //params.q.enqueueWriteBuffer(params.spheres, CL_TRUE, 0, sizeof(float) * 3 * nparticles, pointClouds[name].data());
}

void setScreenSize(int width, int height)
{
    // resize buffers and stuff
}

void computeCameraPosition()
{
    float sinPhi = sin(cam.phi);
    if (sinPhi == 0) {
        cam.phi += EPS_F;
        sinPhi = sin(cam.phi);
    }
    float dirToCamera[3] = {cam.r * sinPhi * sin(cam.theta),
                            cam.r * cos(cam.phi),
                            cam.r * sinPhi * cos(cam.theta)};

    for (int i = 0; i < 3; i++)
    {
        params.cameraPos.s[i] = cam.targetPos[i] + dirToCamera[i];
    }
    float upVec[3] = {0.0f, sinPhi > 0 ? 1.0f : -1.0f, 0.0f};
    float c2w[9];

    // Column 0 = cross(upVec, dirToCamera);
    c2w[0] = upVec[1] * dirToCamera[2] - upVec[2] * dirToCamera[1];
    c2w[3] = upVec[2] * dirToCamera[0] - upVec[0] * dirToCamera[2];
    c2w[6] = upVec[0] * dirToCamera[1] - upVec[1] * dirToCamera[0];
    float norm = sqrt(c2w[0] * c2w[0] + c2w[3] * c2w[3] + c2w[6] * c2w[6]);
    c2w[0] /= norm; c2w[3] /= norm; c2w[6] /= norm;

    // Column 1 = cross(dirToCamera, column 0)
    c2w[1] = dirToCamera[1] * c2w[6] - dirToCamera[2] * c2w[3];
    c2w[4] = dirToCamera[2] * c2w[0] - dirToCamera[0] * c2w[6];
    c2w[7] = dirToCamera[0] * c2w[3] - dirToCamera[1] * c2w[0];
    norm = sqrt(c2w[1] * c2w[1] + c2w[4] * c2w[4] + c2w[7] * c2w[7]);
    c2w[1] /= norm; c2w[4] /= norm; c2w[7] /= norm;

    // Column 2 = normalize(dirToCamera)
    norm = sqrt(dirToCamera[0] * dirToCamera[0] + dirToCamera[1] * dirToCamera[1] + dirToCamera[2] * dirToCamera[2]);
    c2w[2] = dirToCamera[0] / norm;
    c2w[5] = dirToCamera[1] / norm;
    c2w[8] = dirToCamera[2] / norm;

    params.q.enqueueWriteBuffer(params.c2w, CL_TRUE, 0, sizeof(float) * 9, c2w);
}

//https://developer.nvidia.com/blog/thinking-parallel-part-iii-tree-construction-gpu/
unsigned int expandBits(unsigned int v)
{
    v = (v * 0x00010001u) & 0xFF0000FFu;
    v = (v * 0x00000101u) & 0x0F00F00Fu;
    v = (v * 0x00000011u) & 0xC30C30C3u;
    v = (v * 0x00000005u) & 0x49249249u;
    return v;
}

// Calculates a 30-bit Morton code for the
// given 3D point located within the unit cube [0,1].
unsigned int morton3D(float x, float y, float z)
{
    x = min(max(x * 1024.0f, 0.0f), 1023.0f);
    y = min(max(y * 1024.0f, 0.0f), 1023.0f);
    z = min(max(z * 1024.0f, 0.0f), 1023.0f);
    unsigned int xx = expandBits((unsigned int)x);
    unsigned int yy = expandBits((unsigned int)y);
    unsigned int zz = expandBits((unsigned int)z);
    return xx * 4 + yy * 2 + zz;
}


static void glfw_error_callback(int error, const char* desc)
{
    fputs(desc,stderr);
}

static void glfw_cursor_position_callback(GLFWwindow* wind, double x, double y)
{
    if (mouse.pressed)
    {
        float dx = x - mouse.x;
        float dy = y - mouse.y;
        float dPhi = -dy * (M_PI / wind_height);
        float dTheta = -dx * (M_PI / wind_width);

        cam.phi = clamp(cam.phi + dPhi, 0.0f, (float)M_PI);
        cam.theta += dTheta;
        computeCameraPosition();
    }

    mouse.x = x;
    mouse.y = y;
}

static void glfw_mouse_button_callback(GLFWwindow* wind, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
        mouse.pressed = true;
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE)
        mouse.pressed = false;
}

static void glfw_scroll_callback(GLFWwindow* wind, double dx, double dy)
{
    if (fabs(dy) > EPS_F)
    {
        cam.r = clamp(cam.r + (float)dy, MIN_R, MAX_R);
        computeCameraPosition();
    }
}

static void glfw_key_callback(GLFWwindow* wind, int key, int scancode, int action, int mods)
{
    if (action == GLFW_PRESS) {
        int curIndex = selectedIndex;
        if (key == GLFW_KEY_ESCAPE)
            glfwSetWindowShouldClose(wind, GL_TRUE);
        else if (key == GLFW_KEY_0)
            selectedIndex = -1;
        else if (key == GLFW_KEY_1)
            selectedIndex = 0;
        else if (key == GLFW_KEY_2)
            selectedIndex = 1;
        else if (key == GLFW_KEY_3)
            selectedIndex = 2;
        else if (key == GLFW_KEY_4)
            selectedIndex = 3;
        else if (key == GLFW_KEY_5)
            selectedIndex = 4;
        else if (key == GLFW_KEY_6)
            selectedIndex = 5;
        else if (key == GLFW_KEY_7)
            selectedIndex = 6;
        else if (key == GLFW_KEY_8)
            selectedIndex = 7;
        else if (key == GLFW_KEY_9)
            selectedIndex = 8;
        else if (key == GLFW_KEY_Q)
            selectedIndex = 9;
        else if (key == GLFW_KEY_W)
            selectedIndex = 10;
        else if (key == GLFW_KEY_E)
            selectedIndex = 11;
        else if (key == GLFW_KEY_R)
            selectedIndex = 12;
        else if (key == GLFW_KEY_T)
            selectedIndex = 13;    
        else if (key == GLFW_KEY_Y)
            selectedIndex = 14;
        else if (key == GLFW_KEY_U)
            selectedIndex = 15;    
        else if (key == GLFW_KEY_P)
            paused = !paused;  
        
        if (selectedIndex != curIndex)
            loadPointCloud(selectedIndex);   
    }
}

static void glfw_framebuffer_size_callback(GLFWwindow* wind, int width, int height)
{
    glViewport(0,0,width,height);
    setScreenSize(width, height);
    //wind_width = width;
    //wind_height = height;
}

void processTimeStep(float);
void renderFrame(void);

inline float radians(float angle)
{
    return angle * M_PI / 180;
}

inline float degrees(float radians)
{
    return radians * 180 / M_PI;
}

//[l, r)
void bvh(std::vector<float> &spheres, std::vector<int> &permutation, int l, int r)
{
    if (r - l <= 256)
        return;

    float min_x = INFINITY, max_x = -INFINITY, min_y = INFINITY, max_y = -INFINITY, min_z = INFINITY, max_z = -INFINITY;
    for (int i = 0; i < spheres.size(); i += 3) {
        min_x = min(min_x, spheres[3 * i + 0]);
        min_y = min(min_y, spheres[3 * i + 1]);
        min_z = min(min_z, spheres[3 * i + 2]);
        max_x = max(max_x, spheres[3 * i + 0]);
        max_y = max(max_y, spheres[3 * i + 1]);
        max_z = max(max_z, spheres[3 * i + 2]);
    }

    int axis_idx = 0;
    float axis_length = max_x - min_x;
    if (max_y - min_y >= axis_length) {
        axis_length = max_y - min_y;
        axis_idx = 1;
    }
    if (max_z - min_z >= axis_length) {
        axis_length = max_z - min_z;
        axis_idx = 2;
    }

    sort(permutation.begin() + l, permutation.begin() + r, [&](int a, int b) -> bool
        {
            return spheres[3 * a + axis_idx] < spheres[3 * b + axis_idx];
        });
    int m = (l + r) / 2;
    //bvh(spheres, permutation, l, m, (sortAxis + 1) % 3);
    //bvh(spheres, permutation, m, r, (sortAxis + 1) % 3);

    // int axis_idx = rand() % 3;

    bvh(spheres, permutation, l, m);
    bvh(spheres, permutation, m, r);
}

int main()
{
    if (!glfwInit())
        return 255;

          GLFWmonitor* monitor = glfwGetPrimaryMonitor();
    const GLFWvidmode* mode    = glfwGetVideoMode(monitor);


    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    glfwWindowHint(GLFW_RED_BITS    , mode->redBits    );
    glfwWindowHint(GLFW_GREEN_BITS  , mode->greenBits  );
    glfwWindowHint(GLFW_BLUE_BITS   , mode->blueBits   );
    glfwWindowHint(GLFW_REFRESH_RATE, mode->refreshRate);

    // wind_width  = mode->width;
    // wind_height = mode->height;

    GLFWwindow* window;

    glfwSetErrorCallback(glfw_error_callback);

    window = glfwCreateWindow(wind_width,wind_height,"Raymarching",NULL,NULL);
    if (!window) {
        glfwTerminate();
        return 254;
    }

    glfwMakeContextCurrent(window);

    if(!gladLoadGL()) {
        printf("gladLoadGL failed!\n");
        return 253;
    }
    printf("OpenGL %d.%d\n", GLVersion.major, GLVersion.minor);

    // Load point clouds
    pointClouds = read_files("/../assets/point-cloud-300M");
    for (auto &p : pointClouds)
    {
        pointCloudNames.push_back(p.first);
    }

    cl_int errCode;
    try {
        Platform lPlatform = getPlatform();
        // Select the default platform and create a context using this platform and the GPU
#ifdef OS_LNX
        cl_context_properties cps[] = {
            CL_GL_CONTEXT_KHR, (cl_context_properties)glfwGetGLXContext(window),
            CL_GLX_DISPLAY_KHR, (cl_context_properties)glfwGetX11Display(),
            CL_CONTEXT_PLATFORM, (cl_context_properties)lPlatform(),
            0
        };
#endif
#ifdef OS_WIN
        cl_context_properties cps[] = {
            CL_GL_CONTEXT_KHR, (cl_context_properties)glfwGetWGLContext(window),
            CL_WGL_HDC_KHR, (cl_context_properties)GetDC(glfwGetWin32Window(window)),
            CL_CONTEXT_PLATFORM, (cl_context_properties)lPlatform(),
            0
        };
#endif

#ifdef __APPLE__
        CGLContextObj     kCGLContext     = CGLGetCurrentContext();
        CGLShareGroupObj  kCGLShareGroup  = CGLGetShareGroup(kCGLContext);

        cl_context_properties cps[] = {
            CL_CONTEXT_PROPERTY_USE_CGL_SHAREGROUP_APPLE,
            (cl_context_properties) kCGLShareGroup,
            0
        };
#endif
        std::vector<Device> devices;
        lPlatform.getDevices(CL_DEVICE_TYPE_GPU, &devices);
        // Get a list of devices on this platform
        for (unsigned d=0; d<devices.size(); ++d) {
            if (checkExtnAvailability(devices[d],CL_GL_SHARING_EXT)) {
                params.d = devices[d];
                break;
            }
        }
        Context context(params.d, cps);
        // Create a command queue and use the first device
        params.q = CommandQueue(context, params.d);
        
        // Set up program and kernels
        ifstream kernelFile1(ASSETS_DIR"/raytrace.cl");
        string kernelSource1((istreambuf_iterator<char>(kernelFile1)), istreambuf_iterator<char>());
        //ifstream kernelFile2(ASSETS_DIR"/raymarch.cl");
        //string kernelSource2((istreambuf_iterator<char>(kernelFile2)), istreambuf_iterator<char>());
        ifstream kernelFile3(ASSETS_DIR"/calculate.cl");
        string kernelSource3((istreambuf_iterator<char>(kernelFile3)), istreambuf_iterator<char>());
        ifstream kernelFile4(ASSETS_DIR"/integrate.cl");    
        string kernelSource4((istreambuf_iterator<char>(kernelFile4)), istreambuf_iterator<char>());

        Program::Sources sources;
        sources.push_back({ kernelSource1.c_str(), kernelSource1.length() });
        //sources.push_back({ kernelSource2.c_str(), kernelSource2.length() });
        sources.push_back({ kernelSource3.c_str(), kernelSource3.length() });
        sources.push_back({ kernelSource4.c_str(), kernelSource4.length() });
        params.p = Program(context, sources);

        params.p.build(std::vector<Device>(1, params.d));
        params.raytrace = Kernel(params.p, "raytrace");
        //params.raymarch = Kernel(params.p, "raymarch");
        params.kc = Kernel(params.p, "calculate");
        params.ki = Kernel(params.p, "integrate");

        // create opengl stuff
        rparams.prg = initShaders(ASSETS_DIR"/fractal.vert", ASSETS_DIR "/fractal.frag");
        rparams.tex = createTexture2D(wind_width,wind_height);
        GLuint vbo  = createBuffer(12,vertices,GL_STATIC_DRAW);
        GLuint tbo  = createBuffer(8,texcords,GL_STATIC_DRAW);
        GLuint ibo;
        glGenBuffers(1,&ibo);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER,sizeof(uint)*6,indices,GL_STATIC_DRAW);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
        // bind vao
        glGenVertexArrays(1,&rparams.vao);
        glBindVertexArray(rparams.vao);
        // attach vbo
        glBindBuffer(GL_ARRAY_BUFFER,vbo);
        glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,NULL);
        glEnableVertexAttribArray(0);
        // attach tbo
        glBindBuffer(GL_ARRAY_BUFFER,tbo);
        glVertexAttribPointer(1,2,GL_FLOAT,GL_FALSE,0,NULL);
        glEnableVertexAttribArray(1);
        // attach ibo
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,ibo);
        glBindVertexArray(0);
        // create opengl texture reference using opengl texture
        params.tex = ImageGL(context,CL_MEM_READ_WRITE, GL_TEXTURE_2D, 0, rparams.tex,&errCode);
        if (errCode!=CL_SUCCESS) {
            std::cout<<"Failed to create OpenGL texture reference: "<<errCode<<std::endl;
            return 250;
        }

        // Initial camera params
        cam.targetPos[0] = cam.targetPos[1] = cam.targetPos[2] = 0;
        cam.theta = 0;
        cam.phi = M_PI / 2;
        cam.r = 3;
        cam.hFov = 50;
        cam.vFov = 35;

        // Configure camera
        float ar1 = tan(radians(cam.hFov) / 2) / tan(radians(cam.vFov) / 2);
        float ar = static_cast<double>(wind_width) / wind_height;
        if (ar1 < ar)
        {
            // hFov is too small
            cam.hFov = 2 * degrees(atan(tan(radians(cam.vFov) / 2) * ar));
        }
        else if (ar1 > ar)
        {
            // vFov is too small
            cam.vFov = 2 * degrees(atan(tan(radians(cam.hFov) / 2) / ar));
        }

        params.c2w = Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * 9);
        computeCameraPosition();
        
        std::vector<int> seed(wind_width * wind_height);
        for (int i = 0; i < wind_width * wind_height; i++)
            seed[i] = rand();
        params.seed = Buffer(context, CL_MEM_READ_WRITE, sizeof(int) * wind_width * wind_height);
        params.q.enqueueWriteBuffer(params.seed, CL_TRUE, 0, sizeof(int) * wind_width * wind_height, seed.data());

        for (int i = 0; i < nparticles; i++)
            permutation[i] = i;
        params.permutation = Buffer(context, CL_MEM_READ_WRITE, sizeof(int) * nparticles);  
        params.bboxes = Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * 6 * nparticles / LOCAL_WORK_SIZE);

        std::random_device rd;
        std::mt19937 eng(rd());
        std::normal_distribution<> dist(0, 1);
        std::vector<float> vel(3 * nparticles), accel(3 * nparticles), color(3 * nparticles);
        for (int i = 0; i < nparticles; i++)
        {
            accel[3 * i] = accel[3 * i + 1] = accel[3 * i + 2] = 0;
            spheres[3 * i] = dist(eng);
            spheres[3 * i + 1] = dist(eng);
            spheres[3 * i + 2] = dist(eng);
            float norm = sqrt(spheres[3 * i] * spheres[3 * i] + spheres[3 * i + 1] * spheres[3 * i + 1] + spheres[3 * i + 2] * spheres[3 * i + 2]);
            spheres[3 * i] /= norm * 2;
            spheres[3 * i + 1] /= norm * 2;
            spheres[3 * i + 2] /= norm * 2;
            vel[3 * i] = dist(eng);
            vel[3 * i + 1] = dist(eng);
            vel[3 * i + 2] = dist(eng);
            norm = sqrt(vel[3 * i] * vel[3 * i] + vel[3 * i + 1] * vel[3 * i + 1] + vel[3 * i + 2] * vel[3 * i + 2]);
            vel[3 * i] /= norm * 2;
            vel[3 * i + 1] /= norm * 2;
            vel[3 * i + 2] /= norm * 2;
            color[3 * i] = color[3 * i + 1] = color[3 * i + 2] = 1;
        }
        params.spheres = Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * 3 * nparticles);
        params.v = Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * 3 * nparticles);
        params.a = Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * 3 * nparticles);
        params.targetPos = Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * 3 * nparticles);
        params.targetColor = Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * 3 * nparticles);
        params.colors = Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * 3 * nparticles);
        params.q.enqueueWriteBuffer(params.spheres, CL_TRUE, 0, sizeof(float) * 3 * nparticles, spheres.data());
        params.q.enqueueWriteBuffer(params.v, CL_TRUE, 0, sizeof(float) * 3 * nparticles, vel.data());
        params.q.enqueueWriteBuffer(params.a, CL_TRUE, 0, sizeof(float) * 3 * nparticles, accel.data());
        params.q.enqueueWriteBuffer(params.targetPos, CL_TRUE, 0, sizeof(float) * 3 * nparticles, spheres.data());
        params.q.enqueueWriteBuffer(params.targetColor, CL_TRUE, 0, sizeof(float) * 3 * nparticles, color.data());
        params.q.enqueueWriteBuffer(params.colors, CL_TRUE, 0, sizeof(float) * 3 * nparticles, color.data());
    } catch(Error error) {
        std::cout << error.what() << "(" << error.err() << ")" << std::endl;
        std::string val = params.p.getBuildInfo<CL_PROGRAM_BUILD_LOG>(params.d);
        std::cout<<"Log:\n"<<val<<std::endl;
        return 249;
    }

    glfwSetKeyCallback(window, glfw_key_callback);
    glfwSetCursorPosCallback(window, glfw_cursor_position_callback);
    glfwSetMouseButtonCallback(window, glfw_mouse_button_callback);
    glfwSetScrollCallback(window, glfw_scroll_callback);
    glfwSetFramebufferSizeCallback(window,glfw_framebuffer_size_callback);

    mouse.pressed = false;

    double prev_time = 0.0;
    double curr_time = 0.0;
    double time_diff;
    unsigned int counter = 0;

    float previousTime = glfwGetTime();

    while (!glfwWindowShouldClose(window)) {

        curr_time = glfwGetTime();
        time_diff = curr_time - prev_time;
        counter ++;
        if (time_diff >= 1.0 / 30.0) 
        {
            std::string FPS = std::to_string((1.0 / time_diff) * counter).substr(0, 4);
            std::string ms = std::to_string((time_diff / counter) * 1000).substr(0, 4);
            std::string newTitle = "Ray Marching - " + FPS + "FPS / " + ms + "ms - " + (selectedIndex == -1 ? "flocking" : pointCloudNames[selectedIndex]);
            glfwSetWindowTitle(window, newTitle.c_str());
            prev_time = curr_time;
            counter = 0;
        }
        float currentTime = glfwGetTime();
        // process call

        processTimeStep(paused ? 0 : currentTime - previousTime);
        // render call

        auto start_render_time = std::chrono::high_resolution_clock::now();
        renderFrame();
#ifdef PROFILE
        auto end_render_time = std::chrono::high_resolution_clock::now();
        auto render_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_render_time - start_render_time).count();
        std::cout << "Render time: " << render_duration << "μs" << std::endl;
#endif

        // swap front and back buffers
        auto start_swap_buffer_time = std::chrono::high_resolution_clock::now();
        glfwSwapBuffers(window);
#ifdef PROFILE
        auto end_swap_buffer_time = std::chrono::high_resolution_clock::now();
        auto swap_buffer_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_swap_buffer_time - start_swap_buffer_time).count();
        std::cout << "Swap buffer time: " << swap_buffer_duration << "μs" << std::endl;
#endif
        
        // poll for events
        auto start_poll_time = std::chrono::high_resolution_clock::now();
        glfwPollEvents();
#ifdef PROFILE
        auto end_poll_time = std::chrono::high_resolution_clock::now();
        auto poll_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_poll_time - start_poll_time).count();
        std::cout << "Poll time: " << poll_duration << "μs" << std::endl;
#endif

        previousTime = currentTime;
    }

    glfwDestroyWindow(window);

    glfwTerminate();
    return 0;
}

inline unsigned divup(unsigned a, unsigned b)
{
    return (a+b-1)/b;
}

void processTimeStep(float deltaTime)
{
    cl::Event ev;
    try {
        //params.q.finish();
        NDRange local(16);
        NDRange global(16 * divup(nparticles, 16));
        // // set kernel arguments
        params.kc.setArg(0, params.spheres);
        params.kc.setArg(1, params.v);
        params.kc.setArg(2, params.a);
        params.kc.setArg(3, params.targetPos);
        params.kc.setArg(4, selectedIndex);
        local = NDRange(16);
        global = NDRange(16 * divup(nparticles, 16));
        params.q.enqueueNDRangeKernel(params.kc, cl::NullRange, global, local);
        params.ki.setArg(0, params.spheres);
        params.ki.setArg(1, params.v);
        params.ki.setArg(2, params.a);
        params.ki.setArg(3, params.colors);
        params.ki.setArg(4, params.targetColor);
        params.ki.setArg(5, deltaTime);
        
        auto particle_start = std::chrono::high_resolution_clock::now();
        params.q.enqueueNDRangeKernel(params.ki, cl::NullRange, global, local);
        params.q.finish();

        glFinish();
#ifdef PROFILE
        auto particle_end = std::chrono::high_resolution_clock::now();
        auto particle_duration = std::chrono::duration_cast<std::chrono::microseconds>(particle_end - particle_start).count();
        std::cout << "Particle sim time: " << particle_duration << "μs" << std::endl;
#endif

        std::vector<Memory> objs;
        objs.clear();
        objs.push_back(params.tex);
        // flush opengl commands and wait for object acquisition
        cl_int res = params.q.enqueueAcquireGLObjects(&objs,NULL,&ev);
        ev.wait();
        if (res!=CL_SUCCESS) {
            std::cout<<"Failed acquiring GL object: "<<res<<std::endl;
            exit(248);
        }
        float hFov_expr = 2 * tan(0.5 * cam.hFov * M_PI / 180);
        float vFov_expr = 2 * tan(0.5 * cam.vFov * M_PI / 180);

        auto bvh_start = std::chrono::high_resolution_clock::now();
        params.q.enqueueReadBuffer(params.spheres, CL_TRUE, 0, sizeof(float) * 3 * nparticles, spheres.data());
        params.q.finish();

        bvh(spheres, permutation, 0, nparticles);
#ifdef PROFILE
        auto bvh_end = std::chrono::high_resolution_clock::now();
        auto bvh_duration = std::chrono::duration_cast<std::chrono::microseconds>(bvh_end - bvh_start).count();
        std::cout << "BVH construction time: " << bvh_duration << "μs" << std::endl;
#endif

        // float minx = spheres[0], miny = spheres[1], minz = spheres[2];
        // float maxx = spheres[0], maxy = spheres[1], maxz = spheres[2];
        // for (int i = 1; i < nparticles; i++)
        // {
        //     minx = min(minx, spheres[3 * i]);
        //     maxx = max(maxx, spheres[3 * i]);
        //     miny = min(miny, spheres[3 * i + 1]);
        //     maxy = max(maxy, spheres[3 * i + 1]);
        //     minz = min(minz, spheres[3 * i + 2]);
        //     maxz = max(maxz, spheres[3 * i + 2]);
        // }
        // maxx -= minx;
        // maxy -= miny;
        // maxz -= minz;

        // for (int i = 0; i < nparticles; i++)
        // {
        //     codes[i] = morton3D((spheres[3 * i] - minx) / maxx, (spheres[3 * i + 1] - miny) / maxy, (spheres[3 * i + 2] - minz) / maxz);
        // }

        // sort(permutation.begin(), permutation.end(), [](int a, int b) -> bool
        // {
        //     return codes[a] < codes[b];
        // });

        params.q.enqueueWriteBuffer(params.permutation, CL_TRUE, 0, sizeof(int) * nparticles, permutation.data());

        for (int i = 0; i < nparticles / LOCAL_WORK_SIZE; i++)
        {
            bboxes[6 * i + 0] = bboxes[6 * i + 1] = bboxes[6 * i + 2] = INFINITY;
            bboxes[6 * i + 3] = bboxes[6 * i + 4] = bboxes[6 * i + 5] = -INFINITY;
            for (int j = 0; j < LOCAL_WORK_SIZE; j++)
            {
                int index = i * LOCAL_WORK_SIZE + j;
                for (int k = 0; k < 3; k++)
                {
                    bboxes[6 * i + k] = min(bboxes[6 * i + k], spheres[3 * permutation[index] + k] - SPHERE_RADIUS);
                    bboxes[6 * i + 3 + k] = max(bboxes[6 * i + 3 + k], spheres[3 * permutation[index] + k] + SPHERE_RADIUS);
                }
            }
        }
        params.q.enqueueWriteBuffer(params.bboxes, CL_TRUE, 0, sizeof(float) * 6 * nparticles / LOCAL_WORK_SIZE, bboxes.data());

        local = NDRange(LOCAL_WORK_SIZE_X, LOCAL_WORK_SIZE_Y);
        global = NDRange( local[0] * divup(wind_width, local[0]),
                        local[1] * divup(wind_height, local[1]));
        // set kernel arguments
        int temp = -1;
        params.raytrace.setArg(0, params.tex);
        params.raytrace.setArg(1, wind_width);
        params.raytrace.setArg(2, wind_height);
        params.raytrace.setArg(3, hFov_expr);
        params.raytrace.setArg(4, vFov_expr);
        params.raytrace.setArg(5, params.cameraPos);
        params.raytrace.setArg(6, temp);
        params.raytrace.setArg(7, params.c2w);
        params.raytrace.setArg(8, params.spheres);
        params.raytrace.setArg(9, params.permutation);
        params.raytrace.setArg(10, params.bboxes);
        params.raytrace.setArg(11, params.colors);
        params.raytrace.setArg(12, params.seed);
        params.q.enqueueNDRangeKernel(params.raytrace, cl::NullRange, global, local);

        // release opengl object
        res = params.q.enqueueReleaseGLObjects(&objs);
        ev.wait();
        if (res!=CL_SUCCESS) {
            std::cout<<"Failed releasing GL object: "<<res<<std::endl;
            exit(247);
        }
        params.q.finish();
    } catch(Error err) {
        std::cout << err.what() << "(" << err.err() << ")" << std::endl;
    }
}

void renderFrame()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0.2,0.2,0.2,0.0);
    glEnable(GL_DEPTH_TEST);
    // bind shader
    glUseProgram(rparams.prg);
    // get uniform locations
    int mat_loc = glGetUniformLocation(rparams.prg,"matrix");
    int tex_loc = glGetUniformLocation(rparams.prg,"tex");
    // bind texture
    glActiveTexture(GL_TEXTURE0);
    glUniform1i(tex_loc,0);
    glBindTexture(GL_TEXTURE_2D,rparams.tex);
    glGenerateMipmap(GL_TEXTURE_2D);
    // set project matrix
    glUniformMatrix4fv(mat_loc,1,GL_FALSE,matrix);
    // now render stuff
    glBindVertexArray(rparams.vao);
    glDrawElements(GL_TRIANGLES,6,GL_UNSIGNED_INT,0);
    glBindVertexArray(0);
}