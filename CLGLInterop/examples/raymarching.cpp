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

using namespace std;
using namespace cl;

#define EPS_F (0.00001f)
#define MIN_R (0.1f)
#define MAX_R (5.0f)

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

static int wind_width = 640;
static int wind_height= 480;
static int const nparticles = 4096;
static int gJuliaSetIndex = 0;

typedef struct {
    Device d;
    CommandQueue q;
    Program p;
    Kernel k;
    ImageGL tex;
    cl_float3 cameraPos;
    Buffer spheres;
    Buffer c2w;
    Buffer seed;
    size_t dims[3];
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
    float screenDist;
} camera_state;

typedef struct {
    float x;
    float y;
    bool pressed; //left mouse button only
} mouse_state;

process_params params;
render_params rparams;
camera_state cam;
mouse_state mouse;

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

    // cout << "upvec" << endl;
    // cout << upVec[0] << " " << upVec[1] << " " << upVec[2] << endl;

    // cout << "c2w" << endl;
    // for (int i = 0; i < 3; i++)
    // {
    //     for (int j = 0; j < 3; j++)
    //     {
    //         cout << c2w[3 * i + j] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << "Camera pos" << endl;
    // for (int i = 0; i < 3; i++)
    //     cout << params.cameraPos.s[i] << " ";
    // cout << endl;
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
        if (key == GLFW_KEY_ESCAPE)
            glfwSetWindowShouldClose(wind, GL_TRUE);
        else if (key == GLFW_KEY_1)
            gJuliaSetIndex = 0;
        else if (key == GLFW_KEY_2)
            gJuliaSetIndex = 1;
        else if (key == GLFW_KEY_3)
            gJuliaSetIndex = 2;
        else if (key == GLFW_KEY_4)
            gJuliaSetIndex = 3;
        else if (key == GLFW_KEY_5)
            gJuliaSetIndex = 4;
        else if (key == GLFW_KEY_6)
            gJuliaSetIndex = 5;
        else if (key == GLFW_KEY_7)
            gJuliaSetIndex = 6;
        else if (key == GLFW_KEY_8)
            gJuliaSetIndex = 7;
        else if (key == GLFW_KEY_9)
            gJuliaSetIndex = 8;
    }
}

static void glfw_framebuffer_size_callback(GLFWwindow* wind, int width, int height)
{
    glViewport(0,0,width,height);
    setScreenSize(width, height);
    //wind_width = width;
    //wind_height = height;
}

void processTimeStep(void);
void renderFrame(void);

inline float radians(float angle)
{
    return angle * M_PI / 180;
}

inline float degrees(float radians)
{
    return radians * 180 / M_PI;
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
    //window = glfwCreateWindow(wind_width,wind_height,"Julia Sets",monitor,NULL);
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
        
        params.p = getProgram(context, ASSETS_DIR"/raytrace.cl",errCode);

        std::ostringstream options;
        options << "-I " << std::string(ASSETS_DIR);

        params.p.build(std::vector<Device>(1, params.d), options.str().c_str());
        params.k = Kernel(params.p, "raytrace", &errCode);

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
        cam.screenDist = 1;

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
        cam.screenDist = ((double) wind_height) / (2.0 * tan(radians(cam.vFov) / 2));

        params.c2w = Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * 9);
        computeCameraPosition();
        
        std::vector<int> seed(wind_width * wind_height);
        for (int i = 0; i < wind_width * wind_height; i++)
            seed[i] = rand();
        params.seed = Buffer(context, CL_MEM_READ_WRITE, sizeof(int) * wind_width * wind_height);
        params.q.enqueueWriteBuffer(params.seed, CL_TRUE, 0, sizeof(int) * wind_width * wind_height, seed.data());

        std::vector<float> spheres(3 * nparticles);
        for (int i = 0; i < 3 * nparticles; i++)
        {
            spheres[i] = 0;
        }
        spheres[3] = 0.5;
        params.spheres = Buffer(context, CL_MEM_READ_WRITE, sizeof(float) * 3 * nparticles);
        params.q.enqueueWriteBuffer(params.spheres, CL_TRUE, 0, sizeof(float) * 3 * nparticles, spheres.data());

        params.dims[0] = wind_width;
        params.dims[1] = wind_height;
        params.dims[2] = 1;
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

    while (!glfwWindowShouldClose(window)) {
        // process call
        processTimeStep();
        // render call
        renderFrame();
        // swap front and back buffers
        glfwSwapBuffers(window);
        // poll for events
        glfwPollEvents();
    }

    glfwDestroyWindow(window);

    glfwTerminate();
    return 0;
}

inline unsigned divup(unsigned a, unsigned b)
{
    return (a+b-1)/b;
}

void processTimeStep()
{
    cl::Event ev;
    try {
        glFinish();

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
        NDRange local(16, 16);
        NDRange global( local[0] * divup(params.dims[0], local[0]),
                        local[1] * divup(params.dims[1], local[1]));
        // set kernel arguments
        params.k.setArg(0, params.tex);
        params.k.setArg(1, (uint)params.dims[0]);
        params.k.setArg(2, (uint)params.dims[1]);
        params.k.setArg(3, cam.hFov);
        params.k.setArg(4, cam.vFov);
        params.k.setArg(5, params.cameraPos);
        params.k.setArg(6, params.c2w);
        params.k.setArg(7, params.spheres);
        params.k.setArg(8, params.seed);
        params.q.enqueueNDRangeKernel(params.k, cl::NullRange, global, local);
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