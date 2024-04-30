#include <glad/glad.h>

#include <common/OpenCLUtil.h>
#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#else
#include <GL/gl.h>
#endif

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
#include <string>
#include <random>
#include <cmath>
#include <algorithm>

using namespace std;
using namespace cl;

typedef unsigned int uint;

static int const wind_width = 640;
static int const wind_height= 480;
static int const nparticles = 4096;

static const float matrix[16] = {
    1.0f, 0.0f, 0.0f, 0.0f,
    0.0f, 1.0f, 0.0f, 0.0f,
    0.0f, 0.0f, 1.0f, 0.0f,
    0.0f, 0.0f, 0.0f, 1.0f};

typedef struct {
    Device d;
    CommandQueue q;
    Program p;
    Kernel kc; // calculate.cl kernel (calculate acceleration given positions and velocities)
    Kernel ki; // integrate.cl kernel (calculate position and velocity given acceleration)
    Buffer i; // position buffer
    Buffer v; // velocity buffer
    Buffer a; // acceleration buffer
    size_t dims[3];
} process_params;

typedef struct {
    GLuint prg;
    GLuint vao;
    GLuint vbo;
    BufferGL tmp;
} render_params;

process_params params;
render_params rparams;

static void glfw_error_callback(int error, const char* desc)
{
    fputs(desc,stderr);
}

static void glfw_key_callback(GLFWwindow* wind, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
        glfwSetWindowShouldClose(wind, GL_TRUE);
    }
}

static void glfw_framebuffer_size_callback(GLFWwindow* wind, int width, int height)
{
    glViewport(0,0,width,height);
}

void processTimeStep(float);
void renderFrame(void);

int main(void)
{
    GLFWwindow* window;

    if (!glfwInit())
        return 255;


    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    glfwSetErrorCallback(glfw_error_callback);

    window = glfwCreateWindow(wind_width,wind_height,"Random particles",NULL,NULL);
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

        ifstream kernelFile1(ASSETS_DIR"/calculate.cl");
        string kernelSource1((istreambuf_iterator<char>(kernelFile1)), istreambuf_iterator<char>());
        ifstream kernelFile2(ASSETS_DIR"/integrate.cl");
        string kernelSource2((istreambuf_iterator<char>(kernelFile2)), istreambuf_iterator<char>());

        Program::Sources sources;
        sources.push_back({ kernelSource1.c_str(), kernelSource1.length() });
        sources.push_back({ kernelSource2.c_str(), kernelSource2.length() });
        params.p = Program(context, sources);

        //params.p = getProgram(context, ASSETS_DIR"/calculate.cl", errCode);
        params.p.build(std::vector<Device>(1, params.d));
        params.kc = Kernel(params.p, "calculate");
        params.ki = Kernel(params.p, "integrate");
        // create opengl stuff
        rparams.prg = initShaders(ASSETS_DIR"/partsim.vert", ASSETS_DIR"/partsim.frag");

        std::random_device rd;
        std::mt19937 eng(rd());
        std::normal_distribution<> dist(10, 100);
        std::vector<float> data(2*nparticles);
        std::vector<float> vel(3*nparticles);
        std::vector<float> accel(3*nparticles);
        for(int n=0; n<nparticles; ++n) {
            data[2*n+0] = std::fmod(dist(eng), wind_width)/wind_width;
            data[2*n+1] = std::fmod(dist(eng), wind_height)/wind_height;
        }
        for(int n=0; n<nparticles; ++n) {
            vel[3*n+0] = 0; // TODO: add random direction vector later
            vel[3*n+1] = 0;
            vel[3*n+2] = 0;
        }
        for(int n=0; n<nparticles; ++n) {
            accel[3*n+0] = 0; // TODO: add random direction vector later
            accel[3*n+1] = 0;
            accel[3*n+2] = 0;
        }

        rparams.vbo = createBuffer(2*nparticles, data.data(), GL_DYNAMIC_DRAW);
        rparams.tmp = BufferGL(context, CL_MEM_READ_WRITE, rparams.vbo, NULL);
        // bind vao
        glGenVertexArrays(1,&rparams.vao);
        glBindVertexArray(rparams.vao);
        // attach vbo
        glBindBuffer(GL_ARRAY_BUFFER,rparams.vbo);
        glVertexAttribPointer(0,2,GL_FLOAT,GL_FALSE,0,NULL);
        glEnableVertexAttribArray(0);
        glBindVertexArray(0);
        // create opencl input and output buffer
        params.i = Buffer(context, CL_MEM_READ_WRITE, sizeof(float)*2*nparticles);
        params.v = Buffer(context, CL_MEM_READ_WRITE, sizeof(float)*3*nparticles);
        params.a = Buffer(context, CL_MEM_READ_WRITE, sizeof(float)*3*nparticles);
        params.q.enqueueWriteBuffer(params.i, CL_TRUE, 0, sizeof(float)*2*nparticles, data.data());
        params.q.enqueueWriteBuffer(params.v, CL_TRUE, 0, sizeof(float)*3*nparticles, vel.data());
        params.q.enqueueWriteBuffer(params.a, CL_TRUE, 0, sizeof(float)*3*nparticles, accel.data());
        params.dims[0] = nparticles;
        params.dims[1] = 1;
        params.dims[2] = 1;
    } catch(Error error) {
        std::cout << error.what() << "(" << error.err() << ")" << std::endl;
        std::string val = params.p.getBuildInfo<CL_PROGRAM_BUILD_LOG>(params.d);
        std::cout<<"Log:\n"<<val<<std::endl;
        return 249;
    }

    glfwSetKeyCallback(window,glfw_key_callback);
    glfwSetFramebufferSizeCallback(window,glfw_framebuffer_size_callback);

    float previousTime = glfwGetTime();

    double prev_time = 0.0;
    double curr_time = 0.0;
    double time_diff;
    unsigned int counter = 0;

    while (!glfwWindowShouldClose(window)) {
        curr_time = glfwGetTime();
        time_diff = curr_time - prev_time;
        counter ++;
        if (time_diff >= 1.0 / 30.0) 
        {
            std::string FPS = std::to_string((1.0 / time_diff) * counter);
            std::string ms = std::to_string((time_diff / counter) * 1000);
            std::string newTitle = "Boids - " + FPS + "FPS / " + ms + "ms";
            glfwSetWindowTitle(window, newTitle.c_str());
            prev_time = curr_time;
            counter = 0;
        }
        float currentTime = glfwGetTime();
        // process call
        processTimeStep(currentTime - previousTime);
        // render call
        renderFrame();
        // swap front and back buffers
        glfwSwapBuffers(window);
        // poll for events
        glfwPollEvents();
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
        NDRange local(16);
        NDRange global(16 * divup(params.dims[0], 16));
        // set kernel arguments
        params.kc.setArg(0, nparticles);
        params.kc.setArg(1, params.i);
        params.kc.setArg(2, params.v);
        params.kc.setArg(3, params.a);
        params.kc.setArg(4, sizeof(float), &deltaTime);
        params.q.enqueueNDRangeKernel(params.kc, cl::NullRange, global, local);
        params.ki.setArg(0, nparticles);
        params.ki.setArg(1, params.i);
        params.ki.setArg(2, params.v);
        params.ki.setArg(3, params.a);
        params.ki.setArg(4, sizeof(float), &deltaTime);
        params.q.enqueueNDRangeKernel(params.ki, cl::NullRange, global, local);
        glFinish();
        std::vector<Memory> objs;
        objs.clear();
        objs.push_back(rparams.tmp);
        // flush opengl commands and wait for object acquisition
        cl_int res = params.q.enqueueAcquireGLObjects(&objs,NULL,&ev);
        ev.wait();
        if (res!=CL_SUCCESS) {
            std::cout<<"Failed acquiring GL object: "<<res<<std::endl;
            exit(248);
        }
        params.q.enqueueCopyBuffer(params.i, rparams.tmp, 0, 0, 2*nparticles*sizeof(float), NULL, NULL);
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
    static const float pcolor[4] = {1.0, 0.0, 0.0, 1.0};
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0.2,0.2,0.2,1.0);
    glEnable(GL_DEPTH_TEST);
    // bind shader
    glUseProgram(rparams.prg);
    // get uniform locations
    int mat_loc = glGetUniformLocation(rparams.prg,"matrix");
    int col_loc = glGetUniformLocation(rparams.prg,"color");
    // set project matrix
    glUniformMatrix4fv(mat_loc,1,GL_FALSE,matrix);
    glUniform4fv(col_loc,1,pcolor);
    // now render stuff
    glPointSize(4);
    glBindVertexArray(rparams.vao);
    glDrawArrays(GL_POINTS,0,nparticles);
    glBindVertexArray(0);
}
