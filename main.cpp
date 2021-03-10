#include <iostream>
#include <GLFW/glfw3.h>
#include "include/ray_tracing.h"
#include "include/bmp.h"

using namespace glm;
using namespace std;

#define SCR_WIDTH 1920 
#define SCR_HEIGHT 1080 
#define ZOOM 1

int main()
{
    Render renderer(SCR_WIDTH, SCR_HEIGHT);

    // Model box("scenes/cornellbox/cornellbox.obj");
    // Obj obj0(box);
    // renderer.set_camera(glm::dvec3(0, 0, 2.5), glm::dvec3(0, 0, 0), glm::dvec3(0, 1, 0), 60);
    // renderer.adapted_lum = 1.4;

    Model box("scenes/diningroom/cornellbox.obj");
    Obj obj0(box);
    Model room("scenes/diningroom/diningroom.obj");
    Obj obj1(room);
    renderer.set_camera(glm::dvec3(0, 12.72, 31.85), glm::dvec3(0, 12.546, 30.865), glm::dvec3(0, 0.985, -0.174), 45);
    renderer.adapted_lum = 1.5;

    // Model car("scenes/car/car.obj");
    // Obj obj0(car);
    // // renderer.set_camera(glm::dvec3(8.22, -0.61, -9.8), glm::dvec3(7.514, -0.702, -9.097), glm::dvec3(-0.065, 0.996, 0.065), 45);
    // renderer.set_camera(glm::dvec3(5.72, 0.12, 9.55), glm::dvec3(5.085, -0.131, 8.819), glm::dvec3(-0.165, 0.968, -0.189), 45);
    // renderer.adapted_lum = 2;
    // renderer.set_skybox("scenes/car/environment_day.hdr");
    // double k = 148 * std::_Pi / 180;  // 72, 148
    // renderer.skybox_rotate = {{glm::cos(k), 0, glm::sin(k)},
    //                           {0, 1, 0},
    //                           {-glm::sin(k), 0, glm::cos(k)}};

    renderer.add_obj(obj0);
    renderer.add_obj(obj1);
    obj0.set_pose(glm::mat4x4(1.0), 1.0);
    obj0.set_pose({{1, 0, 0, -2},
                   {0, 1, 0, 14},
                   {0, 0, 1, 20},
                   {0, 0, 0, 1}},
                  3.0);

    ////////////////////
    // GLFW
    ////////////////////
    glfwInit();
    GLFWwindow *window = glfwCreateWindow(SCR_WIDTH * ZOOM, SCR_HEIGHT * ZOOM, "x", NULL, NULL);
    glfwMakeContextCurrent(window);
    double glfw_time = glfwGetTime();
    int count = 0;
    int n_samples = 0;
    while (!glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
    {
        n_samples = renderer.render();
        count++;
        glPixelZoom(ZOOM, ZOOM);
        glDrawPixels(SCR_WIDTH, SCR_HEIGHT, GL_RGBA, GL_UNSIGNED_BYTE, renderer.get_framebuffer());
        glfwSwapBuffers(window);
        glfwPollEvents();
        if (count % 300 == 0)
        {
            Bitmap results(SCR_WIDTH, SCR_HEIGHT, renderer.get_framebuffer());
            results.SaveFile("results/" + to_string(n_samples) + ".bmp");
        }
    }
    glfwTerminate();
    return 0;
}
