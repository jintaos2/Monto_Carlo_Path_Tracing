#include <iostream>
#include <GLFW/glfw3.h>
#include "include/ray_tracing.h"
#include "include/bmp.h"

using namespace glm;
using namespace std;

#define SCR_WIDTH 400
#define SCR_HEIGHT 400
#define ZOOM 2

int main()
{

    // string test("f 1/1/1 2/2/2 3/3/3");
    // smatch result;
    // regex r("\\s*f\\s+(\\d+)/(\\d+)/(\\d+)\\s+(\\d+)/(\\d+)/(\\d+)\\s+(\\d+)/(\\d+)/(\\d+)");
    // if (regex_search(test, result, r))
    // {
    //     for (int i = 0; i < result.size(); ++i)
    //     {
    //         cout << result[i] << endl;
    //     }
    // }
    Render renderer(SCR_WIDTH, SCR_HEIGHT);

    Model box("scenes/cornellbox/cornellbox.obj");
    Obj obj0(box);
    renderer.set_camera(vec3(0, 0, 2.5), vec3(0, 0, 0), vec3(0, 1, 0), 60);

    // Model room("scenes/diningroom/diningroom.obj");
    // Obj obj0(room);
    // renderer.set_camera(vec3(0, 12.72,31.85), vec3(0, 12.546, 30.865), vec3(0, 0.985, -0.174), 45);

    // Model car("scenes/car/car.obj");
    // Obj obj0(car);
    // renderer.set_camera(vec3(0, 0, 2.5), vec3(0, 0, 0), vec3(0, 1, 0), 60);

    renderer.add_obj(obj0);
    obj0.set_pose(mat4x4(1.0f), 1.0f);

    renderer.set_background(vec3(0, 0, 0));
    renderer.SAMPLES = 100000000;

    // while(1){}
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
