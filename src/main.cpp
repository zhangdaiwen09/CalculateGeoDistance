

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "GeodesicDistanceCalculator.h"
#include "OBJLoader.h"
#include "Renderer.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <dirent.h>
#include <vector>
#include <string>
#include <iostream>

#include <gl/glut.h> 

// Renderer object to render 3D meshes
Renderer renderer;

// Pointer to a GeodesicDistanceCalculator, initially set to nullptr
GeodesicDistanceCalculator* calculator = nullptr;

// Vec3 objects to represent the start and end points for geodesic distance calculation
Vec3 startPoint, endPoint;

// Flags to indicate whether the start and end points have been selected
bool startSelected = false, endSelected = false;

// Variable to store the calculated geodesic distance
double geodesicDistance = 0.0;

// Matrices for projection, view, and model transformations
glm::mat4 projection, view, model;

// Variables to handle mouse input
bool mousePressed = false;
double lastX, lastY;

// Variables to control camera orientation
float yaw = -90.0f, pitch = 0.0f;
float fov = 45.0f;

// Camera position and orientation vectors
glm::vec3 cameraPos = glm::vec3(0.0f, 0.0f, 3.0f);
glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
glm::vec3 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);

// Function to retrieve OBJ files from a specified folder
std::vector<std::string> getOBJFiles(const std::string& folder) {
    std::vector<std::string> files;
    DIR* dir;
    struct dirent* ent;
    if ((dir = opendir(folder.c_str())) != NULL) {
        while ((ent = readdir(dir)) != NULL) {
            std::string fileName = ent->d_name;
            if (fileName.size() > 4 && fileName.substr(fileName.size() - 4) == ".obj") {
                files.push_back(fileName);  // Add .obj files to the list
            }
        }
        closedir(dir);
    } else {
        perror("Could not open directory");
    }
    return files;
}

// Function to perform a matrix transformation on a point
static void transform_point(GLdouble out[4], const GLdouble m[16], const GLdouble in[4]) {
#define M(row,col) m[col*4+row]
    out[0] =
        M(0, 0) * in[0] + M(0, 1) * in[1] + M(0, 2) * in[2] + M(0, 3) * in[3];
    out[1] =
        M(1, 0) * in[0] + M(1, 1) * in[1] + M(1, 2) * in[2] + M(1, 3) * in[3];
    out[2] =
        M(2, 0) * in[0] + M(2, 1) * in[1] + M(2, 2) * in[2] + M(2, 3) * in[3];
    out[3] =
        M(3, 0) * in[0] + M(3, 1) * in[1] + M(3, 2) * in[2] + M(3, 3) * in[3];
#undef M
}

// Function to project a 3D point onto the screen (gluProject source code)
static GLint gluProject(GLdouble objx, GLdouble objy, GLdouble objz, const GLdouble modelMatrix[16], const GLdouble projMatrix[16], const GLint viewport[4], GLdouble* winx, GLdouble* winy, GLdouble* winz) {
    GLdouble in[4], out[4];
    in[0] = objx;
    in[1] = objy;
    in[2] = objz;
    in[3] = 1.0;
    transform_point(out, modelMatrix, in);  // Apply model view matrix
    transform_point(in, projMatrix, out);   // Apply projection matrix
    if (in[3] == 0.0)
        return GL_FALSE;
    in[0] /= in[3];
    in[1] /= in[3];
    in[2] /= in[3];
    *winx = viewport[0] + (1 + in[0]) * viewport[2] / 2;
    *winy = viewport[1] + (1 + in[1]) * viewport[3] / 2;
    *winz = (1 + in[2]) / 2;
    return GL_TRUE;
}

// Callback function for mouse button events
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);

        float x = xpos, y = ypos, z = 1.0;
        float ox, oy, oz;

        // Screen to world coordinate conversion
        GLint viewport[4];
        GLdouble modelview[16];
        GLdouble projection[16];
        GLfloat winX, winY, winZ;
        GLdouble posX, posY, posZ;

        glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
        glGetDoublev(GL_PROJECTION_MATRIX, projection);
        glGetIntegerv(GL_VIEWPORT, viewport);

        winX = (float)x;
        winY = (float)viewport[3] - (float)y;
        glReadPixels(int(winX), int(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);
        gluUnProject(winX, winY, winZ, modelview, projection, viewport, &posX, &posY, &posZ);

        oy = (float)posY;
        ox = (float)posX;
        oz = (float)posZ;

        glm::vec3 glm_cp2 = glm::vec3(ox, oy, oz);
        Vec3 clickedPoint2 = Vec3(glm_cp2.x + 0.5, glm_cp2.y + 0.5, glm_cp2.z + 0.5);

        if (!startSelected) {
            startPoint = clickedPoint2;
            startSelected = true;
        } else if (!endSelected) {
            endPoint = clickedPoint2;
            endSelected = true;
            if (calculator) {
                std::vector<Vec3> new_verts;
                std::vector<std::array<int, 3>> new_fcs;
                calculator->data(new_verts, new_fcs);
                calculator->exchange(startPoint, new_verts, new_fcs);
                calculator->exchange(endPoint, new_verts, new_fcs);

                GeodesicDistanceCalculator gdc(new_verts, new_fcs);
                gdc.analysis();
                geodesicDistance = gdc.computeGeodesicDistance(startPoint, endPoint);
            }
        } else {
            startSelected = false;
            endSelected = false;
            geodesicDistance = 0.0;
        }
    }
}

// Setup function for ImGui
void setupImGui(GLFWwindow* window) {
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    ImGui::StyleColorsDark();

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 130");
}

// Function to render the ImGui interface
void renderImGui(std::vector<std::string>& objFiles, std::string& selectedFile, bool& fileLoaded) {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    ImGui::Begin("Geodesic Distance Calculator");
    ImGui::Text("OBJ Files in models folder:");
    for (const auto& file : objFiles) {
        if (ImGui::Selectable(file.c_str(), file == selectedFile)) {
            selectedFile = file;
            fileLoaded = false;  // Mark file as not loaded
        }
    }

    ImGui::Text("Start Point: (%.3f, %.3f, %.3f) %s", startPoint[0], startPoint[1], startPoint[2], startSelected ? "ok" : "");
    ImGui::Text("End Point: (%.3f, %.3f, %.3f) %s", endPoint[0], endPoint[1], endPoint[2], endSelected ? "ok" : "");
    
    if (startSelected && endSelected) {
        ImGui::Text("Geodesic Distance: %.3f", geodesicDistance);
    }
    ImGui::End();

    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

// Function to process keyboard and mouse input
void processInput(GLFWwindow* window) {
    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS) {
        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);
        if (!mousePressed) {
            mousePressed = true;
            lastX = xpos;
            lastY = ypos;
        }

        double xoffset = xpos - lastX;
        double yoffset = ypos - lastY;  // Reversed since y-coordinates range from bottom to top
        lastX = xpos;
        lastY = ypos;

        float angleX = (5 * xoffset);
        float angleY = (5 * yoffset);

        // Create rotation matrices and apply them to the model matrix
        model = glm::rotate(model, -glm::radians(angleX), glm::vec3(0.0f, 1.0f, 0.0f));
        model = glm::rotate(model, -glm::radians(angleY), glm::vec3(1.0f, 0.0f, 0.0f));
    } else {
        mousePressed = false;
    }

    // Handle camera movement
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        cameraPos += 0.05f * cameraFront;
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        cameraPos -= 0.05f * cameraFront;
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        cameraPos -= glm::normalize(glm::cross(cameraFront, cameraUp)) * 0.05f;
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        cameraPos += glm::normalize(glm::cross(cameraFront, cameraUp)) * 0.05f;
}

int main() {
    // Initialize GLFW
    if (!glfwInit()) {
        return -1;
    }

    // Create a windowed mode window and its OpenGL context
    GLFWwindow* window = glfwCreateWindow(800, 600, "Geodesic Distance Calculator", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetMouseButtonCallback(window, mouse_button_callback);

    // Initialize GLEW
    if (glewInit() != GLEW_OK) {
        return -1;
    }

    // Set OpenGL state
    glEnable(GL_DEPTH_TEST);

    // Initialize ImGui
    setupImGui(window);

    // Set up projection and view matrices
    projection = glm::perspective(glm::radians(fov), 4.0f / 3.0f, 0.1f, 100.0f);
    view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);

    model = glm::mat4(1.0f);

    // Get OBJ files from the models folder
    std::string modelsFolder = "models";
    std::vector<std::string> objFiles = getOBJFiles(modelsFolder);
    std::string selectedFile;
    bool fileLoaded = false;

    // Main loop
    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Process input
        processInput(window);

        // Update projection and view matrices
        projection = glm::perspective(glm::radians(fov), 4.0f / 3.0f, 0.1f, 100.0f);
        view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);

        // Load and render the selected OBJ file
        if (!fileLoaded && !selectedFile.empty()) {
            OBJLoader objLoader;
            objLoader.load(modelsFolder + "/" + selectedFile);
            renderer.setMesh(objLoader.vertices, objLoader.faces);

            if (calculator) {
                delete calculator;
                calculator = nullptr;
            }

            calculator = new GeodesicDistanceCalculator(objLoader.vertices, objLoader.faces);
            fileLoaded = true;
        }

        // Render the model
        glm::mat4 mvp = projection * view * model;
        glLoadMatrixf(glm::value_ptr(mvp));
        renderer.render();

        // Render ImGui interface
        renderImGui(objFiles, selectedFile, fileLoaded);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // Clean up
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();

    delete calculator;

    return 0;
}
