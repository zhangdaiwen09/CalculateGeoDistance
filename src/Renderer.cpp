#include "Renderer.h"
#include <iostream>

// Set the mesh data for rendering, including vertices and faces
void Renderer::setMesh(const std::vector<Vec3>& verts, const std::vector<std::array<int, 3>>& fcs) {
    vertices = verts;  // Store the vertices
    faces = fcs;  // Store the faces
    calculateFaceNormals();  // Calculate normals for each face
}

// Render the mesh using OpenGL
void Renderer::render() const {
    glPushMatrix();  // Save the current matrix state
 
    glTranslatef(-0.5, -0.5, -0.5);  // Translate the model to center it
    glBegin(GL_TRIANGLES);  // Begin drawing triangles

    // Loop over each face in the mesh
    for (size_t i = 0; i < faces.size(); ++i) {
        const auto& f = faces[i];  // Get the current face (indices of vertices)
        const Vec3& normal = fnormals[i];  // Get the normal for this face

        Vec3 color = fnormals[i];  // Use the normal as the color (for debugging purposes)

        // Calculate color components based on the normal direction (for visualization)
        float r = 0.5 + (color[0] >= 0 ? 1 : -1) * 0.5 * fabs(color[0]);
        float g = 0.5 + (color[1] >= 0 ? 1 : -1) * 0.5 * fabs(color[1]);
        float b = 0.5 + (color[2] >= 0 ? 1 : -1) * 0.5 * fabs(color[2]);

        glColor3f(r, g, b);  // Set the color for the current face
 
        // Set the normal for the current face
        glNormal3f(static_cast<GLfloat>(normal[0]), static_cast<GLfloat>(normal[1]), static_cast<GLfloat>(normal[2]));

        // Loop over the three vertices of the face
        for (int j = 0; j < 3; ++j) {
            const Vec3& vertex = vertices[f[j]];  // Get the vertex by index
            glVertex3f(static_cast<GLfloat>(vertex[0]), static_cast<GLfloat>(vertex[1]), static_cast<GLfloat>(vertex[2]));  // Draw the vertex
        }
    }

    glEnd();  // End drawing triangles
    glTranslatef(0.5, 0.5, 0.5);  // Translate the model back
    glPopMatrix();  // Restore the previous matrix state
}

// Get the vertices of the mesh (const reference to avoid copying)
const std::vector<Vec3>& Renderer::getVertices() const {
    return vertices;
}

// Calculate the normal vectors for each face in the mesh
void Renderer::calculateFaceNormals() {
    fnormals.clear();  // Clear any existing normals
    for (const auto& f : faces) {
        // Calculate the normal for the current face
        Vec3 normal = NormalCalculator::calculateNormal(vertices[f[0]], vertices[f[1]], vertices[f[2]]);
        fnormals.push_back(normal);  // Store the normal
    }

    // Debugging information: print the normal for each face
    for (size_t i = 0; i < fnormals.size(); ++i) {
        const Vec3& normal = fnormals[i];
        std::cout << "Face " << i << " normal: " << normal[0] << ", " << normal[1] << ", " << normal[2] << std::endl;
    }
}
