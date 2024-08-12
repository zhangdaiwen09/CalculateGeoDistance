#ifndef RENDERER_HPP
#define RENDERER_HPP

#include "Vec3.h"
#include "NormalCalculator.h"
#include <vector>
#include <array>
#include <GL/glew.h>

// Class for rendering a 3D mesh, storing vertices, faces, and face normals
class Renderer {
private:
    // Vector to store the vertices of the mesh
    std::vector<Vec3> vertices;

    // Vector to store the faces of the mesh, each face is represented by an array of 3 integers (vertex indices)
    std::vector<std::array<int, 3>> faces;

    // Vector to store the normals of the faces in the mesh
    std::vector<Vec3> fnormals;

public:
    // Function to set the mesh data by providing vertices and faces
    void setMesh(const std::vector<Vec3>& verts, const std::vector<std::array<int, 3>>& fcs);

    // Function to render the mesh (const means this function does not modify any member variables)
    void render() const;

    // Function to get the vertices of the mesh (const reference, so the original vector is not copied)
    const std::vector<Vec3>& getVertices() const; // Newly added method

private:
    // Function to calculate the normals for each face in the mesh
    void calculateFaceNormals();
};


#endif // RENDERER_HPP
