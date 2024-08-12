#include "OBJLoader.h"
#include <fstream>
#include <iostream>
#include <algorithm>

// Clear the vertices and faces stored in the OBJLoader
void OBJLoader::clear() {
    vertices.clear();  // Clear the vertices vector
    faces.clear();     // Clear the faces vector
}

// Load a .obj file and populate the vertices and faces vectors
void OBJLoader::load(const std::string& filename) {
    std::ifstream fin(filename.c_str());  // Open the .obj file for reading
    if (!fin) {
        std::cerr << "Cannot open file " << filename << std::endl;  // Error if file can't be opened
        return;
    }
    OBJLoader::clear();  // Clear any existing data

    char line[1024];  // Buffer to read lines from the file
    while (fin.getline(line, 1024)) {  // Read the file line by line
        if (line[0] == 'v') {  // If the line starts with 'v', it's a vertex
            Vec3 v;
            sscanf_s(line, "v %lf %lf %lf", &v[0], &v[1], &v[2]);  // Parse the vertex coordinates
            vertices.push_back(v);  // Add the vertex to the vertices vector
        }
        else if (line[0] == 'f') {  // If the line starts with 'f', it's a face
            std::array<int, 3> f;
            sscanf_s(line, "f %d %d %d", &f[0], &f[1], &f[2]);  // Parse the face indices
            for (int& index : f) {
                index -= 1;  // Convert 1-based OBJ indices to 0-based indices
            }
            faces.push_back(f);  // Add the face to the faces vector
        }
    }

    fin.close();  // Close the file

    // Print Debug information: print the number of vertices and faces loaded
    std::cout << "Loaded " << vertices.size() << " vertices and " << faces.size() << " faces from " << filename << std::endl;
    
    // Debug information: print each vertex
    for (const auto& v : vertices) {
        std::cout << "Vertex: " << v[0] << ", " << v[1] << ", " << v[2] << std::endl;
    }
    
    // Print Debug information: print each face
    for (const auto& f : faces) {
        std::cout << "Face: " << f[0] << ", " << f[1] << ", " << f[2] << std::endl;
    }
}
