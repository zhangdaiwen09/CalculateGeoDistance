#ifndef OBJLOADER_HPP
#define OBJLOADER_HPP

#include "Vec3.h"
#include <vector>
#include <string>
#include <array>

// Class for loading OBJ files and storing vertices and faces
class OBJLoader {
public:
    // Vector to store the vertices loaded from the OBJ file
    std::vector<Vec3> vertices;

    // Vector to store the faces loaded from the OBJ file (each face is an array of 3 integers representing vertex indices)
    std::vector<std::array<int, 3>> faces;

    // Function to clear the vertices and faces vectors
    void clear();

    // Function to load an OBJ file given its filename
    void load(const std::string& filename);
};


#endif // OBJLOADER_HPP
