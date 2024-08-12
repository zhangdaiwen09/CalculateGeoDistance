#ifndef NORMALCALCULATOR_HPP
#define NORMALCALCULATOR_HPP

#include "Vec3.h"

// Class for calculating normals of a 3D triangle
class NormalCalculator {
public:
    // Static function to calculate the normal vector of a triangle given three vertices
    static Vec3 calculateNormal(const Vec3& v1, const Vec3& v2, const Vec3& v3);
};

#endif // NORMALCALCULATOR_HPP
