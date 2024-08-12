#include "NormalCalculator.h"

// Function to calculate the normal vector of a triangle given three vertices
Vec3 NormalCalculator::calculateNormal(const Vec3& v1, const Vec3& v2, const Vec3& v3) {
    // Calculate the cross product of two edges of the triangle to get the normal vector
    Vec3 normal = cross(v2 - v1, v3 - v1);

    // Normalize the normal vector to make it a unit vector
    unitize(normal);

    // Return the normalized normal vector
    return normal;
}
