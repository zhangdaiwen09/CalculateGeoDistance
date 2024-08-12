#include "Vec3.h"
#include <cmath>

// Default constructor initializes the vector to (0, 0, 0)
Vec3::Vec3() : data{ 0, 0, 0 } {}

// Constructor initializes the vector with specific x, y, z values
Vec3::Vec3(double x, double y, double z) : data{ x, y, z } {}

// Overload of the subscript operator to access vector elements by index (non-const)
double& Vec3::operator[](size_t i) {
    return data[i];
}

// Overload of the subscript operator to access vector elements by index (const)
const double& Vec3::operator[](size_t i) const {
    return data[i];
}

// Calculate the Euclidean norm (magnitude) of the vector
double Vec3::norm() const {
    return std::sqrt(data[0] * data[0] + data[1] * data[1] + data[2] * data[2]);
}

// Overload of the subtraction operator to subtract one vector from another
Vec3 Vec3::operator-(const Vec3& other) const {
    return Vec3(data[0] - other[0], data[1] - other[1], data[2] - other[2]);
}

// Overload of the addition operator to add two vectors together
Vec3 Vec3::operator+(const Vec3& other) const {
    return Vec3(data[0] + other[0], data[1] + other[1], data[2] + other[2]);
}

// Overload of the multiplication operator to scale the vector by a scalar
Vec3 Vec3::operator*(double scalar) const {
    return Vec3(data[0] * scalar, data[1] * scalar, data[2] * scalar);
}

// Overload of the division operator to scale the vector by a scalar
Vec3 Vec3::operator/(double scalar) const {
    return Vec3(data[0] / scalar, data[1] / scalar, data[2] / scalar);
}

// Overload of the multiplication operator to scale the vector by a scalar (scalar on the left)
Vec3 operator*(double scalar, const Vec3& v) {
    return Vec3(v[0] * scalar, v[1] * scalar, v[2] * scalar);
}

// Overload of the equality operator to compare two vectors for equality
bool Vec3::operator==(const Vec3& other) const {
    // Check if the difference between the two vectors is within a small threshold
    return (*this - other).norm() < SMALL_DATASETS;
}

// Overload of the inequality operator to compare two vectors for inequality
bool Vec3::operator!=(const Vec3& other) const {
    return !(*this == other);
}

// Compute the cross product of two vectors
Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3(
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    );
}

// Compute the dot product of two vectors
double dot(const Vec3& a, const Vec3& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Normalize the vector to make it a unit vector
void unitize(Vec3& v) {
    double length = v.norm();
    if (length > 0) {
        v = v / length;
    }
}

// Overload of the stream insertion operator to output the vector to a stream
std::ostream& operator<<(std::ostream& os, const Vec3& vec) {
    os << "(" << vec[0] << ", " << vec[1] << ", " << vec[2] << ")";
    return os;
}
