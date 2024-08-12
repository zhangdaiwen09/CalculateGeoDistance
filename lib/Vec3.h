#ifndef VEC3_H
#define VEC3_H

#include <array>
#include <iostream>

#define SMALL_DATASETS 1e-6;

// Class representing a 3D vector
class Vec3 {
public:
    // Default constructor initializing the vector to (0, 0, 0)
    Vec3();

    // Constructor initializing the vector with specific x, y, z components
    Vec3(double x, double y, double z);

    // Overload of the subscript operator for non-const objects, allows modification of elements
    double& operator[](size_t i);

    // Overload of the subscript operator for const objects, does not allow modification of elements
    const double& operator[](size_t i) const;

    // Function to calculate the Euclidean norm (magnitude) of the vector
    double norm() const;

    // Overload of the subtraction operator to subtract two vectors
    Vec3 operator-(const Vec3& other) const;

    // Overload of the addition operator to add two vectors
    Vec3 operator+(const Vec3& other) const;

    // Overload of the multiplication operator to scale the vector by a scalar
    Vec3 operator*(double scalar) const;

    // Overload of the division operator to scale the vector by a scalar
    Vec3 operator/(double scalar) const;

    // Overload of the equality operator to compare two vectors
    bool operator==(const Vec3& other) const;

    // Overload of the inequality operator to compare two vectors
    bool operator!=(const Vec3& other) const;

    // Friend function to multiply a vector by a scalar (scalar on the left)
    friend Vec3 operator*(double scalar, const Vec3& v);

    // Friend function to compute the cross product of two vectors
    friend Vec3 cross(const Vec3& a, const Vec3& b);

    // Friend function to compute the dot product of two vectors
    friend double dot(const Vec3& a, const Vec3& b);

    // Friend function to normalize (unitize) the vector
    friend void unitize(Vec3& v);

    // Friend function to output the vector to an output stream
    friend std::ostream& operator<<(std::ostream& os, const Vec3& vec);

    // Declare Vec3Hash as a friend struct to allow access to private members
    friend struct Vec3Hash;

private:
    // Array to store the x, y, z components of the vector
    std::array<double, 3> data;
};


#endif // VEC3_H
