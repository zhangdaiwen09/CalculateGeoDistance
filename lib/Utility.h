

#include <functional>
#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <functional>
#include "Vec3.h"

// Custom hash function struct for std::unordered_set to hash pairs of values
struct pair_hash {
    // Template operator() to hash a pair of elements
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1, T2>& p) const {
        // Hash the first element of the pair
        auto hash1 = std::hash<T1>{}(p.first);
        // Hash the second element of the pair
        auto hash2 = std::hash<T2>{}(p.second);
        // Combine the two hashes using XOR and return the result
        return hash1 ^ hash2;
    }
};

// Custom hash function struct for Vec3 objects
struct Vec3Hash {
    // Operator() to hash a Vec3 object
    std::size_t operator()(const Vec3& v) const {
        // Create a hasher for double values
        std::hash<double> hasher;
        // Hash each component of the Vec3 and combine them using XOR
        return hasher(v[0]) ^ hasher(v[1]) ^ hasher(v[2]);
    }
};


#endif // UTILITY_HPP
