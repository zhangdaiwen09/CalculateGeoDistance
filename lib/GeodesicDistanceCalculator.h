
#ifndef GEODESIC_DISTANCE_CALCULATOR_H
#define GEODESIC_DISTANCE_CALCULATOR_H

#include "Vec3.h"
#include <vector>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <Utility.h>

// Define a type alias for a nested unordered_map to represent the adjacency list
typedef std::unordered_map<int, std::unordered_map<int, double>> AdjacencyListMap;

// Class to calculate geodesic distances in a 3D mesh
class GeodesicDistanceCalculator {
public:
    // Constructor that takes vertices and faces as input
    GeodesicDistanceCalculator(const std::vector<Vec3>& verts, const std::vector<std::array<int, 3>>& fcs);

    // Function to compute the geodesic distance between two points
    double computeGeodesicDistance(const Vec3& startPoint, const Vec3& endPoint);

    // Function to get the adjacency list (as a constant reference)
    const AdjacencyListMap& getAdjacencyList() const;

    // Function to compute the geodesic distances for multiple start and end vertices
    void computeGeodesicDistance(const std::vector<int>& startVertices, const std::vector<int>& endVertices);

    // Function to add a virtual point and get its index
    void addVirtualPoint(const Vec3& point, int& virtualIndex);

    // Function to get the geodesic distance between two vertices by their indices
    double getGeodesicDistance(int startVertex, int endVertex);

    // Function to get the total number of vertices
    int getVerticesSize();

    // Function to get the shortest path between two vertices
    void getShortestPath(int startVertex, int endVertex, std::vector<int>& path);

    // Static function to modify vertices and faces; it returns a boolean
    static bool exchange(const Vec3& point, std::vector<Vec3>& verts, std::vector<std::array<int, 3>>& fcs);

    // Function to perform analysis by initializing the adjacency list
    void analysis() {   
        std::cout << "Initializing adjacency list..." << std::endl;  // Print initializing message
        initializeAdjacencyList();  // Initialize the adjacency list
    }

    // Function to output vertices and faces
    void data(std::vector<Vec3>& verts, std::vector<std::array<int, 3>>& fcs)
    {
        verts = vertices;  // Assign the class vertices to the provided verts vector
        fcs = faces;  // Assign the class faces to the provided fcs vector
    }

private:
    // Function to perform initial setup
    void initialize();

    // Function to initialize the adjacency list
    void initializeAdjacencyList();

    // Function to check if two normals are coplanar
    bool areCoplanar(const Vec3& normal1, const Vec3& normal2);

    // Function to check if four points are coplanar
    bool isCoplanar(const Vec3& normal, const Vec3& point1, const Vec3& point2, const Vec3& point3, const Vec3& point4);

    // Function to update the adjacency list with projections for a virtual vertex
    void updateAdjacencyListWithProjections(size_t virtualVertexIndex, const Vec3& point);

    // Function to project a point onto an edge defined by two points
    Vec3 projectPointOntoEdge(const Vec3& point, const Vec3& edgeStart, const Vec3& edgeEnd);

    // Function to calculate the center of a face
    Vec3 calculateFaceCenter(const std::array<int, 3>& face);

    // Function to calculate the midpoint of an edge
    Vec3 calculateEdgeMidpoint(const Vec3& v1, const Vec3& v2);

    // Function to add a point to the adjacency list and return its index
    size_t addPointToAdjacencyList(const Vec3& point);

    // Function to merge and connect vertices
    void mergeAndConnectVertices();

    // Function to remove duplicate points
    void removeDuplicatePoints();

    // Function to rebuild the adjacency list after modifications
    void rebuildAdjacencyList();

    // Function to add a new point to the internal data structures
    void addNewPoint(const Vec3& point);

    // Function to reinitialize the adjacency list when new points are added
    void reinitializeAdjacencyListForNewPoints();

    // Internal storage for vertices
    std::vector<Vec3> vertices;

    // Internal storage for faces
    std::vector<std::array<int, 3>> faces;

    // Storage for geodesic distances
    std::vector<double> geoDistances;

    // The adjacency list representing the mesh connectivity
    std::unordered_map<int, std::unordered_map<int, double>> adjacencyList;
};


#endif // GEODESIC_DISTANCE_CALCULATOR_H
