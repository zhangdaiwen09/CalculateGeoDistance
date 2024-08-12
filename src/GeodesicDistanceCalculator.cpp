

#include "GeodesicDistanceCalculator.h"
#include <iostream>
#include <algorithm>
#include <limits>
#include <cmath>
#include <queue>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
#include <NormalCalculator.h>

// Constructor for GeodesicDistanceCalculator, initializes vertices, faces, and geodesic distances
GeodesicDistanceCalculator::GeodesicDistanceCalculator(const std::vector<Vec3>& verts, const std::vector<std::array<int, 3>>& fcs)
    : vertices(verts), faces(fcs), geoDistances(verts.size(), std::numeric_limits<double>::max()) {
    std::cout << "Initializing GeodesicDistanceCalculator..." << std::endl;
    initialize();  // Calls the initialize function
}

// Initialize the GeodesicDistanceCalculator, primarily by removing duplicate points
void GeodesicDistanceCalculator::initialize() {
    std::cout << "Removing duplicate points..." << std::endl;
    removeDuplicatePoints();  // Remove duplicate vertices
}

// Initialize the adjacency list based on vertices and faces
void GeodesicDistanceCalculator::initializeAdjacencyList() {
    adjacencyList.clear();  // Clear existing adjacency list

    // Compute face normals for each face
    std::vector<Vec3> faceNormals(faces.size());
    for (size_t i = 0; i < faces.size(); ++i) {
        const auto& face = faces[i];
        faceNormals[i] = NormalCalculator::calculateNormal(vertices[face[0]], vertices[face[1]], vertices[face[2]]);
    }

    // Populate the adjacency list by checking coplanarity
    for (size_t i = 0; i < faces.size(); ++i) {
        const auto& face = faces[i];
        std::unordered_set<int> coplanarVertices = { face[0], face[1], face[2] };

        for (size_t j = 0; j < faces.size(); ++j) {
            if (i != j && areCoplanar(faceNormals[i], faceNormals[j])) {
                const auto& otherFace = faces[j];
                if (isCoplanar(faceNormals[i], vertices[face[0]], vertices[otherFace[0]], vertices[otherFace[1]], vertices[otherFace[2]])) {
                    coplanarVertices.insert(otherFace[0]);
                    coplanarVertices.insert(otherFace[1]);
                    coplanarVertices.insert(otherFace[2]);
                }
            }
        }

        std::vector<int> coplanarVerticesVector(coplanarVertices.begin(), coplanarVertices.end());
        for (size_t m = 0; m < coplanarVerticesVector.size(); ++m) {
            for (size_t n = m + 1; n < coplanarVerticesVector.size(); ++n) {
                int v1 = coplanarVerticesVector[m];
                int v2 = coplanarVerticesVector[n];
                if (v1 != v2) {
                    double length = (vertices[v1] - vertices[v2]).norm();
                    adjacencyList[v1][v2] = length;
                    adjacencyList[v2][v1] = length;
                }
            }
        }
    }
    std::cout << "Adjacency list initialized." << std::endl;
}

// Reinitialize the adjacency list to include new points and projections
void GeodesicDistanceCalculator::reinitializeAdjacencyListForNewPoints() {
    std::unordered_set<std::pair<int, int>, pair_hash> addedEdges;

    for (const auto& face : faces) {
        std::vector<int> faceVertices = { face[0], face[1], face[2] };
        for (size_t i = 3; i < vertices.size(); ++i) {
            if (isCoplanar(NormalCalculator::calculateNormal(vertices[face[0]], vertices[face[1]], vertices[face[2]]), vertices[face[0]], vertices[face[1]], vertices[face[2]], vertices[i])) {
                faceVertices.push_back(static_cast<int>(i));
            }
        }

        for (size_t i = 0; i < faceVertices.size(); ++i) {
            for (size_t j = i + 1; j < faceVertices.size(); ++j) {
                int v1 = faceVertices[i];
                int v2 = faceVertices[j];
                if (v1 != v2 && addedEdges.find({ v1, v2 }) == addedEdges.end() && addedEdges.find({ v2, v1 }) == addedEdges.end()) {
                    double length = (vertices[v1] - vertices[v2]).norm();
                    adjacencyList[v1][v2] = length;
                    adjacencyList[v2][v1] = length;
                    addedEdges.insert({ v1, v2 });
                    addedEdges.insert({ v2, v1 });
                }
            }
        }
    }
    std::cout << "Reinitialized adjacency list for new points." << std::endl;
}

// Check if two face normals are coplanar by checking if their cross product is near zero
bool GeodesicDistanceCalculator::areCoplanar(const Vec3& normal1, const Vec3& normal2) {
    return cross(normal1, normal2).norm() < SMALL_DATASETS;
}

// Check if four points are coplanar by measuring their distance from a plane defined by the first three points
bool GeodesicDistanceCalculator::isCoplanar(const Vec3& normal, const Vec3& point1, const Vec3& point2, const Vec3& point3, const Vec3& point4) {
    auto pointToPlaneDistance = [&normal, &point1](const Vec3& point) {
        Vec3 vec = point - point1;
        return std::fabs(dot(vec, normal));
    };

    double threshold = 1e-6;
    return pointToPlaneDistance(point2) < threshold && pointToPlaneDistance(point3) < threshold && pointToPlaneDistance(point4) < threshold;
}

// Update the adjacency list with projections of a point onto the edges of faces
void GeodesicDistanceCalculator::updateAdjacencyListWithProjections(size_t virtualVertexIndex, const Vec3& point) {
    std::unordered_map<Vec3, int, Vec3Hash> uniquePoints;
    for (const auto& face : faces) {
        Vec3 faceCenter = calculateFaceCenter(face);
        if (uniquePoints.find(faceCenter) == uniquePoints.end()) {
            uniquePoints[faceCenter] = static_cast<int>(vertices.size());
            addNewPoint(faceCenter);
        }

        for (int i = 0; i < 3; ++i) {
            Vec3 edgeMidpoint = calculateEdgeMidpoint(vertices[face[i]], vertices[face[(i + 1) % 3]]);
            if (uniquePoints.find(edgeMidpoint) == uniquePoints.end()) {
                uniquePoints[edgeMidpoint] = static_cast<int>(vertices.size());
                addNewPoint(edgeMidpoint);
            }

            Vec3 projection = projectPointOntoEdge(point, vertices[face[i]], vertices[face[(i + 1) % 3]]);
            if (uniquePoints.find(projection) == uniquePoints.end()) {
                uniquePoints[projection] = static_cast<int>(vertices.size());
                addNewPoint(projection);
            }
        }
    }

    std::cout << "Updated adjacency list with projections for virtual vertex index " << virtualVertexIndex << std::endl;
}

// Project a point onto an edge defined by two vertices
Vec3 GeodesicDistanceCalculator::projectPointOntoEdge(const Vec3& point, const Vec3& edgeStart, const Vec3& edgeEnd) {
    Vec3 edge = edgeEnd - edgeStart;
    Vec3 pointToStart = point - edgeStart;
    double t = dot(pointToStart, edge) / dot(edge, edge);
    t = std::max(0.0, std::min(1.0, t));
    return edgeStart + (t * edge);
}

// Calculate the center of a face (centroid) given its vertices
Vec3 GeodesicDistanceCalculator::calculateFaceCenter(const std::array<int, 3>& face) {
    Vec3 v1 = vertices[face[0]];
    Vec3 v2 = vertices[face[1]];
    Vec3 v3 = vertices[face[2]];
    return (v1 + v2 + v3) / 3.0;
}

// Calculate the midpoint of an edge defined by two vertices
Vec3 GeodesicDistanceCalculator::calculateEdgeMidpoint(const Vec3& v1, const Vec3& v2) {
    return (v1 + v2) / 2.0;
}

// Add a new point to the adjacency list, also updating the list with projections
size_t GeodesicDistanceCalculator::addPointToAdjacencyList(const Vec3& point) {
    for (size_t i = 0; i < vertices.size(); ++i) {
        if (vertices[i] == point) { // If point already exists, return its index
            std::cout << "Point already exists in vertices, point: " << point << std::endl;
            return i;
        }
    }
    addNewPoint(point);
    size_t virtualVertexIndex = vertices.size() - 1;
    updateAdjacencyListWithProjections(virtualVertexIndex, point);
    reinitializeAdjacencyListForNewPoints(); // Call the reinitialization function
    return virtualVertexIndex;
}

// Remove duplicate points from the vertices list
void GeodesicDistanceCalculator::removeDuplicatePoints() {
    std::unordered_map<Vec3, int, Vec3Hash> uniqueVertices;
    std::vector<Vec3> newVertices;
    std::vector<double> newGeoDistances;

    for (const auto& vertex : vertices) {
        if (uniqueVertices.find(vertex) == uniqueVertices.end()) {
            uniqueVertices[vertex] = static_cast<int>(newVertices.size());
            newVertices.push_back(vertex);
            newGeoDistances.push_back(std::numeric_limits<double>::max());
        }
    }

    vertices = std::move(newVertices);
    geoDistances = std::move(newGeoDistances);
    std::cout << "Removed duplicate points, new vertex count: " << vertices.size() << std::endl;
}

// Rebuild the adjacency list after changes to the mesh
void GeodesicDistanceCalculator::rebuildAdjacencyList() {
    adjacencyList.clear();
    initializeAdjacencyList();
    std::cout << "Rebuilt adjacency list." << std::endl;
}

// Add a new point to the internal data structures (vertices, geoDistances)
void GeodesicDistanceCalculator::addNewPoint(const Vec3& point) {
    vertices.push_back(point);
    geoDistances.push_back(std::numeric_limits<double>::max());
    size_t virtualVertexIndex = vertices.size();
    std::cout << "Added point to vertices and geoDistances, virtual vertex index: " << virtualVertexIndex << " " << point << std::endl;
}

// Calculate the normal of a triangle defined by three vertices
Vec3 calTriNormal(Vec3 ver1, Vec3 ver2, Vec3 ver3) {
    double temp1[3], temp2[3], normal[3];
    double length = 0.0;
    temp1[0] = ver2[0] - ver1[0];
    temp1[1] = ver2[1] - ver1[1];
    temp1[2] = ver2[2] - ver1[2];
    temp2[0] = ver3[0] - ver2[0];
    temp2[1] = ver3[1] - ver2[1];
    temp2[2] = ver3[2] - ver2[2];

    // Compute the normal vector using the cross product
    normal[0] = temp1[1] * temp2[2] - temp1[2] * temp2[1];
    normal[1] = -(temp1[0] * temp2[2] - temp1[2] * temp2[0]);
    normal[2] = temp1[0] * temp2[1] - temp1[1] * temp2[0];

    // Normalize the normal vector
    length = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
    if (length == 0.0f) { length = 1.0f; }
    normal[0] /= length;
    normal[1] /= length;
    normal[2] /= length;

    Vec3 e_normal(normal[0], normal[1], normal[2]);
    return e_normal;
}

// Calculate the plane equation from three points on the plane
bool ThreePt2PanelEquation(const Vec3& p1, const Vec3& p2, const Vec3& p3, double& a, double& b, double& c, double& d) {   
    a = ((p2[1] - p1[1]) * (p3[2] - p1[2]) - (p2[2] - p1[2]) * (p3[1] - p1[1]));
    b = ((p2[2] - p1[2]) * (p3[0] - p1[0]) - (p2[0] - p1[0]) * (p3[2] - p1[2]));
    c = ((p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0]));
    d = (0 - (a * p1[0] + b * p1[1] + c * p1[2]));

    return true;
}

// Calculate the distance from a point to a plane defined by ax + by + cz + d = 0
double Point2PanelDistance(const Vec3& p, double a, double b, double c, double d) {
    return fabs(a * p[0] + b * p[1] + c * p[2] + d) / sqrt(a * a + b * b + c * c);
}

// Calculate the distance from a point to a plane defined by three points
double Point2PanelDistance(const Vec3& p, const Vec3& p1, const Vec3& p2, const Vec3& p3) {
    double a, b, c, d;
    ThreePt2PanelEquation(p1, p2, p3, a, b, c, d);
    return Point2PanelDistance(p, a, b, c, d);
}

// Check if a point exists in the vertices, if not, modify faces to include it
bool GeodesicDistanceCalculator::exchange(const Vec3& point, std::vector<Vec3>& vertices, std::vector<std::array<int, 3>>& fcs) {
    size_t i = 0;
    for (i = 0; i < vertices.size(); ++i) {
        if (vertices[i] == point) { // If point already exists, return its index
            std::cout << "Point already exists in vertices, point: " << point << std::endl;
            return true;
        }
    }

    std::array<int, 3> face;
    Vec3 p1, p2, p3;
    for (i = 0; i < fcs.size(); ++i) {        
        face = fcs[i];
        p1 = vertices[face[0]];
        p2 = vertices[face[1]];
        p3 = vertices[face[2]];

        double d = Point2PanelDistance(point, p1, p2, p3);

        if (d < 0.01) {  // If point is close enough to the plane of the face
            break;
        }
    }

    if (i < fcs.size()) {
        vertices.push_back(point);
        int index = vertices.size() - 1;
        fcs.erase(fcs.begin() + i);
        fcs.push_back({ index , face[0], face[1] });
        fcs.push_back({ index , face[0], face[2] });
        fcs.push_back({ index , face[1], face[2] });

        std::cout << "face match : " << i << std::endl;

        return true;
    }

    std::cout << "face no match " << std::endl;
    return false;
}

// Compute the geodesic distance between two points by adding them to the adjacency list
double GeodesicDistanceCalculator::computeGeodesicDistance(const Vec3& startPoint, const Vec3& endPoint) {
    size_t startVertex = addPointToAdjacencyList(startPoint);
    size_t endVertex = addPointToAdjacencyList(endPoint);

    std::cout << "Start vertex index: " << startVertex << ", End vertex index: " << endVertex << std::endl;

    computeGeodesicDistance({ static_cast<int>(startVertex) }, { static_cast<int>(endVertex) });

    std::vector<int> path;
    getShortestPath(static_cast<int>(startVertex), static_cast<int>(endVertex), path);

    std::cout << "Computed geodesic distance, path length: " << path.size() << std::endl;
    if (!path.empty()) {
        std::cout << "Shortest path: ";
        for (const int vertex : path) {
            std::cout << vertex << " ";
        }
        std::cout << std::endl;
    } else {
        std::cout << "No path found." << std::endl;
    }

    return getGeodesicDistance(static_cast<int>(startVertex), static_cast<int>(endVertex));
}

// Get the geodesic distance between two vertices, return a large value if indices are invalid
double GeodesicDistanceCalculator::getGeodesicDistance(int startVertex, int endVertex) {
    if (startVertex >= static_cast<int>(geoDistances.size()) || endVertex >= static_cast<int>(geoDistances.size())) {
        return std::numeric_limits<double>::max();
    }
    return geoDistances[endVertex];
}

// Get the number of vertices in the mesh
int GeodesicDistanceCalculator::getVerticesSize() {
    return static_cast<int>(vertices.size());
}

// Compute the geodesic distance from a set of start vertices to a set of end vertices
void GeodesicDistanceCalculator::computeGeodesicDistance(const std::vector<int>& startVertices, const std::vector<int>& endVertices) {
    std::unordered_set<int> visited;
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<>> pq;
    std::unordered_map<int, double> dist;
    std::unordered_map<int, int> prev;

    for (int vertex : startVertices) {
        geoDistances[vertex] = 0.0;
        dist[vertex] = 0.0;
        pq.push({ 0.0, vertex });
    }

    while (!pq.empty()) {
        auto [currentDistance, u] = pq.top();
        pq.pop();

        if (visited.find(u) != visited.end()) {
            continue;
        }
        visited.insert(u);

        if (adjacencyList.find(u) != adjacencyList.end()) {
            for (const auto& [v, weight] : adjacencyList[u]) {
                double alt = dist[u] + weight;
                if (alt < dist[v] || dist.find(v) == dist.end()) {
                    dist[v] = alt;
                    prev[v] = u;
                    pq.push({ alt, v });
                }
            }
        }
    }

    for (const auto& endVertex : endVertices) {
        if (prev.find(endVertex) == prev.end()) {
            std::cout << "No path to vertex " << endVertex << std::endl;
        } else {
            std::cout << "Path to vertex " << endVertex << " found." << std::endl;
        }
    }

    std::cout << "Computed geodesic distances from start vertices." << std::endl;
}

// Add a virtual point to the mesh and assign its index
void GeodesicDistanceCalculator::addVirtualPoint(const Vec3& point, int& virtualIndex) {
    for (size_t i = 0; i < vertices.size(); ++i) {
        if ((vertices[i] - point).norm() < 1e-6) { // If point already exists, return its index
            virtualIndex = static_cast<int>(i);
            std::cout << "Point already exists in vertices, index: " << virtualIndex << std::endl;
            return;
        }
    }
    vertices.push_back(point);
    virtualIndex = static_cast<int>(vertices.size() - 1);
    geoDistances.push_back(std::numeric_limits<double>::max());
    std::cout << "Added point to vertices and geoDistances, virtual vertex index: " << virtualIndex << std::endl;
}

// Merge vertices and connect them based on coplanarity
void GeodesicDistanceCalculator::mergeAndConnectVertices() {
    std::unordered_set<std::pair<int, int>, pair_hash> addedEdges;
    std::vector<int> virtualIndices(vertices.size() - faces.size());
    std::iota(virtualIndices.begin(), virtualIndices.end(), static_cast<int>(faces.size()));

    for (const auto& face : faces) {
        std::vector<int> faceVertices = { face[0], face[1], face[2] };
        for (int virtualIndex : virtualIndices) {
            if (isCoplanar(NormalCalculator::calculateNormal(vertices[face[0]], vertices[face[1]], vertices[face[2]]), vertices[face[0]], vertices[face[1]], vertices[face[2]], vertices[virtualIndex])) {
                faceVertices.push_back(virtualIndex);
            }
        }

        for (size_t i = 0; i < faceVertices.size(); ++i) {
            for (size_t j = i + 1; j < faceVertices.size(); ++j) {
                int v1 = faceVertices[i];
                int v2 = faceVertices[j];
                if (v1 != v2 && addedEdges.find({ v1, v2 }) == addedEdges.end() && addedEdges.find({ v2, v1 }) == addedEdges.end()) {
                    double length = (vertices[v1] - vertices[v2]).norm();
                    adjacencyList[v1][v2] = length;
                    adjacencyList[v2][v1] = length;
                    addedEdges.insert({ v1, v2 });
                    addedEdges.insert({ v2, v1 });
                    std::cout << "Added edge (" << v1 << ", " << v2 << ") with length " << length << std::endl;
                }
            }
        }
    }
    std::cout << "Merged and connected vertices." << std::endl;
}

// Get the shortest path between two vertices
void GeodesicDistanceCalculator::getShortestPath(int startVertex, int endVertex, std::vector<int>& path) {
    std::unordered_map<int, int> previous;
    std::unordered_set<int> visited;
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<>> pq;

    pq.push({ 0.0, startVertex });
    previous[startVertex] = -1;

    while (!pq.empty()) {
        auto [currentDistance, u] = pq.top();
        pq.pop();

        if (u == endVertex) {
            break;
        }

        if (visited.find(u) != visited.end()) {
            continue;
        }
        visited.insert(u);

        if (adjacencyList.find(u) != adjacencyList.end()) {
            for (const auto& [v, weight] : adjacencyList[u]) {
                if (visited.find(v) == visited.end() && geoDistances[u] + weight < geoDistances[v]) {
                    geoDistances[v] = geoDistances[u] + weight;
                    pq.push({ geoDistances[v], v });
                    previous[v] = u;
                }
            }
        }
    }

    if (previous.find(endVertex) != previous.end()) {
        int current = endVertex;
        while (current != -1) {
            path.push_back(current);
            current = previous[current];
        }
        std::reverse(path.begin(), path.end());
    } else {
        std::cout << "No path found in getShortestPath" << std::endl;
    }
}

// Get the adjacency list, return as a const reference
const AdjacencyListMap& GeodesicDistanceCalculator::getAdjacencyList() const {
    return adjacencyList;
}
