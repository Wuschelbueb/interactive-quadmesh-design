#ifndef OPENFLIPPER_DIJKSTRADISTANCE_HH
#define OPENFLIPPER_DIJKSTRADISTANCE_HH

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <OpenFlipper/BasePlugin/BaseInterface.hh>
#include <OpenFlipper/BasePlugin/PluginFunctions.hh>
#include <MeshTools/MeshSelectionT.hh>
#include <iostream>
#include <vector>
#include <float.h>
#include <algorithm>


class DijkstraDistance {
public:
    using Point = ACG::Vec3d;
public:
    DijkstraDistance(TriMesh &trimesh) : trimesh_{trimesh} {}

    ~DijkstraDistance() {
    }

public:

    void getDualGraph(const std::vector<int> &faces);

    void getDijkstraSingularities(std::vector<int> &complementHEdges, std::vector<int> &singularities);

    std::vector<int> createVerticesVector(std::vector<int> &complementHEdges, std::vector<int> &singularities);

    void initVertexProp(std::vector<int> &dualGraphVertices, const bool flag);

    int vertexGetSmallestDist();

    void calculateVDijkstra(const int i);

    void addPathToCutGraph(std::vector<int> &dualGraphVertices, const int i);

    std::vector<int>
    calculateDijkstra(const std::vector<int> HeConstraints, const double refDist, const bool includeBoundary);

    std::vector<int> getHeVectorOfSelection(const std::vector<int> &includedNodes);

    void colorizeEdges(const std::vector<int> &includedHEdges);

    std::vector<int> getHeConstraints();

private:

    int dualGraphGetSmallestDist(const std::vector<int> &faces);

    void includeBoundaryFaces(std::vector<int> &includedFaces, const double refDist);

    std::vector<int> transformHehToFaces(const std::vector<int> &constraintHeh);

    void dijkstraDistBaryCenter(std::vector<int> &includedNodes, const double refDist);

    int getSmallestDistProp(const double refDist);

    void getSelectedVertices(std::vector<int> &constraints);

    void getSelectedEdges(std::vector<int> &constraints);

    void getSelectedHEdges(std::vector<int> &constraints);

    void getSelectedFaces(std::vector<int> &constraints);

    TriMesh &trimesh_;
};


#endif //OPENFLIPPER_DIJKSTRADISTANCE_HH
