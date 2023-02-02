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

/**
 * class which uses the dijkstra algorithm on the mesh with various different conditions.
 */
class DijkstraDistance {
public:
    using Point = ACG::Vec3d;

public:
    DijkstraDistance(TriMesh &trimesh) : trimesh_{trimesh} {}

    ~DijkstraDistance() {}

public:

    /**
     * removes all the custom properties from the mesh. This function only gets used at the start.\n
     * This makes sure that there are no old values which could interfere with new calculations.
    */
    void cleanMeshOfProps();

    /**
     * creates a dualGraph (dual spanning tree) starting from a random face from faces vector
     * additionally divides edges in three categories:\n
     * - outside of selection\n
     * - inside of selection\n
     * - on the border\n
     * @param faces
     */
    void getDualGraph(const std::vector<int> &faces);

    /**
     * uses complementHEdges and singularities to create the cutgraph with the help of Dijkstra.\n
     * adds halfedges to complementHEdges and cutGraphWoBoundary which get used again later.\n
     * @param complementHEdges edges which are not part of the dual spanning tree which got calculated with the getDualGraph
     * @param singularities vertices with a crossfield index != 0
     * @param cutGraphWoBoundary same vector as complementHEdges but without borders/boundaries
     */
    void calcDijkstraWSingularities(std::vector<int> &complementHEdges,
                                    std::vector<int> &singularities, std::vector<int> &cutGraphWoBoundary);

    /**
     * calculate Dijkstra in order to get selection which gets used for quadrangulation
     * @param HeConstraints halfedges which are constrained, dijkstra starts from these
     * @param refDist distance on how far the dijkstra algo needs to spread out
     * @param includeBoundary flag decides if it should include faces which are on the edge of refDist
     * @return
     */
    std::vector<int>
    calculateDijkstra(const std::vector<int> HeConstraints, const double refDist);

    /**
     * returns a vector of halfedges from dijkstra calculation .
     * @param includedFaces is a vector which contains all the faces from dijkstra calculation
     * @return a vector of halfedges
     */
    std::vector<int> getAllHeFromFaces(const std::vector<int> &includedFaces);

    /**
     * colorizes edges which got included in dijkstra calculation
     * @param includedHEdges vector of Hedges in selection
     */
    void colorizeEdges(const std::vector<int> &includedHEdges);

    /**
     * get halfedges which got selected with the selection tool
     * @return a vector with halfedges
     */
    std::vector<int> getHeFromVertex(OpenMesh::VertexHandle selectedVertex, const std::vector<int> &originVertices);

private:

    /**
     * gets the smallest distance (face property) from the origin.\n
     * this is how dijkstra propagates
     * @param faces checks all the faces and their distance
     * @return returns idx of face with shortest distance
     */
    int dualGraphGetSmallestDist(const std::vector<int> &faces);

    /**
     * colors DualGraph. There are 3 colors\n
     * 1 is for the border\n
     * 2 is for the selection inside the border\n
     * 3 is for the selection outside of the border
     */
    void colorDualGraph(const std::vector<int> &faces);

    /**
     * initialize different properties which get used to create the DualGraph
     * @param faces
     */
    void initDualGraphProp(const std::vector<int> &faces);

    /**
     * creates DualGraph with the help of the Dijkstra algorithm
     * @param faces
     */
    void calculateDGDijkstra(const std::vector<int> &faces);

    /**
     * creates vector of vertices. takes complementHEdges and singularities as input.\n
     * extracts from these two vectors the vertices
     * @param complementHEdges vector of halfedges
     * @param singularities vector of vertices
     * @return vector of vertices
     */
    std::vector<int> createVerticesVector(std::vector<int> &complementHEdges, std::vector<int> &singularities);

    /**
     * initializes vertex properties which then get used to calculate the cutgraph with the help of dijsktra.\n
     * @param cutGraphVertices vector of vertices
     * @param first_it this flag is only true for the first use of function
     */
    void initProperties(std::vector<int> &cutGraphVertices, const bool first_it);

    /**
     * get smallest distance of vertex selection. is used as a helper function of dijkstra.\n
     * @return return vertex with smallest distance
     */
    int vertexGetSmallestDist();

    /**
     * directed dijkstra calculation with vertices. distance between vertices is always 1.\n
     * has always the goal of reaching the singularity.\n
     * @param i is a index of a singularity
     */
    void calculateVDijkstra(const int i);

    /**
     * when dijkstra stops path to the singularity gets added to complementHEdges and cutGraphWoboundary.\n
     * vertex of singularity gets added to cutGraphVertcies for the next dijkstra iteration.
     * @param cutGraphVertices vector which gets used for dijkstra iteration
     * @param complementHEdges vector of halfedges
     * @param i singularity
     * @param cutGraphWoBoundary vector of halfedges
     */
    void addPathToCutGraph(std::vector<int> &cutGraphVertices, std::vector<int> &complementHEdges, const int i,
                           std::vector<int> &cutGraphWoBoundary);

    /**
     * transforms constraintHeh vector to a face vector. faces which are assigned to halfedges\n
     * initializes a lot of properties.\n
     * @param constraintHeh vector of halfedges
     * @return vector of faces
     */
    std::vector<int> transformHehToFaces(const std::vector<int> &constraintHeh);

    /**
     * calculates Dijkstra. uses distance from one face (barycenter) to another.\n
     * @param includedFaces vector which contains faces which get added by dijkstra
     * @param refDist reference distance
     */
    void dijkstraDistBaryCenter(std::vector<int> &includedFaces, const double refDist);

    /**
     * is a helper function of the Dijkstra calculation. gets the face with the smallest distance.\n
     * the distance can't be bigger than the reference distance.
     * @param refDist reference distance
     * @return index of face
     */
    int getSmallestDistProp(const double refDist);

    TriMesh &trimesh_;
};


#endif //OPENFLIPPER_DIJKSTRADISTANCE_HH
