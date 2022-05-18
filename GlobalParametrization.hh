//
// Created by wuschelbueb on 05.05.22.
//

#ifndef OPENFLIPPER_GLOBALPARAMETRIZATION_H
#define OPENFLIPPER_GLOBALPARAMETRIZATION_H

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <MeshTools/MeshSelectionT.hh>
#include <ACG/Geometry/Types/PlaneT.hh>
#include <ACG/Scenegraph/TransformNode.hh>
#include <ACG/Utils/ColorCoder.hh>
#include <CoMISo/Solver/ConstrainedSolver.hh>
#include <gmm/gmm.h>
#include <functional>
#include <iostream>
#include <vector>
#include <float.h>
#include <cmath>
#include "DijkstraDistance.hh"

class GlobalParametrization {
public:
    using Point = OpenMesh::Vec3d;
    /// Handy typedefs sparse gmm vector matrix types
    typedef gmm::col_matrix<gmm::wsvector<double>> CMatrixType;
    typedef gmm::row_matrix<gmm::wsvector<double>> RMatrixType;
    typedef gmm::wsvector<double> CVectorType;
public:

    /*
     * is the same as:>
     * Crossfield(TriMesh &trimesh, std::vector<int> &heInRange) {
     *      trimesh_ = trimesh;
     *      heInRange_ = heInRange;
     * }
     */
    GlobalParametrization(TriMesh &trimesh)
            : trimesh_{trimesh} {
    }

    ~GlobalParametrization() {
    }

    void getGlobalParam();

private:

    void meshCutting();

    void getFaceVec(std::vector<int> &faces);

    int createVertexPosParamDom(std::vector<int> &faces);

    void setUpLocFaceCoordSys(const std::vector<int> &faces);

    void createEdgesAndLocalVUi(const OpenMesh::FaceHandle fh, int &counter, std::vector<Point> &edges);

    void createBasisTfMtx(const OpenMesh::FaceHandle fh, gmm::col_matrix<std::vector<double> > &_C,
                          const std::vector<Point> &edges);

    void connectSingularityToCutGraph(std::vector<int> complementHEdges, DijkstraDistance &dualGraph);

    gmm::col_matrix<std::vector<double>> createCMatrix();

    std::vector<double> getRhs(const std::vector<int> &faces, const int rhsSize);

    CMatrixType getHessian(const std::vector<int> &faces, const int rhsSize);

    void getRhsEntryForVertex(const OpenMesh::FaceHandle fh, const Point CrossFieldAxis, const bool flagUorV,
                              std::vector<double> &_rhs);

    int mapLocCoordToGlobCoordSys(const OpenMesh::FaceHandle fh, const OpenMesh::VertexHandle vh);

    void getDiaEntriesHessian(const OpenMesh::FaceHandle fh, CMatrixType &_H);

    void getEntriesHessian(const OpenMesh::FaceHandle fh, CMatrixType &_H);

    void colorCompEdges(const std::vector<int> &complementEdges);

    void getComplementMeshSel(std::vector<int> &complementEdges);

    void removeOpenPaths(std::vector<int> &complementHEdges);

    void removeEdgeFromGraph(const int i, std::vector<int> &complementHEdges);

    void tagEdgesFromDualSpanningTree();

    TriMesh &trimesh_;

};


#endif //OPENFLIPPER_GLOBALPARAMETRIZATION_H
