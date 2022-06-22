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

    std::vector<int> getFaceVec();

    int createVertexPosParamDomain(std::vector<int> &faces);

    void checkCGandSetPos(OpenMesh::VertexHandle fv_it, int &countVertices);

    void getPositionInnerNode(OpenMesh::VertexHandle &fv_it, int &countVertices);

    void setUpLocFaceCoordSys(const std::vector<int> &faces);

    void createEdgesAndLocalVUi(const OpenMesh::FaceHandle fh, std::vector<Point> &edges);

    void createBasisTfMtx(const OpenMesh::FaceHandle fh, gmm::col_matrix<std::vector<double> > &_C,
                          const std::vector<Point> &edges);

    std::vector<int> getSingularities();

    gmm::col_matrix<std::vector<double>> createCMatrix();

    std::vector<double> getRhs(const std::vector<int> &faces, const int rhsSizePartOne, const int rhsSizePartTwo);

    CMatrixType getHessian(const int rhsSizePartOne, const int rhsSizePartTwo);

    gmm::row_matrix<gmm::wsvector<double>>
    getConstraints(const int nbVerticesUaV, std::vector<int> &cutGraphWoBoundary,
                   std::vector<int> &singularities);

    int getRowSizeAndStatus();

    void setZeroPoint(std::vector<int> &singularities, gmm::row_matrix<gmm::wsvector<double>> &_constraints);

    void getConstraintsMatrix(int &jkStartCounter, std::vector<int> &cutGraphWoBoundary,
                              gmm::row_matrix<gmm::wsvector<double>> &_constraints);

    void setConRows(int &counter, int &jkStartCounter, const int diff, OpenMesh::SmartHalfedgeHandle &he,
                    gmm::row_matrix<gmm::wsvector<double>> &_constraints);

    int getPositionConstraintRow(OpenMesh::SmartVertexHandle &vh, int cutGraphZone);

    void getRhsEntryForVertex(const OpenMesh::FaceHandle fh, const Point CrossFieldAxis, const bool flagUorV,
                              std::vector<double> &_rhs);

    int mapLocCoordToGlobCoordSys(const OpenMesh::FaceHandle fh, const OpenMesh::VertexHandle vh);

    void getDiaEntriesHessian(const OpenMesh::FaceHandle fh, CMatrixType &_H);

    void getEntriesHessian(const OpenMesh::SmartHalfedgeHandle he, CMatrixType &_H);

    void colorCompHEdges(const std::vector<int> &complementEdges);

    std::vector<int> getComplementMeshSel();

    void removeOpenPaths(std::vector<int> &complementHEdges);

    void removeRedundantEdges(std::vector<int> &complementHEdges);

    void removeEdgeFromGraph(const int i, std::vector<int> &complementHEdges);

    void tagEdgesFromDualSpanningTree();

    void checkIfFaceInSelection(OpenMesh::FaceHandle &face);

    void checkIfEBetweenTriangleInDualGraph(OpenMesh::FaceHandle &face, OpenMesh::FaceHandle &fh_pred);

    void tagEdgeIfInDualGraph(TriMesh::FaceHalfedgeIter &fhe_pred_it, OpenMesh::HalfedgeHandle &oheh);

    void createSectorsCutGraph(std::vector<int> &singularities);

    void initVectorStartSec(const int i, std::vector<OpenMesh::HalfedgeHandle> &startOfSectors);

    void propagateForSectors(int &sector, const std::vector<OpenMesh::HalfedgeHandle> &startOfSectors,
                             std::vector<int> &singularities);

    void propagation(OpenMesh::HalfedgeHandle &heh, int &sector, const std::vector<int> &singularities);

    bool checkIfLeaf(const OpenMesh::VertexHandle &heToVertex);

    void fixRotationsCrossBoundaryComp(std::vector<int> &complementHEdges, std::vector<int> &singularities,
                                       std::vector<int> &faces);

    Point rotPointWithRotMatrix(const OpenMesh::FaceHandle fh, const Point vec, const double theta);

    std::pair<int, int> getPJ(TriMesh::FaceHalfedgeIter &fhe_it, OpenMesh::SmartHalfedgeHandle &ohe);

    void updatePJandCrossfield(std::pair<int, int> &pj, OpenMesh::SmartHalfedgeHandle &ohe);

    void updateStack(TriMesh::FaceHalfedgeIter &he, std::queue<OpenMesh::FaceHandle> &stack);

    void setFaceStatusToFalse();

    std::vector<int>
    getIdxToRound(int nbVerticesUaV, int jkValues, std::vector<int> &singularities, std::vector<int> &onlyBoundaries);

    TriMesh &trimesh_;
};


#endif //OPENFLIPPER_GLOBALPARAMETRIZATION_H
