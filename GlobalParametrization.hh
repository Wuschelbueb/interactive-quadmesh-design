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
#include <fstream>
#include <iostream>
#include <vector>
#include <float.h>
#include <cmath>
#include "DijkstraDistance.hh"

/**
 * compute global parametrization. prepare mesh for mapping on a disk-shaped param domain.\n
 */
class GlobalParametrization {
public:
    using Point = OpenMesh::Vec3d;
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
    GlobalParametrization(TriMesh &trimesh, const double &hValue)
            : trimesh_{trimesh}, hVal{hValue} {
    }

    ~GlobalParametrization() {
    }

    /**
     * get global parametrization computed.\n
     */
    void getGlobalParam();

private:

    TriMesh &trimesh_;
    const double &hVal;

    /**
     * get faces part of the selection.\n
     * @return vector with indices of faces
     */
    std::vector<int> getFaceVec();

    /**
     * create for each vertex of a face an entry in the global param domain.\n
     * @param faces vector of faces in selection
     * @return int, nb of uv values.
     */
    int createVertexPosParamDomain(std::vector<int> &faces);

    /**
     * check if fv_it is on cutgraph, and assign position accordingly.\n
     * @param fv_it vertex handle
     * @param countVertices int, nb of uv values.
     */
    void checkCGandSetPos(OpenMesh::VertexHandle fv_it, int &countVertices);

    /**
     * if vertex appears multiple time on global param domain count how many times and assign positions accordingly.\n
     * @param fv_it vertex handle
     * @param countVertices int, nb of uv values.
     */
    void getPositionInnerNode(OpenMesh::VertexHandle &fv_it, int &countVertices);

    /**
     * create a local coordsystem for each face.\n
     * calls createMatrix3v3, createEdgesAndLocalVUi and createBasisTfMtx.\n
     * @param faces vector of faces of selection
     */
    void setUpLocFaceCoordSys(const std::vector<int> &faces);

    /**
     * takes an empty vector and returns it with three edges.\n
     * edge zero and one are  adjacent edges of triangle and edge two is normal of these two.\n
     * @param fh face handle
     * @param edges empty vector
     */
    void createEdgesAndLocalVUi(const OpenMesh::FaceHandle fh, std::vector<Point> &edges);

    /**
     * create basis transformation matrix for each face.\n
     * takes edges vector, CMatrix and creates based on equations from jupyternotebook the basis transformation matrix.\n
     * @param fh face handle
     * @param _abc matrix from createMatrix3v3
     * @param edges vector with three edges
     */
    void createBasisTfMtx(const OpenMesh::FaceHandle fh, gmm::col_matrix<std::vector<double> > &_abc,
                          const std::vector<Point> &edges);

    /**
     * extract singularities from crossfield index.\n
     * @return vector with singularities
     */
    std::vector<int> getSingularities(std::vector<int> &faces);

    /**
     * is a 3x3 matrix. u_x calculations get expressed with this matrix.\n
     * @return return matrix
     */
    static gmm::col_matrix<std::vector<double>> createMatrix3v3();

    /**
     * get rhs vector.\n
     * calls getRhsEntryForVertex.\n
     * @param faces vector of faces in selection
     * @param rhsSizePartOne nb of uv values; based on number of vertices on disk topology.
     * @param rhsSizePartTwo nb of jk Values; based on the number of halfedges in cutgraph.
     * @return vector with rhs values.
     */
    std::vector<double> getRhs(const int rhsSizePartOne, const int rhsSizePartTwo);

    /**
     * get entry for each vertex.\n
     * dotproduct of crossfieldAxis and basisTransformationMatrix times some other minor components.\n
     * calls mapLocCoordToGlobCoordSys to map local vertices to global param.\n
     * @param he halfedge handle
     * @param CrossFieldAxis is Point of U or V Axis of local coord system
     * @param flagUorV flag differs between U or V
     * @param _rhs vector. add entries
     */
    void getRhsEntryForVertex(const OpenMesh::SmartHalfedgeHandle he, const Point CrossFieldAxis,
                              const bool flagUorV,
                              std::vector<double> &_rhs);

    /**
     * get Hessian matrix.\n
     * iterates through all halfedges of the mesh and if they are in face selection calls getDiaEntriesHessian and getEntriesHessian.\n
     * @param rhsSizePartOne nb of uv values; based on number of vertices on disk topology.
     * @param rhsSizePartTwo nb of jk Values; based on the number of halfedges in cutgraph.
     * @return returns hessian sparse column matrix.
     */
    GlobalParametrization::CMatrixType getHessian(const int rhsSizePartOne, const int rhsSizePartTwo,
                                                  std::vector<int> &cutGraphWoBoundary,
                                                  const std::vector<int> &faces);

    /**
     * responsible to fill up diagonal of hessian matrix.\n
     * dotproduct of basisTransformationMatrix[i] and basisTransformationMatrix[i] times some other minor components.\n
     * calls helper function mapLocCoordToGlobCoordSys.\n
     * @param he halfedge handle
     * @param _hessian hessian sparse column matrix.\n
     */
    void getDiaEntriesHessian(const OpenMesh::SmartHalfedgeHandle he, CMatrixType &_hessian);

    /**
     * fills the rest of the hessian matrix up.\n
     * dotproduct of basisTransformationMatrix[i] and basisTransformationMatrix[j] times some other minor components.\n
     * calls helper function mapLocCoordToGlobCoordSys.\n
     * @param he halfedge
     * @param _hessian essian sparse column matrix.\n
     */
    void getEntriesHessian(const OpenMesh::SmartHalfedgeHandle he, CMatrixType &_hessian);

    /**
     * creates empty constraint row matrix matrix.\n
     * calls helper functions getRowSizeAndSetStatus, setZeroPointConstraint and getConstraintsMatrix to fill matrix up.\n
     * @param nbVerticesUaV first half of column size.
     * @param cutGraphWoBoundary second half of column size.
     * @param singularities vector with vertices
     * @return constraints row matrix.
     */
    gmm::row_matrix<gmm::wsvector<double>>
    getConstraints(const int nbVerticesUaV, std::vector<int> &cutGraphWoBoundary,
                   std::vector<int> &onlyBoundaries, std::vector<int> &singularities,
                   std::vector<int> &faces);

    /**
     * counts the amount of vertices, including if appearing multiple times.\n
     * helper function of getConstraints.\n
     * @return number of vertices.
     */
    int getRowSizeAndSetStatus();

    /**
     * set a unique point, preferably a singularity which is the origin {0,0} of the global param system.\n
     * helper function of getConstraints.\n
     * @param singularities vector of vertices
     * @param _constraints row matrix
     */
    void setZeroPointConstraint(std::vector<int> &singularities,
                                gmm::row_matrix<gmm::wsvector<double>> &_constraints,
                                std::vector<int> &faces);

    /**
     * iterate through cutGraphWoBoundary to fill up rows of constraint matrix.\n
     * calls setConRows to fill rows.\n
     * helper function of getConstraints.\n
     * @param jkStartCounter
     * @param cutGraphWoBoundary vector with halfedges
     * @param _constraints row matrix
     */
    void getConstraintsMatrix(int &jkStartCounter, std::vector<int> &cutGraphWoBoundary,
                              gmm::row_matrix<gmm::wsvector<double>> &_constraints);

    /**
     * based on periodjump set up row with constraints.\n
     * there are for different cases for pj: 0, 90/-270, 180/-180 and 270/-90 degrees.\n
     * calls getPositionConstraintRow.\n
     * @param counter responsible for row in constraint matrix. gets increased
     * @param jkStartCounter where jk columns of constraint matrix start
     * @param diff periodjump difference
     * @param he halfedge
     * @param _constraints row matrix
     */
    void setConRows(int &counter, int &jkStartCounter, const int diff, OpenMesh::SmartHalfedgeHandle &he,
                    OpenMesh::SmartVertexHandle &vh,
                    gmm::row_matrix<gmm::wsvector<double>> &_constraints);

    /**
     * gets position of uv value, i.e. which column of the row the values need to be assigned.\n
     * gets all adjacent sectors, except 0. sorts them and returns position which matches sector of face.\n
     * @param vh vertex handle
     * @param cutGraphZone sector of face
     * @return position
     */
    int getPositionConstraintRow(OpenMesh::SmartVertexHandle &vh, int cutGraphZone);

    /**
     * sets up the mapping of the local coord system to the global param one.\n
     * each face has a localVUi saved as the local coord system. retrieve position of vertex.\n
     * @param fh face handle
     * @param vh vertex handle
     * @return position of vertex
     */
    int mapLocCoordToGlobCoordSys(const OpenMesh::FaceHandle fh, const OpenMesh::VertexHandle vh);

    /**
     * colors halfedges with true of false. if in cutgraph == true.\n
     * @param complementEdges vector of halfedges
     */
    void colorCompHEdges(const std::vector<int> &complementEdges);

    /**
     * gets the cutgraph, i.e. the primal of non-spanning tree edges.\n
     * add to cutgraph if edges aren't tagged, i.e. not part of dual spanning tree graph.\n
     * @return vector of halfedges
     */
    std::vector<int> getComplementMeshSel();

    /**
     * removes leaf paths, i.e. vertices in cutgraph with valence 1
     * @param complementHEdges vector of halfedges
     */
    void removeOpenPaths(std::vector<int> &complementHEdges);

    /**
     * remove border edges which aren't part of dual graph.
     * @param complementHEdges vector of halfedges
     */
    void removeRedundantEdges(std::vector<int> &complementHEdges);

    /**
     * helper function of removeOpenPaths. removes an edge which is in cutgraph and adjoins a vertex with valence 1.\n
     * @param i halfedge
     * @param complementHEdges vector of halfedges
     */
    void removeEdgeFromGraph(const int i, std::vector<int> &complementHEdges);

    /**
     * helper function of getComplementMeshSel.\n
     * iterates through faces and calls checkIfFaceInSelection.\n
     */
    void tagEdgesFromDualSpanningTree();

    /**
     * calls checkIfEBetweenTriangleInDualGraph.\n
     * if face is in selection and the dual graph pred is valid.\n
     * @param face face handle
     */
    void checkIfFaceInSelection(OpenMesh::FaceHandle &face);

    /**
     * calls checkIfEBetweenTriangleInDualGraph.\n
     * iterates through both triangles (face, fh_pred)
     * @param face face handle
     * @param fh_pred parent face of dual spanning tree face
     */
    void checkIfEBetweenTriangleInDualGraph(OpenMesh::FaceHandle &face, OpenMesh::FaceHandle &fh_pred);

    /**
     * tag both halfedges if they are opposite of each other.\n
     * @param fhe_pred_it halfedge of fh_pred
     * @param oheh opposite halfedge of face's halfedge
     */
    void tagEdgeIfInDualGraph(TriMesh::FaceHalfedgeIter &fhe_pred_it, OpenMesh::HalfedgeHandle &oheh);

    /**
     * create sectors on cutgraph in order to know how many times a vertex appears on cutgraph.\n
     * iterates through singularities and checks if they are leaf vertices. if so assign sectors from there.\n
     * calls initVectorStartSec and propagateForSectors.\n
     * @param singularities vertices with crossfield idx != 0
     */
    void createSectorsOnCutGraph(std::vector<int> &singularities);

    /**
     * checks if vertex is singularity on cutgraph.\n
     * @param heToVertex vertex handle
     * @return true if leaf else false
     */
    bool checkIfLeaf(const OpenMesh::VertexHandle &heToVertex);

    /**
     * comb the mesh. propagate through mesh via pj and rotate local coord system with the help of pj.\n
     * ingore pj on cutgraph but rotate the rest.\n
     * calls updateStack.
     * @param faces vector of faces in selection
     */
    void fixRotationsCrossBoundaryComp(std::vector<int> &faces);

    /**
     * rotate vec with Rodrigues rotation formula.\n
     * is used when a vector needs to be rotated a certain angle in space, given an axis and angle of rotation.\n
     * @param fh face handle
     * @param vec Vec3d
     * @param theta angle
     * @return Vec3d
     */
    Point rotPointWithRotMatrix(const OpenMesh::FaceHandle fh, const Point vec, const double theta);

    /**
     * adds adjacent face to stack if the following conditions are met:\n
     * halfedge and opposite halfedge are not part of the cutgraph, adjacent face is part of selection and isn't already tagged.\n
     * calls getPJ and updatePJandCrossfield.\n
     * @param he halfedge handle
     * @param stack queue of faces
     */
    void updateStack(TriMesh::FaceHalfedgeIter &he, std::queue<OpenMesh::FaceHandle> &stack);

    /**
     * checks both halfedges between the two faces and returns the one with the pj.\n
     * @param fhe_it halfedge 1
     * @param ohe halfedge 2
     * @return pair of <halfedge idx, periodjump>
     */
    std::pair<int, int> getPJ(TriMesh::FaceHalfedgeIter &fhe_it, OpenMesh::SmartHalfedgeHandle &ohe);

    /**
     * update the periodjump and rotate crossfield accordingly.\n
     * @param pj pair of <halfedge idx, periodjump>
     * @param ohe halfede of ajdacient face
     */
    void updatePJandCrossfield(std::pair<int, int> &pj, OpenMesh::SmartHalfedgeHandle &ohe);

    /**
     * set status of faces to false since release_face_status() doesn't work.\n
     */
    void setFaceStatusToFalse();

    /**
     * round indices based on this vector.\n
     * there are boundary-, feature edge- and singularity constraints which need to be rounded.\n
     * additionally the j and k values need to be rounded aswell.\n
     * @param nbVerticesUaV number of vertices in global param system
     * @param jkValues number of j and k values
     * @param singularities vector of vertices
     * @param onlyBoundaries vector of halfedges
     * @return vector of indices which are to be rounded.
     */
    std::vector<int>
    getIdxToRound(int nbVerticesUaV, int jkValues, std::vector<int> &singularities, std::vector<int> &onlyBoundaries);

    /**
     * extract params from solution to assign vertices with coordinates of disk topology.\n
     * calls initPropForSolVector, getSolFromVerticesWMoreOneApp and saveSolToVertices.\n
     * @param _x solution vector
     * @param faces vector of faces in selection
     * @param cutGraphWoBoundary vector of halfedges
     * @param nbVerticesUaV start of jkValues
     */
    void
    saveSolAsCoord(std::vector<double> &_x, std::vector<int> &faces);

    /**
     * set up status for use in getSolFromVerticesWMoreOneApp and getSolFromVerticesWOneApp.\n
     */
    void initPropForSolVector();

    /**
     * gets Coordinate for vertex with one apperance on disk topology.\n
     * @param fh halfedge handle
     * @param _x solution vector
     */
    void saveSolToVertices(OpenMesh::SmartHalfedgeHandle he, std::vector<double> &_x);

    void colorCompBoundaries(std::vector<int> &onlyBoundaries);

    void
    setFeatureLineConstraint(gmm::row_matrix<gmm::wsvector<double>> &_constraints, std::vector<int> &onlyBoundaries,
                             const int startingPoint);

    void propagation(const int &singularity, OpenMesh::SmartHalfedgeHandle he);
};


#endif //OPENFLIPPER_GLOBALPARAMETRIZATION_H
