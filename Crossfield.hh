//
// Created by wuschelbueb on 15.09.21.
//

#ifndef OPENFLIPPER_CROSSFIELD_HH
#define OPENFLIPPER_CROSSFIELD_HH

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

/**
 * class is used to create crossfield of each triangle included in selection
 */
class Crossfield {
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
    Crossfield(TriMesh &trimesh, std::vector<int> &heInRange, std::vector<int> &heConstraints, ACG::Vec3d refVector = {0,0,0})
            : trimesh_{trimesh}, heInRange_{heInRange}, heConstraints_{heConstraints}, refVector_{refVector} {
    }

    ~Crossfield() {
    }

    /**
     * get the crossfield
     */
    void getCrossfield();

private:

    /**
     * create crossfield with solution from CoMiSo solver.\n
     * for each face which is in selection there is a crossfield.\n
     * @param faces ector of faces of selection
     */
    void createCrossfields(const std::vector<int> &faces);

    /**
     * gets and creates constraint matrix
     * @param heKappa map with halfedges and their assigned kappa
     * @param faces vector of faces of selection
     * @return constraint matrix
     */
    gmm::row_matrix<gmm::wsvector<double>>
    getConstraintMatrix(const std::map<int, double> &heKappa, const std::vector<int> &faces);

    /**
     * create the theta constraints.\n
     * assign in n_col+1 the theta angle and at the position of the according face a one.\n
     * increase counter for next iteration
     * @param n_col size of columns
     * @param counter decides which row is worked on
     * @param faceConstraints vector of face constraints
     * @param _constraints row matrix
     */
    void getThetaConstraints(const int n_col, int &counter, std::vector<int> &faceConstraints,
                             gmm::row_matrix<gmm::wsvector<double>> &_constraints);

    /**
     * create periodjump constraints.\n
     * if two faces have the same origin and one is the predecessor of the other add an entry in the constraint matrix.\n
     * @param heKappa map with halfedges and their assigned kappa
     * @param counter decides which row is worked on
     * @param pj_start at which position (column) periodjump start
     * @param _constraints row matrix
     */
    void getPJConstraints(const std::map<int, double> &heKappa, int &counter, const int pj_start,
                          gmm::row_matrix<gmm::wsvector<double>> &_constraints);

    /**
     * gets the number of periodjump constraints a mesh has.\n
     * i.e. the number of faces which aren't their own origin.\n
     * @param faces vector of faces of selection
     * @return the number of periodjump constraints
     */
    int getAmountPJConstraints(const std::vector<int> &faces, const std::map<int, double> &heKappa);

    /**
     * gets the indices which need to be rounded which are all periodjump entries.\n     *
     * @param heKappa map of halfedges with their assigned kappa
     * @param faces vector of faces of selection
     * @return amount of indicies to be rounded
     */
    std::vector<int> getIdxToRound(const std::map<int, double> &heKappa, int pj_start);

    /**
     * creates righthandside (b). Ax = b.\n
     * needed for the constrained solver. is split in two parts.\n
     * @param heKappa map of halfedges with their assigned kappa
     * @param faces vector of faces of selection
     * @return vector of righthandside
     */
    std::vector<double> getRHS(const std::map<int, double> &heKappa, const std::vector<int> &faces);

    /**
     * first part of creating rhs vector.\n
     * handles all the angles (thetas) for the righthandside.\n
     * each face has a certain amount of kappas attached to it. add them up.\n
     * @param faces vector of faces of selection
     * @param _rhs vector of righthandside
     * @param heKappa map of halfedges with their assigned kappa
     */
    void getRhsFirstHalf(const std::vector<int> &faces, std::vector<double> &_rhs,
                         const std::map<int, double> &heKappa);

    /**
     * helper function of the getRhsFirstHalf.\n
     * iterates through each halfedge of triangle. and adds up the attached kappas\n
     * @param i index of a face
     * @param _rhs vector of righthandside
     * @param heKappa map of halfedges with their assigned kappa
     */
    void getSum(const int i, std::vector<double> &_rhs, const std::map<int, double> &heKappa);

    /**
     * extracts the kappa of a halfedge
     * @param fh_it halfedge
     * @param heKappa map of halfedges with their assigned kappa
     * @return kappa
     */
    double extractKappa(OpenMesh::HalfedgeHandle fh_it, const std::map<int, double> &heKappa);

    /**
     * adds kappa of halfedge to the according position of rhs.\n
     * @param _rhs vector of righthandside
     * @param heKappa map of halfedges with their assigned kappa
     * @param facesPlusOne position where on rhs to start
     */
    void getRhsSecondHalf(std::vector<double> &_rhs,
                          const std::map<int, double> &heKappa, const int facesPlusOne);

    /**
     * get hessian matrix for the CoMiSo solver.\n
     * there are multiple helper function to create the hessian matrix.\n
     * @param faces vector of faces of selection
     * @param heKappa map of halfedges with their assigned kappa
     * @return column matrix
     */
    CMatrixType getHessianMatrix(const std::vector<int> &faces, const std::map<int, double> &heKappa);

    /**
     * multiplies entry of face fh in Hessian matrix by the amount of adjacent faces in the faces vector.\n
     * this is only relevant for the diagonal of the hessian matrix.\n
     * helper function of hessian matrix.\n
     * @param fh face handle
     * @param faces vector of faces of selection
     * @return factor which gets multiplied
     */
    int getFactor(const OpenMesh::FaceHandle fh, const std::vector<int> &faces);

    /**
     * each face gets a unique number (property) to make sure that the position in the hessian is unique aswell.\n
     * sometimes not all faces are in selection, that is why we need a numbering.\n
     * helper function of hessian matrix.\n
     * @param heKappa map of halfedges with their assigned kappa
     */
    void setPositionInHessianForFaces(const std::map<int, double> &heKappa);

    /**
     * creates a map of halfedges and their assigned kappas.\n
     * @param faces vector containing faces of selection
     * @return map which contains halfedges and their assigned kappa
     */
    std::map<int, double> getMapHeKappa(const std::vector<int> &faces);

    /**
     * checks if the neighbouring face is aswell part of the face selection.\n
     * helper function of getMapHeKappa.\n
     * @param fh facehandle
     * @param fh_neigh neighbouring face handle
     * @param heKappa map of halfedges with their assigned kappa
     */
    void getStatusNeigh(const OpenMesh::FaceHandle fh, const OpenMesh::FaceHandle fh_neigh,
                        std::map<int, double> &heKappa);

    /**
     * adds kappa from getKappa function to the heKappa map if the opposing halfedge wasn't already added.\n
     * helper function of getMapHeKappa.\n
     * @param commonEdge pair of two opposing halfedges
     * @param kappa angle kappa
     * @param heKappa map of halfedges with the assigned kappa
     */
    void addKappaHeToMap(const std::pair<int, int> commonEdge, const double kappa, std::map<int, double> &heKappa);

    /**
     * returns common edge, i.e. the two opposing halfedges between two neighbouring triangles.\n
     * helper function of getStatusNeigh.\n
     * @param fh face handle
     * @param fh_neigh face handle neighbour
     * @return pair of two opposing halfedges
     */
    std::pair<int, int>
    getCommonEdgeBetweenTriangles(const OpenMesh::FaceHandle fh, const OpenMesh::FaceHandle fh_neigh);

    /**
     * calculates the kappa between two reference edges, i.e. two neighbouring faces.\n
     * helper function of getMapHeKappa.\n
     * @param refEdgeMain reference Edge of face 1
     * @param refEdgeNeigh reference Edge of face 2
     * @param commonEdge opposing halfedges which connect those two faces
     * @return angle kappa
     */
    double getKappa(const int refEdgeMain, const int refEdgeNeigh, const std::pair<int, int> commonEdge);

    /**
     * initialize point (uVectorFieldRot1,2,3,4) properties with {0,0,0}
     */
    void setVecFieldProp();

    /**
     * calculates the constraint angle with atan2 and sets x-Axis of local Vectorfield.
     * @param faces vector containing faces of selection
     */
    void getConstraintAngleAndVecField(const std::vector<int> &faces);

    /**
     * gets crossfield Index for each vertex. used to find singularities.
     * has helper functions
     * @param faces vector containing faces of selection
     * @param heKappa map of halfedges with their assigned kappa
     * @param _x solution vector from solver
     */
    void getCrossFieldIdx(const std::vector<int> &faces, const std::map<int, double> &heKappa,
                          const std::vector<double> &_x);

    /**
     * rotate crossfield of each triangle with theta from solution vector _x\n
     * save as rotated as crossfield property
     * @param faces included in selection
     * @param _x solution vector from CoMiSo constrained solver
     */
    void setRotThetaOfVectorField(const std::vector<int> &faces, const std::vector<double> _x);

    /**
     * color faces depending if they are in selection or not
     * @param faces vector containing faces in selection
     */
    void colorFaces(const std::vector<int> &faces);

    /**
     * color halfedges depending if they are in selection, constraints or none of these two
     */
    void colorHEdges();

    /**
     * get the vector containing faces from selection. additionally set the reference He for each triangle in selection.\n
     * reference halfedge is a property stored in a face.
     * @return vector containing faces of selection
     */
    std::vector<int> getFacesVecWithRefHeProp();

    /**
     * add faces from halfedges to faces vector and use these halfedges as reference halfedge.\n
     * check if halfedge is boundary, if yes break.
     * @param i halfedge index
     * @param faces vector containing faces of selection
     */
    void setRefHeToFace(const int i, std::vector<int> &faces);

    /**
     * basically like modulo, but with 2*pi.\n
     * return angle is always between (-pi, pi]
     * @param kappa angle
     * @return shorted angle
     */
    double shortenKappa(const double kappa);

    /**
     * rotate vec with Rodrigues rotation formula.\n
     * is used when a vector needs to be rotated a certain angle in space, given an axis and angle of rotation.\n
     * @param fh face handle
     * @param vec Vec3d
     * @param theta angle
     * @return Vec3d
     */
    Point rotPointWithRotMatrix(const OpenMesh::FaceHandle fh, const Point vec, const double theta = 1.0);

    /**
     * get one constrained per face even if there are multiple constrained halfedges.\n
     * @return vector with constrained halfedges
     */
    std::vector<int> getConstrainedHe();

    /**
     * set periodjump property for each halfedge. extract periodjump from solution, _x vector.\n
     * @param heKappa map of halfedges with their assigned kappa
     * @param _x vector from solver containing solution
     * @param faceSize position where periodjump start
     */
    void setPJProp(std::map<int, double> &heKappa, std::vector<double> &_x, const int faceSize);

    TriMesh &trimesh_;
    /**
     * halfedges which were included of dijkstra calc.
     */
    std::vector<int> &heInRange_;
    /**
     * halfedges which were selected by selection tool.
     */
    std::vector<int> &heConstraints_;
    ACG::Vec3d refVector_;


    void getCrossFldIdxConstraints(const std::map<int, double> &heKappa, int &counter, const int pj_start,
                                   gmm::row_matrix<gmm::wsvector<double>> &_constraints,
                                   const std::vector<int> &faces);
};


#endif //OPENFLIPPER_CROSSFIELD_HH
