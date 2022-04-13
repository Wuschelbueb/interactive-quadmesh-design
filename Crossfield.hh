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
    Crossfield(TriMesh &trimesh, std::vector<int> &heInRange, std::vector<int> &heConstraints)
            : trimesh_{trimesh}, heInRange_{heInRange}, heConstraints_{heConstraints} {
    }

    ~Crossfield() {
    }

    void getCrossfield();

    void getEnergy();


    /// test functions using given linear system from quadratic problem A, x, b: x^T * A * x + b^T x
    double computeEnergy(const CMatrixType &A, const CVectorType &x, const CVectorType b) const {
        int m = gmm::mat_nrows(A);
        int n = gmm::mat_ncols(A);
        assert(m == n);

        std::vector<double> y(n);
        gmm::mult(A, x, y);    // y = A*x
        double p = gmm::vect_sp(x, y); // x^T * y = x^T * A * x
        double q = gmm::vect_sp(x, b); // scalar product

        return (p + q);
    }

    /// check that angle from first to second equals negative of second to first
    inline bool testKappa(const int refEdgeMain, const int refEdgeNeigh, const std::pair<int, int> commonEdge) {
        const double tol = 1e-7;
        double kappa01 = getKappa(refEdgeMain, refEdgeNeigh, commonEdge);
        double kappa10 = getKappa(refEdgeNeigh, refEdgeMain, {commonEdge.second, commonEdge.first});
        double adiff = std::abs(kappa01 + kappa10);
        bool ok = adiff < tol;
        if (!ok) {
            std::cerr << "kappa is not skew symmetric: |" << kappa01 << " + " << kappa10 << "| = " << adiff << "!= 0"
                      << std::endl;
        }
        return ok;
    }

    /// check constraints are satisfied
    inline bool check_constraints(const RMatrixType &C, const CVectorType &x, double &max_violation) const {
        const double tol = 1e-6;
        int m = gmm::mat_nrows(C);
        max_violation = 0.;
        for (int r = 0; r < m; ++r) {
            gmm::wsvector<double> number = C.row(r);
            number.size();
            x.size();
            double val = gmm::vect_sp(C.row(r), x);
            assert(!std::isnan(val));
            double nval = std::abs(val);
            if (nval > max_violation) {
                max_violation = val;
            }
        }
        return (max_violation < tol);
    }

private:

    void createCrossfields(const std::vector<int> &faces);

    double
    getEnergy(const CMatrixType &A, const std::vector<double> &x, const std::vector<double> &b);

    gmm::row_matrix<gmm::wsvector<double>>
    getConstraintMatrix(const std::map<int, double> &heKappa, const std::vector<int> &faces);

    void getThetaConstraints(const int n_col, int &counter, std::vector<int> &faceConstraints,
                             gmm::row_matrix<gmm::wsvector<double>> &_constraints);

    void getPJConstraints(const std::map<int, double> &heKappa, int &counter, const int pj_start,
                          const int &noOriginConst,
                          gmm::row_matrix<gmm::wsvector<double>> &_constraints);

    int getAmountPJConstraints(const std::vector<int> &faces);

    std::vector<int> getIdxToRound(const std::map<int, double> &heKappa, const std::vector<int> &faces);

    std::vector<double> getRHS(const std::map<int, double> &heKappa, const std::vector<int> &faces);

    void getRhsFirstHalf(const std::vector<int> &faces, std::vector<double> &_rhs,
                         const std::map<int, double> &heKappa);

    void getSum(const int i, std::vector<double> &_rhs, const std::map<int, double> &heKappa);

    double setSum(OpenMesh::FaceHandle fh, OpenMesh::HalfedgeHandle fh_it, const std::map<int, double> &heKappa);

    void getRhsSecondHalf(std::vector<double> &_rhs,
                          const std::map<int, double> &heKappa, const int facesPlusOne);

    RMatrixType getMatrixA(const std::vector<int> &faces, const std::map<int, double> &heKappa);

    CMatrixType getHessianMatrix(const std::vector<int> &faces, const std::map<int, double> &heKappa);

    int getFactor(const OpenMesh::FaceHandle fh);

    void setPositionInHessianForFaces(const std::map<int, double> &heKappa);

    std::map<int, double> getMapHeKappa(const std::vector<int> &faces);

    void getStatusNeigh(const OpenMesh::FaceHandle fh, const OpenMesh::FaceHandle fh_neigh,
                        std::map<int, double> &heKappa);

    void addKappaHeToMap(const std::pair<int, int> commonEdge, const double kappa, std::map<int, double> &heKappa);

    std::pair<int, int>
    getCommonEdgeBetweenTriangles(const OpenMesh::FaceHandle fh, const OpenMesh::FaceHandle fh_neigh);

    double getKappa(const int refEdgeMain, const int refEdgeNeigh, const std::pair<int, int> commonEdge);

    void setVecFieldProp();

    void getConstraintAngleAndVecField(const std::vector<int> &faces);

    void setRotThetaOfVectorField(const std::vector<int> &faces, const std::vector<double> _x);

    gmm::dense_matrix<double> getRotMatrix(const double theta);

    void colorFaces(const std::vector<int> &faces);

    void colorHEdges();

    std::vector<int> getFacesVecWithRefHeProp();

    void setRefHeToFace(const int i, std::vector<int> &faces);

    double shortenKappa(const double kappa);

    Point multPointWithRotMatrix(const Point rotVec, const double angle = 1.0);

    std::vector<int> getFaceConstraints();


    TriMesh &trimesh_;
    std::vector<int> &heInRange_;
    std::vector<int> &heConstraints_;
};


#endif //OPENFLIPPER_CROSSFIELD_HH
