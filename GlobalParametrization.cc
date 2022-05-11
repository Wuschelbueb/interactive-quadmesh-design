//
// Created by wuschelbueb on 05.05.22.
//

#include "GlobalParametrization.hh"

void GlobalParametrization::getGlobalParam() {
    std::vector<int> faces;
    getFaceVec(faces);
    // *2 because we have a V and a U param
    int rhsSize = createVertexPosParamDom(faces) * 2;
    setUpLocFaceCoordSys(faces);
    std::vector<double> _x(rhsSize);
    std::vector<double> _rhs = getRhs(faces, rhsSize);
    CMatrixType _H = getHessian(faces, rhsSize);
    std::cout << "ehllo world\n";
}

void GlobalParametrization::getFaceVec(std::vector<int> &faces) {
    auto FaceSel = OpenMesh::FProp<bool>(trimesh_, "FaceSel");
    for (TriMesh::FaceIter f_it = trimesh_.faces_begin(); f_it != trimesh_.faces_end(); ++f_it) {
        if (FaceSel[*f_it] = true) {
            faces.push_back(f_it->idx());
        }
    }
}

int GlobalParametrization::createVertexPosParamDom(std::vector<int> &faces) {
    trimesh_.release_vertex_status();
    trimesh_.request_vertex_status();
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    int countVertices = 0;
    for (int i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        for (TriMesh::FaceVertexIter fv_it = trimesh_.fv_iter(fh); fv_it.is_valid(); ++fv_it) {
            if (!trimesh_.status(*fv_it).tagged()) {
                vertexPosUi[*fv_it] = countVertices++;
                vertexPosVi[*fv_it] = countVertices++;
                trimesh_.status(*fv_it).set_tagged(true);
            }
        }
    }
    return countVertices;
}

void GlobalParametrization::setUpLocFaceCoordSys(const std::vector<int> &faces) {
    gmm::col_matrix<std::vector<double> > _C = createCMatrix();
    std::vector<Point> edges(3);
    for (int i: faces) {
        int counter = 0;
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        createEdgesAndLocalVUi(fh, counter, edges);
        createBasisTfMtx(fh, _C, edges);
    }
}

void
GlobalParametrization::createEdgesAndLocalVUi(const OpenMesh::FaceHandle fh, int &counter, std::vector<Point> &edges) {
    auto localVUi = OpenMesh::FProp<std::vector<int>>(trimesh_, "localVUi");
    for (TriMesh::FaceHalfedgeIter fh_it = trimesh_.fh_iter(fh); fh_it.is_valid(); ++fh_it) {
        if (counter == 0) {
            Point temp = trimesh_.calc_edge_vector(*fh_it).normalize();
            localVUi[fh].push_back(trimesh_.to_vertex_handle(*fh_it).idx());
            localVUi[fh].push_back(trimesh_.from_vertex_handle(*fh_it).idx());
            // minus because common vertex is vec[0]
            edges[0] = -temp;
        } else if (counter == 1) {
            Point temp = trimesh_.calc_edge_vector(*fh_it).normalize();
            localVUi[fh].push_back(trimesh_.to_vertex_handle(*fh_it).idx());
            edges[1] = temp;
        }
        counter++;
    }
    // calculate normal of triangle
    edges[2] = (edges[0] % edges[1]).normalize();
}

void GlobalParametrization::createBasisTfMtx(const OpenMesh::FaceHandle fh, gmm::col_matrix<std::vector<double> > &_C,
                                             const std::vector<Point> &edges) {
    auto basisTransformationMtx = OpenMesh::FProp<gmm::row_matrix<std::vector<double>>>(trimesh_,
                                                                                        "basisTransformationMtx");
    gmm::row_matrix<std::vector<double> > _D(3, 3);
    gmm::row_matrix<std::vector<double> > _basis(3, 3);
    // build basis vector out of edges
    _basis(0, 0) = edges[0][0];
    _basis(1, 0) = edges[0][1];
    _basis(2, 0) = edges[0][2];
    _basis(0, 1) = edges[1][0];
    _basis(1, 1) = edges[1][1];
    _basis(2, 1) = edges[1][2];
    _basis(0, 2) = edges[2][0];
    _basis(1, 2) = edges[2][1];
    _basis(2, 2) = edges[2][2];
    gmm::lu_inverse(_basis);
    gmm::mult(_basis, _C, _D);
    basisTransformationMtx[fh] = _D;
}

gmm::col_matrix<std::vector<double>> GlobalParametrization::createCMatrix() {
    gmm::col_matrix<std::vector<double> > _C(3, 3);
    _C(0, 0) = -1;
    _C(1, 0) = 1;
    _C(2, 0) = 0;
    _C(0, 1) = -1;
    _C(1, 1) = 0;
    _C(2, 1) = 1;
    _C(0, 2) = 0;
    _C(1, 2) = 0;
    _C(2, 2) = 0;
    return _C;
}

std::vector<double> GlobalParametrization::getRhs(const std::vector<int> &faces, const int rhsSize) {
    auto uVectorFieldRotOne = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotOne");
    auto uVectorFieldRotTwo = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotTwo");
    std::vector<double> _rhs(rhsSize);

    // add vertex numbering in order to know where to start
    // use rhsSize and lhs formulas with dkm -> dkm is entry in D matrix
    // utk is entry of vectorRotOne
    for (int i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        Point uCrossField = uVectorFieldRotOne[fh];
        Point vCrossField = uVectorFieldRotTwo[fh];
        getRhsEntryForVertex(fh, uCrossField, true, _rhs);
        getRhsEntryForVertex(fh, vCrossField, false, _rhs);
    }
    return _rhs;
}

void
GlobalParametrization::getRhsEntryForVertex(const OpenMesh::FaceHandle fh, const Point CrossFieldAxis,
                                            const bool flagUorV,
                                            std::vector<double> &_rhs) {
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    auto localVUi = OpenMesh::FProp<std::vector<int>>(trimesh_, "localVUi");
    auto basisTransformationMtx = OpenMesh::FProp<gmm::row_matrix<std::vector<double>>>(trimesh_,
                                                                                        "basisTransformationMtx");
    std::vector<int> vertices = localVUi[fh];
    gmm::row_matrix<std::vector<double> > _D = basisTransformationMtx[fh];
    double weight = 1, area = 0, h = 1, col = 0;
    for (TriMesh::FaceVertexIter fv_it = trimesh_.fv_iter(fh); fv_it.is_valid(); ++fv_it) {
        //map U coordinates to vertices
        if (fv_it->idx() == vertices[0]) {
            col = 0;
        } else if (fv_it->idx() == vertices[1]) {
            col = 1;
        } else if (fv_it->idx() == vertices[2]) {
            col = 2;
        }
        //get sum of ut*dkm
        double sum = 0;
        for (int j = 0; j < 3; ++j) {
            sum += CrossFieldAxis[j] * basisTransformationMtx[fh][j][col];
        }

        //check if U or V position is needed
        int pos;
        if (flagUorV) {
            pos = vertexPosUi[*fv_it];
        } else {
            pos = vertexPosVi[*fv_it];
        }

        _rhs[pos] += 2 * weight * area * h * sum;
    }
}

GlobalParametrization::CMatrixType GlobalParametrization::getHessian(const std::vector<int> &faces, const int rhsSize) {

}