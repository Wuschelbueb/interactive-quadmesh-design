//
// Created by wuschelbueb on 05.05.22.
//

#include "GlobalParametrization.hh"

void GlobalParametrization::getGlobalParam() {
    std::vector<int> faces;
    getFaceVec(faces);
    // *2 because we have a V and a U param
    int nbVerticesUaV = createVertexPosParamDom(faces) * 2;
    setUpLocFaceCoordSys(faces);
    DijkstraDistance dualGraph(trimesh_);
    dualGraph.getDualGraph(faces);
    std::vector<double> _x(nbVerticesUaV);
    std::vector<double> _rhs = getRhs(faces, nbVerticesUaV);
    CMatrixType _H = getHessian(faces, nbVerticesUaV);
    std::vector<int> complementHEdges;
    getComplementMeshSel(complementHEdges);
    removeOpenPaths(complementHEdges);
    colorCompEdges(complementHEdges);
    connectSingularityToCutGraph(complementHEdges, dualGraph);
    std::cout << "ehllo world\n";
}

void
GlobalParametrization::connectSingularityToCutGraph(std::vector<int> complementHEdges, DijkstraDistance &dualGraph) {
    auto crossFieldIdx = OpenMesh::VProp<double>(trimesh_, "crossFieldIdx");
    std::vector<int> singularities;
    for (auto vh: trimesh_.vertices()) {
        if (crossFieldIdx[vh] < -1E-1) {
            singularities.push_back(vh.idx());
        } else if (crossFieldIdx[vh] > 1E-1) {
            singularities.push_back(vh.idx());
        }
    }
    dualGraph.getDijkstraSingularities(complementHEdges, singularities);
}

void GlobalParametrization::removeOpenPaths(std::vector<int> &complementHEdges) {
    trimesh_.release_halfedge_status();
    trimesh_.request_halfedge_status();
    int currentIter = complementHEdges.size(), prevIter = INT_MAX;
    while (currentIter != prevIter) {
        for (int i: complementHEdges) {
            removeEdgeFromGraph(i, complementHEdges);
        }
        prevIter = currentIter;
        currentIter = complementHEdges.size();
    }
}

void GlobalParametrization::removeEdgeFromGraph(const int i, std::vector<int> &complementHEdges) {
    auto crossFieldIdx = OpenMesh::VProp<double>(trimesh_, "crossFieldIdx");
    OpenMesh::HalfedgeHandle he = trimesh_.halfedge_handle(i);
    OpenMesh::HalfedgeHandle ohe = trimesh_.opposite_halfedge_handle(he);
    OpenMesh::VertexHandle vh = trimesh_.to_vertex_handle(he);
    int counter = 0;
    for (TriMesh::VertexOHalfedgeIter voh_it = trimesh_.voh_iter(vh); voh_it.is_valid(); ++voh_it) {
        if ((std::find(complementHEdges.begin(), complementHEdges.end(), voh_it->idx()) != complementHEdges.end())) {
            counter++;
        }
    }
    //valence check
    if (counter == 1) {
        complementHEdges.erase(std::remove(complementHEdges.begin(), complementHEdges.end(), he.idx()),
                               complementHEdges.end());
        complementHEdges.erase(std::remove(complementHEdges.begin(), complementHEdges.end(), ohe.idx()),
                               complementHEdges.end());
    }
}

void GlobalParametrization::colorCompEdges(const std::vector<int> &complementEdges) {
    auto colorBlub = OpenMesh::EProp<int>(trimesh_, "colorBlub");
    for (auto eh: trimesh_.edges()) {
        colorBlub[eh] = 0;
    }
    for (int i: complementEdges) {
        OpenMesh::HalfedgeHandle he = trimesh_.halfedge_handle(i);
        OpenMesh::EdgeHandle eh = trimesh_.edge_handle(he);
        colorBlub[eh] = 1;
    }
}

void GlobalParametrization::getComplementMeshSel(std::vector<int> &complementHEdges) {
    trimesh_.release_halfedge_status();
    trimesh_.request_halfedge_status();
    tagEdgesFromDualSpanningTree();
    for (auto he: trimesh_.halfedges()) {
        if (!trimesh_.status(he).tagged()) {
            complementHEdges.push_back(he.idx());
        }
    }
}

void GlobalParametrization::tagEdgesFromDualSpanningTree() {
    auto dualGraphPred = OpenMesh::FProp<int>(trimesh_, "dualGraphPred");
    auto FaceSel = OpenMesh::FProp<bool>(trimesh_, "FaceSel");
    for (auto face: trimesh_.faces()) {
        if (dualGraphPred[face] != INT_MAX && dualGraphPred[face] != face.idx() && FaceSel[face] == true) {
            OpenMesh::FaceHandle fh_pre = trimesh_.face_handle(dualGraphPred[face]);
            for (TriMesh::FaceHalfedgeIter fhe_it = trimesh_.fh_iter(face); fhe_it.is_valid(); ++fhe_it) {
                OpenMesh::HalfedgeHandle oheh = trimesh_.opposite_halfedge_handle(*fhe_it);
                for (TriMesh::FaceHalfedgeIter fhe_pre_it = trimesh_.fh_iter(
                        fh_pre); fhe_pre_it.is_valid(); ++fhe_pre_it) {
                    if (fhe_pre_it->idx() == oheh.idx()) {
                        trimesh_.status(oheh).set_tagged(true);
                        trimesh_.status(*fhe_it).set_tagged(true);
                    }
                }
            }
        }
    }
}

void GlobalParametrization::getFaceVec(std::vector<int> &faces) {
    auto FaceSel = OpenMesh::FProp<bool>(trimesh_, "FaceSel");
    for (auto face: trimesh_.faces()) {
        if (FaceSel[face] == true) {
            faces.push_back(face.idx());
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

std::vector<double> GlobalParametrization::getRhs(const std::vector<int> &faces, const int nbVerticesUaV) {
    auto uVectorFieldRotOne = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotOne");
    auto uVectorFieldRotTwo = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotTwo");
    std::vector<double> _rhs(nbVerticesUaV);

    // add vertex numbering in order to know where to start
    // use nbVerticesUaV and lhs formulas with dkm -> dkm is entry in D matrix
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
    auto basisTransformationMtx = OpenMesh::FProp<gmm::row_matrix<std::vector<double>>>(trimesh_,
                                                                                        "basisTransformationMtx");
    gmm::row_matrix<std::vector<double> > _D = basisTransformationMtx[fh];
    double weight = 1, area = trimesh_.calc_face_area(fh), h = 1, col = 0;
    for (TriMesh::FaceVertexIter fv_it = trimesh_.fv_iter(fh); fv_it.is_valid(); ++fv_it) {
        //map U coordinates to vertices
        col = mapLocCoordToGlobCoordSys(fh, *fv_it);
        //get sum of ut*dkm
        double sum = 0;
        for (int j = 0; j < 3; ++j) {
            sum += CrossFieldAxis[j] * basisTransformationMtx[fh][j][col];
        }
        //check if U or V position is needed
        int pos = (flagUorV) ? vertexPosUi[*fv_it] : vertexPosVi[*fv_it];
        _rhs[pos] += 2 * weight * area * h * sum;
    }
}

GlobalParametrization::CMatrixType
GlobalParametrization::getHessian(const std::vector<int> &faces, const int nbVerticesUaV) {
    CMatrixType _H(nbVerticesUaV, nbVerticesUaV);
    for (int i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        getDiaEntriesHessian(fh, _H);
//        getEntriesHessian(fh, _H);
    }
    return _H;
}

//map U coordinates to vertices
int GlobalParametrization::mapLocCoordToGlobCoordSys(const OpenMesh::FaceHandle fh, const OpenMesh::VertexHandle vh) {
    auto localVUi = OpenMesh::FProp<std::vector<int>>(trimesh_, "localVUi");
    double column = 0;
    std::vector<int> vertices = localVUi[fh];
    if (vh.idx() == vertices[0]) {
        column = 0;
    } else if (vh.idx() == vertices[1]) {
        column = 1;
    } else if (vh.idx() == vertices[2]) {
        column = 2;
    }
    return column;
}

void GlobalParametrization::getDiaEntriesHessian(const OpenMesh::FaceHandle fh, CMatrixType &_H) {
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    auto basisTransformationMtx = OpenMesh::FProp<gmm::row_matrix<std::vector<double>>>(trimesh_,
                                                                                        "basisTransformationMtx");
    double weight = 1, area = trimesh_.calc_face_area(fh), h = 1, col;
    for (TriMesh::FaceVertexIter fv_it = trimesh_.fv_iter(fh); fv_it.is_valid(); ++fv_it) {
        col = mapLocCoordToGlobCoordSys(fh, *fv_it);
        double sum = 0;
        for (int j = 0; j < 3; ++j) {
            sum += basisTransformationMtx[fh][j][col] * basisTransformationMtx[fh][j][col];
        }
        _H(vertexPosUi[*fv_it], vertexPosUi[*fv_it]) += 2 * weight * area * pow(h, 2) * sum;
        _H(vertexPosVi[*fv_it], vertexPosVi[*fv_it]) += 2 * weight * area * pow(h, 2) * sum;
    }
}

void GlobalParametrization::getEntriesHessian(const OpenMesh::FaceHandle fh, CMatrixType &_H) {
    //doesn't work yet
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    auto basisTransformationMtx = OpenMesh::FProp<gmm::row_matrix<std::vector<double>>>(trimesh_,
                                                                                        "basisTransformationMtx");
    double weight = 1, area, h = 1, col1, col2;
    for (TriMesh::FaceVertexIter vh_i = trimesh_.fv_iter(fh); vh_i.is_valid(); ++vh_i) {
        for (TriMesh::VertexOHalfedgeIter voh_it = trimesh_.voh_iter(*vh_i); voh_it.is_valid(); ++voh_it) {
            OpenMesh::FaceHandle fh2 = trimesh_.face_handle(*voh_it);
            area = trimesh_.calc_face_area(fh2);
            OpenMesh::VertexHandle vh_j = trimesh_.to_vertex_handle(*voh_it);
            col1 = mapLocCoordToGlobCoordSys(fh2, *vh_i);
            col2 = mapLocCoordToGlobCoordSys(fh2, vh_j);
            //get sum of ut*dkm
            double sum = 0;
            for (int j = 0; j < 3; ++j) {
                sum += basisTransformationMtx[fh2][j][col1] * basisTransformationMtx[fh2][j][col2];
            }
            _H(vertexPosUi[*vh_i], vertexPosUi[vh_j]) += 2 * weight * area * pow(h, 2) * sum;
            _H(vertexPosUi[vh_j], vertexPosUi[*vh_i]) += 2 * weight * area * pow(h, 2) * sum;
            _H(vertexPosVi[*vh_i], vertexPosVi[vh_j]) += 2 * weight * area * pow(h, 2) * sum;
            _H(vertexPosVi[vh_j], vertexPosVi[*vh_i]) += 2 * weight * area * pow(h, 2) * sum;
        }
    }
}