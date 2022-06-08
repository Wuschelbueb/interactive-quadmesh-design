//
// Created by wuschelbueb on 05.05.22.
//

#include "GlobalParametrization.hh"

void GlobalParametrization::getGlobalParam() {
    //set up vectors
    std::vector<int> faces = getFaceVec();
    std::vector<int> singularities = getSingularities();
    setUpLocFaceCoordSys(faces);
    DijkstraDistance dualGraph(trimesh_);
    dualGraph.getDualGraph(faces);

    std::vector<int> complementHEdges = getComplementMeshSel();
    removeOpenPaths(complementHEdges);
    removeRedundantEdges(complementHEdges);
    //just for visualization
    colorCompHEdges(complementHEdges);
    std::vector<int> onlyBoundaries = complementHEdges;
    std::vector<int> cutGraphWoBoundary;
    //add path from boundary to singularity to cut graph
    dualGraph.completeDijkstraWSingularities(complementHEdges, singularities, cutGraphWoBoundary);
    createSectorsCutGraph(singularities);
    fixRotationsCrossBoundaryComp(complementHEdges, singularities, faces);

    int nbVerticesUaV = createVertexPosParamDomain(faces);
    std::vector<double> _x(nbVerticesUaV);
    std::vector<double> _rhs = getRhs(faces, nbVerticesUaV);
    CMatrixType _H = getHessian(faces, nbVerticesUaV);
//    RMatrixType _constraints = getConstraints(faces);
    std::cout << "end of code\n";
}

std::vector<int> GlobalParametrization::getSingularities() {
    auto crossFieldIdx = OpenMesh::VProp<double>(trimesh_, "crossFieldIdx");
    std::vector<int> singularities;
    for (auto vh: trimesh_.vertices()) {
        if (crossFieldIdx[vh] < -1E-1) {
            singularities.push_back(vh.idx());
        } else if (crossFieldIdx[vh] > 1E-1) {
            singularities.push_back(vh.idx());
        }
    }
    return singularities;
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

void GlobalParametrization::removeRedundantEdges(std::vector<int> &complementHEdges) {
    auto borderDualG = OpenMesh::EProp<int>(trimesh_, "borderDualG");
    auto faceSel = OpenMesh::FProp<bool>(trimesh_, "faceSel");
    for (auto eh: trimesh_.edges()) {
        auto heh = make_smart(trimesh_.halfedge_handle(eh, 1), trimesh_);
        auto oheh = make_smart(trimesh_.halfedge_handle(eh, 0), trimesh_);
        auto it = std::find(complementHEdges.begin(), complementHEdges.end(), heh.idx());
        auto it2 = std::find(complementHEdges.begin(), complementHEdges.end(), oheh.idx());
        if ((it != complementHEdges.end() || it2 != complementHEdges.end()) && borderDualG[eh] != 1) {
            complementHEdges.erase(it);
            complementHEdges.erase(it2);
        }
        if (heh.face().is_valid() && !faceSel[heh.face()] && it != complementHEdges.end() && borderDualG[eh] == 1) {
            complementHEdges.erase(it);
        }
        if (oheh.face().is_valid() && !faceSel[oheh.face()] && it2 != complementHEdges.end() && borderDualG[eh] == 1) {
            complementHEdges.erase(it2);
        }
    }
}

void GlobalParametrization::colorCompHEdges(const std::vector<int> &complementEdges) {
    auto cutGraphHe = OpenMesh::HProp<bool>(trimesh_, "cutGraphHe");
    for (auto heh: trimesh_.halfedges()) {
        cutGraphHe[heh] = false;
    }
    for (int i: complementEdges) {
        OpenMesh::HalfedgeHandle he = trimesh_.halfedge_handle(i);
        cutGraphHe[he] = true;
    }
}

std::vector<int> GlobalParametrization::getComplementMeshSel() {
    std::vector<int> complementHEdges;
    trimesh_.release_halfedge_status();
    tagEdgesFromDualSpanningTree();
    for (auto heh: trimesh_.halfedges()) {
        if (!trimesh_.status(heh).tagged()) {
            complementHEdges.push_back(heh.idx());
        }
    }
    return complementHEdges;
}

void GlobalParametrization::tagEdgesFromDualSpanningTree() {
    for (auto face: trimesh_.faces()) {
        checkIfFaceInSelection(face);
    }
}

void GlobalParametrization::checkIfFaceInSelection(OpenMesh::FaceHandle &face) {
    auto dualGraphPred = OpenMesh::FProp<int>(trimesh_, "dualGraphPred");
    auto faceSel = OpenMesh::FProp<bool>(trimesh_, "faceSel");
    if (dualGraphPred[face] != INT_MAX && dualGraphPred[face] != face.idx() && faceSel[face] == true) {
        OpenMesh::FaceHandle fh_pre = trimesh_.face_handle(dualGraphPred[face]);
        checkIfEBetweenTriangleInDualGraph(face, fh_pre);

    }
}

void
GlobalParametrization::checkIfEBetweenTriangleInDualGraph(OpenMesh::FaceHandle &face, OpenMesh::FaceHandle &fh_pred) {
    for (TriMesh::FaceHalfedgeIter fhe_it = trimesh_.fh_iter(face); fhe_it.is_valid(); ++fhe_it) {
        OpenMesh::HalfedgeHandle oheh = trimesh_.opposite_halfedge_handle(*fhe_it);
        for (TriMesh::FaceHalfedgeIter fhe_pred_it = trimesh_.fh_iter(
                fh_pred); fhe_pred_it.is_valid(); ++fhe_pred_it) {
            tagEdgeIfInDualGraph(fhe_pred_it, oheh);
        }
    }
}

void
GlobalParametrization::tagEdgeIfInDualGraph(TriMesh::FaceHalfedgeIter &fhe_pred_it, OpenMesh::HalfedgeHandle &oheh) {
    trimesh_.request_halfedge_status();
    auto eh = trimesh_.edge_handle(*fhe_pred_it);
    if ((fhe_pred_it->idx() == oheh.idx())) {
        auto heh = trimesh_.opposite_halfedge_handle(oheh);
        trimesh_.status(oheh).set_tagged(true);
        trimesh_.status(heh).set_tagged(true);
    }
}

std::vector<int> GlobalParametrization::getFaceVec() {
    std::vector<int> faces;
    auto faceSel = OpenMesh::FProp<bool>(trimesh_, "faceSel");
    for (auto face: trimesh_.faces()) {
        if (faceSel[face] == true) {
            faces.push_back(face.idx());
        }
    }
    return faces;
}

int GlobalParametrization::createVertexPosParamDomain(std::vector<int> &faces) {
    trimesh_.release_vertex_status();
    trimesh_.request_vertex_status();
    int countVertices = 0;
    for (int i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        for (TriMesh::FaceVertexIter fv_it = trimesh_.fv_iter(fh); fv_it.is_valid(); ++fv_it) {
            if (!trimesh_.status(*fv_it).tagged()) {
                checkCGandSetPos(*fv_it, countVertices);
            }
        }
    }
    return countVertices;
}

void GlobalParametrization::checkCGandSetPos(OpenMesh::VertexHandle fv_it, int &countVertices) {
    auto vertexDist = OpenMesh::VProp<double>(trimesh_, "vertexDist");
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    if (vertexDist[fv_it] == 0) {
        //check if on cutgraph
        if (checkIfLeaf(fv_it)) { // leaf
            //add increase number normally
            vertexPosUi[fv_it] = countVertices++;
            vertexPosVi[fv_it] = countVertices++;
            trimesh_.status(fv_it).set_tagged(true);
        } else { //inner node
            getPositionInnerNode(fv_it, countVertices);
        }
    } else {
        vertexPosUi[fv_it] = countVertices++;
        vertexPosVi[fv_it] = countVertices++;
        trimesh_.status(fv_it).set_tagged(true);
    }
}

void GlobalParametrization::getPositionInnerNode(OpenMesh::VertexHandle &fv_it, int &countVertices) {
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    std::vector<int> adjacentSectors;
    //if on inner cutgraph add additional space for matrix
    for (TriMesh::VertexFaceIter vf_it = trimesh_.vf_iter(fv_it); vf_it.is_valid(); ++vf_it) {
        auto it = find(adjacentSectors.begin(), adjacentSectors.end(), cutGraphFZone[*vf_it]);
        if (it == adjacentSectors.end() && cutGraphFZone[*vf_it] != 0) {
            adjacentSectors.push_back(cutGraphFZone[*vf_it]);
        }
    }
    //this is in case there are borders which are in the cutgraph
    int appearance = ((int) adjacentSectors.size() == 0) ? 1 : (int) adjacentSectors.size();
    vertexAppearanceCG[fv_it] = appearance;
    vertexPosUi[fv_it] = countVertices;
    countVertices += appearance;
    vertexPosVi[fv_it] = countVertices;
    countVertices += appearance;
    trimesh_.status(fv_it).set_tagged(true);
}

void GlobalParametrization::setUpLocFaceCoordSys(const std::vector<int> &faces) {
    gmm::col_matrix<std::vector<double>> _C = createCMatrix();
    std::vector<Point> edges(3);
    for (int i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        createEdgesAndLocalVUi(fh, edges);
        createBasisTfMtx(fh, _C, edges);
    }
}

void
GlobalParametrization::createEdgesAndLocalVUi(const OpenMesh::FaceHandle fh,
                                              std::vector<Point> &edges) {
    auto localVUi = OpenMesh::FProp<std::vector<int >>(trimesh_, "localVUi");
    int counter = 0;
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

void
GlobalParametrization::createBasisTfMtx(const OpenMesh::FaceHandle fh, gmm::col_matrix<std::vector<double>> &_C,
                                        const std::vector<Point> &edges) {
    auto basisTransformationMtx = OpenMesh::FProp<gmm::row_matrix<std::vector<double>>>(trimesh_,
                                                                                        "basisTransformationMtx");
    gmm::row_matrix<std::vector<double>> _basisInv(3, 3);
    gmm::row_matrix<std::vector<double>> _D(3, 3);
    // build inverse of basis vector out of edges
    _basisInv(0, 0) = edges[1][1] * edges[2][2] - edges[1][2] * edges[2][1];
    _basisInv(1, 0) = edges[1][2] * edges[2][0] - edges[1][0] * edges[2][2];
    _basisInv(2, 0) = edges[1][0] * edges[2][1] - edges[1][1] * edges[2][0];
    _basisInv(0, 1) = edges[0][2] * edges[2][1] - edges[0][1] * edges[2][2];
    _basisInv(1, 1) = edges[0][0] * edges[2][2] - edges[0][2] * edges[2][0];
    _basisInv(2, 1) = edges[0][1] * edges[2][0] - edges[0][0] * edges[2][1];
    _basisInv(0, 2) = edges[0][1] * edges[1][2] - edges[0][2] * edges[1][1];
    _basisInv(1, 2) = edges[0][2] * edges[1][0] - edges[0][0] * edges[1][2];
    _basisInv(2, 2) = edges[0][0] * edges[1][1] - edges[0][1] * edges[1][0];

    double denominator = 1 / (edges[0][0] * edges[1][1] * edges[2][2] - edges[0][0] * edges[1][2] * edges[2][1] -
                              edges[0][1] * edges[1][0] * edges[2][2] + edges[0][1] * edges[1][2] * edges[2][0] +
                              edges[0][2] * edges[1][0] * edges[2][1] - edges[0][2] * edges[1][1] * edges[2][0]);
    gmm::scale(_basisInv, denominator);
    gmm::mult(_basisInv, _C, _D);
    basisTransformationMtx[fh] = _D;
}

gmm::col_matrix<std::vector<double>> GlobalParametrization::createCMatrix() {
    gmm::col_matrix<std::vector<double>> _C(3, 3);
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
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    auto basisTransformationMtx = OpenMesh::FProp<gmm::row_matrix<std::vector<double>>>(trimesh_,
                                                                                        "basisTransformationMtx");
    gmm::row_matrix<std::vector<double>> D = basisTransformationMtx[fh];
    double weight = 1, area = trimesh_.calc_face_area(fh), h = 1, col;
    for (TriMesh::FaceVertexIter fv_it = trimesh_.fv_iter(fh); fv_it.is_valid(); ++fv_it) {
        int appearance = vertexAppearanceCG[*fv_it];
        //map U coordinates to vertices
        col = mapLocCoordToGlobCoordSys(fh, *fv_it);
        //get sum of ut*dkm
        double sum = 0;
        for (int j = 0; j < 3; ++j) {
            sum += CrossFieldAxis[j] * D(j, col);
        }
        //check if U or V position is needed
        int pos = (flagUorV) ? vertexPosUi[*fv_it] : vertexPosVi[*fv_it];
        for (int i = 0; i < appearance; ++i) {
            _rhs[pos + i] += 2 * weight * area * h * sum;
        }
    }
}

//map U coordinates to vertices
int
GlobalParametrization::mapLocCoordToGlobCoordSys(const OpenMesh::FaceHandle fh, const OpenMesh::VertexHandle vh) {
    auto localVUi = OpenMesh::FProp<std::vector<int >>(trimesh_, "localVUi");
    auto faceSel = OpenMesh::FProp<bool>(trimesh_, "faceSel");
    double column;
//    std::vector<int> vertices = localVUi[fh];
    if (vh.idx() == localVUi[fh][0]) {
        column = 0;
    } else if (vh.idx() == localVUi[fh][1]) {
        column = 1;
    } else if (vh.idx() == localVUi[fh][2]) {
        column = 2;
    }
    return column;
}

GlobalParametrization::CMatrixType
GlobalParametrization::getHessian(const std::vector<int> &faces, const int nbVerticesUaV) {
    CMatrixType _H(nbVerticesUaV, nbVerticesUaV);
    for (int i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        getDiaEntriesHessian(fh, _H);
        getEntriesHessian(fh, _H);
    }
    return _H;
}

void GlobalParametrization::getDiaEntriesHessian(const OpenMesh::FaceHandle fh, CMatrixType &_H) {
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    auto basisTransformationMtx = OpenMesh::FProp<gmm::row_matrix<std::vector<double>>>(trimesh_,
                                                                                        "basisTransformationMtx");
    double weight = 1, area = trimesh_.calc_face_area(fh), h = 1;
    for (TriMesh::FaceVertexIter fv_it = trimesh_.fv_iter(fh); fv_it.is_valid(); ++fv_it) {
        int occurence = vertexAppearanceCG[*fv_it];
        gmm::row_matrix<std::vector<double>> D = basisTransformationMtx[fh];
        double sum = 0, col = mapLocCoordToGlobCoordSys(fh, *fv_it);
        for (int j = 0; j < 3; ++j) {
            sum += D(j, col) * D(j, col);
        }
        for (int i = 0; i < occurence; ++i) {
            _H(vertexPosUi[*fv_it] + i, vertexPosUi[*fv_it] + i) += 2 * weight * area * pow(h, 2) * sum;
            _H(vertexPosVi[*fv_it] + i, vertexPosVi[*fv_it] + i) += 2 * weight * area * pow(h, 2) * sum;
        }
    }
}

void GlobalParametrization::getEntriesHessian(const OpenMesh::FaceHandle fh, CMatrixType &_H) {
    //todo: might be different for same vertices but on different sides of the param domain
    //since the boundary cuts off some sectors - 90% sure that its different
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    auto basisTransformationMtx = OpenMesh::FProp<gmm::row_matrix<std::vector<double>>>(trimesh_,
                                                                                        "basisTransformationMtx");
    for (TriMesh::FaceVertexIter vh_i = trimesh_.fv_iter(fh); vh_i.is_valid(); ++vh_i) {
        checkBoundaryAndUniqueFace(*vh_i, fh, _H);
    }
}

void
GlobalParametrization::checkBoundaryAndUniqueFace(const OpenMesh::VertexHandle &vh_i, const OpenMesh::FaceHandle &fh,
                                                  CMatrixType &_H) {
    auto faceSel = OpenMesh::FProp<bool>(trimesh_, "faceSel");
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    auto basisTransformationMtx = OpenMesh::FProp<gmm::row_matrix<std::vector<double>>>(trimesh_,
                                                                                        "basisTransformationMtx");
    for (TriMesh::VertexOHalfedgeIter voh_it = trimesh_.voh_iter(vh_i); voh_it.is_valid(); ++voh_it) {
        OpenMesh::FaceHandle fh_neigh = voh_it->face();
        if (!trimesh_.is_boundary(*voh_it) && faceSel[fh_neigh]) {
            OpenMesh::VertexHandle vh_j = voh_it->to();
            gmm::row_matrix<std::vector<double>> D = basisTransformationMtx[fh_neigh];
            double area = trimesh_.calc_face_area(fh_neigh), weight = 1, h = 1;
            double col1 = mapLocCoordToGlobCoordSys(fh_neigh, vh_i);
            double col2 = mapLocCoordToGlobCoordSys(fh_neigh, vh_j);
            //get sum of ut*dkm
            double sum = 0;
            for (int j = 0; j < 3; ++j) {
                sum += D(j, col1) * D(j, col2);
            }
            _H(vertexPosUi[vh_i], vertexPosUi[vh_j]) += 2 * weight * area * pow(h, 2) * sum;
            _H(vertexPosUi[vh_j], vertexPosUi[vh_i]) += 2 * weight * area * pow(h, 2) * sum;
            _H(vertexPosVi[vh_i], vertexPosVi[vh_j]) += 2 * weight * area * pow(h, 2) * sum;
            _H(vertexPosVi[vh_j], vertexPosVi[vh_i]) += 2 * weight * area * pow(h, 2) * sum;
        }
    }
}

void GlobalParametrization::createSectorsCutGraph(std::vector<int> &singularities) {
    trimesh_.release_halfedge_status();
    trimesh_.request_halfedge_status();
    int sector = 1;
    for (int singularity: singularities) {
        auto vh = trimesh_.vertex_handle(singularity);
        if (checkIfLeaf(vh)) {
            std::vector<OpenMesh::HalfedgeHandle> startOfSectors;
            initVectorStartSec(singularity, startOfSectors);
            propagateForSectors(sector, startOfSectors, singularities);
        }
    }
}

void GlobalParametrization::initVectorStartSec(const int singularity,
                                               std::vector<OpenMesh::HalfedgeHandle> &startOfSectors) {
    auto cutGraphHe = OpenMesh::HProp<bool>(trimesh_, "cutGraphHe");
    auto vh = trimesh_.vertex_handle(singularity);
    //initialize values with first cycle
    for (TriMesh::VertexOHalfedgeIter voh_it = trimesh_.voh_iter(vh); voh_it.is_valid(); ++voh_it) {
        if (cutGraphHe[*voh_it] && cutGraphHe[voh_it->opp()]) {
            startOfSectors.push_back(*voh_it);
        }
    }
}

void GlobalParametrization::propagateForSectors(int &sector,
                                                const std::vector<OpenMesh::HalfedgeHandle> &startOfSectors,
                                                std::vector<int> &singularities) {
    for (auto heh: startOfSectors) {
        propagation(heh, sector, singularities);
    }
}

void
GlobalParametrization::propagation(OpenMesh::HalfedgeHandle &heh, int &sector, const std::vector<int> &singularities) {
    auto cutGraphHe = OpenMesh::HProp<bool>(trimesh_, "cutGraphHe");
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");
    OpenMesh::SmartHalfedgeHandle newOutgoHe = make_smart(heh, trimesh_), incHeOpp = newOutgoHe.opp();
    auto heToVertex = newOutgoHe.to();
    bool next_sing_found = false;
    cutGraphFZone[newOutgoHe.face()] = sector;
    int idx = -1;
    while (!next_sing_found) {
        std::vector<OpenMesh::HalfedgeHandle> outGoingHe;
        for (TriMesh::VertexOHalfedgeIter voh_it = trimesh_.voh_iter(heToVertex); voh_it.is_valid(); ++voh_it) {
            if (cutGraphHe[*voh_it]) {
                outGoingHe.push_back(*voh_it);
            }
        }
        //get pos of incHeOpp in outGoingHe.
        if (cutGraphHe[incHeOpp]) {
            auto it = find(outGoingHe.begin(), outGoingHe.end(), incHeOpp);
            if (it != outGoingHe.end()) {
                idx = it - outGoingHe.begin();
            }
            //get element after incHeOpp in outGoingHe
            if (idx == ((int) outGoingHe.size() - 1)) {
                //set incHeOpp = element
                newOutgoHe = make_smart(outGoingHe[0], trimesh_);
            } else {
                newOutgoHe = make_smart(outGoingHe[idx + 1], trimesh_);
            }
        } else {
            for (auto he: outGoingHe) {
                auto temp = make_smart(he, trimesh_);
                if (!cutGraphHe[temp.opp()]) {
                    newOutgoHe = temp;
                }
            }
        }

        //set heToVertex = element.to()
        heToVertex = newOutgoHe.to();
        //set cutGraphFZone = sector
        cutGraphFZone[newOutgoHe.face()] = sector;
        //check if heToVertex is singularity
        bool leafCheck = checkIfLeaf(heToVertex);
        if (leafCheck &&
            (std::find(singularities.begin(), singularities.end(), heToVertex.idx()) != singularities.end())) {
            next_sing_found = true;
        }
        //make it outgoing for the next cycle
        incHeOpp = newOutgoHe.opp();
    }
    sector++;
}

bool GlobalParametrization::checkIfLeaf(const OpenMesh::VertexHandle &heToVertex) {
    auto cutGraphHe = OpenMesh::HProp<bool>(trimesh_, "cutGraphHe");
    int counter = 0;
    for (TriMesh::VertexOHalfedgeIter voh_it = trimesh_.voh_iter(heToVertex); voh_it.is_valid(); ++voh_it) {
        if (cutGraphHe[*voh_it] && cutGraphHe[voh_it->opp()]) {
            counter++;
        }
    }
    if (counter == 1) {
        return true;
    } else {
        return false;
    }
}

void GlobalParametrization::fixRotationsCrossBoundaryComp(std::vector<int> &complementHEdges,
                                                          std::vector<int> &singularities, std::vector<int> &faces) {
    auto currentPJ = OpenMesh::FProp<int>(trimesh_, "currentPJ");
    trimesh_.request_face_status();
    trimesh_.release_halfedge_status();
    trimesh_.request_halfedge_status();
    auto temp = trimesh_.face_handle(faces[0]);
    std::queue<OpenMesh::FaceHandle> stack;
    setFaceStatusToFalse();
    currentPJ[temp] = 0;
    stack.push(temp);
    trimesh_.status(temp).set_tagged2(true);
    while (!stack.empty()) {
        for (TriMesh::FaceHalfedgeIter fhe_it = trimesh_.fh_iter(stack.front()); fhe_it.is_valid(); ++fhe_it) {
            updateStack(fhe_it, stack);
        }
        stack.pop();
    }
    trimesh_.release_face_status();
}

void GlobalParametrization::updateStack(OpenMesh::PolyConnectivity::FaceHalfedgeIter &he,
                                        std::queue<OpenMesh::FaceHandle> &stack) {
    auto cutGraphHe = OpenMesh::HProp<bool>(trimesh_, "cutGraphHe");
    auto faceSel = OpenMesh::FProp<bool>(trimesh_, "faceSel");
    auto ohe = he->opp();
    if (!(cutGraphHe[*he] || cutGraphHe[ohe]) && ohe.face().is_valid() && faceSel[ohe.face()] &&
        !trimesh_.status(ohe.face()).tagged2()) {
        std::pair<int, int> pj = getPJ(he, ohe);
        updatePJandCrossfield(pj, ohe);
        trimesh_.status(ohe.face()).set_tagged2(true);
        stack.push(ohe.face());
    }
}

std::pair<int, int>
GlobalParametrization::getPJ(TriMesh::FaceHalfedgeIter &fhe_it, OpenMesh::SmartHalfedgeHandle &ohe) {
    auto cutGraphHe = OpenMesh::HProp<bool>(trimesh_, "cutGraphHe");
    auto faceSel = OpenMesh::FProp<bool>(trimesh_, "faceSel");
    auto periodJump = OpenMesh::HProp<int>(trimesh_, "periodJump");
    std::pair<int, int> pj = {0, 0};
    if ((periodJump[*fhe_it] != INT_MAX && periodJump[*fhe_it] != 0)) {
        pj = {fhe_it->idx(), periodJump[*fhe_it]};
    } else if ((periodJump[ohe] != INT_MAX && periodJump[ohe] != 0)) {
        pj = {ohe.idx(), periodJump[ohe]};
    }
    return pj;
}

void GlobalParametrization::updatePJandCrossfield(std::pair<int, int> &pj, OpenMesh::SmartHalfedgeHandle &ohe) {
    ;
    auto currentPJ = OpenMesh::FProp<int>(trimesh_, "currentPJ");
    auto uVectorFieldRotOne = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotOne");
    auto uVectorFieldRotTwo = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotTwo");
    if (pj.first == ohe.opp().idx()) {
        int pjVal = -pj.second + currentPJ[ohe.opp().face()];
        uVectorFieldRotOne[ohe.face()] = rotPointWithRotMatrix(ohe.face(), uVectorFieldRotOne[ohe.face()],
                                                               M_PI / 2 *
                                                               pjVal);
        uVectorFieldRotTwo[ohe.face()] = rotPointWithRotMatrix(ohe.face(), uVectorFieldRotTwo[ohe.face()],
                                                               M_PI / 2 *
                                                               pjVal);
        currentPJ[ohe.face()] = pjVal;
    } else {
        int pjVal = pj.second + currentPJ[ohe.opp().face()];
        uVectorFieldRotOne[ohe.face()] = rotPointWithRotMatrix(ohe.face(), uVectorFieldRotOne[ohe.face()],
                                                               M_PI / 2 *
                                                               pjVal);
        uVectorFieldRotTwo[ohe.face()] = rotPointWithRotMatrix(ohe.face(), uVectorFieldRotTwo[ohe.face()],
                                                               M_PI / 2 *
                                                               pjVal);
        currentPJ[ohe.face()] = pjVal;
    }
}

GlobalParametrization::Point
GlobalParametrization::rotPointWithRotMatrix(const OpenMesh::FaceHandle fh, const Point vec, const double theta) {
    Point f_normal = trimesh_.calc_face_normal(fh).normalize();
    //Rodrigues rotation formula
    Point rotVec = vec * cos(theta) - (f_normal % vec) * sin(theta) + f_normal * (vec | f_normal) * (1 - cos(theta));
    return rotVec;
}

gmm::row_matrix<gmm::wsvector<double>> GlobalParametrization::getConstraints(const std::vector<int> &faces) {
    //1. boundary constraints, feature edge (edges which should be preserved)
    //2. cutedge constaints
    //3. cross boundary compatibility constraints
    gmm::row_matrix<gmm::wsvector<double>> _constraints(5, 5);
    return _constraints;
}

void GlobalParametrization::setFaceStatusToFalse() {

    for (auto face: trimesh_.faces()) {
        trimesh_.status(face).set_tagged(false);
    }
}
