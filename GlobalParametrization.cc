//
// Created by wuschelbueb on 05.05.22.
//

#include "GlobalParametrization.hh"
#include <fstream>

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
//    //just for visualization
    std::vector<int> onlyBoundaries = complementHEdges;
    std::vector<int> cutGraphWoBoundary;
//    add path from boundary to singularity to cut graph
    dualGraph.calcDijkstraWSingularities(complementHEdges, singularities, cutGraphWoBoundary);
    createSectorsOnCutGraph(singularities);
    fixRotationsCrossBoundaryComp(faces);
//    to visualize cutGraph with or without edges
    colorCompHEdges(complementHEdges);
    int nbVerticesUaV = createVertexPosParamDomain(faces);
    int jkValues = cutGraphWoBoundary.size();
    std::cout << "nbVerticesUaV: " << nbVerticesUaV << " and jkValues: " << jkValues << std::endl;

    std::vector<double> _x(nbVerticesUaV + jkValues, 0.0);
    std::vector<double> _rhs = getRhs(nbVerticesUaV, jkValues);
    std::vector<int> _idx_to_round = getIdxToRound(nbVerticesUaV, jkValues, singularities, onlyBoundaries);
    CMatrixType _hessian = getHessian(nbVerticesUaV, jkValues, cutGraphWoBoundary);
    RMatrixType _constraints = getConstraints(nbVerticesUaV, cutGraphWoBoundary, singularities);
    std::ofstream hessian("/home/wuschelbueb/Desktop/hessian_matrix.txt");
    hessian << _hessian;
    hessian.close();
    std::ofstream rhs("/home/wuschelbueb/Desktop/rhs.txt");
    for (auto i: _rhs) {
        rhs << i << std::endl;
    }
    rhs.close();
    std::ofstream idx_to_round("/home/wuschelbueb/Desktop/idx_to_round.txt");
    for (auto i: _idx_to_round) {
        idx_to_round << i << std::endl;
    }
    idx_to_round.close();
    std::ofstream cons("/home/wuschelbueb/Desktop/constraints.txt");
    cons << _constraints << std::endl;
    cons.close();
    std::cout << "Calculation of Global Parametrization\n" << "Dimensions of parameters:\n" << "Constraints Rows:\t"
              << gmm::mat_nrows(_constraints)
              << " and columns: " << gmm::mat_ncols(_constraints) << std::endl
              << "A Rows:\t\t\t\t" << gmm::mat_nrows(_hessian) << " and columns: " << gmm::mat_ncols(_hessian)
              << std::endl
              << "Size of _x:\t\t\t" << _x.size() << std::endl << "Size of _rhs:\t\t" << _rhs.size() << std::endl
              << "Size of idx:\t\t" << _idx_to_round.size() << std::endl;
    COMISO::ConstrainedSolver csolver;
    csolver.misolver().set_iter_full(false);
    csolver.misolver().set_local_iters(50000);
    csolver.misolver().set_cg_iters(20);
    csolver.misolver().set_local_error(1e-3);
    csolver.misolver().set_cg_error(1e-3);
    csolver.misolver().set_multiple_rounding();
    csolver.solve(_constraints, _hessian, _x, _rhs, _idx_to_round);
    saveSolAsCoord(_x, faces, cutGraphWoBoundary, nbVerticesUaV);
}

std::vector<int> GlobalParametrization::getSingularities() {
    auto crossFieldIdx = OpenMesh::VProp<double>(trimesh_, "crossFieldIdx");
    std::vector<int> singularities;
    for (auto vh: trimesh_.vertices()) {
        if (crossFieldIdx[vh] < -1E-1 || crossFieldIdx[vh] > 1E-1) {
            singularities.push_back(vh.idx());
        }
    }
    return singularities;
}

void GlobalParametrization::removeOpenPaths(std::vector<int> &complementHEdges) {
    trimesh_.request_halfedge_status();
    int currentIter = complementHEdges.size(), prevIter = INT_MAX;
    while (currentIter != prevIter) {
        for (int i: complementHEdges) {
            removeEdgeFromGraph(i, complementHEdges);
        }
        prevIter = currentIter;
        currentIter = complementHEdges.size();
    }
    trimesh_.release_halfedge_status();
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
        if (oheh.face().is_valid() && !faceSel[oheh.face()] && it2 != complementHEdges.end() &&
            borderDualG[eh] == 1) {
            complementHEdges.erase(it2);
        }
        if (heh.is_boundary() && it != complementHEdges.end() && borderDualG[eh] == 1) {
            complementHEdges.erase(it);
        } else if (oheh.is_boundary() && it2 != complementHEdges.end() && borderDualG[eh] == 1) {
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
    trimesh_.request_halfedge_status();
    for (auto heh: trimesh_.halfedges()) {
        trimesh_.status(heh).set_tagged(false);
    }
    std::vector<int> complementHEdges;
    tagEdgesFromDualSpanningTree();
    for (auto heh: trimesh_.halfedges()) {
        if (!trimesh_.status(heh).tagged()) {
            complementHEdges.push_back(heh.idx());
        }
    }
    trimesh_.release_halfedge_status();
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
    if (dualGraphPred[face] != INT_MAX && dualGraphPred[face] != face.idx() && faceSel[face]) {
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
        if (faceSel[face]) {
            faces.push_back(face.idx());
        }
    }
    return faces;
}

int GlobalParametrization::createVertexPosParamDomain(std::vector<int> &faces) {
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
    trimesh_.release_vertex_status();
    return countVertices;
}

void GlobalParametrization::checkCGandSetPos(OpenMesh::VertexHandle fv_it, int &countVertices) {
    auto vertexDist = OpenMesh::VProp<double>(trimesh_, "vertexDist");
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    if (vertexDist[fv_it] == 0) {
        getPositionInnerNode(fv_it, countVertices);
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
    if (checkIfLeaf(fv_it)) {
        appearance = 1;
    }
    vertexAppearanceCG[fv_it] = appearance;
    vertexPosUi[fv_it] = countVertices;
    vertexPosVi[fv_it] = (countVertices + appearance);
    countVertices += (appearance * 2);
    trimesh_.status(fv_it).set_tagged(true);
}

void GlobalParametrization::setUpLocFaceCoordSys(const std::vector<int> &faces) {
    gmm::col_matrix<std::vector<double>> _abc = createSomeMatrix();
    std::vector<Point> edges(3);
    for (int i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        createEdgesAndLocalVUi(fh, edges);
        createBasisTfMtx(fh, _abc, edges);
    }
}

void
GlobalParametrization::createEdgesAndLocalVUi(const OpenMesh::FaceHandle fh,
                                              std::vector<Point> &edges) {
    auto localVUi = OpenMesh::FProp<std::vector<int>>(trimesh_, "localVUi");
    int counter = 0;
    for (TriMesh::FaceHalfedgeIter fh_it = trimesh_.fh_iter(fh); fh_it.is_valid(); ++fh_it) {
        if (counter == 0) {
            Point temp = trimesh_.calc_edge_vector(*fh_it);
            localVUi[fh].push_back(trimesh_.to_vertex_handle(*fh_it).idx());
//            Point x = trimesh_.point(trimesh_.to_vertex_handle(*fh_it));
//            Point y = trimesh_.point(trimesh_.from_vertex_handle(*fh_it));
            localVUi[fh].push_back(trimesh_.from_vertex_handle(*fh_it).idx());
            // minus because common vertex is vec[0]
            edges[0] = -temp;
        } else if (counter == 1) {
            Point temp = trimesh_.calc_edge_vector(*fh_it);
            localVUi[fh].push_back(trimesh_.to_vertex_handle(*fh_it).idx());
//            Point z = trimesh_.point(trimesh_.to_vertex_handle(*fh_it));
            edges[1] = temp;
        }
        counter++;
    }
    // calculate normal of triangle
    edges[2] = (edges[0] % edges[1]);
}

void
GlobalParametrization::createBasisTfMtx(const OpenMesh::FaceHandle fh, gmm::col_matrix<std::vector<double>> &_abc,
                                        const std::vector<Point> &edges) {
    auto basisTransformationMtx = OpenMesh::FProp<gmm::row_matrix<std::vector<double>>>(trimesh_,
                                                                                        "basisTransformationMtx");
    gmm::row_matrix<std::vector<double>> _basisInv(3, 3);
    gmm::row_matrix<std::vector<double>> _transformationMatrix(3, 3);
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
    gmm::mult(_basisInv, _abc, _transformationMatrix);
//    std::cout << "D(" << fh.idx() << "): " << _transformationMatrix;
    basisTransformationMtx[fh] = _transformationMatrix;
}

gmm::col_matrix<std::vector<double>> GlobalParametrization::createSomeMatrix() {
    gmm::col_matrix<std::vector<double>> _abc(3, 3);
    _abc(0, 0) = -1;
    _abc(1, 0) = -1;
    _abc(2, 0) = 0;
    _abc(0, 1) = 1;
    _abc(1, 1) = 0;
    _abc(2, 1) = 0;
    _abc(0, 2) = 0;
    _abc(1, 2) = 1;
    _abc(2, 2) = 0;
    return _abc;
}

std::vector<double>
GlobalParametrization::getRhs(const int rhsSizePartOne, const int rhsSizePartTwo) {
    auto faceSel = OpenMesh::FProp<bool>(trimesh_, "faceSel");
    auto uVectorFieldRotOne = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotOne");
    auto uVectorFieldRotTwo = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotTwo");
    int size = rhsSizePartOne + rhsSizePartTwo;
    std::vector<double> _rhs(size, 0);

    // use nbVerticesUaV and lhs formulas with dkm -> dkm is entry in D matrix
    // utk is entry of vectorRotOne
    for (auto he: trimesh_.halfedges()) {
        if (!trimesh_.is_boundary(he) && faceSel[he.face()]) {
            Point uCrossField = uVectorFieldRotOne[he.face()];
            Point vCrossField = uVectorFieldRotTwo[he.face()];
            getRhsEntryForVertex(he, uCrossField, true, _rhs);
            getRhsEntryForVertex(he, vCrossField, false, _rhs);
        }
    }
    //set entries smaller than 1E-10 to zero
    for (size_t i = 0; i < _rhs.size(); ++i) {
        if (_rhs[i] < 1E-10 && _rhs[i] > -1E-10) {
            _rhs[i] = 0;
        }
    }
    return _rhs;
}

void
GlobalParametrization::getRhsEntryForVertex(const OpenMesh::SmartHalfedgeHandle he, const Point CrossFieldAxis,
                                            const bool flagUorV,
                                            std::vector<double> &_rhs) {
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    auto basisTransformationMtx = OpenMesh::FProp<gmm::row_matrix<std::vector<double>>>(trimesh_,
                                                                                        "basisTransformationMtx");
    OpenMesh::FaceHandle adjFace = he.face();
    OpenMesh::SmartVertexHandle vh = he.to();
    int posAppearance = 0, appearance = vertexAppearanceCG[vh];
    gmm::row_matrix<std::vector<double>> transformationMatrix = basisTransformationMtx[adjFace];
    double area = trimesh_.calc_face_area(adjFace), col = mapLocCoordToGlobCoordSys(adjFace, vh), h = 1;

    if (appearance > 1) {
        posAppearance = getPositionConstraintRow(vh, cutGraphFZone[adjFace]);
    }
    //get sum of ut*dkm
    double sum = 0;
    for (int j = 0; j < 3; ++j) {
        sum += CrossFieldAxis[j] * basisTransformationMtx[adjFace](j, col);
    }
    //check if U or V position is needed
    int UVPos = (flagUorV) ? vertexPosUi[vh] : vertexPosVi[vh];
//        std::cout << "position of vertex " << fv_it->idx() << " is " << UVPos + posAppearance << std::endl;
    _rhs[UVPos + posAppearance] += 2 * area * h * sum;
}

//map U coordinates to vertices
int
GlobalParametrization::mapLocCoordToGlobCoordSys(const OpenMesh::FaceHandle fh, const OpenMesh::VertexHandle vh) {
    auto localVUi = OpenMesh::FProp<std::vector<int>>(trimesh_, "localVUi");
    auto faceSel = OpenMesh::FProp<bool>(trimesh_, "faceSel");
    int column;
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
GlobalParametrization::getHessian(const int rhsSizePartOne, const int rhsSizePartTwo,
                                  std::vector<int> &cutGraphWoBoundary) {
//    rhsSizePartOne == u,v values of each vertex
//    rhsSizePartTwo == j&k values
    auto faceSel = OpenMesh::FProp<bool>(trimesh_, "faceSel");
    trimesh_.request_halfedge_status();
    trimesh_.request_face_status();
    int size = rhsSizePartOne + rhsSizePartTwo;
    CMatrixType _hessian(size, size);
    //initialize halfedges and faces
    for (auto he: trimesh_.halfedges()) {
        if (!trimesh_.is_boundary(he) && faceSel[he.face()]) {
            getDiaEntriesHessian(he, _hessian);
            getEntriesHessian(he, _hessian);
        }
    }
//    for (int i = rhsSizePartOne; i < rhsSizePartOne + rhsSizePartTwo; ++i) {
//        _hessian(i, i) = 1;
//    }
    trimesh_.release_halfedge_status();
    trimesh_.release_face_status();
    return _hessian;
}

void GlobalParametrization::getDiaEntriesHessian(const OpenMesh::SmartHalfedgeHandle he, CMatrixType &_hessian) {
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    auto basisTransformationMtx = OpenMesh::FProp<gmm::row_matrix<std::vector<double>>>(trimesh_,
                                                                                        "basisTransformationMtx");
    OpenMesh::FaceHandle adjFace = he.face();
    int posOne = 0, idx;
    double area = trimesh_.calc_face_area(adjFace), weight = 1, h = 1;
    OpenMesh::SmartVertexHandle vh = he.to();
    int appearance = vertexAppearanceCG[vh];
    gmm::row_matrix<std::vector<double>> transformationMatrix = basisTransformationMtx[adjFace];
    double col = mapLocCoordToGlobCoordSys(adjFace, vh);
    if (appearance > 1) {
        posOne = getPositionConstraintRow(vh, cutGraphFZone[adjFace]);
    }
    //get sum of ut*dkm
    double sum = 0;
    for (int j = 0; j < 3; ++j) {
        sum += transformationMatrix(j, col) * transformationMatrix(j, col);
    }
    _hessian(vertexPosUi[vh] + posOne, vertexPosUi[vh] + posOne) += 2 * weight * area * pow(h, 2) * sum;
    _hessian(vertexPosVi[vh] + posOne, vertexPosVi[vh] + posOne) += 2 * weight * area * pow(h, 2) * sum;
}

void GlobalParametrization::getEntriesHessian(const OpenMesh::SmartHalfedgeHandle he, CMatrixType &_hessian) {
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    auto basisTransformationMtx = OpenMesh::FProp<gmm::row_matrix<std::vector<double>>>(trimesh_,
                                                                                        "basisTransformationMtx");
    OpenMesh::FaceHandle adjFace = he.face();
//            std::cout << "Face " << voh_it->face().idx() << "\n start V " << vh_i.idx() << "\t end V " << voh_it->to()
//                      << std::endl;
    int posOne = 0, posTwo = 0, idx;
    OpenMesh::SmartVertexHandle vh_i = he.from();
    OpenMesh::SmartVertexHandle vh_j = he.to();
    int appearanceI = vertexAppearanceCG[vh_i];
    int appearanceJ = vertexAppearanceCG[vh_j];
    gmm::row_matrix<std::vector<double>> D = basisTransformationMtx[adjFace];
    double area = trimesh_.calc_face_area(adjFace), weight = 1, h = 1;
    double col1 = mapLocCoordToGlobCoordSys(adjFace, vh_i);
    double col2 = mapLocCoordToGlobCoordSys(adjFace, vh_j);
//    std::cout << "face: " << he.face().idx() << "\nvertex 1: " << vh_i.idx() << "\tvertex 2: " << vh_j.idx()
//              << "\ncol 1: " << col1 << "\ncol 2: " << col2 << "\nappearance i: "
//              << appearanceI << "\nappearance j: " << appearanceJ << std::endl;
    if (appearanceI > 1) {
        posOne = getPositionConstraintRow(vh_i, cutGraphFZone[he.face()]);
    }
    if (appearanceJ > 1) {
        posTwo = getPositionConstraintRow(vh_j, cutGraphFZone[he.face()]);
    }
    //get sum of ut*dkm
    double sum = 0;
//            std::cout << "area " << area << " col 1: " << col1 << "\tcol 2: " << col2 << "\nsum: ";
    for (int j = 0; j < 3; ++j) {
        sum += D(j, col1) * D(j, col2);
//                std::cout << sum << "\t";
    }
//            std::cout << "\n";
    // u position
    _hessian(vertexPosUi[vh_i] + posOne, vertexPosUi[vh_j] + posTwo) += 2 * weight * area * pow(h, 2) * sum;
    _hessian(vertexPosUi[vh_j] + posTwo, vertexPosUi[vh_i] + posOne) += 2 * weight * area * pow(h, 2) * sum;

    // v position
    _hessian(vertexPosVi[vh_i] + posOne, vertexPosVi[vh_j] + posTwo) += 2 * weight * area * pow(h, 2) * sum;
    _hessian(vertexPosVi[vh_j] + posTwo, vertexPosVi[vh_i] + posOne) += 2 * weight * area * pow(h, 2) * sum;
//    std::cout << "hessian U val (" << vertexPosUi[vh_i] + posOne << "," << vertexPosUi[vh_j] + posTwo << ") "
//              << _hessian(vertexPosUi[vh_i] + posOne, vertexPosUi[vh_j] + posTwo) << std::endl;
}

void GlobalParametrization::createSectorsOnCutGraph(std::vector<int> &singularities) {
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
    trimesh_.release_halfedge_status();
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
    OpenMesh::SmartHalfedgeHandle newOutgoHe = make_smart(heh, trimesh_);
    OpenMesh::SmartHalfedgeHandle newOutgoHeOpp = newOutgoHe.opp();
    auto toVhOfNewOutgoHe = newOutgoHe.to();
    bool next_sing_found = false;
    cutGraphFZone[newOutgoHe.face()] = sector;
    while (!next_sing_found) {
        std::vector<OpenMesh::HalfedgeHandle> vecOfOutgoHe = getListOfOutgoHe(toVhOfNewOutgoHe);
        getNextHe(newOutgoHe, newOutgoHeOpp, vecOfOutgoHe);

        // special case; is necessary to avoid wrong coloring of cutgraphs which are only 1 edge long
        if (newOutgoHe == newOutgoHeOpp) {
            cutGraphFZone[newOutgoHe.face()] = ++sector;
            break;
        }
        //set toVhOfNewOutgoHe = element.to()
        toVhOfNewOutgoHe = newOutgoHe.to();
        //make it outgoing for the next cycle
        newOutgoHeOpp = newOutgoHe.opp();

        colorFaces(newOutgoHe, sector);
        //check if toVhOfNewOutgoHe is singularity
        if (checkIfLeaf(toVhOfNewOutgoHe) &&
            (std::find(singularities.begin(), singularities.end(), toVhOfNewOutgoHe.idx()) != singularities.end())) {
            next_sing_found = true;
        }
    }
    sector++;
}

std::vector<OpenMesh::HalfedgeHandle>
GlobalParametrization::getListOfOutgoHe(OpenMesh::SmartVertexHandle toVhOfNewOutgoHe) {
    auto cutGraphHe = OpenMesh::HProp<bool>(trimesh_, "cutGraphHe");
    std::vector<OpenMesh::HalfedgeHandle> vecOfOutgoHe;
    //fill vector with halfedges who are in cutgraph
    for (TriMesh::VertexOHalfedgeIter voh_it = trimesh_.voh_iter(toVhOfNewOutgoHe); voh_it.is_valid(); ++voh_it) {
        if (cutGraphHe[*voh_it]) {
            vecOfOutgoHe.push_back(*voh_it);
        }
    }
    return vecOfOutgoHe;
}

void GlobalParametrization::getNextHe(OpenMesh::SmartHalfedgeHandle &newOutgoHe,
                                      OpenMesh::SmartHalfedgeHandle &newOutgoHeOpp,
                                      std::vector<OpenMesh::HalfedgeHandle> &vecOfOutgoHe) {
    auto cutGraphHe = OpenMesh::HProp<bool>(trimesh_, "cutGraphHe");
    int idx = -1;
//get pos of newOutgoHeOpp in vecOfOutgoHe.
    if (cutGraphHe[newOutgoHeOpp]) {
        auto it = find(vecOfOutgoHe.begin(), vecOfOutgoHe.end(), newOutgoHeOpp);
        if (it != vecOfOutgoHe.end()) {
            idx = it - vecOfOutgoHe.begin();
        }
        //get next element in vecOfOutgoHe
        if (idx == ((int) vecOfOutgoHe.size() - 1)) {
            //set newOutgoHeOpp = element
            newOutgoHe = make_smart(vecOfOutgoHe[0], trimesh_);
        } else {
            newOutgoHe = make_smart(vecOfOutgoHe[idx + 1], trimesh_);
        }
    } else {
        for (auto he: vecOfOutgoHe) {
            auto temp = make_smart(he, trimesh_);
            if (!cutGraphHe[temp.opp()]) {
                newOutgoHe = temp;
            }
        }
    }
}

void GlobalParametrization::colorFaces(OpenMesh::SmartHalfedgeHandle newOutgoHe, int sector) {
    auto cutGraphHe = OpenMesh::HProp<bool>(trimesh_, "cutGraphHe");
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");
    auto faceOppOfHeAfterToVertex = newOutgoHe.next().opp().face();
    auto faceOppOfHeBeforeToVertex = newOutgoHe.prev().opp().face();
    //set cutGraphFZone = sector
    cutGraphFZone[newOutgoHe.face()] = sector;
    if (!cutGraphHe[newOutgoHe.next()] && faceOppOfHeAfterToVertex.is_valid()) {
        cutGraphFZone[faceOppOfHeAfterToVertex] = sector;
    }
    if (!cutGraphHe[newOutgoHe.prev()] && faceOppOfHeBeforeToVertex.is_valid()) {
        cutGraphFZone[faceOppOfHeBeforeToVertex] = sector;
    }
}

bool GlobalParametrization::checkIfLeaf(const OpenMesh::VertexHandle &heToVertex) {
    auto cutGraphHe = OpenMesh::HProp<bool>(trimesh_, "cutGraphHe");
    int counter = 0;
    for (TriMesh::VertexOHalfedgeIter voh_it = trimesh_.voh_iter(heToVertex); voh_it.is_valid(); ++voh_it) {
        //&& cutGraphHe[voh_it->opp()] -> prob not necessary was to avoid
        if (cutGraphHe[*voh_it]) {
            counter++;
        }
    }
    if (counter == 1) {
        return true;
    } else {
        return false;
    }
}

void GlobalParametrization::fixRotationsCrossBoundaryComp(std::vector<int> &faces) {
    auto currentPJ = OpenMesh::FProp<int>(trimesh_, "currentPJ");
    trimesh_.request_face_status();
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
    trimesh_.release_halfedge_status();
    trimesh_.release_face_status();
}

void GlobalParametrization::updateStack(OpenMesh::PolyConnectivity::FaceHalfedgeIter &he,
                                        std::queue<OpenMesh::FaceHandle> &stack) {
    auto cutGraphHe = OpenMesh::HProp<bool>(trimesh_, "cutGraphHe");
    auto faceSel = OpenMesh::FProp<bool>(trimesh_, "faceSel");
    auto ohe = he->opp();
    if (!(cutGraphHe[*he] || cutGraphHe[ohe]) && ohe.face().is_valid() && faceSel[ohe.face()] &&
        !trimesh_.status(ohe.face()).tagged2()) {
        //std::pair<halfedge idx, pj>
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
    auto currentPJ = OpenMesh::FProp<int>(trimesh_, "currentPJ");
    auto uVectorFieldRotOne = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotOne");
    auto uVectorFieldRotTwo = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotTwo");
    if (pj.first == ohe.opp().idx()) {
        int pjVal = -pj.second + currentPJ[ohe.opp().face()];
        uVectorFieldRotOne[ohe.face()] =
                rotPointWithRotMatrix(ohe.face(), uVectorFieldRotOne[ohe.face()], M_PI / 2 * pjVal);
        uVectorFieldRotTwo[ohe.face()] =
                rotPointWithRotMatrix(ohe.face(), uVectorFieldRotTwo[ohe.face()], M_PI / 2 * pjVal);
        currentPJ[ohe.face()] = pjVal;
    } else {
        int pjVal = pj.second + currentPJ[ohe.opp().face()];
        uVectorFieldRotOne[ohe.face()] =
                rotPointWithRotMatrix(ohe.face(), uVectorFieldRotOne[ohe.face()], M_PI / 2 * pjVal);
        uVectorFieldRotTwo[ohe.face()] =
                rotPointWithRotMatrix(ohe.face(), uVectorFieldRotTwo[ohe.face()], M_PI / 2 * pjVal);
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

gmm::row_matrix<gmm::wsvector<double>>
GlobalParametrization::getConstraints(const int nbVerticesUaV, std::vector<int> &cutGraphWoBoundary,
                                      std::vector<int> &singularities) {
    trimesh_.request_vertex_status();
    trimesh_.request_halfedge_status();
    int row_size = getRowSizeAndSetStatus();
    int col_size = nbVerticesUaV + cutGraphWoBoundary.size(); //vertex U a
    // nd V positions + j & k values
    int singularity = 2, jkStartCounter = nbVerticesUaV;
//    gmm::row_matrix<gmm::wsvector<double>> _constraints(singularity, col_size + 1);
    gmm::row_matrix<gmm::wsvector<double>> _constraints(row_size + singularity, col_size + 1);
//    std::cout << "constraints size " << gmm::mat_nrows(_constraints) << "x" << gmm::mat_ncols(_constraints)
//              << " and jk couting starts at " << jkStartCounter << std::endl;
    setZeroPointConstraint(singularities, _constraints);
    getConstraintsMatrix(jkStartCounter, cutGraphWoBoundary, _constraints);
    trimesh_.release_vertex_status();
    return _constraints;
}

int GlobalParametrization::getRowSizeAndSetStatus() {
    int row_size = 0;
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    for (auto he: trimesh_.halfedges()) {
        trimesh_.status(he).set_tagged(false);
    }
    for (auto vh: trimesh_.vertices()) {
        trimesh_.status(vh).set_tagged(true);
        if (vertexAppearanceCG[vh] > 1) {
            row_size += vertexAppearanceCG[vh] * 2; //u&v values
            trimesh_.status(vh).set_tagged(false);
        }
    }
    return row_size;
}

void GlobalParametrization::setZeroPointConstraint(std::vector<int> &singularities,
                                                   gmm::row_matrix<gmm::wsvector<double>> &_constraints) {
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    // set point (0,0) of coord system from a leaf singularity
    for (auto i: singularities) {
        auto vh = trimesh_.vertex_handle(i);
        if (checkIfLeaf(vh)) {
//            std::cout << "origin singularity: " << vh.idx() << std::endl;
            _constraints(0, vertexPosUi[vh]) = 1;
            _constraints(1, vertexPosVi[vh]) = 1;
            break;
        }
    }
    if (singularities.empty()) {
        for (auto vh: trimesh_.vertices()) {
//            std::cout << "origin non-singularity: " << vh.idx() << std::endl;
            _constraints(0, vertexPosUi[vh]) = 1;
            _constraints(1, vertexPosVi[vh]) = 1;
            break;
        }
    }
}

void
GlobalParametrization::getConstraintsMatrix(int &jkStartCounter, std::vector<int> &cutGraphWoBoundary,
                                            gmm::row_matrix<gmm::wsvector<double>> &_constraints) {
    auto currentPJ = OpenMesh::FProp<int>(trimesh_, "currentPJ");
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    int counter = 2; // because singularity constraint is at pos 0&1

    //jkstartcounter == nbOfVerticesUaV
    for (auto it = cutGraphWoBoundary.begin(); it != cutGraphWoBoundary.end(); it++) {
        auto he1 = make_smart(trimesh_.halfedge_handle(*it++), trimesh_);
        auto he2 = make_smart(trimesh_.halfedge_handle(*it), trimesh_);
//        std::cout << "halfedge " << he1.idx() << " (" << he1.tagged() << ") with vertex " << he1.to().idx() << " ("
//                  << he1.to().tagged() << ")\nhalfedge " << he2.idx() << " (" << he2.tagged() << ") with vertex "
//                  << he2.to().idx() << " (" << he2.to().tagged() << ")\n";
        try {
            if (he1.opp().idx() != he2.idx()) {
                throw 404;
            }
        } catch (int x) {
            std::cerr << "getConstraintMatrix: halfedges have to be opposite, something went wrong\n";
        }
        int diff = (currentPJ[he1.face()] + currentPJ[he2.face()]) % 4;
        if (!he1.tagged() && !he1.to().tagged()) {
//            std::cout << "enter he1 condition\n";
            auto vh1 = he1.to();
            setConRows(counter, jkStartCounter, diff, he1, vh1, _constraints);
            trimesh_.status(he1).set_tagged(true);
        }
        if (!he2.tagged() && !he2.to().tagged()) {
//            auto vh2 = he2.to();
//            std::cout << "enter he2 condition\n";
            auto vh2 = he1.from();
            setConRows(counter, jkStartCounter, diff, he1, vh2, _constraints);
            trimesh_.status(he2).set_tagged(true);
        }
        jkStartCounter += 2;
    }
}

void
GlobalParametrization::setConRows(int &counter, int &jkStartCounter, const int diff,
                                  OpenMesh::SmartHalfedgeHandle &he, OpenMesh::SmartVertexHandle &vh,
                                  gmm::row_matrix<gmm::wsvector<double>> &_constraints) {
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");
    int posOne = getPositionConstraintRow(vh, cutGraphFZone[he.face()]);
    int posTwo = getPositionConstraintRow(vh, cutGraphFZone[he.opp().face()]);
//    std::cout << "vertex " << vh.idx() << " with position U/V " << vertexPosUi[vh] << " / " << vertexPosVi[vh]
//              << "\nhalfedge " << he.idx() << "\nFace " << he.face().idx() << "\tposOne " << posOne << "\nFace Opp "
//              << he.opp().face().idx() << "\tposTwo " << posTwo << std::endl;
    switch (diff) {
        case 0: //0
            // u position of point p
            _constraints(counter, vertexPosUi[vh] + posTwo) = -1; //u_p
            _constraints(counter, vertexPosUi[vh] + posOne) = 1; // u'_p
            _constraints(counter++, jkStartCounter) = -1; //j
            // v position of point p
            _constraints(counter, vertexPosVi[vh] + posTwo) = -1;
            _constraints(counter, vertexPosVi[vh] + posOne) = 1;
            _constraints(counter++, jkStartCounter + 1) = -1;
            break;
        case -3://-270 || 90
        case 1:
            // u position of point p
            _constraints(counter, vertexPosVi[vh] + posTwo) = 1; //-v_p
            _constraints(counter, vertexPosUi[vh] + posOne) = 1; //u'_p
            _constraints(counter++, jkStartCounter) = -1; //j
            // v position of point p
            _constraints(counter, vertexPosUi[vh] + posTwo) = -1;
            _constraints(counter, vertexPosVi[vh] + posOne) = 1;
            _constraints(counter++, jkStartCounter + 1) = -1;
            break;
        case -2://-180 || 180
        case 2:
            // u position of point p
            _constraints(counter, vertexPosUi[vh] + posTwo) = 1; //-u_p
            _constraints(counter, vertexPosUi[vh] + posOne) = 1; // u'_p
            _constraints(counter++, jkStartCounter) = -1; //j
            // v position of point p
            _constraints(counter, vertexPosVi[vh] + posTwo) = 1;
            _constraints(counter, vertexPosVi[vh] + posOne) = 1;
            _constraints(counter++, jkStartCounter + 1) = -1;
            break;
        case -1://-90 || 270
        case 3:
            // u position of point p
            _constraints(counter, vertexPosVi[vh] + posTwo) = -1; //v_p
            _constraints(counter, vertexPosUi[vh] + posOne) = 1; //u'_p
            _constraints(counter++, jkStartCounter) = -1; //j
            // v position of point p
            _constraints(counter, vertexPosUi[vh] + posTwo) = 1;
            _constraints(counter, vertexPosVi[vh] + posOne) = 1;
            _constraints(counter++, jkStartCounter + 1) = -1;
            break;
        default:
            break;
    }
}

int GlobalParametrization::getPositionConstraintRow(OpenMesh::SmartVertexHandle &vh, int cutGraphZone) {
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    std::vector<int> sectors;
    for (TriMesh::VertexFaceIter vf_it = trimesh_.vf_begin(vh); vf_it.is_valid(); ++vf_it) {
        int zone = cutGraphFZone[*vf_it];
        if (zone != 0 && (std::find(sectors.begin(), sectors.end(), zone) == sectors.end())) {
            sectors.push_back(zone);
        }
    }
    try {
        if ((int) sectors.size() != vertexAppearanceCG[vh]) {
            throw 404;
        }
    } catch (int x) {
        std::cerr << "getPositionConstraintRow: vector " << vh.idx()
                  << " doesn't have the same size as vertexAppearance property v: "
                  << sectors.size()
                  << " != p: " << vertexAppearanceCG[vh] << "\n";
    }
    std::sort(sectors.begin(), sectors.end());
    auto it1 = find(sectors.begin(), sectors.end(), cutGraphZone);
    int pos = it1 - sectors.begin();
//    std::cout << "for vertex " << vh.idx() << " with cutzoneGraph " << cutGraphZone << " we get position " << pos
//              << std::endl << "sector vector\n";
//    for (auto j: sectors) {
//        std::cout << j << "\t";
//    }
//    std::cout << std::endl;
    return pos;
}

void GlobalParametrization::setFaceStatusToFalse() {
    for (auto face: trimesh_.faces()) {
        trimesh_.status(face).set_tagged(false);
    }
}

std::vector<int> GlobalParametrization::getIdxToRound(int nbVerticesUaV, int jkValues, std::vector<int> &singularities,
                                                      std::vector<int> &onlyBoundaries) {
    trimesh_.request_vertex_status();
    for (auto vh: trimesh_.vertices()) {
        trimesh_.status(vh).set_tagged(false);
    }
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    std::vector<int> _idx_to_round;
    //boundary constraints
    for (auto i: onlyBoundaries) {
        auto he = make_smart(trimesh_.halfedge_handle(i), trimesh_);
        if (!he.to().tagged()) {
            _idx_to_round.push_back(vertexPosUi[he.to()]);
            _idx_to_round.push_back(vertexPosVi[he.to()]);
            trimesh_.status(he.to()).set_tagged(true);
            std::cout << "idx to round, boundary " << he.to().idx() << std::endl;
        }
        if (!he.from().tagged()) {
            _idx_to_round.push_back(vertexPosUi[he.from()]);
            _idx_to_round.push_back(vertexPosVi[he.from()]);
            trimesh_.status(he.from()).set_tagged(true);
            std::cout << "idx to round, boundary " << he.to().idx() << std::endl;
        }

    }
    //feature edge constraint; not yet implemented

    //singularity constraint
    for (auto i: singularities) {
        auto vh = make_smart(trimesh_.vertex_handle(i), trimesh_);
        if (!vh.tagged()) {
            for (int j = 0; j < vertexAppearanceCG[vh]; ++j) {
                _idx_to_round.push_back(vertexPosUi[vh] + j);
                _idx_to_round.push_back(vertexPosVi[vh] + j);
            }
            std::cout << "idx to round, singularity " << i << std::endl;
//            _idx_to_round.push_back(vertexPosUi[vh]);
//            _idx_to_round.push_back(vertexPosVi[vh]);
            trimesh_.status(vh).set_tagged(true);
        }
    }
    //jk values are integers as well
    for (int i = nbVerticesUaV; i < nbVerticesUaV + jkValues; ++i) {
        std::cout << "idx to round, jk Values " << i << std::endl;
        _idx_to_round.push_back(i);
    }
    trimesh_.release_vertex_status();
    return _idx_to_round;
}

void GlobalParametrization::saveSolAsCoord(std::vector<double> &_x, std::vector<int> &faces,
                                           std::vector<int> &cutGraphWoBoundary, const int nbVerticesUaV) {
    int jkStartValue = nbVerticesUaV;
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    auto solCoordSysUV = OpenMesh::VProp<std::vector<OpenMesh::Vec2d>>(trimesh_, "solCoordSysUV");
    auto faceSel = OpenMesh::FProp<bool>(trimesh_, "faceSel");

    for (auto he: trimesh_.halfedges()) {
        if (!trimesh_.is_boundary(he) && faceSel[he.face()]) {
            if (vertexAppearanceCG[he.to()] < 2) {
                solCoordSysUV[he.to()] = std::vector<OpenMesh::Vec2d>(vertexAppearanceCG[he.to()], {0, 0});
                getSolFromVerticesWOneApp(he, _x);
            } else {
                solCoordSysUV[he.to()] = std::vector<OpenMesh::Vec2d>(vertexAppearanceCG[he.to()], {0, 0});
            }
        }
    }
    for (auto it = cutGraphWoBoundary.begin(); it != cutGraphWoBoundary.end(); it++) {
        auto he1 = make_smart(trimesh_.halfedge_handle(*it++), trimesh_);
        auto he2 = make_smart(trimesh_.halfedge_handle(*it), trimesh_);
//        std::cout << "halfedge " << he1.idx() << " (" << he1.tagged() << ") with vertex " << he1.to().idx() << " ("
//                  << he1.to().tagged() << ")\nhalfedge " << he2.idx() << " (" << he2.tagged() << ") with vertex "
//                  << he2.to().idx() << " (" << he2.to().tagged() << ")\n";
        try {
            if (he1.opp().idx() != he2.idx()) {
                throw 404;
            }
        } catch (int x) {
            std::cerr << "getConstraintMatrix: halfedges have to be opposite, something went wrong\n";
        }
        getSolFromVerticesWMoreOneApp(he1, _x, jkStartValue);
        getSolFromVerticesWMoreOneApp(he2, _x, jkStartValue);
        jkStartValue += 2;
    }

    std::ofstream xVector("/home/wuschelbueb/Desktop/xVectorGlobalParam.txt");
    for (auto vh: trimesh_.vertices()) {
        xVector << "vertex (" << solCoordSysUV[vh].size() << ") " << vh.idx() << " with coord:\n";
        for (size_t i = 0; i < solCoordSysUV[vh].size(); ++i) {
            xVector << "x[" << i << "]: " << solCoordSysUV[vh][i][0] << "\ny[" << i << "]: " << solCoordSysUV[vh][i][1]
                    << std::endl;
        }
    }
    xVector << "\n";
    for (auto it = _x.begin(); it != _x.end(); ++it) {
        xVector << *it << std::endl;
    }
    xVector.close();
}

void GlobalParametrization::getSolFromVerticesWMoreOneApp(OpenMesh::SmartHalfedgeHandle he, std::vector<double> &_x,
                                                          int &jkStartValues) {
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    auto solCoordSysUV = OpenMesh::VProp<std::vector<OpenMesh::Vec2d>>(trimesh_, "solCoordSysUV");
    if (vertexAppearanceCG[he.to()] > 1) {
        auto vh = he.to();
        int posOne = getPositionConstraintRow(vh, cutGraphFZone[he.face()]);
        int posU = vertexPosUi[vh] + posOne;
        int posV = vertexPosVi[vh] + posOne;
        //todo not add jkvalues
        double u = _x[posU] + _x[jkStartValues];
        double v = _x[posV] + _x[jkStartValues + 1];
//        std::cout << "vertex " << vh.idx() << "\twith PosU " << posU - posOne << "\twith posOne: " << posOne
//                  << "\nwith jk: " << _x[jkStartValues - 2] << "\tand: " << _x[jkStartValues - 1] << std::endl;
        solCoordSysUV[vh][posOne] = {u, v};
    }
}

void GlobalParametrization::getSolFromVerticesWOneApp(OpenMesh::SmartHalfedgeHandle he, std::vector<double> &_x) {
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    auto solCoordSysUV = OpenMesh::VProp<std::vector<OpenMesh::Vec2d>>(trimesh_, "solCoordSysUV");
    int posU = vertexPosUi[he.to()];
    int posV = vertexPosVi[he.to()];
    double u = _x[posU];
    double v = _x[posV];
    OpenMesh::Vec2d p = {u, v};
    solCoordSysUV[he.to()][0] = p;
}



