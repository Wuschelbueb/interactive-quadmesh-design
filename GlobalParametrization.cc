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
    dualGraph.completeDijkstraWSingularities(complementHEdges, singularities, cutGraphWoBoundary);
    createSectorsCutGraph(singularities);
    fixRotationsCrossBoundaryComp(complementHEdges, singularities, faces);
//    to visualize cutGraph with or without edges
    colorCompHEdges(complementHEdges);
    int nbVerticesUaV = createVertexPosParamDomain(faces);
    int jkValues = cutGraphWoBoundary.size();
    std::cout << "nbVerticesUaV: " << nbVerticesUaV << " and jkValues: " << jkValues << std::endl;
//
    std::vector<double> _x(nbVerticesUaV + jkValues, 0.0);
    std::vector<double> _rhs = getRhs(faces, nbVerticesUaV, jkValues);
    std::vector<int> _idx_to_round = getIdxToRound(nbVerticesUaV, jkValues, singularities, onlyBoundaries);
    CMatrixType _H = getHessian(nbVerticesUaV, jkValues);
    RMatrixType _constraints = getConstraints(nbVerticesUaV, cutGraphWoBoundary, singularities);
//    std::ofstream hessian("/home/wuschelbueb/Desktop/hessian_matrix.txt");
//    hessian << _H;
//    hessian.close();
//    std::ofstream rhs("/home/wuschelbueb/Desktop/rhs.txt");
//    for (auto i: _rhs) {
//        rhs << "element: " << i << std::endl;
//    }
//    rhs.close();
//    std::ofstream idx_to_round("/home/wuschelbueb/Desktop/idx_to_round.txt");
//    for (auto i: _idx_to_round) {
//        idx_to_round << i << std::endl;
//    }
//    idx_to_round.close();
//    std::ofstream cons("/home/wuschelbueb/Desktop/constraints.txt");
//    cons << _constraints << std::endl;
//    cons.close();
    std::cout << "Calculation of Global Parametrization\n" << "Dimensions of parameters:\n" << "Constraints Rows:\t"
              << gmm::mat_nrows(_constraints)
              << " and columns: " << gmm::mat_ncols(_constraints) << std::endl
              << "A Rows:\t\t\t\t" << gmm::mat_nrows(_H) << " and columns: " << gmm::mat_ncols(_H) << std::endl
              << "Size of _x:\t\t\t" << _x.size() << std::endl << "Size of _rhs:\t\t" << _rhs.size() << std::endl
              << "Size of idx:\t\t" << _idx_to_round.size() << std::endl;
    COMISO::ConstrainedSolver csolver;
    csolver.misolver().set_iter_full(false);
    csolver.misolver().set_local_iters(50000);
    csolver.misolver().set_cg_iters(20);
    csolver.misolver().set_local_error(1e-3);
    csolver.misolver().set_cg_error(1e-3);
    csolver.misolver().set_multiple_rounding();
    csolver.solve(_constraints, _H, _x, _rhs, _idx_to_round);
    saveSolAsCoord(_x, faces, singularities, cutGraphWoBoundary, nbVerticesUaV, jkValues);
//    std::ofstream xVector("/home/wuschelbueb/Desktop/xVectorGlobalParam.txt");
//    for (auto it = _x.begin(); it != _x.end(); ++it) {
//        xVector << "Pos U: " << *it++ << " Pos V: " << *it << std::endl;
//    }
//    xVector.close();
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
    auto localVUi = OpenMesh::FProp<std::vector<int >>
            (trimesh_, "localVUi");
    int counter = 0;
    for (TriMesh::FaceHalfedgeIter fh_it = trimesh_.fh_iter(fh); fh_it.is_valid(); ++fh_it) {
        if (counter == 0) {
            Point temp = trimesh_.calc_edge_vector(*fh_it);
            localVUi[fh].push_back(trimesh_.to_vertex_handle(*fh_it).idx());
            Point x = trimesh_.point(trimesh_.to_vertex_handle(*fh_it));
            Point y = trimesh_.point(trimesh_.from_vertex_handle(*fh_it));
            localVUi[fh].push_back(trimesh_.from_vertex_handle(*fh_it).idx());
            // minus because common vertex is vec[0]

            edges[0] = -temp;
        } else if (counter == 1) {
            Point temp = trimesh_.calc_edge_vector(*fh_it);
            localVUi[fh].push_back(trimesh_.to_vertex_handle(*fh_it).idx());
            Point z = trimesh_.point(trimesh_.to_vertex_handle(*fh_it));
            edges[1] = temp;
        }
        counter++;
    }
    // calculate normal of triangle
    edges[2] = (edges[0] % edges[1]);
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
//    std::cout << "D(" << fh.idx() << "): " << _D;
    basisTransformationMtx[fh] = _D;
}

gmm::col_matrix<std::vector<double>> GlobalParametrization::createCMatrix() {
    gmm::col_matrix<std::vector<double>> _C(3, 3);
    _C(0, 0) = -1;
    _C(1, 0) = -1;
    _C(2, 0) = 0;
    _C(0, 1) = 1;
    _C(1, 1) = 0;
    _C(2, 1) = 0;
    _C(0, 2) = 0;
    _C(1, 2) = 1;
    _C(2, 2) = 0;
    return _C;
}

std::vector<double>
GlobalParametrization::getRhs(const std::vector<int> &faces, const int rhsSizePartOne, const int rhsSizePartTwo) {
    auto uVectorFieldRotOne = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotOne");
    auto uVectorFieldRotTwo = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotTwo");
    int size = rhsSizePartTwo + rhsSizePartOne;
    std::vector<double> _rhs(size, 0);

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
    //set entries smaller than 1E-10 to zero
    for (std::size_t i = 0; i < _rhs.size(); ++i) {
        if (_rhs[i] < 1E-10 && _rhs[i] > -1E-10) {
            _rhs[i] = 0.0;
        }
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
    double weight = 1, area = trimesh_.calc_face_area(fh), h = 1, col;
//    std::cout << "face " << fh.idx() << " area: " << area << std::endl;
    for (TriMesh::FaceVertexIter fv_it = trimesh_.fv_iter(fh); fv_it.is_valid(); ++fv_it) {
        int appearance = vertexAppearanceCG[*fv_it];
        //map U coordinates to vertices
        col = mapLocCoordToGlobCoordSys(fh, *fv_it);
//        std::cout << "rhs calc, vertex: " << fv_it->idx() << " and face: " << fh.idx() << ") with col: " << col
//                  << " and app: " << appearance << std::endl;
        //get sum of ut*dkm
        double sum = 0;
        for (int j = 0; j < 3; ++j) {
//            std::cout << "sum += " << CrossFieldAxis[j] << " * " << _D(j, col) << "\t" << CrossFieldAxis[j] * _D(j, col)
//                      << std::endl;
            sum += CrossFieldAxis[j] * basisTransformationMtx[fh](j, col);
        }
        //check if U or V position is needed
        int pos = (flagUorV) ? vertexPosUi[*fv_it] : vertexPosVi[*fv_it];
//        std::cout << "position of vertex " << fv_it->idx() << " is " << pos << std::endl;
        for (int i = 0; i < appearance; ++i) {
            double totSum = 2 * weight * area * h * sum;
//            std::cout << "total sum: " << 2 * weight * area * h * sum << std::endl;
            _rhs[pos + i] += totSum;
        }
    }
}

//map U coordinates to vertices
int
GlobalParametrization::mapLocCoordToGlobCoordSys(const OpenMesh::FaceHandle fh, const OpenMesh::VertexHandle vh) {
    auto localVUi = OpenMesh::FProp<std::vector<int>>(trimesh_, "localVUi");
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
GlobalParametrization::getHessian(const int rhsSizePartOne, const int rhsSizePartTwo) {
    auto faceSel = OpenMesh::FProp<bool>(trimesh_, "faceSel");
    trimesh_.request_halfedge_status();
    trimesh_.request_face_status();
    int size = rhsSizePartOne + rhsSizePartTwo;
    CMatrixType _H(size, size);
    for (auto he: trimesh_.halfedges()) {
        trimesh_.status(he).set_tagged(false);
        if (he.face().is_valid()) {
            trimesh_.status(he.face()).set_tagged(false);
        }
    }
    for (auto he: trimesh_.halfedges()) {
        if (!trimesh_.is_boundary(he) && faceSel[he.face()]) {
            if (!he.face().tagged()) {
                getDiaEntriesHessian(he.face(), _H);
            }
            if (!he.tagged()) {
                getEntriesHessian(he, _H);
            }
        }
    }
    trimesh_.release_halfedge_status();
    trimesh_.release_face_status();
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
        int occurrence = vertexAppearanceCG[*fv_it];
        gmm::row_matrix<std::vector<double>> D = basisTransformationMtx[fh];
        double sum = 0, col = mapLocCoordToGlobCoordSys(fh, *fv_it);
//        std::cout << "face: " << fh.idx() << "\nvertex: " << fv_it->idx() << "\ncol: " << col << "\nsum: ";
        for (int j = 0; j < 3; ++j) {
            sum += D(j, col) * D(j, col);
            //            std::cout << sum << "\t";
        }
//        std::cout << std::endl;
        for (int i = 0; i < occurrence; ++i) {
            _H(vertexPosUi[*fv_it] + i, vertexPosUi[*fv_it] + i) += 2 * weight * area * pow(h, 2) * sum;
            _H(vertexPosVi[*fv_it] + i, vertexPosVi[*fv_it] + i) += 2 * weight * area * pow(h, 2) * sum;
//            std::cout << "hessian U val (" << vertexPosUi[*fv_it] + i << "," << vertexPosUi[*fv_it] + i << ") "
//                      << _H(vertexPosUi[*fv_it] + i, vertexPosUi[*fv_it] + i) << " and V val " <<
//                      _H(vertexPosVi[*fv_it] + i, vertexPosVi[*fv_it] + i) << std::endl;
        }
        trimesh_.status(fh).set_tagged(true);
    }

}

void GlobalParametrization::getEntriesHessian(const OpenMesh::SmartHalfedgeHandle he, CMatrixType &_H) {
    auto faceSel = OpenMesh::FProp<bool>(trimesh_, "faceSel");
    auto vertexDist = OpenMesh::VProp<double>(trimesh_, "vertexDist");
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    auto basisTransformationMtx = OpenMesh::FProp<gmm::row_matrix<std::vector<double>>>(trimesh_,
                                                                                        "basisTransformationMtx");
    OpenMesh::FaceHandle adjFace = he.face();
//            std::cout << "Face " << voh_it->face().idx() << "\n start V " << vh_i.idx() << "\t end V " << voh_it->to()
//                      << std::endl;
    trimesh_.status(he).set_tagged(true);
    OpenMesh::VertexHandle vh_i = he.from();
    OpenMesh::VertexHandle vh_j = he.to();
    gmm::row_matrix<std::vector<double>> D = basisTransformationMtx[adjFace];
    double area = trimesh_.calc_face_area(adjFace), weight = 1, h = 1;
    double col1 = mapLocCoordToGlobCoordSys(adjFace, vh_i);
    double col2 = mapLocCoordToGlobCoordSys(adjFace, vh_j);
    //get sum of ut*dkm
    double sum = 0;
//            std::cout << "area " << area << " col 1: " << col1 << "\tcol 2: " << col2 << "\nsum: ";
    for (int j = 0; j < 3; ++j) {
        sum += D(j, col1) * D(j, col2);
//                std::cout << sum << "\t";
    }
//            std::cout << "\n";
    // u position
    _H(vertexPosUi[vh_i], vertexPosUi[vh_j]) += 2 * weight * area * pow(h, 2) * sum;
    _H(vertexPosUi[vh_j], vertexPosUi[vh_i]) += 2 * weight * area * pow(h, 2) * sum;
    // v position
    _H(vertexPosVi[vh_i], vertexPosVi[vh_j]) += 2 * weight * area * pow(h, 2) * sum;
    _H(vertexPosVi[vh_j], vertexPosVi[vh_i]) += 2 * weight * area * pow(h, 2) * sum;
//            std::cout << "hessian U val (" << vertexPosUi[vh_i] << "," << vertexPosUi[vh_j] << ") "
//                      << _H(vertexPosUi[vh_i], vertexPosUi[vh_j]) << std::endl;
}

void GlobalParametrization::createSectorsCutGraph(std::vector<int> &singularities) {
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
    int row_size = getRowSizeAndStatus();
    int col_size = nbVerticesUaV + cutGraphWoBoundary.size(); //j & k values + vertex U and V positions
    int singularity = 2, jkStartCounter = nbVerticesUaV;
    gmm::row_matrix<gmm::wsvector<double>> _constraints(row_size + singularity, col_size + 1);
    setZeroPoint(singularities, _constraints);
    getConstraintsMatrix(jkStartCounter, cutGraphWoBoundary, _constraints);
    trimesh_.release_vertex_status();
    return _constraints;
}

int GlobalParametrization::getRowSizeAndStatus() {
    int row_size = 0;
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    for (auto vh: trimesh_.vertices()) {
        trimesh_.status(vh).set_tagged(true);
        if (vertexAppearanceCG[vh] > 1) {
            row_size += vertexAppearanceCG[vh] * 2; // u and v value
            trimesh_.status(vh).set_tagged(false);
        }
    }
    return row_size;
}

void GlobalParametrization::setZeroPoint(std::vector<int> &singularities,
                                         gmm::row_matrix<gmm::wsvector<double>> &_constraints) {
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    // set point (0,0) of coord system from a leaf singularity
    for (auto i: singularities) {
        auto vh = trimesh_.vertex_handle(i);
        if (checkIfLeaf(vh)) {
            _constraints(0, vertexPosUi[vh]) = 1;
            _constraints(1, vertexPosVi[vh]) = 1;
            break;
        }
    }
    if (singularities.empty()) {
        for (auto vh: trimesh_.vertices()) {
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
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    int counter = 2; // because singularity is at pos 0&1
    for (auto it = cutGraphWoBoundary.begin(); it != cutGraphWoBoundary.end(); it++) {
        auto he1 = make_smart(trimesh_.halfedge_handle(*it++), trimesh_);
        auto he2 = make_smart(trimesh_.halfedge_handle(*it), trimesh_);
        try {
            if (he1.opp().idx() != he2.idx()) {
                throw 404;
            }
        } catch (int x) {
            std::cerr << "getConstraintMatrix: halfedges have to be opposite, something went wrong\n";
        }
        int diff = (currentPJ[he1.face()] + currentPJ[he2.face()]) % 4;
        if (!he1.to().tagged()) {
            setConRows(counter, jkStartCounter, diff, he1, _constraints);
        }
        if (!he2.to().tagged()) {
            setConRows(counter, jkStartCounter, diff, he2, _constraints);
        }
        jkStartCounter += 2;
    }
}

void
GlobalParametrization::setConRows(int &counter, int &jkStartCounter, const int diff,
                                  OpenMesh::SmartHalfedgeHandle &he,
                                  gmm::row_matrix<gmm::wsvector<double>> &_constraints) {
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");
    auto vh = he.to();
    int posOne = getPositionConstraintRow(vh, cutGraphFZone[he.face()]);
    int posTwo = getPositionConstraintRow(vh, cutGraphFZone[he.opp().face()]);
    switch (diff) {
        case 0: //0
            // u position of point p
            _constraints(counter, vertexPosUi[vh] + posOne) = -1; //u_p
            _constraints(counter, vertexPosUi[vh] + posTwo) = 1; // u'_p
            _constraints(counter++, jkStartCounter) = -1; //j
            // v position of point p
            _constraints(counter, vertexPosVi[vh] + posOne) = -1;
            _constraints(counter, vertexPosVi[vh] + posTwo) = 1;
            _constraints(counter++, jkStartCounter + 1) = -1;
            break;
        case -3://-270 || 90
        case 1:
            // u position of point p
            _constraints(counter, vertexPosVi[vh] + posOne) = 1; //-v_p
            _constraints(counter, vertexPosUi[vh] + posTwo) = 1; //u'_p
            _constraints(counter++, jkStartCounter) = -1; //j
            // v position of point p
            _constraints(counter, vertexPosUi[vh] + posOne) = -1;
            _constraints(counter, vertexPosVi[vh] + posTwo) = 1;
            _constraints(counter++, jkStartCounter + 1) = -1;
            break;
        case -2://-180 || 180
        case 2:
            // u position of point p
            _constraints(counter, vertexPosUi[vh] + posOne) = 1; //-u_p
            _constraints(counter, vertexPosUi[vh] + posTwo) = 1; // u'_p
            _constraints(counter++, jkStartCounter) = -1; //j
            // v position of point p
            _constraints(counter, vertexPosVi[vh] + posOne) = 1;
            _constraints(counter, vertexPosVi[vh] + posTwo) = 1;
            _constraints(counter++, jkStartCounter + 1) = -1;
            break;
        case -1://-90 || 270
        case 3:
            // u position of point p
            _constraints(counter, vertexPosVi[vh] + posOne) = -1; //v_p
            _constraints(counter, vertexPosUi[vh] + posTwo) = 1; //u'_p
            _constraints(counter++, jkStartCounter) = -1; //j
            // v position of point p
            _constraints(counter, vertexPosUi[vh] + posOne) = 1;
            _constraints(counter, vertexPosVi[vh] + posTwo) = 1;
            _constraints(counter++, jkStartCounter + 1) = -1;
            break;
    }
}

int GlobalParametrization::getPositionConstraintRow(OpenMesh::SmartVertexHandle &vh, int cutGraphZone) {
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    std::vector<int> sectors;
    for (auto it = vh.outgoing_halfedges().begin(); it.is_valid(); ++it) {
        int zone = cutGraphFZone[it->face()];
        if (zone != 0 && (std::find(sectors.begin(), sectors.end(), zone) == sectors.end())) {
            sectors.push_back(zone);
        }
    }
    try {
        if ((int) sectors.size() != vertexAppearanceCG[vh]) {
            throw 404;
        }
    } catch (int x) {
        std::cerr << "getPositionConstraintRow: vector doesn't have the same size as vertexAppearance property v: "
                  << sectors.size()
                  << " != p: " << vertexAppearanceCG[vh] << "\n";
    }
    std::sort(sectors.begin(), sectors.end());
    auto it1 = find(sectors.begin(), sectors.end(), cutGraphZone);
    int pos = it1 - sectors.begin();
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
        }
        if (!he.from().tagged()) {
            _idx_to_round.push_back(vertexPosUi[he.from()]);
            _idx_to_round.push_back(vertexPosVi[he.from()]);
            trimesh_.status(he.from()).set_tagged(true);
        }

    }
    //feature edge constraint
    //not yet implemented

    //singularity constraint
    for (auto i: singularities) {
        auto vh = make_smart(trimesh_.vertex_handle(i), trimesh_);
        if (!vh.tagged()) {
            _idx_to_round.push_back(vertexPosUi[vh]);
            _idx_to_round.push_back(vertexPosVi[vh]);
            trimesh_.status(vh).set_tagged(true);
        }
    }
    //jk values are integers as well
    for (int i = nbVerticesUaV; i < nbVerticesUaV + jkValues; ++i) {
        _idx_to_round.push_back(i);
    }
    trimesh_.release_vertex_status();
    return _idx_to_round;
}

void GlobalParametrization::saveSolAsCoord(std::vector<double> &_x, std::vector<int> &faces,
                                           std::vector<int> &singularities, std::vector<int> &cutGraphWoBoundary,
                                           int nbVerticesUaV, int jkValues) {
    trimesh_.request_vertex_status();
    auto solCoordSysUV = OpenMesh::VProp<std::vector<Point>>(trimesh_, "solCoordSysUV");
    //init status
    initPropForSolVector();

    // this is how i get j,k values. check getConstraintMatrix function in order to deduce ordering
    for (auto it = cutGraphWoBoundary.begin(); it != cutGraphWoBoundary.end(); ++it) {
        auto he = make_smart(trimesh_.halfedge_handle(*it), trimesh_);
        getSolFromVerticesWMoreOneApp(he, _x, jkValues);
    }

    // then we can run through the faces and ignore the tagged vertices
    for (auto f: faces) {
        auto fh = trimesh_.face_handle(f);
        getSolFromVerticesWOneApp(fh, _x);
    }
    //print solution
//    for (auto vh: trimesh_.vertices()) {
//        trimesh_.status(vh).set_tagged(false);
//    }
//    std::ofstream solutionVector("/home/wuschelbueb/Desktop/solutionVector.txt");
//    for (auto f: faces) {
//        auto fh = trimesh_.face_handle(f);
//        for (TriMesh::FaceVertexIter fv_it = trimesh_.fv_iter(fh); fv_it.is_valid(); ++fv_it) {
//            if (!fv_it->tagged()) {
//                auto he1 = make_smart(trimesh_.halfedge_handle(*fv_it), trimesh_);
//                auto vh = he1.to();
//                solutionVector << "vertex: " << vh.idx() << std::endl;
//                for (size_t i = 0; i < solCoordSysUV[vh].size(); ++i) {
//                    solutionVector << "Position " << i << " with point " << solCoordSysUV[vh][i] << std::endl;
//                }
//                trimesh_.status(*fv_it).set_tagged(true);
//            }
//
//        }
//    }
//    solutionVector.close();
    trimesh_.release_vertex_status();
}

void GlobalParametrization::getSolFromVerticesWMoreOneApp(OpenMesh::SmartHalfedgeHandle he, std::vector<double> &_x,
                                                          int &jkValues) {
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");
    auto solCoordSysUV = OpenMesh::VProp<std::vector<Point>>(trimesh_, "solCoordSysUV");
    if (!he.to().tagged2()) {
        // should more or less look like that.
        auto vh = he.to();
        int posTwo = getPositionConstraintRow(vh, cutGraphFZone[he.opp().face()]);
        int posU = vertexPosUi[vh] + posTwo;
        int posV = vertexPosVi[vh] + posTwo;
        double u = _x[posU] + _x[jkValues++];
        double v = _x[posV] + _x[jkValues++];
        Point p = {u, v, 0};
        solCoordSysUV[vh].push_back(p);
    }
}

void GlobalParametrization::initPropForSolVector() {
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    auto solCoordSysUV = OpenMesh::VProp<std::vector<Point>>(trimesh_, "solCoordSysUV");
    for (auto vh: trimesh_.vertices()) {
        trimesh_.status(vh).set_tagged(false);
        // needed to handle cutGraphWoBoundary
        trimesh_.status(vh).set_tagged2(true);
        if (vertexAppearanceCG[vh] > 1) {
            trimesh_.status(vh).set_tagged2(false);
            trimesh_.status(vh).set_tagged(true);
        }
        solCoordSysUV[vh];
    }
}

void GlobalParametrization::getSolFromVerticesWOneApp(OpenMesh::FaceHandle fh, std::vector<double> &_x) {
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    auto solCoordSysUV = OpenMesh::VProp<std::vector<Point>>(trimesh_, "solCoordSysUV");
    for (TriMesh::FaceVertexIter fv_it = trimesh_.fv_iter(fh); fv_it.is_valid(); ++fv_it) {
        if (!fv_it->tagged()) {
            int posU = vertexPosUi[*fv_it];
            int posV = vertexPosVi[*fv_it];
            double u = _x[posU];
            double v = _x[posV];
            Point p = {u, v, 0};
            solCoordSysUV[*fv_it].push_back(p);
            trimesh_.status(*fv_it).set_tagged(true);
        }
    }
}


