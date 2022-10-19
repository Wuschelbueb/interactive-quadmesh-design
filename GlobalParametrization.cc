//
// Created by wuschelbueb on 05.05.22.
//

#include "GlobalParametrization.hh"

void GlobalParametrization::getGlobalParam() {
    //set up vectors
    std::vector<int> faces = getFaceVec();
    std::vector<int> singularities = getSingularities(faces);
    setUpLocFaceCoordSys(faces);
    DijkstraDistance dualGraph(trimesh_);
    dualGraph.getDualGraph(faces);
    std::vector<int> complementHEdges = getComplementMeshSel();
    removeOpenPaths(complementHEdges);
    removeRedundantEdges(complementHEdges);
    std::vector<int> onlyBoundaries = complementHEdges;
    std::vector<int> cutGraphWoBoundary;
//    add path from boundary to singularity to cut graph
    dualGraph.calcDijkstraWSingularities(complementHEdges, singularities, cutGraphWoBoundary);
//    to visualize cutGraph with or without edges
    colorCompBoundaries(onlyBoundaries);
    colorCompHEdges(complementHEdges);
    createSectorsOnCutGraph(singularities);
    fixRotationsCrossBoundaryComp(faces);
    int nbVerticesUaV = createVertexPosParamDomain(faces);
    int jkValues = (int) cutGraphWoBoundary.size();
    std::cout << "nbVerticesUaV: " << nbVerticesUaV << " and jkValues: " << jkValues << std::endl;
    std::vector<double> _x(nbVerticesUaV + jkValues, 0.0);
    std::vector<double> _rhs = getRhs(nbVerticesUaV, jkValues);
    std::vector<int> _idx_to_round = getIdxToRound(nbVerticesUaV, jkValues, singularities, onlyBoundaries);
    CMatrixType _hessian = getHessian(nbVerticesUaV, jkValues, cutGraphWoBoundary, faces);
    RMatrixType _constraints = getConstraints(nbVerticesUaV, cutGraphWoBoundary, onlyBoundaries, singularities,
                                              faces);
//    std::ofstream hessian("/home/wuschelbueb/Desktop/data/hessian_matrix.txt");
//    hessian << _hessian;
//    hessian.close();
//    std::ofstream rhs("/home/wuschelbueb/Desktop/data/rhs.txt");
//    for (auto i: _rhs) {
//        rhs << i << std::endl;
//    }
//    rhs.close();
//    std::ofstream idx_to_round("/home/wuschelbueb/Desktop/data/idx_to_round.txt");
//    for (auto i: _idx_to_round) {
//        idx_to_round << i << std::endl;
//    }
//    idx_to_round.close();
//    std::ofstream cons("/home/wuschelbueb/Desktop/data/constraints.txt");
//    cons << _constraints << std::endl;
//    cons.close();
//    std::ofstream vertices("/home/wuschelbueb/Desktop/data/vertices.txt");
//    for (auto vh: trimesh_.vertices()) {
//        vertices << "vertex " << vh.idx() << "= [" << trimesh_.calc_centroid(vh) << "]" << std::endl;
//    }
//    vertices.close();
    std::cout << "Calculation of Global Parametrization\n" << "Dimensions of parameters:\n" << "Constraints Rows:\t"
              << gmm::mat_nrows(_constraints)
              << " and columns: " << gmm::mat_ncols(_constraints) << std::endl
              << "A Rows:\t\t\t\t" << gmm::mat_nrows(_hessian) << " and columns: " << gmm::mat_ncols(_hessian)
              << std::endl
              << "Size of _x:\t\t\t" << _x.size() << std::endl << "Size of _rhs:\t\t" << _rhs.size() << std::endl
              << "Size of idx:\t\t" << _idx_to_round.size() << std::endl;
    COMISO::ConstrainedSolver csolver;
    csolver.misolver().set_iter_full(false);
    csolver.misolver().set_local_iters(0);
    csolver.misolver().set_cg_iters(20);
    csolver.solve(_constraints, _hessian, _x, _rhs, _idx_to_round);
    saveSolAsCoord(_x, faces);
}

std::vector<int> GlobalParametrization::getSingularities(std::vector<int> &faces) {
    auto crossFieldIdx = OpenMesh::VProp<double>(trimesh_, "crossFieldIdx");
    std::vector<int> singularities;
    for (const auto &i: faces) {
        auto fh = trimesh_.face_handle(i);
        for (auto fv_it = trimesh_.fv_iter(fh); fv_it.is_valid(); ++fv_it) {
            // try with != 0
            if (crossFieldIdx[*fv_it] < -1E-1 || crossFieldIdx[*fv_it] > 1E-1) {
                singularities.push_back(fv_it->idx());
            }
        }
    }
//    for (auto vh: trimesh_.vertices()) {
//        // try with != 0
//        if (crossFieldIdx[vh] < -1E-1 || crossFieldIdx[vh] > 1E-1) {
//            singularities.push_back(vh.idx());
//        }
//    }
    return singularities;
}

void GlobalParametrization::removeOpenPaths(std::vector<int> &complementHEdges) {
    trimesh_.request_halfedge_status();
    int currentIter = complementHEdges.size(), prevIter = INT_MAX;
    while (currentIter != prevIter) {
        for (const int &i: complementHEdges) {
            removeEdgeFromGraph(i, complementHEdges);
        }
        prevIter = currentIter;
        currentIter = complementHEdges.size();
    }
    trimesh_.release_halfedge_status();
}

void GlobalParametrization::removeEdgeFromGraph(const int i, std::vector<int> &complementHEdges) {
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
    for (const auto &eh: trimesh_.edges()) {
        auto heh = make_smart(trimesh_.halfedge_handle(eh, 1), trimesh_);
        auto oheh = make_smart(trimesh_.halfedge_handle(eh, 0), trimesh_);
        auto iter1 = std::find(complementHEdges.begin(), complementHEdges.end(), heh.idx());
        auto iter2 = std::find(complementHEdges.begin(), complementHEdges.end(), oheh.idx());
        // 1 = border of selection, 2 = dual graph, 3 = non-dual graph
        bool checkHeCompForIter1 = iter1 != complementHEdges.end();
        bool checkHeCompForIter2 = iter2 != complementHEdges.end();
        bool borderFalse = borderDualG[eh] != 1;
        bool borderTrue = borderDualG[eh] == 1;
        bool checkFaceSelection = heh.face().is_valid() && !faceSel[heh.face()];
        bool checkOppFaceSelection = oheh.face().is_valid() && !faceSel[oheh.face()];
        bool heIsBoundary = heh.is_boundary();
        bool oheIsBoundary = oheh.is_boundary();
        if (borderFalse) {
            if (checkHeCompForIter1 || checkHeCompForIter2) {
                complementHEdges.erase(iter1);
                complementHEdges.erase(iter2);
            }
        } else if (borderTrue) {
            if (checkHeCompForIter1) {
                if (checkFaceSelection || heIsBoundary) {
                    complementHEdges.erase(iter1);
                }
            }
            if (checkHeCompForIter2) {
                if (checkOppFaceSelection || oheIsBoundary) {
                    complementHEdges.erase(iter2);
                }
            }
        }
    }
}

void GlobalParametrization::colorCompHEdges(const std::vector<int> &complementEdges) {
    auto cutGraphHe = OpenMesh::HProp<bool>(trimesh_, "cutGraphHe");
    for (const auto &heh: trimesh_.halfedges()) {
        cutGraphHe[heh] = false;
    }
    for (const int &i: complementEdges) {
        OpenMesh::HalfedgeHandle he = trimesh_.halfedge_handle(i);
        cutGraphHe[he] = true;
    }
}

void GlobalParametrization::colorCompBoundaries(std::vector<int> &onlyBoundaries) {
    auto boundaryHe = OpenMesh::HProp<bool>(trimesh_, "boundaryHe");
    for (const auto &heh: trimesh_.halfedges()) {
        boundaryHe[heh] = false;
    }
    for (const int &i: onlyBoundaries) {
        OpenMesh::HalfedgeHandle he = trimesh_.halfedge_handle(i);
        boundaryHe[he] = true;
    }
}

std::vector<int> GlobalParametrization::getComplementMeshSel() {
    trimesh_.request_halfedge_status();
    for (auto heh: trimesh_.halfedges()) {
        trimesh_.status(heh).set_tagged(false);
    }
    std::vector<int> complementHEdges;
    tagEdgesFromDualSpanningTree();
    for (const auto &heh: trimesh_.halfedges()) {
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
    for (const int &i: faces) {
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
    gmm::col_matrix<std::vector<double>> _abc = createMatrix3v3();
    std::vector<Point> edges(3);
    for (const int &i: faces) {
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

gmm::col_matrix<std::vector<double>> GlobalParametrization::createMatrix3v3() {
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
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");
    auto uVectorFieldRotOne = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotOne");
    auto uVectorFieldRotTwo = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotTwo");
    int size = rhsSizePartOne + rhsSizePartTwo;
    std::vector<double> _rhs(size, 0);
//    std::ofstream vectorFieldRot("/home/wuschelbueb/Desktop/data/vectorFieldRot.txt");
    // use nbVerticesUaV and lhs formulas with dkm -> dkm is entry in D matrix
    // utk is entry of vectorRotOne
    for (const auto &he: trimesh_.halfedges()) {
        if (!trimesh_.is_boundary(he) && faceSel[he.face()]) {
            Point uCrossField = uVectorFieldRotOne[he.face()];
            Point vCrossField = uVectorFieldRotTwo[he.face()];
//            vectorFieldRot << "Face " << he.face().idx() << " = [[" << uCrossField << "],[" << vCrossField << "]]\n";
            getRhsEntryForVertex(he, uCrossField, true, _rhs);
            getRhsEntryForVertex(he, vCrossField, false, _rhs);
        }
    }
    //set entries smaller than 1E-10 to ze ro
    for (auto &i: _rhs) {
        if (i < 1E-10 && i > -1E-10) {
            i = 0;
        }
    }
//    vectorFieldRot.close();
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
    double area = trimesh_.calc_face_area(adjFace), col = mapLocCoordToGlobCoordSys(adjFace, vh);

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
    _rhs[UVPos + posAppearance] += 2 * area * hVal * sum;
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
                                  std::vector<int> &cutGraphWoBoundary,
                                  const std::vector<int> &faces) {
    auto faceSel = OpenMesh::FProp<bool>(trimesh_, "faceSel");
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");
    trimesh_.request_halfedge_status();
    trimesh_.request_face_status();
    int size = rhsSizePartOne + rhsSizePartTwo;
    CMatrixType _hessian(size, size);
    //initialize halfedges and faces
    for (const int &i: faces) {
        auto fh = trimesh_.face_handle(i);
        for (auto fh_it = trimesh_.fh_iter(fh); fh_it.is_valid(); ++fh_it) {
            if (!trimesh_.is_boundary(*fh_it)) {
                getDiaEntriesHessian(*fh_it, _hessian);
                getEntriesHessian(*fh_it, _hessian);
            }
        }

    }
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
    double area = trimesh_.calc_face_area(adjFace), weight = 1;
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
    _hessian(vertexPosUi[vh] + posOne, vertexPosUi[vh] + posOne) += 2 * weight * area * pow(hVal, 2) * sum;
    _hessian(vertexPosVi[vh] + posOne, vertexPosVi[vh] + posOne) += 2 * weight * area * pow(hVal, 2) * sum;
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
    double area = trimesh_.calc_face_area(adjFace), weight = 1;
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
    _hessian(vertexPosUi[vh_i] + posOne, vertexPosUi[vh_j] + posTwo) += 2 * weight * area * pow(hVal, 2) * sum;
    _hessian(vertexPosUi[vh_j] + posTwo, vertexPosUi[vh_i] + posOne) += 2 * weight * area * pow(hVal, 2) * sum;

    // v position
    _hessian(vertexPosVi[vh_i] + posOne, vertexPosVi[vh_j] + posTwo) += 2 * weight * area * pow(hVal, 2) * sum;
    _hessian(vertexPosVi[vh_j] + posTwo, vertexPosVi[vh_i] + posOne) += 2 * weight * area * pow(hVal, 2) * sum;
//    std::cout << "hessian U val (" << vertexPosUi[vh_i] + posOne << "," << vertexPosUi[vh_j] + posTwo << ") "
//              << _hessian(vertexPosUi[vh_i] + posOne, vertexPosUi[vh_j] + posTwo) << std::endl;
}

void GlobalParametrization::createSectorsOnCutGraph(std::vector<int> &singularities) {
    auto cutGraphHe = OpenMesh::HProp<bool>(trimesh_, "cutGraphHe");
    bool singTest = false;
    OpenMesh::SmartHalfedgeHandle ohe;
    for (const int &singularity: singularities) {
        auto vh = trimesh_.vertex_handle(singularity);
        if (checkIfLeaf(vh)) {
            for (TriMesh::VertexOHalfedgeIter voh_it = trimesh_.voh_iter(vh); voh_it.is_valid(); ++voh_it) {
                if (cutGraphHe[*voh_it] && cutGraphHe[voh_it->opp()]) {
                    singTest = true;
                    ohe = *voh_it;
                    break;
                }
            }
            if (singTest) {
                propagation(singularity, ohe);
                break;
            }
        }
    }
}

void GlobalParametrization::propagation(const int &singularity, OpenMesh::SmartHalfedgeHandle he) {
    auto cutGraphHe = OpenMesh::HProp<bool>(trimesh_, "cutGraphHe");
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");
    auto vh = trimesh_.vertex_handle(singularity);
    int color = 1;
    cutGraphFZone[he.face()] = color;
    OpenMesh::SmartHalfedgeHandle prevCutHe = he, currentCutHe = he, firstHe = he;
    bool done = false, specialCase = false;
    while (!(done && specialCase)) {
        he = he.next();
        cutGraphFZone[he.face()] = color;
        // if new he is on cutgraph continue this path
        if (cutGraphHe[he]) {
            prevCutHe = currentCutHe;
            currentCutHe = he;
            // special case 1:
            // increase color if path between singularity and singularity/border is only one edge long
            if (prevCutHe.opp() == currentCutHe) {
                color++;
            }
//            std::cout << "he is " << he.idx() << " part of Cut\n"
//                      << "prevCut " << prevCutHe << "\ncurrentCut "
//                      << currentCutHe << std::endl;
            // else check go to opposite he, in order to cycle through
        } else if (!he.opp().is_boundary()) {
            he = he.opp();
//            std::cout << "he is " << he.idx() << " bound" << std::endl;
        }
//        std::cout << "new he is " << he.idx() << std::endl;
        if (firstHe == he) {
            cutGraphFZone[he.face()] = 1;
            cutGraphFZone[he.prev().opp().face()] = 1;
//            std::cout << "reached end\n";
            done = true;
        }
        if (done && color == 2 && cutGraphHe[prevCutHe.opp()] && !cutGraphHe[currentCutHe.opp()]) {
            cutGraphFZone[he.next().opp().face()] = color;
            specialCase = true;
        }
        if (done && color > 2) {
            specialCase = true;
        }
    }
}


bool GlobalParametrization::checkIfLeaf(const OpenMesh::VertexHandle &heToVertex) {
    auto cutGraphHe = OpenMesh::HProp<bool>(trimesh_, "cutGraphHe");
    int counter = 0;
    for (TriMesh::VertexOHalfedgeIter voh_it = trimesh_.voh_iter(heToVertex); voh_it.is_valid(); ++voh_it) {
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
    int pjVal = 0;
    if (pj.first == ohe.opp().idx()) {
        pjVal = -pj.second + currentPJ[ohe.opp().face()];
    } else {
        pjVal = pj.second + currentPJ[ohe.opp().face()];
    }
    uVectorFieldRotOne[ohe.face()] =
            rotPointWithRotMatrix(ohe.face(), uVectorFieldRotOne[ohe.face()], M_PI / 2 * pjVal);
    uVectorFieldRotTwo[ohe.face()] =
            rotPointWithRotMatrix(ohe.face(), uVectorFieldRotTwo[ohe.face()], M_PI / 2 * pjVal);
    currentPJ[ohe.face()] = pjVal;
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
                                      std::vector<int> &onlyBoundaries, std::vector<int> &singularities,
                                      std::vector<int> &faces) {
    // nd V positions + j & k values
    int singularity = 2, jkStartCounter = nbVerticesUaV;
    int partTwo = onlyBoundaries.size(), partOne = cutGraphWoBoundary.size() * 2 + singularity;
    int row_size = partOne + partTwo;
    int col_size = nbVerticesUaV + cutGraphWoBoundary.size();
    gmm::row_matrix<gmm::wsvector<double>> _constraints(row_size, col_size + 1);
    setZeroPointConstraint(singularities, _constraints, faces);
    getConstraintsMatrix(jkStartCounter, cutGraphWoBoundary, _constraints);
//    setFeatureLineConstraint(_constraints, onlyBoundaries, partOne);
    return _constraints;
}

void GlobalParametrization::setZeroPointConstraint(std::vector<int> &singularities,
                                                   gmm::row_matrix<gmm::wsvector<double>> &_constraints,
                                                   std::vector<int> &faces) {
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    // set point (0,0) of coord system from a leaf singularity
    for (const auto &i: singularities) {
        auto vh = trimesh_.vertex_handle(i);
        if (checkIfLeaf(vh)) {
            _constraints(0, vertexPosUi[vh]) = 1;
            _constraints(1, vertexPosVi[vh]) = 1;
            break;
        }
    }
    if (singularities.empty()) {
        auto fh = trimesh_.face_handle(faces[0]);
        auto fv_it = trimesh_.fv_iter(fh);
        _constraints(0, vertexPosUi[*fv_it]) = 1;
        _constraints(1, vertexPosVi[*fv_it]) = 1;
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
    auto periodJump = OpenMesh::HProp<int>(trimesh_, "periodJump");
    int counter = 2; // because singularity constraint is at pos 0&1
//    std::ofstream abcd("/home/wuschelbueb/data/Desktop/abcd.txt");
//    abcd.close();

    //jkstartcounter == nbOfVerticesUaV
    for (auto it = cutGraphWoBoundary.begin(); it != cutGraphWoBoundary.end(); it++) {
        auto he1 = make_smart(trimesh_.halfedge_handle(*it++), trimesh_);
        auto he2 = make_smart(trimesh_.halfedge_handle(*it), trimesh_);
//        abcd << "halfedge " << he1.idx() << " with vertex " << he1.to().idx() << "\nhalfedge " << he2.idx()
//             << " with vertex " << he2.to().idx() << "\n";
        try {
            if (he1.opp().idx() != he2.idx()) {
                throw 404;
            }
        } catch (int x) {
            std::cerr << "getConstraintMatrix: halfedges have to be opposite, something went wrong\n";
        }
        int pj = 0;
        if ((periodJump[he1] != INT_MAX && periodJump[he1] != 0)) {
            pj = periodJump[he1];
        } else if ((periodJump[he2] != INT_MAX && periodJump[he2] != 0)) {
            pj = -periodJump[he2];
        }
        auto vh1 = he1.to();
        auto vh2 = he1.from();
        // diff = currentPJ difference + pj over cutgraph edge
        int diff = (currentPJ[he2.face()] - currentPJ[he1.face()] + pj) % 4;
        if ((periodJump[he1] != INT_MAX && periodJump[he1] != 0)) {
            periodJump[he1] = diff;
        } else if ((periodJump[he2] != INT_MAX && periodJump[he2] != 0)) {
            periodJump[he2] = diff;
        }

        setConRows(counter, jkStartCounter, diff, he1, vh1, _constraints);
        setConRows(counter, jkStartCounter, diff, he1, vh2, _constraints);
        jkStartCounter += 2;
    }
//    abcd.close();
}

void
GlobalParametrization::setConRows(int &counter, int &jkStartCounter, const int diff, OpenMesh::SmartHalfedgeHandle &he,
                                  OpenMesh::SmartVertexHandle &vh,
                                  gmm::row_matrix<gmm::wsvector<double>> &_constraints) {
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");
    bool leafcheck = checkIfLeaf(vh);
    int posOne = 0, posTwo = 0;
    if (!leafcheck) {
        posOne = getPositionConstraintRow(vh, cutGraphFZone[he.face()]);
        posTwo = getPositionConstraintRow(vh, cutGraphFZone[he.opp().face()]);
    }

    switch (diff) {
        case 0: //0
            if (!leafcheck) {
                // u position of point p
                _constraints(counter, vertexPosUi[vh] + posTwo) += -1; //u_p
                _constraints(counter, vertexPosUi[vh] + posOne) += 1; // u'_p
                _constraints(counter++, jkStartCounter) = -1; //j
                // v position of point p
                _constraints(counter, vertexPosVi[vh] + posTwo) += -1;
                _constraints(counter, vertexPosVi[vh] + posOne) += 1;
                _constraints(counter++, jkStartCounter + 1) = -1; //k
            } else {
                _constraints(counter++, jkStartCounter) = -1; //j
                _constraints(counter++, jkStartCounter + 1) = -1; //k
            }
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
            if (!leafcheck) {
                // u position of point p
                _constraints(counter, vertexPosUi[vh] + posTwo) = 1; //-u_p
                _constraints(counter, vertexPosUi[vh] + posOne) = 1; // u'_p
                _constraints(counter++, jkStartCounter) = -1; //j
                // v position of point p
                _constraints(counter, vertexPosVi[vh] + posTwo) = 1;
                _constraints(counter, vertexPosVi[vh] + posOne) = 1;
                _constraints(counter++, jkStartCounter + 1) = -1;
            } else {
                _constraints(counter, vertexPosUi[vh] + posTwo) = 2; //-u_p
                _constraints(counter++, jkStartCounter) = -1; //j
                _constraints(counter, vertexPosVi[vh] + posOne) = 2;
                _constraints(counter++, jkStartCounter + 1) = -1; //k
            }
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
//    std::string leaf = (leafcheck) ? " leaf" : " nonLeaf";
//    abcd << "\tvertex " << vh.idx() << leaf << " origin from he " << he.idx() << " (adj. Face "
//         << he.face().idx() << ")\n\t\tposition U / V: " << vertexPosUi[vh] + posOne << "_(" << vertexPosUi[vh] << "+"
//         << posOne << ") / " << vertexPosVi[vh] + posOne << "_(" << vertexPosVi[vh] << "+"
//         << posOne << ")\n\tvertex " << vh.idx() << leaf << " origin from he " << he.opp().idx() << " (adj. Face "
//         << he.opp().face().idx() << ")\n\t\tposition U / V: "
//         << vertexPosUi[vh] + posTwo << "_(" << vertexPosUi[vh] << "+"
//         << posTwo << ") / " << vertexPosVi[vh] + posTwo << "_(" << vertexPosVi[vh] << "+"
//         << posTwo << ")" << std::endl << "\tu row: " << _constraints[counter - 2] << std::endl << "\tv row: "
//         << _constraints[counter - 1] << std::endl << std::endl;
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
//    std::cout << "for vertex " << vh.idx() << " with cutGraphFZone " << cutGraphZone << " we get position " << pos
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
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    std::vector<int> _idx_to_round;

    //singularity constraint
    for (const auto &i: singularities) {
        auto vh = make_smart(trimesh_.vertex_handle(i), trimesh_);
        if (!vh.tagged()) {
            for (int j = 0; j < vertexAppearanceCG[vh]; ++j) {
                _idx_to_round.push_back(vertexPosUi[vh] + j);
                _idx_to_round.push_back(vertexPosVi[vh] + j);
            }
            trimesh_.status(vh).set_tagged(true);
        }
    }

    //boundary constraints
    for (const auto &i: onlyBoundaries) {
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
    //feature edge constraint; not yet implemented

    //jk values are integers as well
    for (int i = nbVerticesUaV; i < nbVerticesUaV + jkValues; ++i) {
        _idx_to_round.push_back(i);
    }
    trimesh_.release_vertex_status();
    return _idx_to_round;
}

void GlobalParametrization::saveSolAsCoord(std::vector<double> &_x, std::vector<int> &faces) {
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");
    auto solCoordSysUV = OpenMesh::VProp<std::vector<OpenMesh::Vec2d>>(trimesh_, "solCoordSysUV");
    auto faceSel = OpenMesh::FProp<bool>(trimesh_, "faceSel");
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    trimesh_.request_vertex_status();
    initPropForSolVector();
//    std::ofstream rest("/home/wuschelbueb/Desktop/data/rest.txt");
//    rest.close();

    for (const auto &i: faces) {
        auto fh = trimesh_.face_handle(i);
        for (auto fh_it = trimesh_.fh_iter(fh); fh_it.is_valid(); ++fh_it) {
            if (!trimesh_.is_boundary(*fh_it)) {
                saveSolToVertices(*fh_it, _x);
            }
        }
    }

    //solution printed to file
//    int nFaces = trimesh_.n_faces(), nFacesIter = 0;
//    std::ofstream xVector("/home/wuschelbueb/Desktop/data/xVectorGlobalParam.txt");
//
//    for (auto vh: trimesh_.vertices()) {
//        xVector << "vertex (" << solCoordSysUV[vh].size() << ") " << vh.idx() << " with coord:\n";
//        if (!solCoordSysUV[vh].empty()) {
//            for (size_t i = 0; i < solCoordSysUV[vh].size(); ++i) {
//                xVector << "x[" << i << "]: " << solCoordSysUV[vh][i][0] << "\ny[" << i << "]: "
//                        << solCoordSysUV[vh][i][1] << std::endl;
//            }
//        }
//    }
//
//    xVector << "\n[";
//    for (auto fh: trimesh_.faces()) {
//        if (faceSel[fh]) {
//            int counter = 0;
//            xVector << "[";
//            for (TriMesh::FaceHalfedgeIter fhe_it = trimesh_.fh_iter(fh); fhe_it.is_valid(); ++fhe_it) {
//                OpenMesh::Vec2d p;
//                if (vertexAppearanceCG[fhe_it->to()] > 1) {
//                    auto vh = fhe_it->to();
//                    int posOne = getPositionConstraintRow(vh, cutGraphFZone[fhe_it->face()]);
//                    p = solCoordSysUV[vh][posOne];
//                } else {
//                    p = solCoordSysUV[fhe_it->to()][0];
//                }
//                if (counter != 2) {
//                    xVector << "[" << p[0] << "," << p[1] << "],";
//                } else {
//                    xVector << "[" << p[0] << "," << p[1] << "]";
//                }
//                counter++;
//            }
//            if (nFacesIter != nFaces - 1) {
//                xVector << "],";
//            } else {
//                xVector << "]";
//            }
//        }
//        nFacesIter++;
//    }
//    xVector << "]";
//    xVector.close();
    trimesh_.release_vertex_status();
}

void GlobalParametrization::initPropForSolVector() {
    auto faceSel = OpenMesh::FProp<bool>(trimesh_, "faceSel");
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    auto solCoordSysUV = OpenMesh::VProp<std::vector<OpenMesh::Vec2d>>(trimesh_, "solCoordSysUV");
    for (const auto &he: trimesh_.halfedges()) {
        if (!trimesh_.is_boundary(he) && faceSel[he.face()]) {
            solCoordSysUV[he.to()] = std::vector<OpenMesh::Vec2d>(vertexAppearanceCG[he.to()], {0, 0});
            trimesh_.status(he.to()).set_tagged(false);
        } else if (!trimesh_.is_boundary(he) && !faceSel[he.face()]) {
            solCoordSysUV[he.to()] = std::vector<OpenMesh::Vec2d>(vertexAppearanceCG[he.to()], {0, 0});
            trimesh_.status(he.to()).set_tagged(true);
        }
    }
}

void GlobalParametrization::saveSolToVertices(OpenMesh::SmartHalfedgeHandle he, std::vector<double> &_x) {
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
        double u = _x[posU];
        double v = _x[posV];
        OpenMesh::Vec2d p = {u, v};
        solCoordSysUV[vh][posOne] = p;
//        rest << "vertex " << vh.idx() << " origin from he " << he.idx() << " (adj. Face "
//             << he.face().idx() << ")\nposition U / V: " << vertexPosUi[vh] + posOne << " / "
//             << vertexPosVi[vh] + posOne << std::endl << std::endl;
    } else {
        int posU = vertexPosUi[he.to()];
        int posV = vertexPosVi[he.to()];
        double u = _x[posU];
        double v = _x[posV];
        OpenMesh::Vec2d p = {u, v};
        solCoordSysUV[he.to()][0] = p;
        trimesh_.status(he.to()).set_tagged(true);
//        rest << "vertex " << he.to().idx() << " origin from he " << he.idx() << " (adj. Face "
//             << he.face().idx() << ")\nposition U / V: " << vertexPosUi[he.to()] << " / "
//             << vertexPosVi[he.to()] << std::endl << std::endl;
    }
}

void GlobalParametrization::setFeatureLineConstraint(gmm::row_matrix<gmm::wsvector<double>> &_constraints,
                                                     std::vector<int> &onlyBoundaries, const int startingPoint) {
    auto currentPJ = OpenMesh::FProp<int>(trimesh_, "currentPJ");
    auto periodJump = OpenMesh::HProp<int>(trimesh_, "periodJump");
    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    int counter = 0;
    for (const int& i: onlyBoundaries) {
        auto heh = make_smart(trimesh_.halfedge_handle(i), trimesh_);
        int posFrom = 0, posTo = 0, rotation = currentPJ[heh.face()] % 4;
        if (vertexAppearanceCG[heh.to()] > 1) {
            auto vh = heh.to();
            posTo = getPositionConstraintRow(vh, cutGraphFZone[heh.face()]);
        }
        if (vertexAppearanceCG[heh.from()] > 1) {
            auto vh = heh.from();
            posFrom = getPositionConstraintRow(vh, cutGraphFZone[heh.face()]);
        }
        switch (rotation) {
            case 0:
            case 2:
            case -2:
                _constraints(startingPoint + counter, vertexPosVi[heh.from()] + posFrom) = 1;
                _constraints(startingPoint + counter++, vertexPosVi[heh.to()] + posTo) = -1;
                break;
            case 1:
            case -1:
            case -3:
            case 3:
                _constraints(startingPoint + counter, vertexPosUi[heh.from()] + posFrom) = 1;
                _constraints(startingPoint + counter++, vertexPosUi[heh.to()] + posTo) = -1;
                break;
            default:
                break;
        }
    }
}



