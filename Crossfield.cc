//
// Created by wuschelbueb on 15.09.21.
//

#include "Crossfield.hh"

void Crossfield::getCrossfield() {
    std::vector<int> faces = getFacesVecWithRefHeProp();
    setVecFieldProp();
    getConstraintAngleAndVecField(faces);
    std::map<int, double> heKappa = getMapHeKappa(faces);

    // indexes the faces from 0 to n
    // this is needed in case a face index is higher than the matrix size
    setPositionInHessianForFaces(heKappa);
    CMatrixType _H = getHessianMatrix(faces, heKappa);
    RMatrixType _constraints = getConstraintMatrix(heKappa, faces);
    std::vector<double> _x(heKappa.size() + faces.size(), 0.0);
    std::vector<double> _rhs = getRHS(heKappa, faces);
    std::vector<int> _idx_to_round = getIdxToRound(heKappa, faces.size());
//    std::cout << "Matrices before solver: " << std::endl;
//    std::cout << "H: " << _H << std::endl;
//
//    for (std::size_t i = 0, max = _rhs.size(); i != max; ++i) {
//        std::cout << "rhs[" << i << "] = " << _rhs[i] << std::endl;
//    }
//    std::cout << std::endl;
//    for (std::size_t i = 0, max = _idx_to_round.size(); i != max; ++i) {
//        std::cout << "idx_to_round[" << i << "] = " << _idx_to_round[i] << std::endl;
//    }
//    std::cout << std::endl;
//    for (std::size_t i = 0, max = _x.size(); i != max; ++i) {
//        std::cout << "_x[" << i << "] = " << _x[i] << std::endl;
//    }
//    std::cout << std::endl;
//    std::cout << "constraints: " << _constraints << std::endl;
    COMISO::ConstrainedSolver csolver;
    csolver.misolver().set_iter_full(false);
    csolver.misolver().set_local_iters(50000);
    csolver.misolver().set_cg_iters(20);
    csolver.misolver().set_local_error(1e-3);
    csolver.misolver().set_cg_error(1e-3);
    csolver.misolver().set_multiple_rounding();
//    std::cout << "Calculation of Smooth Crossfield\n" << "Dimensions of parameters:\n" << "Constraints Rows:\t"
//              << gmm::mat_nrows(_constraints)
//              << " and columns: " << gmm::mat_ncols(_constraints) << std::endl
//              << "A Rows:\t\t\t\t" << gmm::mat_nrows(_H) << " and columns: " << gmm::mat_ncols(_H) << std::endl
//              << "Size of_constraints _x:\t\t\t" << _x.size() << std::endl << "Size of _rhs:\t\t" << _rhs.size()
//              << std::endl
//              << "Size of heKappa:\t" << heKappa.size() << std::endl << "Size of idx:\t\t" << _idx_to_round.size()
//              << std::endl;
    csolver.solve(_constraints, _H, _x, _rhs, _idx_to_round);
    for (auto &it: _x) {
        if (std::abs(it) < 1E-10) {
            it = 0;
        }
    }
    setRotThetaOfVectorField(faces, _x);
    createCrossfields(faces);
    setPJProp(heKappa, _x, faces.size());
    getCrossFieldIdx(faces, heKappa, _x);

//    std::cout << "Energy: " << energy << std::endl;
//        std::cout << "Matrices after solver: " << std::endl;
//    std::cout << "H: " << _H << std::endl;
//    for (std::size_t i = 0, max = _rhs.size(); i != max; ++i) {
//        std::cout << "_rhs[" << i << "] = " << _rhs[i] << std::endl;
//    }
//    std::cout << std::endl;
//    std::ofstream xVector("/home/wuschelbueb/Desktop/xVector.txt");
//    for (std::size_t i = 0, max = _x.size(); i != max; ++i) {
//        xVector << "_x[" << i << "] = " << _x[i] << std::endl;
//    }
//    xVector.close();
//    std::cout << std::endl;
//    std::cout << "constraints after: " << _constraints << std::endl;
}

void Crossfield::createCrossfields(const std::vector<int> &faces) {
    auto uVectorFieldRotOne = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotOne");
    auto uVectorFieldRotTwo = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotTwo");
    for (const int &i: faces) {
        std::vector<Point> PiRotVec;
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        Point rotVec = uVectorFieldRotOne[fh];
        for (int j = 1; j < 4; ++j) {
            Point temp = rotPointWithRotMatrix(fh, rotVec, M_PI / 2 * j);
            PiRotVec.push_back(temp);
        }
        uVectorFieldRotTwo[fh] = PiRotVec[0];
    }
}

gmm::row_matrix<gmm::wsvector<double>>
Crossfield::getConstraintMatrix(const std::map<int, double> &heKappa, const std::vector<int> &faces) {
    trimesh_.request_vertex_status();
    int cNplusOne = 1, counter = 0, pj_start = faces.size(), counterCxFldConstraints = 0;
    std::vector<int> faceConstraints = getConstrainedHe();
    int pjConstraints = getAmountPJConstraints(faces, heKappa);
    int n_row = faceConstraints.size() + pjConstraints;
    int n_col = heKappa.size() + faces.size();
    gmm::row_matrix<gmm::wsvector<double>> _constraints(n_row, n_col + cNplusOne);
    getThetaConstraints(n_col, counter, faceConstraints, _constraints);
    getPJConstraints(heKappa, counter, pj_start, _constraints);
    getCrossFldIdxConstraints(heKappa, counter, pj_start, _constraints, faces);
    for (auto vh: trimesh_.vertices()) {
        trimesh_.status(vh).set_tagged(false);
    }
    try {
        if (counter != n_row) {
            throw counter;
        }
    } catch (int x) {
        std::cerr << "getConstraintMatrix: amount of rows has to be the: " << n_row << " and not: " << counter << "\n";
    }
    gmm::clean(_constraints, 1E-10);
    trimesh_.release_vertex_status();
    return _constraints;
}

std::vector<int> Crossfield::getConstrainedHe() {
    trimesh_.request_face_status();
    std::vector<int> faceConst, heConst;
    for (auto fh: trimesh_.faces()) {
        trimesh_.status(fh).set_tagged(false);
    }
    for (const int &i: heConstraints_) {
        auto heh = make_smart(trimesh_.halfedge_handle(i), trimesh_);
        //(std::find(faceConst.begin(), faceConst.end(), heh.face().idx()) == faceConst.end())
        if (heh.face().is_valid() && !heh.face().tagged()) {
            trimesh_.status(heh.face()).set_tagged(true);
            heConst.push_back(heh.idx());
        }
    }
    trimesh_.release_face_status();
    return heConst;
}

int Crossfield::getAmountPJConstraints(const std::vector<int> &faces, const std::map<int, double> &heKappa) {
    auto origin_constraint = OpenMesh::FProp<int>(trimesh_, "origin_constraint");
    auto predecessor_face = OpenMesh::FProp<int>(trimesh_, "predecessor_face");
    for (auto vh: trimesh_.vertices()) {
        trimesh_.status(vh).set_tagged(false);
    }
    int counter = 0;
    for (const int &i: faces) {
        auto fh = trimesh_.face_handle(i);
        if (fh.idx() != origin_constraint[fh]) {
            counter++;
        }
    }
    for (const auto &idx: heKappa) {
        OpenMesh::SmartHalfedgeHandle he = make_smart(trimesh_.halfedge_handle(idx.first), trimesh_);
        if (!trimesh_.is_boundary(he.to()) && !trimesh_.status(he.to()).tagged()) {
            trimesh_.status(he.to()).set_tagged(true);
            counter++;
        }
    }
    return counter;
}

void Crossfield::getThetaConstraints(const int n_col, int &counter, std::vector<int> &faceConstraints,
                                     gmm::row_matrix<gmm::wsvector<double>> &_constraints) {
    auto constraint_angle = OpenMesh::FProp<double>(trimesh_, "constraint_angle");
    auto positionHessianMatrix = OpenMesh::FProp<int>(trimesh_, "positionHessianMatrix");
    auto referenceHeIdx = OpenMesh::FProp<int>(trimesh_, "referenceHeIdx");
    for (const int &i: faceConstraints) {
        OpenMesh::HalfedgeHandle hehConst = trimesh_.halfedge_handle(i);
        OpenMesh::FaceHandle fh = trimesh_.face_handle(hehConst);
        double constraint = constraint_angle[fh];
        int position = positionHessianMatrix[fh];
        _constraints(counter, position) = 1.0;
        _constraints(counter, n_col) = constraint;
        counter++;
    }
}

void Crossfield::getPJConstraints(const std::map<int, double> &heKappa, int &counter, const int pj_start,
                                  gmm::row_matrix<gmm::wsvector<double>> &_constraints) {
    std::vector<int> faces;
    auto origin_constraint = OpenMesh::FProp<int>(trimesh_, "origin_constraint");
    auto predecessor_face = OpenMesh::FProp<int>(trimesh_, "predecessor_face");
    int iteration = 0, temp = 0;
    for (const auto &j: heKappa) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(j.first);
        OpenMesh::FaceHandle fh1 = trimesh_.face_handle(heh);
        OpenMesh::FaceHandle fh2 = trimesh_.face_handle(trimesh_.opposite_halfedge_handle(heh));
        int fh1_origin = origin_constraint[fh1];
        int fh2_origin = origin_constraint[fh2];
        int pre_fh1 = predecessor_face[fh1];
        int pre_fh2 = predecessor_face[fh2];
        if ((fh1_origin == fh2_origin) && (pre_fh1 == fh2.idx() || pre_fh2 == fh1.idx())) {
            _constraints(counter++, pj_start + iteration++) = 1;
        } else {
            iteration++;
        }
    }
//    std::cout << "there are " << temp << " skips of " << heKappa.size() << " which means there are "
//              << heKappa.size() - temp << " pj constraints\n";
}

void Crossfield::getCrossFldIdxConstraints(const std::map<int, double> &heKappa, int &counter, const int pj_start,
                                           gmm::row_matrix<gmm::wsvector<double>> &_constraints,
                                           const std::vector<int> &faces) {
    auto heColor = OpenMesh::HProp<int>(trimesh_, "heColor");
    auto vertexColor = OpenMesh::VProp<int>(trimesh_, "vertexColor");
    for (auto vh: trimesh_.vertices()) {
        trimesh_.status(vh).set_tagged(false);
    }
    for (const auto &fhIdx: faces) {
        auto fh = trimesh_.face_handle(fhIdx);
        for (TriMesh::FaceVertexIter fv_it = trimesh_.fv_iter(fh); fv_it.is_valid(); ++fv_it) {
            for (TriMesh::VertexOHalfedgeIter vohe_it = trimesh_.voh_iter(*fv_it); vohe_it.is_valid(); ++vohe_it) {
                if (heColor[*vohe_it] == 1) {
                    trimesh_.status(*fv_it).set_tagged(true);
                    vertexColor[*fv_it] = 1;
                }
            }
        }
    }
    //calculate -4*I_o(v_i)
    for (const auto &it: heKappa) {
        double intValBaseIdx = 0, sumKappa = 0, angleDefect = 0;
        OpenMesh::SmartHalfedgeHandle he = make_smart(trimesh_.halfedge_handle(it.first), trimesh_);
        if (trimesh_.status(he.to()).tagged() || trimesh_.is_boundary(he.to())) {
            continue;
        }
//        std::cout << "START, he idx: " << he.to().idx() << std::endl;
        for (TriMesh::VertexOHalfedgeIter voh_it = trimesh_.voh_iter(he.to()); voh_it.is_valid(); ++voh_it) {
            vertexColor[he.to()] = 2;
            int position = 0;
            //check if voh part of hekappa
            auto itHe = heKappa.find(voh_it->idx());
            //check if voh.opp part of hekappa
            auto itOpp = heKappa.find(voh_it->opp().idx());
            angleDefect += trimesh_.calc_sector_angle(voh_it->prev());
            if (itHe != heKappa.end()) {
                // if so get position in heKappa
                position = distance(heKappa.begin(), itHe);
//                std::cout << "he idx " << itHe->first << " with value: " << itHe->second
//                          << " and position: " << position << std::endl;
                sumKappa -= itHe->second;
                // add to constraint matrix
                _constraints(counter, pj_start + position) = -1;
            }
            if (itOpp != heKappa.end()) {
                position = distance(heKappa.begin(), itOpp);
//                std::cout << "heOpp idx " << itOpp->first << " with value: " << itOpp->second
//                          << " and position: " << position << std::endl;
                sumKappa += itOpp->second;
                _constraints(counter, pj_start + position) = 1;
            }
        }
        trimesh_.status(he.to()).set_tagged(true);
        intValBaseIdx = 4 * ((2 * M_PI - angleDefect) + sumKappa) / (2 * M_PI);
//        std::cout << "vertex Value base index: " << intValBaseIdx << std::endl;
        _constraints(counter++, faces.size() + heKappa.size()) = intValBaseIdx;
    }
}

std::vector<int> Crossfield::getIdxToRound(const std::map<int, double> &heKappa, int pj_start) {
    std::vector<int> _idx_to_round;
    for (const auto &i: heKappa) {
        _idx_to_round.push_back(pj_start++);
    }
    return _idx_to_round;
}

std::vector<double> Crossfield::getRHS(const std::map<int, double> &heKappa, const std::vector<int> &faces) {
    std::vector<double> _rhs(faces.size() + heKappa.size());
    int facesPlusOne = faces.size();
    getRhsFirstHalf(faces, _rhs, heKappa);
    getRhsSecondHalf(_rhs, heKappa, facesPlusOne);
    // scalar multiplication with vector, b*-1 = -b
    gmm::scale(_rhs, -1.0);
    return _rhs;
}

void Crossfield::getRhsFirstHalf(const std::vector<int> &faces, std::vector<double> &_rhs,
                                 const std::map<int, double> &heKappa) {
    for (const int &i: faces) {
        getSum(i, _rhs, heKappa);
    }
}

void Crossfield::getSum(const int i, std::vector<double> &_rhs,
                        const std::map<int, double> &heKappa) {
    auto positionHessianMatrix = OpenMesh::FProp<int>(trimesh_, "positionHessianMatrix");
    double sum = 0;
    OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
    for (TriMesh::FaceHalfedgeIter fh_it = trimesh_.fh_iter(fh); fh_it.is_valid(); ++fh_it) {
        // this condition is necessary because an opposite face doesn't exist if the opposite heh is boundary
        sum += extractKappa(*fh_it, heKappa);
    }
    int position = positionHessianMatrix[fh];
//    std::cout << "face index:\t\t" << fh.idx() << " \nposition in hessian matrix: " << position << std::endl;
    _rhs[position] = sum;
}

double
Crossfield::extractKappa(OpenMesh::HalfedgeHandle fh_it, const std::map<int, double> &heKappa) {
    double sum = 0;
    OpenMesh::HalfedgeHandle opposite_heh = trimesh_.opposite_halfedge_handle(fh_it);
    if (!trimesh_.is_boundary(opposite_heh)) {
        auto it = heKappa.find(fh_it.idx());
        auto it2 = heKappa.find(opposite_heh.idx());
        if (it != heKappa.end()) {
            sum += (it->second * 2);
//            std::cout << "Kappa value main triangle: " << it->second << "\nsum main:\t" << sum << std::endl;
        }
        if (it2 != heKappa.end()) {
            sum -= (it2->second * 2);
//            std::cout << "Kappa value neighbour triangle: " << it2->second << "\nsum neigh:\t" << sum << std::endl;
        }
        try {
            if (it != heKappa.end() && it2 != heKappa.end()) {
                throw 404;
            }
        } catch (int x) {
            std::cerr << "extractKappa: Opposite Halfedges can't both be in heKappa vector!\n";
        }
    }
    return sum;
}

void
Crossfield::getRhsSecondHalf(std::vector<double> &_rhs, const std::map<int, double> &heKappa, const int facesPlusOne) {
    int counter = facesPlusOne;
    for (const auto &i: heKappa) {
        _rhs[counter++] = (i.second * M_PI);
    }
}

Crossfield::CMatrixType
Crossfield::getHessianMatrix(const std::vector<int> &faces, const std::map<int, double> &heKappa) {
    auto positionHessianMatrix = OpenMesh::FProp<int>(trimesh_, "positionHessianMatrix");
    int counter = 0, iteration = 0, pj_start = faces.size(), n = heKappa.size() + faces.size();
    gmm::col_matrix<gmm::wsvector<double>> _H(n, n);

    // fills up sparse column matrix _H
    for (const auto &i: heKappa) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i.first);
        OpenMesh::FaceHandle fh1 = trimesh_.face_handle(heh);
        OpenMesh::FaceHandle fh2 = trimesh_.face_handle(trimesh_.opposite_halfedge_handle(heh));
        int factorI = getFactor(fh1, faces);
        int factorJ = getFactor(fh2, faces);
        int pos_i = positionHessianMatrix[fh1];
        int pos_j = positionHessianMatrix[fh2];
//        std::cout << "face: " << fh1.idx() << " (" << pos_i << ") with factor: " << factorI << " and face: "
//                  << fh2.idx()
//                  << "(" << pos_j << ") with factor: " << factorJ << std::endl;
        // |  2*w_ij   -2*w_ij     pi*w_ij |
        // | -2*w_ij    2*w_ij    -pi*w_ij |
        // | pi*w_ij  -pi*w_ij pi^2/2*w_ij |
        _H(pos_i, pos_i) = 2.0 * factorI;
        _H(pos_j, pos_j) = 2.0 * factorJ;
        _H(pos_i, pos_j) = -2.0;
        _H(pos_j, pos_i) = -2.0;
        _H(pos_i, pj_start + iteration) = M_PI;
        _H(pj_start + iteration, pos_i) = M_PI;
        _H(pos_j, pj_start + iteration) = -M_PI;
        _H(pj_start + iteration, pos_j) = -M_PI;
        _H(pj_start + iteration, pj_start + iteration) = ((pow(M_PI, 2) / 2));
        iteration++;
    }
    gmm::clean(_H, 1E-10);
    return _H;
}

int Crossfield::getFactor(const OpenMesh::FaceHandle fh, const std::vector<int> &faces) {
    int factor = 0;
    for (TriMesh::FaceFaceIter ff_it = trimesh_.ff_iter(fh); ff_it.is_valid(); ++ff_it) {
        if ((std::find(faces.begin(), faces.end(), ff_it->idx()) != faces.end()))
            factor++;
    }
    return factor;
}

void Crossfield::setPositionInHessianForFaces(const std::map<int, double> &heKappa) {
    auto positionHessianMatrix = OpenMesh::FProp<int>(trimesh_, "positionHessianMatrix");
    trimesh_.request_face_status();
    for (auto fh: trimesh_.faces()) {
        trimesh_.status(fh).set_tagged(false);
    }
    int counter = 0;
    for (const auto &i: heKappa) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i.first);
        OpenMesh::FaceHandle fh1 = trimesh_.face_handle(heh);
        OpenMesh::FaceHandle fh2 = trimesh_.face_handle(trimesh_.opposite_halfedge_handle(heh));
        if (!trimesh_.status(fh1).tagged()) {
            positionHessianMatrix[fh1] = counter;
            trimesh_.status(fh1).set_tagged(true);
//            std::cout << "face: " << fh1.idx() << " with pos: " << positionHessianMatrix[fh1] << " and counter: " << counter << std::endl;
            counter++;
        }
        if (!trimesh_.status(fh2).tagged()) {
            positionHessianMatrix[fh2] = counter;
            trimesh_.status(fh2).set_tagged(true);
//            std::cout << "face: " << fh2.idx() << " with pos: " << positionHessianMatrix[fh2] << " and counter: " << counter << std::endl;
            counter++;
        }
    }
    trimesh_.release_face_status();
}

std::map<int, double> Crossfield::getMapHeKappa(const std::vector<int> &faces) {
    std::map<int, double> heKappa;
    for (const int &i: faces) {
//        std::cout << "face: " << i << std::endl;
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        for (TriMesh::FaceFaceIter ff_it = trimesh_.ff_iter(fh); ff_it.is_valid(); ++ff_it) {
//            std::cout << "\tgetStatusNeigh call with facehandle: " << fh.idx() << " and neighbour: " << ff_it->idx()
//                      << std::endl;
            getStatusNeigh(fh, *ff_it, heKappa);
        }
    }
    try {
        if (heKappa.empty()) {
            throw 404;
        }
    } catch (int x) {
        std::cerr << "getMapHeKappa: heKappa Map can't be empty!\n";
    }
    return heKappa;
}

void Crossfield::getStatusNeigh(const OpenMesh::FaceHandle fh, const OpenMesh::FaceHandle fh_neigh,
                                std::map<int, double> &heKappa) {
    auto faceSel = OpenMesh::FProp<bool>(trimesh_, "faceSel");
    auto referenceHeIdx = OpenMesh::FProp<int>(trimesh_, "referenceHeIdx");
    auto constraint_angle = OpenMesh::FProp<double>(trimesh_, "constraint_angle");
    double kappa = DBL_MAX;
    std::pair<int, int> commonEdge = {INT_MAX, INT_MAX};
//    std::cout << "\t\tin getStatusNeigh function" << std::endl;
//    neighbour needs to be in faces
//    trimesh_.status(fh_neigh).tagged()
    if (faceSel[fh_neigh]) {
        // get index of ref edges
        int refEdgeMain = referenceHeIdx[fh];
        int refEdgeNeigh = referenceHeIdx[fh_neigh];
        // get the common edge between the triangles
        commonEdge = getCommonEdgeBetweenTriangles(fh, fh_neigh);
//        std::cout << "\t\t\tcall getKappa, refEdgeMain: " << refEdgeMain << " and refEdgeNeigh: " << refEdgeNeigh
//        << std::endl;
        kappa = getKappa(refEdgeMain, refEdgeNeigh, commonEdge);
        addKappaHeToMap(commonEdge, kappa, heKappa);
    }
}

std::pair<int, int>
Crossfield::getCommonEdgeBetweenTriangles(const OpenMesh::FaceHandle fh, const OpenMesh::FaceHandle fh_neigh) {
    std::pair<int, int> commonEdge = {INT_MAX, INT_MAX};
    for (TriMesh::FaceHalfedgeIter fhe_it_main = trimesh_.fh_iter(fh); fhe_it_main.is_valid(); ++fhe_it_main) {
        if (fhe_it_main->opp().face() == fh_neigh) {
            commonEdge.first = fhe_it_main->idx();
            commonEdge.second = fhe_it_main->opp().idx();
//            std::cout << "\t\t\t\tcommonedge: " << commonEdge.first << " and " << commonEdge.second
//                      << "\n\t\t\t\twhile opposites are " << trimesh_.opposite_halfedge_handle(*fhe_it_main).idx()
//                      << " and "
//                      << trimesh_.opposite_halfedge_handle(oheh).idx() << std::endl;
        }
    }
    try {
        if (commonEdge.first == INT_MAX || commonEdge.second == INT_MAX) {
            throw 404;
        }
    } catch (int n) {
        std::cerr << "getCommonEdgeBetweenTriangles: commonEdge can't have the value INT_MAX!\n";
    }
    return commonEdge;
}

double
Crossfield::getKappa(const int refEdgeMain, const int refEdgeNeigh, const std::pair<int, int> commonEdge) {
    // refEdge shares common edge
    double alpha = 0.0, beta = alpha, kappa = alpha;
    // counter starts at one because there is always at least a PI to simulate the turn form one triangle to the other
    int tempHe = refEdgeMain, counter = 1;
    while (tempHe != commonEdge.first) {
//        std::cout << tempHe << std::endl;
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(tempHe);
        alpha += trimesh_.calc_sector_angle(heh);
        tempHe = trimesh_.next_halfedge_handle(heh).idx();
        counter++;
    }
    tempHe = commonEdge.second;
    while (tempHe != refEdgeNeigh) {
//        std::cout << tempHe << std::endl;
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(tempHe);
        beta += trimesh_.calc_sector_angle(heh);
        tempHe = trimesh_.next_halfedge_handle(heh).idx();
        counter++;
    }
    kappa = counter * M_PI - alpha - beta;
    // make sure kappa is in range of [-pi, pi)
    double kappa_final = shortenKappa(kappa);
//    std::cout << "\t\t\t\tkappa: " << kappa << " shorted: " << kappa_final << std::endl;
    return kappa_final;
}

void Crossfield::addKappaHeToMap(const std::pair<int, int> commonEdge, const double kappa,
                                 std::map<int, double> &heKappa) {
    // if the opposite he isn't already in the list, add this he to the list
    if (heKappa.find(commonEdge.second) == heKappa.end()) {
        heKappa[commonEdge.first] = kappa;
//        std::cout << "\t\t\t\tkappa and halfedge were added to map\n";
    } else {
//        std::cout << "\t\t\t\tkappa and halfedge weren't added to map\n";
    }
}

std::vector<int> Crossfield::getFacesVecWithRefHeProp() {

    // request to change the status
    trimesh_.request_halfedge_status();
    trimesh_.request_face_status();
    for (auto fh: trimesh_.faces()) {
        trimesh_.status(fh).set_tagged(false);
    }
    std::vector<int> faces;
    // add faces of constraint halfedges to faces vector and add reference halfedge
    for (int i: heConstraints_) {
        setRefHeToFace(i, faces);
    }

    // add faces of non-constraint halfedges to faces vector
    for (int i: heInRange_) {
        setRefHeToFace(i, faces);
    }
    colorHEdges();
    colorFaces(faces);
    trimesh_.release_halfedge_status();
    trimesh_.release_face_status();
    return faces;
}

void Crossfield::setRefHeToFace(const int i, std::vector<int> &faces) {
    auto referenceHeIdx = OpenMesh::FProp<int>(trimesh_, "referenceHeIdx");
    OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i);
    if (!trimesh_.is_boundary(heh)) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(heh);
        if (!trimesh_.status(fh).tagged()) {
            // add reference edge to face
            referenceHeIdx[fh] = heh.idx();
            trimesh_.status(heh).set_tagged(true);
            trimesh_.status(fh).set_tagged(true);
            faces.push_back(fh.idx());
        }
    }
}


void Crossfield::setVecFieldProp() {
    auto uVectorField = OpenMesh::FProp<Point>(trimesh_, "uVectorField");
    auto uVectorFieldRotOne = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotOne");
    auto uVectorFieldRotTwo = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotTwo");
    Point zero_point = {0, 0, 0};
    for (TriMesh::FaceIter f_it = trimesh_.faces_begin(); f_it != trimesh_.faces_end(); ++f_it) {
        uVectorField[*f_it] = zero_point;
        uVectorFieldRotOne[*f_it] = zero_point;
        uVectorFieldRotTwo[*f_it] = zero_point;
    }
}

void Crossfield::getConstraintAngleAndVecField(const std::vector<int> &faces) {
    auto referenceHeIdx = OpenMesh::FProp<int>(trimesh_, "referenceHeIdx");
    auto constraint_angle = OpenMesh::FProp<double>(trimesh_, "constraint_angle");
    auto uVectorField = OpenMesh::FProp<Point>(trimesh_, "uVectorField");
    for (const int &i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(referenceHeIdx[fh]);
        Point u = trimesh_.calc_edge_vector(heh);
        Point v = refVector_;
        double alpha = std::atan2(v[1], v[0]) - std::atan2(u[1], u[0]);
        if (alpha < 1E-10 && alpha > -1E-10) {
            alpha = 0;
        }
        if (alpha < 0) {
            alpha += 2 * M_PI;
        }
        uVectorField[fh] = u;
        constraint_angle[fh] = alpha;
//        std::cout << "theta_c (" << i << ") is: " << alpha * 180 / M_PI << " degrees." << std::endl;
    }
}

void Crossfield::setRotThetaOfVectorField(const std::vector<int> &faces, const std::vector<double> _x) {
    auto positionHessianMatrix = OpenMesh::FProp<int>(trimesh_, "positionHessianMatrix");
    auto uVectorField = OpenMesh::FProp<Point>(trimesh_, "uVectorField");
    auto uVectorFieldRotOne = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotOne");
    for (const int &i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        int position = positionHessianMatrix[fh];
        Point p_x = uVectorField[fh];
        double theta = shortenKappa(_x[position]);
        Point temp = rotPointWithRotMatrix(fh, p_x, theta);
        uVectorFieldRotOne[fh] = temp.normalize();
//        std::cout << "face (" << i << ") with position " << position << " and value (" << _x[position]
//                  << "): with angle: " << theta << "(rad) and "
//                  << theta * 180 / M_PI << "(deg)"
//                  << std::endl;
    }
}

void Crossfield::getCrossFieldIdx(const std::vector<int> &faces, const std::map<int, double> &heKappa,
                                  const std::vector<double> &_x) {
    trimesh_.request_vertex_status();
    auto periodJump = OpenMesh::HProp<int>(trimesh_, "periodJump");
    auto vertexColor = OpenMesh::VProp<int>(trimesh_, "vertexColor");
    auto crossFieldIdx = OpenMesh::VProp<double>(trimesh_, "crossFieldIdx");
    auto positionHessianMatrix = OpenMesh::FProp<int>(trimesh_, "positionHessianMatrix");
    for (auto vh: trimesh_.vertices()) {
        crossFieldIdx[vh] = 0;
        trimesh_.status(vh).set_tagged(false);
    }
//    std::ofstream xVector("/home/wuschelbueb/Desktop/newGetCxIdx.txt");
    for (const auto &heIdx: heKappa) {
        auto he = make_smart(trimesh_.halfedge_handle(heIdx.first), trimesh_);
        if (he.to().tagged() || he.to().is_boundary() || vertexColor[he.to()] == 1) {
            continue;
        }
//        xVector << "vertex idx: " << he.to().idx() << std::endl;
        double sumKappa = 0.0, angleDefect = 0.0, sumPJ = 0.0, intValBaseIdx = 0.0;
        for (TriMesh::VertexOHalfedgeIter vohe_it = trimesh_.voh_iter(he.to()); vohe_it.is_valid(); ++vohe_it) {
            auto it = heKappa.find(vohe_it->idx());
            auto it2 = heKappa.find(vohe_it->opp().idx());
            angleDefect += trimesh_.calc_sector_angle(vohe_it->prev());
//            xVector << "angleDefect individual " << trimesh_.calc_sector_angle(vohe_it->prev()) << std::endl;
            if (it != heKappa.end()) {
                //with position and faceSize we can extract p_ij value of _x vector (solution)
                sumPJ -= periodJump[*vohe_it];
                sumKappa -= it->second;
//                xVector << "he idx: " << it->first << " with kappa: " << -it->second <<
//                        "\npj: " << -periodJump[*vohe_it] << std::endl;
            } else if (it2 != heKappa.end()) {
                sumPJ += periodJump[vohe_it->opp()];
                sumKappa += it2->second;
//                xVector << "he idx: " << it2->first << " with kappa: " << it2->second
//                        << "\npj: " << periodJump[vohe_it->opp()] << std::endl;
            }
        }
        intValBaseIdx = ((2 * M_PI - angleDefect) + sumKappa) / (2 * M_PI);
        double crossFldIdx = intValBaseIdx + (sumPJ / 4);
        //clean values close to zero
        if (std::abs(crossFldIdx) < 1E-10) {
            crossFldIdx = 0;
        }
//        xVector << "Crossfield Index (" << he.idx() << "): " << crossFldIdx
//                << "\nsum Kappa: " << sumKappa
//                << "\nsum PJ: " << sumPJ
//                << "\nangle defect: " << angleDefect
//                << "\nintValeBaseIdx: " << intValBaseIdx << std::endl;
        crossFieldIdx[he.to()] = crossFldIdx;
        trimesh_.status(he.to()).set_tagged(true);
    }
//    xVector.close();
    trimesh_.release_vertex_status();
}

void Crossfield::colorFaces(const std::vector<int> &faces) {
    auto faceSel = OpenMesh::FProp<bool>(trimesh_, "faceSel");
    for (TriMesh::FaceIter f_it = trimesh_.faces_begin(); f_it != trimesh_.faces_end(); ++f_it) {
        faceSel[*f_it] = false;
    }
    for (const int &i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        faceSel[fh] = true;
    }
}

void Crossfield::colorHEdges() {
    auto heColor = OpenMesh::HProp<int>(trimesh_, "heColor");
    for (TriMesh::HalfedgeIter he_it = trimesh_.halfedges_begin(); he_it != trimesh_.halfedges_end(); ++he_it) {
        heColor(*he_it) = 1;
    }
    for (const int &i: heInRange_) {
        heColor(trimesh_.halfedge_handle(i)) = 2;
    }
    for (const int &i: heConstraints_) {
        heColor(trimesh_.halfedge_handle(i)) = 3;
    }
}

double Crossfield::shortenKappa(const double kappa) {
    double temp = kappa;
    while (temp > M_PI) {
        temp -= 2 * M_PI;
    }
    while (temp <= -M_PI) {
        temp += 2 * M_PI;
    }
//    std::cout << "kappa is: " << kappa << " and the contracted kappa is: " << temp << std::endl;
    return temp;
}

Crossfield::Point
Crossfield::rotPointWithRotMatrix(const OpenMesh::FaceHandle fh, const Point vec, const double theta) {
    Point f_normal = trimesh_.calc_face_normal(fh).normalize();
    Point rotVec = vec * cos(theta) - (f_normal % vec) * sin(theta) + f_normal * (vec | f_normal) * (1 - cos(theta));
    return rotVec;
}

void Crossfield::setPJProp(std::map<int, double> &heKappa, std::vector<double> &_x, const int faceSize) {
    auto periodJump = OpenMesh::HProp<int>(trimesh_, "periodJump");
    int position = 0;
    for (const auto &it: heKappa) {
        auto he = trimesh_.halfedge_handle(it.first);
        periodJump[he] = std::round(_x[faceSize + position++]);
    }
}
