//
// Created by wuschelbueb on 15.09.21.
//

#include "Crossfield.hh"

void Crossfield::getCrossfield() {
    std::vector<int> faces = getFacesVecWithRefHeProp();
    setVecFieldProp();
    getConstraintAngleAndVecField(faces);
    std::map<int, double> heKappa = getMapHeKappa(faces);
    CMatrixType _H = getHessianMatrix(faces, heKappa);
    RMatrixType _constraints = getConstraintMatrix(heKappa, faces);
    std::vector<double> _x(heKappa.size() + faces.size(), 0.0);
    std::vector<double> _rhs = getRHS(heKappa, faces);
    std::vector<int> _idx_to_round = getIdxToRound(heKappa, faces);
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
    std::cout << "Dimensions of parameters:\n" << "Constraints Rows:\t" << gmm::mat_nrows(_constraints)
              << " and columns: " << gmm::mat_ncols(_constraints) << std::endl
              << "A Rows:\t\t\t\t" << gmm::mat_nrows(_H) << " and columns: " << gmm::mat_ncols(_H) << std::endl
              << "Size of _x:\t\t\t" << _x.size() << std::endl << "Size of _rhs:\t\t" << _rhs.size() << std::endl
              << "Size of heKappa:\t" << heKappa.size() << std::endl << "Size of idx:\t\t" << _idx_to_round.size()
              << std::endl;
    csolver.solve(_constraints, _H, _x, _rhs, _idx_to_round);
    setRotThetaOfVectorField(faces, _x);
    createCrossfields(faces);
    RMatrixType _A = getMatrixA(faces, heKappa);
    std::vector<double> _b = getVectorb(heKappa);
    double energy = getEnergy(_A, _x, _b);
    getCrossFieldIdx(faces, heKappa, _x);
//    std::cout << "Energy: " << energy << std::endl;
//        std::cout << "Matrices after solver: " << std::endl;
//    std::cout << "H: " << _H << std::endl;
//    for (std::size_t i = 0, max = _rhs.size(); i != max; ++i) {
//        std::cout << "_rhs[" << i << "] = " << _rhs[i] << std::endl;
//    }
//    std::cout << std::endl;
//    for (std::size_t i = 0, max = _x.size(); i != max; ++i) {
//        std::cout << "_x[" << i << "] = " << _x[i] << std::endl;
//    }
//    std::cout << std::endl;
//    std::cout << "constraints after: " << _constraints << std::endl;
}

void Crossfield::createCrossfields(const std::vector<int> &faces) {
    auto uVectorFieldRotOne = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotOne");
    auto uVectorFieldRotTwo = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotTwo");
    auto uVectorFieldRotThree = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotThree");
    auto uVectorFieldRotFour = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotFour");
    for (int i: faces) {
        std::vector<Point> PiRotVec;
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        Point rotVec = uVectorFieldRotOne[fh];
        for (int j = 1; j < 4; ++j) {
            Point temp = rotPointWithRotMatrix(fh, rotVec, M_PI / 2 * j);
            PiRotVec.push_back(temp);
        }
        uVectorFieldRotTwo[fh] = PiRotVec[0];
        uVectorFieldRotThree[fh] = PiRotVec[1];
        uVectorFieldRotFour[fh] = PiRotVec[2];
    }
}

double Crossfield::getEnergy(const RMatrixType &_A, const std::vector<double> &_x, const std::vector<double> &_b) {
    std::vector<double> _energy(_b.size());
    double finalEnergy = 0;
    gmm::mult(_A, _x, _b, _energy);
    gmm::clean(_energy, 1E-10);
    return gmm::vect_sp(_energy, _energy);
}

Crossfield::RMatrixType Crossfield::getMatrixA(const std::vector<int> &faces, const std::map<int, double> &heKappa) {
    auto positionHessianMatrix = OpenMesh::FProp<int>(trimesh_, "positionHessianMatrix");
    int columns = faces.size() + heKappa.size(), rows = heKappa.size(), counter = 0, pj_start = faces.size();
    RMatrixType _A(rows, columns);
    for (auto i: heKappa) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i.first);
        OpenMesh::FaceHandle fh1 = trimesh_.face_handle(heh);
        OpenMesh::FaceHandle fh2 = trimesh_.face_handle(trimesh_.opposite_halfedge_handle(heh));
        int pos_i = positionHessianMatrix[fh1];
        int pos_j = positionHessianMatrix[fh2];
        _A(counter, pos_i) = 1;
        _A(counter, pos_j) = -1;
        _A(counter, pj_start + counter) = (M_PI / 2);
        counter++;
    }
    return _A;
}

std::vector<double> Crossfield::getVectorb(const std::map<int, double> &heKappa) {
    std::vector<double> _b;
    for (auto i: heKappa) {
        _b.push_back(i.second);
    }
    return _b;
}

gmm::row_matrix<gmm::wsvector<double>>
Crossfield::getConstraintMatrix(const std::map<int, double> &heKappa, const std::vector<int> &faces) {
    int cNplusOne = 1, counter = 0, pj_start = faces.size();
    std::vector<int> faceConstraints = getFaceConstraints();
    int pjConstraints = getAmountPJConstraints(faces);
    int n_row = faceConstraints.size() + pjConstraints;
    int n_col = heKappa.size() + faces.size();
    gmm::row_matrix<gmm::wsvector<double>> _constraints(n_row, n_col + cNplusOne);
    getThetaConstraints(n_col, counter, faceConstraints, _constraints);
    getPJConstraints(heKappa, counter, pj_start, pjConstraints, _constraints);
    try {
        if (counter != n_row) {
            throw counter;
        }
    } catch (int x) {
        std::cerr << "getConstraintMatrix: amount of rows has to be the: " << n_row << " and not: " << counter << "\n";
    }
    gmm::clean(_constraints, 1E-10);
    return _constraints;
}

std::vector<int> Crossfield::getFaceConstraints() {
    std::vector<int> faceConst;
    for (int i: heConstraints_) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i);
        if ((std::find(faceConst.begin(), faceConst.end(), heh.idx()) == faceConst.end()) &&
            !trimesh_.is_boundary(heh)) {
            faceConst.push_back(heh.idx());
        }
    };
    return faceConst;
}

int Crossfield::getAmountPJConstraints(const std::vector<int> &faces) {
    auto origin_constraint = OpenMesh::FProp<int>(trimesh_, "origin_constraint");
    auto predecessor_face = OpenMesh::FProp<int>(trimesh_, "predecessor_face");
    int counter = 0;
    for (int i: faces) {
        auto fh = trimesh_.face_handle(i);
        if (fh.idx() != origin_constraint[fh]) {
            counter += 1;
        }
    }
    return counter;
}

void Crossfield::getThetaConstraints(const int n_col, int &counter, std::vector<int> &faceConstraints,
                                     gmm::row_matrix<gmm::wsvector<double>> &_constraints) {
    auto constraint_angle = OpenMesh::FProp<double>(trimesh_, "constraint_angle");
    auto positionHessianMatrix = OpenMesh::FProp<int>(trimesh_, "positionHessianMatrix");
    auto referenceHeIdx = OpenMesh::FProp<int>(trimesh_, "referenceHeIdx");
    for (int i: faceConstraints) {
        OpenMesh::HalfedgeHandle hehConst = trimesh_.halfedge_handle(i);
        OpenMesh::FaceHandle fh = trimesh_.face_handle(hehConst);
        double constraint = constraint_angle[fh];
        int position = positionHessianMatrix[fh];
        _constraints(counter, position) = 1.0;
        _constraints(counter, n_col) = 1.0 * constraint;
        counter++;
    }
}

void Crossfield::getPJConstraints(const std::map<int, double> &heKappa, int &counter, const int pj_start,
                                  const int &noOriginConst,
                                  gmm::row_matrix<gmm::wsvector<double>> &_constraints) {
    std::vector<int> faces;
    auto origin_constraint = OpenMesh::FProp<int>(trimesh_, "origin_constraint");
    auto predecessor_face = OpenMesh::FProp<int>(trimesh_, "predecessor_face");
    int iteration = 0, temp = 0;
    for (auto j: heKappa) {
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

std::vector<int> Crossfield::getIdxToRound(const std::map<int, double> &heKappa, const std::vector<int> &faces) {
    std::vector<int> _idx_to_round;
    int pj_start = faces.size();
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
    for (int i: faces) {
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
        sum += setSum(fh, *fh_it, heKappa);
    }
    int position = positionHessianMatrix[fh];
//    std::cout << "face index:\t\t" << fh.idx() << " \nposition in hessian matrix: " << position << std::endl;
    _rhs[position] = sum;
}

double
Crossfield::setSum(OpenMesh::FaceHandle fh, OpenMesh::HalfedgeHandle fh_it, const std::map<int, double> &heKappa) {
    double sum = 0;
    OpenMesh::HalfedgeHandle opposite_heh = trimesh_.opposite_halfedge_handle(fh_it);
    if (!trimesh_.is_boundary(opposite_heh)) {
        OpenMesh::FaceHandle opposite_fh = trimesh_.face_handle(trimesh_.opposite_halfedge_handle(fh_it));
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
            std::cerr << "setSum: Opposite Halfedges can't both be in heKappa vector!\n";
        }
    }
    return sum;
}

void
Crossfield::getRhsSecondHalf(std::vector<double> &_rhs, const std::map<int, double> &heKappa, const int facesPlusOne) {
    int counter = facesPlusOne;
    for (const auto &i: heKappa) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i.first);
        _rhs[counter++] = (i.second * M_PI);
    }
}

Crossfield::CMatrixType
Crossfield::getHessianMatrix(const std::vector<int> &faces, const std::map<int, double> &heKappa) {
    trimesh_.release_halfedge_status();
    trimesh_.release_face_status();
    // request to change the status
    trimesh_.request_halfedge_status();
    trimesh_.request_face_status();
    auto positionHessianMatrix = OpenMesh::FProp<int>(trimesh_, "positionHessianMatrix");
    int counter = 0, iteration = 0, pj_start = faces.size(), n = heKappa.size() + faces.size();
    gmm::col_matrix<gmm::wsvector<double>> _H(n, n);
    // indexes the faces from 0 to n
    // this is needed in case a face index is higher than the matrix size
    setPositionInHessianForFaces(heKappa);

    // fills up sparse column matrix _H
    for (auto i: heKappa) {
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
}

std::map<int, double> Crossfield::getMapHeKappa(const std::vector<int> &faces) {
    std::map<int, double> heKappa;
    for (int i: faces) {
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
    auto referenceHeIdx = OpenMesh::FProp<int>(trimesh_, "referenceHeIdx");
    auto constraint_angle = OpenMesh::FProp<double>(trimesh_, "constraint_angle");
    double kappa = DBL_MAX;
    std::pair<int, int> commonEdge = {INT_MAX, INT_MAX};
//    std::cout << "\t\tin getStatusNeigh function" << std::endl;
    // neighbour needs to be in faces
    if (trimesh_.status(fh_neigh).tagged()) {
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
        OpenMesh::HalfedgeHandle oheh = trimesh_.opposite_halfedge_handle(*fhe_it_main);
        OpenMesh::FaceHandle fh_temp = trimesh_.face_handle(oheh);
        if (fh_temp == fh_neigh) {
            commonEdge.first = fhe_it_main->idx();
            commonEdge.second = oheh.idx();
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
    // status needs to be released before using, in case it still has saved some stuff from before
    trimesh_.release_halfedge_status();
    trimesh_.release_face_status();
    // request to change the status
    trimesh_.request_halfedge_status();
    trimesh_.request_face_status();
    std::vector<int> faces;
    int temp = INT_MAX;

    // assign constraint edges to faces as property
    for (int i: heConstraints_) {
        setRefHeToFace(i, faces);
    }

    // assign non-constraint edges to faces as property
    for (int i: heInRange_) {
        setRefHeToFace(i, faces);
    }
    colorHEdges();
    colorFaces(faces);
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
            // both halfedges need to be tagged, else it is possible opposite faces share an edge
            trimesh_.status(heh).set_tagged(true);
//            trimesh_.status(trimesh_.opposite_halfedge_handle(heh)).set_tagged(true);
            trimesh_.status(fh).set_tagged(true);
            faces.push_back(fh.idx());
        }
    }
}


void Crossfield::setVecFieldProp() {
    auto uVectorField = OpenMesh::FProp<Point>(trimesh_, "uVectorField");
    auto uVectorFieldRotOne = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotOne");
    auto uVectorFieldRotTwo = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotTwo");
    auto uVectorFieldRotThree = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotThree");
    auto uVectorFieldRotFour = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotFour");
    Point zero_point = {0, 0, 0};
    for (TriMesh::FaceIter f_it = trimesh_.faces_begin(); f_it != trimesh_.faces_end(); ++f_it) {
        uVectorField[*f_it] = zero_point;
        uVectorFieldRotOne[*f_it] = zero_point;
        uVectorFieldRotTwo[*f_it] = zero_point;
        uVectorFieldRotThree[*f_it] = zero_point;
        uVectorFieldRotFour[*f_it] = zero_point;
    }
}

void Crossfield::getConstraintAngleAndVecField(const std::vector<int> &faces) {
    auto referenceHeIdx = OpenMesh::FProp<int>(trimesh_, "referenceHeIdx");
    auto constraint_angle = OpenMesh::FProp<double>(trimesh_, "constraint_angle");
    auto uVectorField = OpenMesh::FProp<Point>(trimesh_, "uVectorField");
    for (int i: faces) {
        double alpha = 0.0;
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(referenceHeIdx[fh]);
        Point u = trimesh_.calc_edge_vector(heh).normalize();
        Point v = u % trimesh_.calc_face_normal(fh);
        v = v.normalize();
        double xComponent = OpenMesh::dot(u, u);
        double yComponent = OpenMesh::dot(u, v);
        alpha = std::atan2(yComponent, xComponent);
        constraint_angle[fh] = alpha;
//        std::cout << "theta_c (" << i << ") is: " << alpha * 180 / M_PI << " degrees." << std::endl;
        uVectorField[fh] = u;
    }
}

void Crossfield::setRotThetaOfVectorField(const std::vector<int> &faces, const std::vector<double> _x) {
    auto positionHessianMatrix = OpenMesh::FProp<int>(trimesh_, "positionHessianMatrix");
    auto uVectorField = OpenMesh::FProp<Point>(trimesh_, "uVectorField");
    auto uVectorFieldRotOne = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotOne");
    for (int i: faces) {
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
    //TODO works for a cube but nothing else
    trimesh_.release_vertex_status();
    // request to change the status
    trimesh_.request_vertex_status();
    auto crossFieldIdx = OpenMesh::VProp<double>(trimesh_, "crossFieldIdx");
    for (TriMesh::VertexIter v_it = trimesh_.vertices_begin(); v_it != trimesh_.vertices_end(); ++v_it) {
        crossFieldIdx[*v_it] = 0;
    }
    for (int i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        double intValBaseIdx = 0.0;
        for (TriMesh::FaceVertexIter v_it = trimesh_.fv_iter(fh); v_it.is_valid(); ++v_it) {
            if (!trimesh_.status(*v_it).tagged() && !trimesh_.is_boundary(*v_it)) {
                setCrossFieldIdx(v_it, faces.size(), heKappa, _x);
            }
        }
    }
}

void
Crossfield::setCrossFieldIdx(TriMesh::FaceVertexIter &fv_it, const int faceSize, const std::map<int, double> &heKappa,
                             const std::vector<double> &_x) {
    auto crossFieldIdx = OpenMesh::VProp<double>(trimesh_, "crossFieldIdx");
    bool checkOHe = true;
    double sumKappa = 0.0, angleDefect = 0.0, sumPJ = 0.0, intValBaseIdx = 0.0;
    for (TriMesh::VertexOHalfedgeIter vohe_it = trimesh_.voh_iter(*fv_it); vohe_it.is_valid(); ++vohe_it) {
        int position = 0;
        OpenMesh::HalfedgeHandle oppHe = trimesh_.opposite_halfedge_handle(*vohe_it);
        auto it = heKappa.find(vohe_it->idx());
        auto it2 = heKappa.find(oppHe.idx());
        angleDefect += trimesh_.calc_sector_angle(oppHe);
        if (it != heKappa.end()) {
            //return position of it in heKappa
            position = std::distance(std::begin(heKappa), it);
            //with position and faceSize we can extract p_ij value of _x vector (solution)
            sumPJ += _x[faceSize + position];
            sumKappa += it->second;
        } else if (it2 != heKappa.end()) {
            position = std::distance(std::begin(heKappa), it2);
            sumPJ -= _x[faceSize + position];
            sumKappa -= it2->second;
        }
    }
    intValBaseIdx = (1 / (2 * M_PI)) * ((2 * M_PI - angleDefect) + sumKappa);
    crossFieldIdx[*fv_it] = intValBaseIdx + (sumPJ / 4);
    trimesh_.status(*fv_it).set_tagged(true);
//    std::cout << "\nVertex: " << fv_it->idx() << " with\nangle defect: " << (2 * M_PI - angleDefect) / (2 * M_PI)
//              << "\nsumKappa: " << sumKappa / (2 * M_PI)
//              << "\nInteger Value base Index: " << intValBaseIdx << "\nCrossfield Idx is: " << crossFieldIdx[*fv_it]
//              << std::endl;
}

void Crossfield::colorFaces(const std::vector<int> &faces) {
    auto FaceSel = OpenMesh::FProp<int>(trimesh_, "FaceSel");
    for (TriMesh::FaceIter f_it = trimesh_.faces_begin(); f_it != trimesh_.faces_end(); ++f_it) {
        FaceSel[*f_it] = 0;
    }
    for (int i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        FaceSel[fh] = 1;
    }
}

void Crossfield::colorHEdges() {
    auto heColor = OpenMesh::HProp<int>(trimesh_, "heColor");
    for (TriMesh::HalfedgeIter he_it = trimesh_.halfedges_begin(); he_it != trimesh_.halfedges_end(); ++he_it) {
        heColor(*he_it) = 1;
    }
    for (int i: heInRange_) {
        heColor(trimesh_.halfedge_handle(i)) = 2;
    }
    for (int i: heConstraints_) {
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
    //Rodrigues rotation formula
    Point rotVec = vec * cos(theta) - (f_normal % vec) * sin(theta) + f_normal * (vec | f_normal) * (1 - cos(theta));
    return rotVec;
}