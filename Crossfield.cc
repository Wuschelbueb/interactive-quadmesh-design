//
// Created by wuschelbueb on 15.09.21.
//

#include "Crossfield.hh"

void Crossfield::getCrossfield() {
    createCrossfields();
}

void Crossfield::createCrossfields() {
//    removeProperties();
    std::vector<int> constrainedHalfEdges = getConstraints();
    std::vector<int> faces = getReferenceEdge(constrainedHalfEdges);
    setlocalCoordFrame(faces);
    std::map<int, double> heKappa = getMapHeKappa(faces);
    CMatrixType _H = getHessianMatrix(faces, heKappa);
    RMatrixType _constraints = getConstraintMatrix(heKappa, constrainedHalfEdges, faces);
    std::vector<double> _x(heKappa.size() + faces.size(), 0.0);
    std::vector<double> _rhs = getRHS(heKappa, faces);
    std::vector<int> _idx_to_round = getIdxToRound(heKappa, faces);
    std::vector<double> _rhsOld = _rhs;
    std::vector<int> _rhsNew;
    std::vector<int> _c_elim;
    std::cout << "Matrices before solver: " << std::endl;
    std::cout << "A: " << _H << std::endl;

    for (std::size_t i = 0, max = _rhs.size(); i != max; ++i) {
        std::cout << "rhs[" << i << "] = " << _rhs[i] << std::endl;
    }
    std::cout << std::endl;
    for (std::size_t i = 0, max = _idx_to_round.size(); i != max; ++i) {
        std::cout << "idx_to_round[" << i << "] = " << _idx_to_round[i] << std::endl;
    }
    std::cout << std::endl;
    for (std::size_t i = 0, max = _x.size(); i != max; ++i) {
        std::cout << "_x[" << i << "] = " << _x[i] << std::endl;
    }
    std::cout << std::endl;
    std::cout << "constraints: " << _constraints << std::endl;

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
    csolver.solve_const(_constraints, _H, _x, _rhs, _idx_to_round);
    std::cout << "Matrices after solver: " << std::endl;
    rotateLocalCoordFrame(faces, _x);
//    std::cout << "A: " << _H << std::endl;
//    for (std::size_t i = 0, max = _rhs.size(); i != max; ++i) {
//        std::cout << "_rhs[" << i << "] = " << _rhs[i] << std::endl;
//    }
//    std::cout << std::endl;
    for (std::size_t i = 0, max = _x.size(); i != max; ++i) {
        std::cout << "_x[" << i << "] = " << _x[i] << std::endl;
    }
//    std::cout << std::endl;
//    std::cout << "constraints after: " << _constraints << std::endl;
//    std::cout << std::endl;
    double energy_after = getEnergy(_H, _x, _rhs);
//    double energy_after = getEnergy(heKappa, faces, _rhs, _x);
    std::cout << "E_smooth = " << energy_after << std::endl;

    // test examples
    const double tol = 1e-6;
    double zero = 0.;
    int n = gmm::mat_nrows(_H);
    CVectorType y(n), b(n);
    gmm::copy(_x, y);
    gmm::copy(_rhs, b);
    double ea = computeEnergy(_H, y, b);
    bool eb_ok = std::abs(ea - energy_after) < tol;
//    assert(eb_ok);
    if (!eb_ok) {
        std::cerr << "the energies do not coincide: " << energy_after
                  << " x^T A x + b^t x = " << ea << std::endl;
    }
//    bool const_check = check_constraints(_constraints, y, zero);
//    std::cout << "constraint check: " << const_check << std::endl;
}

double Crossfield::getEnergy(const CMatrixType &_H, const std::vector<double> &_x, const std::vector<double> &_rhs) {
    int n = gmm::mat_nrows(_H);
    CVectorType temp(n);
    gmm::mult(_H, _x, temp); // temp = A*x
    double p = gmm::vect_sp(_x, temp); // x^T * temp = x^T * A * x
    double q = gmm::vect_sp(_x, _rhs); // scalar product
    return (p + q);
}

gmm::row_matrix<gmm::wsvector<double>>
Crossfield::getConstraintMatrix(const std::map<int, double> &heKappa, const std::vector<int> &constraintHalfEdges,
                                const std::vector<int> &faces) {
    auto constraint_angle = OpenMesh::FProp<double>(trimesh_, "constraint_angle");
    auto positionHessianMatrix = OpenMesh::FProp<int>(trimesh_, "positionHessianMatrix");
    auto referenceHeIdx = OpenMesh::FProp<int>(trimesh_, "referenceHeIdx");
    //cNplusOne = angle between local coordinate x-axis and constraint
    int cNplusOne = 1, counter = 0, pj_start = faces.size();
    // this results in have the rows of the constrainedhalfedges
    std::vector<int> constraintEdges = getConstraintEdges(
            constraintHalfEdges), noOriginConst = getNumberOriginConstraints(faces);
    int n_row = constraintEdges.size() + noOriginConst.size(), n_col = heKappa.size() + faces.size();
    gmm::row_matrix<gmm::wsvector<double>> _constraints(n_row, n_col + cNplusOne);
    for (int i: constraintEdges) {
        // what if edge is between faces, do both faces get a 1 in their respective position in the row of the matrix
        OpenMesh::EdgeHandle eh = trimesh_.edge_handle(i);
        OpenMesh::HalfedgeHandle hehConst = trimesh_.halfedge_handle(eh, 0);
        if (trimesh_.is_boundary(hehConst)) {
            hehConst = trimesh_.opposite_halfedge_handle(hehConst);
        }
        OpenMesh::FaceHandle fh = trimesh_.face_handle(hehConst);
        double constraint = constraint_angle[fh];
        int position = positionHessianMatrix[fh];
        _constraints(counter, position) = 1.0;
        _constraints(counter, n_col) = 1.0 * constraint;
        counter++;
    }
    dualSpanningTreeConstraint(heKappa, counter, pj_start, noOriginConst, _constraints);
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

void Crossfield::dualSpanningTreeConstraint(const std::map<int, double> &heKappa, int &counter, const int pj_start,
                                            const std::vector<int> &noOriginConst,
                                            gmm::row_matrix<gmm::wsvector<double>> &_constraints) {
    auto origin_constraint = OpenMesh::FProp<int>(trimesh_, "origin_constraint");
    for (int i: noOriginConst) {
        int iteration = 0;
        for (auto j: heKappa) {
            OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(j.first);
            OpenMesh::FaceHandle fh1 = trimesh_.face_handle(heh);
            OpenMesh::FaceHandle fh2 = trimesh_.face_handle(trimesh_.opposite_halfedge_handle(heh));
            int fh1_origin = origin_constraint[fh1];
            int fh2_origin = origin_constraint[fh2];
//            std::cout << "Halfedge: " << heh.idx() << " with faces: " << fh1.idx() << " (originconst: "
//                      << origin_constraint[fh1] << ") and " << fh2.idx() << " (originconst: " << origin_constraint[fh2]
//                      << ")\n";
            if (fh1_origin == fh2_origin && fh1_origin == i) {
                _constraints(counter, pj_start + iteration++) = 1;
            } else {
                _constraints(counter, pj_start + iteration++) = 0;
            }
        }
        counter++;
    }
}

std::vector<int> Crossfield::getNumberOriginConstraints(const std::vector<int> &faces) {
    auto origin_constraint = OpenMesh::FProp<int>(trimesh_, "origin_constraint");
    std::vector<int> noOrigins;
    for (int i: faces) {
        auto fh = trimesh_.face_handle(i);
        if (std::find(noOrigins.begin(), noOrigins.end(), origin_constraint[fh]) == noOrigins.end()) {
            noOrigins.push_back(origin_constraint[fh]);
        }
    }
    return noOrigins;
}

std::vector<int> Crossfield::getConstraintEdges(const std::vector<int> &constrainedHEdges) {
    std::vector<int> edges;
    for (int i: constrainedHEdges) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i);
        OpenMesh::EdgeHandle eh = trimesh_.edge_handle(heh);
        if (std::find(edges.begin(), edges.end(), eh.idx()) == edges.end()) {
            edges.push_back(eh.idx());
        }
    };
    return edges;
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
    double totalArea = getTotalArea(faces);
    int facesPlusOne = faces.size();

    getRhsFirstHalf(totalArea, faces, _rhs, heKappa);
    getRhsSecondHalf(totalArea, _rhs, heKappa, facesPlusOne);

    // scalar multiplication with vector, b*-1 = -b
//    gmm::scale(_rhs, -1.0);
    return _rhs;
}

void Crossfield::getRhsFirstHalf(double const totalArea, const std::vector<int> &faces, std::vector<double> &_rhs,
                                 const std::map<int, double> &heKappa) {
    for (int i: faces) {
        getSum(totalArea, i, _rhs, heKappa);
    }
}

void Crossfield::getSum(double const totalArea, const int i, std::vector<double> &_rhs,
                        const std::map<int, double> &heKappa) {
    auto positionHessianMatrix = OpenMesh::FProp<int>(trimesh_, "positionHessianMatrix");
    double sum = 0;
    OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
    for (TriMesh::FaceHalfedgeIter fh_it = trimesh_.fh_iter(fh); fh_it.is_valid(); ++fh_it) {
        // this condition is necessary because an opposite face doesn't exist if the opposite heh is boundary
        sum += setSum(fh, *fh_it, heKappa, totalArea);
    }
    int position = positionHessianMatrix[fh];
    _rhs[position] = sum;
}

double
Crossfield::setSum(OpenMesh::FaceHandle fh, OpenMesh::HalfedgeHandle fh_it, const std::map<int, double> &heKappa,
                   double const totalArea) {
    double sum = 0;
    OpenMesh::HalfedgeHandle opposite_heh = trimesh_.opposite_halfedge_handle(fh_it);
    if (!trimesh_.is_boundary(opposite_heh)) {
        OpenMesh::FaceHandle opposite_fh = trimesh_.face_handle(trimesh_.opposite_halfedge_handle(fh_it));
        double edge_weight = 1;
//        edge_weight = ((trimesh_.calc_face_area(fh) / 3) + (trimesh_.calc_face_area(opposite_fh) / 3)) / totalArea;
        auto it = heKappa.find(fh_it.idx());
        auto it2 = heKappa.find(opposite_heh.idx());
        // it doesn't matter that edge_weight doesn't have a value, since this condition is never true,
        // if edge_weight is valueless
        if (it != heKappa.end()) {
            sum += (it->second * edge_weight * 2);
//            std::cout << "Kappa value main triangle: " << it->second << "\nsum main:\t" << sum << std::endl;
        }
        if (it2 != heKappa.end()) {
            sum -= (it2->second * edge_weight * 2);
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

void Crossfield::getRhsSecondHalf(double const totalArea, std::vector<double> &_rhs,
                                  const std::map<int, double> &heKappa, const int facesPlusOne) {
    int counter = facesPlusOne;
    for (const auto &i: heKappa) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i.first);
        OpenMesh::FaceHandle fh = trimesh_.face_handle(heh);
        OpenMesh::FaceHandle opposite_fh = trimesh_.face_handle(trimesh_.opposite_halfedge_handle(heh));
        double edge_weight = 1;
//        edge_weight = ((trimesh_.calc_face_area(fh) / 3) + (trimesh_.calc_face_area(opposite_fh) / 3)) / totalArea;
        _rhs[counter++] = (i.second * edge_weight);
    }
}

gmm::col_matrix<gmm::wsvector<double>>
Crossfield::getHessianMatrix(const std::vector<int> &faces, const std::map<int, double> &heKappa) {
    trimesh_.release_halfedge_status();
    trimesh_.release_face_status();
    // request to change the status
    trimesh_.request_halfedge_status();
    trimesh_.request_face_status();
    auto positionHessianMatrix = OpenMesh::FProp<int>(trimesh_, "positionHessianMatrix");
    int counter = 0, iteration = 0, pj_start = faces.size(), n = heKappa.size() + faces.size();
    double edge_weight = 1, totalArea = getTotalArea(faces);
    gmm::col_matrix<gmm::wsvector<double>> _H(n, n);
    // indexes the faces from 0 to n
    // this is needed in case a face index is higher than the matrix size
    setPositioninHessianForFaces(heKappa);

    // fills up sparse column matrix _H
    for (auto i: heKappa) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i.first);
        OpenMesh::FaceHandle fh1 = trimesh_.face_handle(heh);
        OpenMesh::FaceHandle fh2 = trimesh_.face_handle(trimesh_.opposite_halfedge_handle(heh));
        int factorI = getFactor(fh1);
        int factorJ = getFactor(fh2);
        std::cout << "face: " << fh1.idx() << " with factor: " << factorI << " and face: " << fh2.idx()
                  << " with factor: " << factorJ << std::endl;
//        edge_weight = ((trimesh_.calc_face_area(fh1) / 3) + (trimesh_.calc_face_area(fh2) / 3)) / totalArea;
        int pos_i = positionHessianMatrix[fh1];
        int pos_j = positionHessianMatrix[fh2];
        // |  2*w_ij   -2*w_ij     pi*w_ij |
        // | -2*w_ij    2*w_ij    -pi*w_ij |
        // | pi*w_ij  -pi*w_ij pi^2/2*w_ij |
        _H(pos_i, pos_i) = 2.0 * factorI * edge_weight;
        _H(pos_j, pos_j) = 2.0 * factorJ * edge_weight;
        _H(pos_i, pos_j) = -2.0 * edge_weight;
        _H(pos_j, pos_i) = -2.0 * edge_weight;
        _H(pos_i, pj_start + iteration) = M_PI * edge_weight;
        _H(pj_start + iteration, pos_i) = M_PI * edge_weight;
        _H(pos_j, pj_start + iteration) = -M_PI * edge_weight;
        _H(pj_start + iteration, pos_j) = -M_PI * edge_weight;
        _H(pj_start + iteration, pj_start + iteration) = (pow(M_PI, 2) / 2) * edge_weight;
        iteration++;
    }
    gmm::clean(_H, 1E-10);

    return _H;
}

int Crossfield::getFactor(const OpenMesh::FaceHandle fh) {
    int factor = 0;
    for (TriMesh::FaceFaceIter ff_it = trimesh_.ff_iter(fh); ff_it.is_valid(); ++ff_it) {
        factor++;
    }
    return factor;
}

void Crossfield::setPositioninHessianForFaces(const std::map<int, double> &heKappa) {
    auto positionHessianMatrix = OpenMesh::FProp<int>(trimesh_, "positionHessianMatrix");
    int counter = 0;
    for (const auto &i: heKappa) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i.first);
        OpenMesh::FaceHandle fh1 = trimesh_.face_handle(heh);
        OpenMesh::FaceHandle fh2 = trimesh_.face_handle(trimesh_.opposite_halfedge_handle(heh));
        if (!trimesh_.status(fh1).tagged()) {
            positionHessianMatrix[fh1] = counter;
            trimesh_.status(fh1).set_tagged(true);
            counter++;
        }
        if (!trimesh_.status(fh2).tagged()) {
            positionHessianMatrix[fh2] = counter;
            trimesh_.status(fh2).set_tagged(true);
            counter++;
        }
    }
}

std::map<int, double> Crossfield::getMapHeKappa(const std::vector<int> &faces) {
    std::map<int, double> heKappa;
    for (int i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        for (TriMesh::FaceFaceIter ff_it = trimesh_.ff_iter(fh); ff_it.is_valid(); ++ff_it) {
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
    // neighbour needs to be in faces
    if (trimesh_.status(fh_neigh).tagged()) {
        // get index of ref edges
        int refEdgeMain = referenceHeIdx[fh];
        int refEdgeNeigh = referenceHeIdx[fh_neigh];
        // get the common edge between the triangles
        commonEdge = getCommonEdgeBetweenTriangles(fh, fh_neigh, refEdgeMain, refEdgeNeigh);
        kappa = getKappa(refEdgeMain, refEdgeNeigh, commonEdge);

        // example test kappa
        bool kappa_ok = testKappa(refEdgeMain, refEdgeNeigh, commonEdge);
//        assert(kappa_ok);
        addKappaHeToMap(commonEdge, kappa, heKappa);
    }
}

std::pair<int, int>
Crossfield::getCommonEdgeBetweenTriangles(const OpenMesh::FaceHandle fh, const OpenMesh::FaceHandle fh_neigh,
                                          const int &refEdgeMain, const int &refEdgeNeigh) {
    std::pair<int, int> commonEdge = {INT_MAX, INT_MAX};
    for (TriMesh::FaceHalfedgeIter fhe_it_main = trimesh_.fh_iter(fh); fhe_it_main.is_valid(); ++fhe_it_main) {
        for (TriMesh::FaceHalfedgeIter fhe_it_neigh = trimesh_.fh_iter(
                fh_neigh); fhe_it_neigh.is_valid(); ++fhe_it_neigh) {
            if (trimesh_.opposite_halfedge_handle(*fhe_it_main).idx() == fhe_it_neigh->idx()) {
                commonEdge.first = fhe_it_main->idx();
                commonEdge.second = fhe_it_neigh->idx();
                // as soon as the common edge is found, the loop is obsolete
                goto END_LOOP;
            }
        }
    }
    END_LOOP:
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
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(tempHe);
        alpha += trimesh_.calc_sector_angle(heh);
        tempHe = trimesh_.next_halfedge_handle(heh).idx();
        counter++;
    }
    tempHe = commonEdge.second;
    while (tempHe != refEdgeNeigh) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(tempHe);
        beta += trimesh_.calc_sector_angle(heh);
        tempHe = trimesh_.next_halfedge_handle(heh).idx();
        counter++;
    }
    kappa = counter * M_PI - alpha - beta;
    // make sure kappa is in range of [-pi, pi)
    double kappa_final = shortenKappa(kappa);
//    std::cout << "kappa: " << kappa_final << std::endl;
    return kappa_final;
}

void Crossfield::addKappaHeToMap(const std::pair<int, int> commonEdge, const double kappa,
                                 std::map<int, double> &heKappa) {
    // if the opposite he isn't already in the list, add this he to the list
    if (heKappa.find(commonEdge.second) == heKappa.end()) {
        heKappa[commonEdge.first] = kappa;
    }
}

std::vector<int> Crossfield::getReferenceEdge(const std::vector<int> &constrainedHEdges) {
    // status needs to be released before using, in case it still has saved some stuff from before
    trimesh_.release_halfedge_status();
    trimesh_.release_face_status();
    // request to change the status
    trimesh_.request_halfedge_status();
    trimesh_.request_face_status();
    std::vector<int> faces;
    int temp = INT_MAX;

    // assign constraint edges to their faces
    for (int i: constrainedHEdges) {
        setRefHeWithConstraint(i, temp, faces);
    }

    // assign reference edges in range to faces
    for (int i: heInRange_) {
        setRefHeWithoutConstraint(i, temp, faces);

    }
    colorHEdges(constrainedHEdges);
    colorFaces(faces);
    return faces;
}

void Crossfield::setRefHeWithConstraint(const int i, int &temp, std::vector<int> &faces) {
    auto referenceHeIdx = OpenMesh::FProp<int>(trimesh_, "referenceHeIdx");
    OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i);
    if (trimesh_.is_boundary(heh)) {
        heh = trimesh_.opposite_halfedge_handle(heh);
    }
    OpenMesh::FaceHandle fh = trimesh_.face_handle(heh);
    if (!trimesh_.status(heh).tagged() && !trimesh_.status(fh).tagged()) {
        // add reference edge to face
        referenceHeIdx[fh] = heh.idx();
        temp = heh.idx();
        // set face and Hedge as used
        trimesh_.status(heh).set_tagged(true);
        trimesh_.status(trimesh_.opposite_halfedge_handle(heh)).set_tagged(true);
        trimesh_.status(fh).set_tagged(true);
        faces.push_back(fh.idx());
    }
}

void Crossfield::setRefHeWithoutConstraint(const int i, int &temp, std::vector<int> &faces) {
    auto referenceHeIdx = OpenMesh::FProp<int>(trimesh_, "referenceHeIdx");
    auto heh_color = OpenMesh::HProp<int>(trimesh_, "heh_color");
    OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i);
    if (trimesh_.is_boundary(heh))
        heh = trimesh_.opposite_halfedge_handle(heh);
    heh_color(heh) = 2;
    OpenMesh::FaceHandle fh = trimesh_.face_handle(heh);
    // avoids including faces which only have one halfedge in heInRange
    OpenMesh::HalfedgeHandle nheh = trimesh_.next_halfedge_handle(heh);
    // halfedges or faces which are already tagged, won't be used anymore, no duplicates are allowed
    if ((std::find(heInRange_.begin(), heInRange_.end(), nheh.idx()) != heInRange_.end()) &&
        (!trimesh_.status(heh).tagged() && !trimesh_.status(fh).tagged())) {
        referenceHeIdx[fh] = heh.idx();
        // both halfedges need to be tagged, else it is possible opposite faces share an edge
        trimesh_.status(heh).set_tagged(true);
        trimesh_.status(trimesh_.opposite_halfedge_handle(heh)).set_tagged(true);
        trimesh_.status(fh).set_tagged(true);
        faces.push_back(fh.idx());
    }
}


// convert array of edges to faces and fill constraints with them
std::vector<int> Crossfield::getConstraints() {
    std::vector<int> constraints;
    getSelectedVertices(constraints);
    getSelectedEdges(constraints);
    getSelectedHEdges(constraints);
    getSelectedFaces(constraints);
    return constraints;
}

void Crossfield::getSelectedVertices(std::vector<int> &constraints) {
    std::vector<int> selection = MeshSelection::getVertexSelection(&trimesh_);
    // if vertex not empty checks neighbouring vertices. if they are also selected, add edge between
    for (int i: selection) {
        OpenMesh::VertexHandle vh = trimesh_.vertex_handle(i);
        for (int j: selection) {
            OpenMesh::VertexHandle vhj = trimesh_.vertex_handle(j);
            OpenMesh::HalfedgeHandle heh1 = trimesh_.find_halfedge(vhj, vh);
            OpenMesh::HalfedgeHandle heh2 = trimesh_.opposite_halfedge_handle(heh1);
            // avoids duplicates with std::find
            if (std::find(constraints.begin(), constraints.end(), heh1.idx()) ==
                constraints.end()) {
                constraints.push_back(heh1.idx());
            }
            if (std::find(constraints.begin(), constraints.end(), heh2.idx()) ==
                constraints.end()) {
                constraints.push_back(heh2.idx());
            }
        }
    }
}

void Crossfield::getSelectedEdges(std::vector<int> &constraints) {
    std::vector<int> selection = MeshSelection::getEdgeSelection(&trimesh_);
    for (int i: selection) {
        // avoids duplicates with std::find
        OpenMesh::EdgeHandle eh = trimesh_.edge_handle(i);
        OpenMesh::HalfedgeHandle heh1 = trimesh_.halfedge_handle(eh, 0);
        OpenMesh::HalfedgeHandle heh2 = trimesh_.halfedge_handle(eh, 1);
        if (std::find(constraints.begin(), constraints.end(), heh1.idx()) == constraints.end()) {
            constraints.push_back(heh1.idx());
        }
        if (std::find(constraints.begin(), constraints.end(), heh2.idx()) == constraints.end()) {
            constraints.push_back(heh2.idx());
        }
    }
}

void Crossfield::getSelectedHEdges(std::vector<int> &constraints) {
    std::vector<int> selection = MeshSelection::getHalfedgeSelection(&trimesh_);
    for (int i: selection) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i);
        // avoids duplicates with std::find
        if (std::find(constraints.begin(), constraints.end(), heh.idx()) == constraints.end())
            constraints.push_back(heh.idx());
    }
}

void Crossfield::getSelectedFaces(std::vector<int> &constraints) {
    std::vector<int> selection = MeshSelection::getFaceSelection(&trimesh_);
    for (int i: selection) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        for (TriMesh::FaceHalfedgeIter fh_it = trimesh_.fh_iter(fh); fh_it.is_valid(); ++fh_it) {
            // avoids duplicates with std::find
            if (std::find(constraints.begin(), constraints.end(), fh_it->idx()) == constraints.end())
                constraints.push_back(fh_it->idx());
        }
    }
}

void Crossfield::setlocalCoordFrame(const std::vector<int> &faces) {
    auto referenceHeIdx = OpenMesh::FProp<int>(trimesh_, "referenceHeIdx");
    auto constraint_angle = OpenMesh::FProp<double>(trimesh_, "constraint_angle");
    auto uVectorField = OpenMesh::FProp<Point>(trimesh_, "uVectorField");
    auto vVectorField = OpenMesh::FProp<Point>(trimesh_, "vVectorField");
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
        vVectorField[fh] = v;
    }
    std::cout << std::endl;
}

void Crossfield::rotateLocalCoordFrame(const std::vector<int> &faces, const std::vector<double> _x) {
    auto positionHessianMatrix = OpenMesh::FProp<int>(trimesh_, "positionHessianMatrix");
    auto uVectorField = OpenMesh::FProp<Point>(trimesh_, "uVectorField");
    auto vVectorField = OpenMesh::FProp<Point>(trimesh_, "vVectorField");
    auto uVectorFieldRot = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRot");
    auto vVectorFieldRot = OpenMesh::FProp<Point>(trimesh_, "vVectorFieldRot");
    for (int i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        int position = positionHessianMatrix[fh];
        Point p_x = uVectorField[fh];
        Point p_y = vVectorField[fh];
        double radians = shortenKappa(_x[position]);
//        std::cout << "face (" << i << ") with (" << _x[position] << "): with angle: " << radians << "(rad) and "
//                  << radians * 180 / M_PI << "(deg)"
//                  << std::endl;
        uVectorFieldRot[fh] = p_x * cos(radians) - p_y * sin(radians);
        vVectorFieldRot[fh] = p_x * sin(radians) + p_y * cos(radians);
    }
}

double Crossfield::getTotalArea(const std::vector<int> &faces) {
    double totalarea = 0.0;
    for (int i: faces)
        totalarea += trimesh_.calc_face_area(trimesh_.face_handle(i));
    return totalarea;
}


void Crossfield::colorFaces(const std::vector<int> &faces) {
    auto face_color = OpenMesh::FProp<int>(trimesh_, "face_color");
    for (auto f_it = trimesh_.faces_begin(); f_it != trimesh_.faces_end(); ++f_it) {
        face_color[*f_it] = 0;
    }
    for (int i: faces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        face_color[fh] = 1;
    }
}

void Crossfield::colorHEdges(const std::vector<int> &constrainedEdges) {
    auto heh_color = OpenMesh::HProp<int>(trimesh_, "heh_color");
    for (int i: heInRange_) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i);
        heh_color(heh) = 1;
    }
    for (int i: constrainedEdges) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i);
        heh_color(heh) = 2;
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
