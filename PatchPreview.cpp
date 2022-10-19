//
// Created by wuschelbueb on 06.10.22.
//

#include "PatchPreview.hh"

void PatchPreview::getCurvature() {
    auto minCurvature = OpenMesh::VProp<OpenMesh::Vec3d>(trimesh_, "minCurvature");
    auto maxCurvature = OpenMesh::VProp<OpenMesh::Vec3d>(trimesh_, "maxCurvature");
    auto normal_approx = OpenMesh::VProp<OpenMesh::Vec3d>(trimesh_, "normal_approx");
    for (auto vh: trimesh_.vertices()) {
        minCurvature[vh] = {0, 0, 0};
        maxCurvature[vh] = {0, 0, 0};
        normal_approx[vh] = {0, 0, 0};
    }
    //for loop over all vertices is only here for testing
//    OpenMesh::VertexHandle vh = getSelectedVertex();
//    OpenMesh::VertexHandle vh = trimesh_.vertex_handle(6);
    for (auto vh: trimesh_.vertices()) {
        if (!trimesh_.is_boundary(vh)) {
            std::vector<OpenMesh::FaceHandle> twoRingFaces;
            get2RingFaceNeighbours(twoRingFaces, vh);
            std::vector<OpenMesh::HalfedgeHandle> twoRingHe = getTwoRingHe(twoRingFaces);
            double totalArea = calculateArea(twoRingFaces);
            Eigen::Matrix3d discreteTensor = Eigen::Matrix3d::Zero();
            for (OpenMesh::HalfedgeHandle he: twoRingHe) {
                // adapt so we can cast trimesh edge vector directly to eigen vector
                OpenMesh::Vec3d edgeVec = trimesh_.calc_edge_vector(he);
                Eigen::Vector3d edgeVecEigen = {edgeVec[0], edgeVec[1], edgeVec[2]};
                double dihedral_angle = trimesh_.calc_dihedral_angle(he);
                discreteTensor += dihedral_angle * edgeVecEigen * edgeVecEigen.normalized().transpose();
            }
            discreteTensor /= totalArea;
            Eigen::EigenSolver<Eigen::Matrix3d> es(discreteTensor);
            std::vector<std::pair<int, double>> IdxAndEigenVal;
            //range-for loop
            for (auto &entry: es.eigenvalues()) {
                auto idx = &entry - &es.eigenvalues()[0];
                IdxAndEigenVal.emplace_back(idx, std::abs(entry.real()));
            }
            //sort by second element of pairs
            std::sort(IdxAndEigenVal.begin(), IdxAndEigenVal.end(),
                      [](const std::pair<int, double> &left, const std::pair<int, double> &right) {
                          return left.second < right.second;
                      });
            // use eigenvalues to extract position of normal, kmin, kmax
            int approx_normal = IdxAndEigenVal[0].first;
            int kmin = IdxAndEigenVal[1].first;
            int kmax = IdxAndEigenVal[2].first;
            minCurvature[vh] = {std::real(es.eigenvectors()(0, kmax)), std::real(es.eigenvectors()(1, kmax)),
                                std::real(es.eigenvectors()(2, kmax))};
            normal_approx[vh] = {std::real(es.eigenvectors()(0, approx_normal)), std::real(es.eigenvectors()(1, approx_normal)),
                                 std::real(es.eigenvectors()(2, approx_normal))};
            maxCurvature[vh] = {std::real(es.eigenvectors()(0, kmin)), std::real(es.eigenvectors()(1, kmin)),
                                std::real(es.eigenvectors()(2, kmin))};
        }
    }

}

void PatchPreview::get2RingFaceNeighbours(std::vector<OpenMesh::FaceHandle> &twoRingFaces, OpenMesh::VertexHandle vh) {
    trimesh_.release_face_status();
    trimesh_.request_face_status();
    auto twoRing = OpenMesh::FProp<bool>(trimesh_, "twoRing");
    for (auto fh: trimesh_.faces()) {
        twoRing[fh] = false;
    }
    std::vector<OpenMesh::SmartHalfedgeHandle> oneRingHe;
    OpenMesh::SmartVertexHandle centerVh = make_smart(vh, trimesh_);

    // add first ring neighbors
    for (TriMesh::VertexOHalfedgeIter voh_it = trimesh_.voh_begin(centerVh); voh_it.is_valid(); ++voh_it) {
        oneRingHe.push_back(*voh_it);
        twoRingFaces.push_back(voh_it->face());
        twoRing[voh_it->face()] = true;
        trimesh_.status(voh_it->face()).set_tagged(true);

    }
    // add second ring neighbors
    for (OpenMesh::SmartHalfedgeHandle elementOneRing: oneRingHe) {
        for (TriMesh::VertexOHalfedgeIter voh_it = trimesh_.voh_begin(
                elementOneRing.to()); voh_it.is_valid(); ++voh_it) {
            if (!trimesh_.status(voh_it->face()).tagged()) {
                twoRingFaces.push_back(voh_it->face());
                twoRing[voh_it->face()] = true;
                trimesh_.status(voh_it->face()).set_tagged(true);
            }
        }
    }
}

double PatchPreview::calculateArea(std::vector<OpenMesh::FaceHandle> twoRingfaces) {
    double totalArea = 0;
    for (OpenMesh::FaceHandle face: twoRingfaces) {
        totalArea += trimesh_.calc_face_area(face);
    }
    return totalArea;
}

std::vector<OpenMesh::HalfedgeHandle> PatchPreview::getTwoRingHe(std::vector<OpenMesh::FaceHandle> twoRingFaces) {
    trimesh_.release_halfedge_status();
    trimesh_.request_halfedge_status();
    auto twoRingHe = OpenMesh::HProp<bool>(trimesh_, "twoRingHe");
    for (auto he: trimesh_.halfedges()) {
        twoRingHe[he] = false;
    }
    std::vector<OpenMesh::HalfedgeHandle> twoRingUniqueHe;
    for (OpenMesh::FaceHandle face: twoRingFaces) {
        for (TriMesh::FaceHalfedgeIter fhe_it = trimesh_.fh_begin(face); fhe_it.is_valid(); ++fhe_it) {
            if (!trimesh_.status(*fhe_it).tagged() && !trimesh_.status(fhe_it->opp()).tagged() &&
                trimesh_.status(fhe_it->opp().face()).tagged()) {
                trimesh_.status(*fhe_it).set_tagged(true);
                trimesh_.status(fhe_it->opp()).set_tagged(true);
                twoRingUniqueHe.push_back(*fhe_it);
                twoRingHe[*fhe_it] = true;
            }
        }
    }
    return twoRingUniqueHe;
}

OpenMesh::VertexHandle PatchPreview::getSelectedVertex() {
    OpenMesh::VertexHandle centerVh;
    std::vector<int> vertexSelection;
    try {
        vertexSelection = MeshSelection::getVertexSelection(&trimesh_);
        if (vertexSelection.empty()) {
            throw std::invalid_argument("No Vertex selected! Select One!");
        } else if (vertexSelection.size() > 1) {
            throw std::invalid_argument("Select ONLY 1 vertex!");
        }
    }
    catch (const std::invalid_argument &e) {
        std::cerr << e.what() << std::endl;
    }
    return centerVh;
}

