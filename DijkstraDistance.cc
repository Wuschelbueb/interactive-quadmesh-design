#include "DijkstraDistance.hh"

std::vector<int>
DijkstraDistance::calculateDijkstra(const std::vector<int> HeConstraints, const double refDist,
                                    const bool includeBoundary) {
    auto distanceBaryCenter = OpenMesh::FProp<double>(trimesh_, "distanceBaryCenter");
    std::vector<int> includedNodes;
    std::vector<int> constraintFaces = transformHehToFaces(HeConstraints);
    dijkstraDistBaryCenter(includedNodes, refDist);
    if (includeBoundary) {
        includeBoundaryFaces(includedNodes, refDist);
    }
    return includedNodes;
}

void DijkstraDistance::dijkstraDistBaryCenter(std::vector<int> &includedNodes, const double refDist) {
    trimesh_.release_face_status();
    trimesh_.request_face_status();
    auto distanceBaryCenter = OpenMesh::FProp<double>(trimesh_, "distanceBaryCenter");
    auto origin_constraint = OpenMesh::FProp<int>(trimesh_, "origin_constraint");
    auto predecessor_face = OpenMesh::FProp<int>(trimesh_, "predecessor_face");
    while (true) {
        double distance = 0.0;
        int faceIdx = getSmallestDistProp(refDist);
        // if all vertices are visited the algo stops
        if (faceIdx == INT_MAX)
            break;
        OpenMesh::FaceHandle fh = trimesh_.face_handle(faceIdx);
        trimesh_.status(fh).set_tagged(true);
        includedNodes.push_back(fh.idx());
        for (TriMesh::FaceFaceIter ff_it = trimesh_.ff_iter(fh); ff_it.is_valid(); ++ff_it) {
            Point origin = trimesh_.calc_face_centroid(fh);
            Point neighbour = trimesh_.calc_face_centroid(*ff_it);
            Point difference = neighbour - origin;
            distance = distanceBaryCenter[fh] + difference.norm();
            if (!trimesh_.status(*ff_it).tagged()
                && distance < distanceBaryCenter[*ff_it]) {
                distanceBaryCenter[*ff_it] = distance;
                origin_constraint[*ff_it] = origin_constraint[fh];
                predecessor_face[*ff_it] = fh.idx();
            }
        }
    }
}

std::vector<int> DijkstraDistance::transformHehToFaces(const std::vector<int> &constraintHeh) {
    auto origin_constraint = OpenMesh::FProp<int>(trimesh_, "origin_constraint");
    auto distanceBaryCenter = OpenMesh::FProp<double>(trimesh_, "distanceBaryCenter");
    auto predecessor_face = OpenMesh::FProp<int>(trimesh_, "predecessor_face");
    auto positionHessianMatrix = OpenMesh::FProp<int>(trimesh_, "positionHessianMatrix");
    double minDistance = DBL_MAX, zeroDist = 0.0;
    std::vector<int> constraintFaces;
    for (auto fh: trimesh_.faces()) {
        distanceBaryCenter[fh] = minDistance;
        positionHessianMatrix[fh] = -1;
    }
    for (int i: constraintHeh) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i);
        if (!trimesh_.is_boundary(heh)) {
            OpenMesh::FaceHandle fh = trimesh_.face_handle(heh);
            constraintFaces.push_back(fh.idx());
            origin_constraint[fh] = fh.idx();
            distanceBaryCenter[fh] = zeroDist;
            predecessor_face[fh] = fh.idx();
        }
    }
    return constraintFaces;
}

void DijkstraDistance::includeBoundaryFaces(std::vector<int> &includedFaces, const double refDist) {
    auto distanceBaryCenter = OpenMesh::FProp<double>(trimesh_, "distanceBaryCenter");
    std::vector<int> tempFaces;
    for (int i: includedFaces) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        for (TriMesh::FaceFaceIter ff_it = trimesh_.ff_iter(fh); ff_it.is_valid(); ++ff_it) {
            if (distanceBaryCenter[*ff_it] != DBL_MAX
                && distanceBaryCenter[*ff_it] >= refDist
                && (std::find(includedFaces.begin(), includedFaces.end(), ff_it->idx()) == includedFaces.end())) {
                includedFaces.push_back(ff_it->idx());
            }
        }
    }
}

std::vector<int>
DijkstraDistance::getHeVectorOfSelection(const std::vector<int> &includedNodes) {
    std::vector<int> includedHe;
    for (int i: includedNodes) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        for (TriMesh::FaceHalfedgeIter fh_it = trimesh_.fh_iter(fh); fh_it.is_valid(); ++fh_it) {
            includedHe.push_back(fh_it->idx());
        }

    }
    return includedHe;
}

void DijkstraDistance::colorizeEdges(const std::vector<int> &includedHEdges) {
    // request color change
    trimesh_.request_edge_colors();
    // define colors
    TriMesh::Color green = {0, 1, 0, 1};
    TriMesh::Color white = {1, 1, 1, 1};
    // colorize all edges white
    for (OpenMesh::EdgeHandle eh: trimesh_.edges()) {
        trimesh_.set_color(eh, white);
    }

    // colorize edges where vertices have a smaller distance than the refDist blue
    // and edges where the refDist is bigger but some vertices of the face are smaller than refDist green
    for (int i: includedHEdges) {
        OpenMesh::HalfedgeHandle ehh = trimesh_.halfedge_handle(i);
        trimesh_.set_color(trimesh_.edge_handle(ehh), green);
    }
}

// check every vertex and return vertex with the smallest "distance" property which is still unvisited
int DijkstraDistance::getSmallestDistProp(const double refDist) {
    auto distanceBaryCenter = OpenMesh::FProp<double>(trimesh_, "distanceBaryCenter");
    double minDistance = DBL_MAX;
    int idx = INT_MAX;
    for (auto fh: trimesh_.faces()) {
        if (!trimesh_.status(fh).tagged() && distanceBaryCenter[fh] < refDist && distanceBaryCenter[fh] < minDistance) {
            minDistance = distanceBaryCenter[fh];
            idx = fh.idx();
        }
    }
    return idx;
}

std::vector<int> DijkstraDistance::getHeConstraints() {
    std::vector<int> constraints;
    getSelectedEdges(constraints);
    getSelectedHEdges(constraints);
    getSelectedFaces(constraints);
    return constraints;
}

void DijkstraDistance::getSelectedEdges(std::vector<int> &constraints) {
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

void DijkstraDistance::getSelectedHEdges(std::vector<int> &constraints) {
    std::vector<int> selection = MeshSelection::getHalfedgeSelection(&trimesh_);
    for (int i: selection) {
        OpenMesh::HalfedgeHandle heh = trimesh_.halfedge_handle(i);
        // avoids duplicates with std::find
        if (std::find(constraints.begin(), constraints.end(), heh.idx()) == constraints.end())
            constraints.push_back(heh.idx());
    }
}

void DijkstraDistance::getSelectedFaces(std::vector<int> &constraints) {
    std::vector<int> selection = MeshSelection::getFaceSelection(&trimesh_);
    for (int i: selection) {
        OpenMesh::FaceHandle fh = trimesh_.face_handle(i);
        for (TriMesh::FaceHalfedgeIter fh_it = trimesh_.fh_iter(fh); fh_it.is_valid(); ++fh_it) {
            // avoids duplicates with std::find
            if (std::find(constraints.begin(), constraints.end(), fh_it->idx()) == constraints.end()) {
                constraints.push_back(fh_it->idx());
            }
        }
    }
}