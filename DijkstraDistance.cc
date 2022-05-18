#include "DijkstraDistance.hh"

void DijkstraDistance::getDijkstraSingularities(std::vector<int> &complementHEdges, std::vector<int> &singularities) {
    std::vector<int> dualGraphVertices = createVerticesVector(complementHEdges, singularities);
    initVertexProp(dualGraphVertices, true);
    for (int i: singularities) {
        calculateVDijkstra(i);
        addPathToCutGraph(dualGraphVertices, i);
        initVertexProp(dualGraphVertices, false);
    }
}

std::vector<int>
DijkstraDistance::createVerticesVector(std::vector<int> &complementHEdges, std::vector<int> &singularities) {
    std::vector<int> dualGraphVertices;
    for (int i: complementHEdges) {
        auto he = trimesh_.halfedge_handle(i);
        auto vht = trimesh_.to_vertex_handle(he);
        auto vhf = trimesh_.from_vertex_handle(he);
        if (std::find(dualGraphVertices.begin(), dualGraphVertices.end(), vhf.idx()) == dualGraphVertices.end()) {
            dualGraphVertices.push_back(vhf.idx());
        }
        if (std::find(dualGraphVertices.begin(), dualGraphVertices.end(), vht.idx()) == dualGraphVertices.end()) {
            dualGraphVertices.push_back(vht.idx());
        }
    }
    if (complementHEdges.empty()) {
        dualGraphVertices.push_back(singularities[0]);
        singularities.erase(singularities.begin());
    }
    return dualGraphVertices;
}

void DijkstraDistance::initVertexProp(std::vector<int> &dualGraphVertices, const bool flag) {
    auto vertexDist = OpenMesh::VProp<double>(trimesh_, "vertexDist");
    auto vertexOrigin = OpenMesh::VProp<int>(trimesh_, "vertexOrigin");
    auto vertexPredecessor = OpenMesh::VProp<int>(trimesh_, "vertexPredecessor");
    int max = INT_MAX, zeroDist = 0.0;
    for (auto vh: trimesh_.vertices()) {
        vertexDist[vh] = max;
        if (flag) {
            vertexOrigin[vh] = max;
            vertexPredecessor[vh] = max;
        }
    }
    for (auto i: dualGraphVertices) {
        auto vh = trimesh_.vertex_handle(i);
        vertexDist[vh] = zeroDist;
        if (flag) {
            vertexOrigin[vh] = vh.idx();
            vertexPredecessor[vh] = vh.idx();
        }
    }
}

int DijkstraDistance::vertexGetSmallestDist() {
    auto vertexDist = OpenMesh::VProp<double>(trimesh_, "vertexDist");
    double minDistance = DBL_MAX;
    int idx = INT_MAX;
    for (auto vh: trimesh_.vertices()) {
        if (!trimesh_.status(vh).tagged() && vertexDist[vh] < minDistance) {
            minDistance = vertexDist[vh];
            idx = vh.idx();
        }
    }
    return idx;
}

void DijkstraDistance::calculateVDijkstra(const int i) {
    trimesh_.release_vertex_status();
    trimesh_.request_vertex_status();
    auto vertexDist = OpenMesh::VProp<double>(trimesh_, "vertexDist");
    auto vertexOrigin = OpenMesh::VProp<int>(trimesh_, "vertexOrigin");
    auto vertexPredecessor = OpenMesh::VProp<int>(trimesh_, "vertexPredecessor");
    auto vhSingularity = trimesh_.vertex_handle(i);
    while (true) {
        int vertexIdx = vertexGetSmallestDist();
        // if destination vertex is visited
        if (trimesh_.status(vhSingularity).tagged()) {
            break;
        }
        auto vh_og = trimesh_.vertex_handle(vertexIdx);
        trimesh_.status(vh_og).set_tagged(true);
        for (TriMesh::VertexOHalfedgeIter vhe_it = trimesh_.voh_iter(vh_og); vhe_it.is_valid(); ++vhe_it) {
            auto vh_next = trimesh_.to_vertex_handle(*vhe_it);
            double distance = vertexDist[vh_og] + 1.0;
            if (!trimesh_.status(vh_next).tagged()
                && distance < vertexDist[vh_next]) {
                vertexDist[vh_next] = distance;
                vertexOrigin[vh_next] = vertexOrigin[vh_og];
                vertexPredecessor[vh_next] = vh_og.idx();
            }
        }
    }
}

void DijkstraDistance::addPathToCutGraph(std::vector<int> &dualGraphVertices, const int i) {
    auto vertexPredecessor = OpenMesh::VProp<int>(trimesh_, "vertexPredecessor");
    auto vertexOrigin = OpenMesh::VProp<int>(trimesh_, "vertexOrigin");
    auto vh = trimesh_.vertex_handle(i);
    bool flag = true;
    while (flag) {
        if (std::find(dualGraphVertices.begin(), dualGraphVertices.end(), vh.idx()) == dualGraphVertices.end()) {
            dualGraphVertices.push_back(vh.idx());
            vh = trimesh_.vertex_handle(vertexPredecessor[vh]);
        } else {
            flag = false;
        }
    }
}

void DijkstraDistance::getDualGraph(const std::vector<int> &faces) {
    trimesh_.release_face_status();
    trimesh_.request_face_status();
    double initValue = INT_MAX, zeroDist = 0.0;
    auto dualGraphDist = OpenMesh::FProp<double>(trimesh_, "dualGraphDist");
    auto dualGraphOrigin = OpenMesh::FProp<int>(trimesh_, "dualGraphOrigin");
    auto dualGraphPred = OpenMesh::FProp<int>(trimesh_, "dualGraphPred");
    for (auto it = std::begin(faces), first = it, end = std::end(faces); it != end; ++it) {
        auto fh = trimesh_.face_handle(*it);
        dualGraphDist[fh] = initValue;
        dualGraphOrigin[fh] = initValue;
        dualGraphPred[fh] = initValue;
        if (it == first) {
            dualGraphDist[fh] = zeroDist;
            dualGraphOrigin[fh] = *it;
            dualGraphPred[fh] = *it;
        }
    }
    while (true) {
        double distance = 0.0;
        int faceIdx = dualGraphGetSmallestDist(faces);
        // if all vertices are visited the algo stops
        if (faceIdx == INT_MAX)
            break;
        OpenMesh::FaceHandle fh = trimesh_.face_handle(faceIdx);
        trimesh_.status(fh).set_tagged(true);
        for (TriMesh::FaceFaceIter ff_it = trimesh_.ff_iter(fh); ff_it.is_valid(); ++ff_it) {
            Point origin = trimesh_.calc_face_centroid(fh);
            Point neighbour = trimesh_.calc_face_centroid(*ff_it);
            Point difference = neighbour - origin;
            distance = dualGraphDist[fh] + difference.norm();
            if (!trimesh_.status(*ff_it).tagged()
                && distance < dualGraphDist[*ff_it]) {
                dualGraphDist[*ff_it] = distance;
                dualGraphOrigin[*ff_it] = dualGraphOrigin[fh];
                dualGraphPred[*ff_it] = fh.idx();
            }
        }
    }
}

int DijkstraDistance::dualGraphGetSmallestDist(const std::vector<int> &faces) {
    auto dualGraphDist = OpenMesh::FProp<double>(trimesh_, "dualGraphDist");
    double minDistance = DBL_MAX;
    int idx = INT_MAX;
    for (int i: faces) {
        auto fh = trimesh_.face_handle(i);
        if (!trimesh_.status(fh).tagged() && dualGraphDist[fh] < minDistance) {
            minDistance = dualGraphDist[fh];
            idx = i;
        }
    }
    return idx;
}

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
    auto periodJump = OpenMesh::HProp<int>(trimesh_, "periodJump");
    double minDistance = INT_MAX, zeroDist = 0.0;
    std::vector<int> constraintFaces;

    for (auto he: trimesh_.halfedges()) {
        periodJump[he] = minDistance;
    }

    for (auto fh: trimesh_.faces()) {
        origin_constraint[fh] = minDistance;
        distanceBaryCenter[fh] = minDistance;
        predecessor_face[fh] = minDistance;
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