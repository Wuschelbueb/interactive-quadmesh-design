#include "DijkstraDistance.hh"

void
DijkstraDistance::calcDijkstraWSingularities(std::vector<int> &complementHEdges, std::vector<int> &singularities,
                                             std::vector<int> &cutGraphWoBoundary) {
    trimesh_.request_vertex_status();
    for (auto vh: trimesh_.vertices()) {
        trimesh_.status(vh).set_tagged(false);
        trimesh_.status(vh).set_tagged2(false);
    }
    std::vector<int> cutGraphVertices = createVerticesVector(complementHEdges, singularities);
    initProperties(cutGraphVertices, true);
    for (int i: singularities) {
        if (!trimesh_.status(trimesh_.vertex_handle(i)).tagged2()) {
            calculateVDijkstra(i);
            addPathToCutGraph(cutGraphVertices, complementHEdges, i, cutGraphWoBoundary);
            initProperties(cutGraphVertices, false);
        }
    }
    trimesh_.release_vertex_status();

}

std::vector<int>
DijkstraDistance::createVerticesVector(std::vector<int> &complementHEdges, std::vector<int> &singularities) {
    std::vector<int> dualGraphVertices;
    //if both singularities and complementHEdges are empty crash
    for (int i: complementHEdges) {
        auto heh = trimesh_.halfedge_handle(i);
        auto vht = trimesh_.to_vertex_handle(heh);
        auto vhf = trimesh_.from_vertex_handle(heh);
        if (std::find(dualGraphVertices.begin(), dualGraphVertices.end(), vhf.idx()) == dualGraphVertices.end()) {
            dualGraphVertices.push_back(vhf.idx());
        }
        if (std::find(dualGraphVertices.begin(), dualGraphVertices.end(), vht.idx()) == dualGraphVertices.end()) {
            dualGraphVertices.push_back(vht.idx());
        }
    }
    if (complementHEdges.empty()) {
        dualGraphVertices.push_back(singularities[0]);
        trimesh_.status(trimesh_.vertex_handle(singularities[0])).set_tagged2(true);
    }
    return dualGraphVertices;
}

void DijkstraDistance::initProperties(std::vector<int> &cutGraphVertices, const bool first_it) {
    auto vertexDist = OpenMesh::VProp<double>(trimesh_, "vertexDist");
    auto vertexOrigin = OpenMesh::VProp<int>(trimesh_, "vertexOrigin");
    auto vertexPredecessor = OpenMesh::VProp<int>(trimesh_, "vertexPredecessor");
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    int max = INT_MAX, zeroDist = 0.0, standardAppearance = 0;
    for (auto vh: trimesh_.vertices()) {
        vertexDist[vh] = max;
        trimesh_.status(vh).set_tagged(false);
        if (first_it) {
            vertexOrigin[vh] = max;
            vertexPredecessor[vh] = max;
            vertexAppearanceCG[vh] = standardAppearance;
        }
    }
    for (auto i: cutGraphVertices) {
        auto vh = trimesh_.vertex_handle(i);
        vertexDist[vh] = zeroDist;
        if (first_it) {
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
            double distance = vertexDist[vh_og] + 1.0;
            if (!trimesh_.status(vhe_it->to()).tagged()
                && distance < vertexDist[vhe_it->to()]) {
                vertexDist[vhe_it->to()] = distance;
                vertexOrigin[vhe_it->to()] = vertexOrigin[vh_og];
                vertexPredecessor[vhe_it->to()] = vh_og.idx();
            }
        }
    }
}

void DijkstraDistance::addPathToCutGraph(std::vector<int> &cutGraphVertices, std::vector<int> &complementHEdges,
                                         const int i, std::vector<int> &cutGraphWoBoundary) {
    auto vertexPredecessor = OpenMesh::VProp<int>(trimesh_, "vertexPredecessor");
    auto vertexOrigin = OpenMesh::VProp<int>(trimesh_, "vertexOrigin");
    auto cutGraphHe = OpenMesh::HProp<bool>(trimesh_, "cutGraphHe");;
    auto vh = trimesh_.vertex_handle(i);
    bool flag = true;
    while (flag) {
        if (std::find(cutGraphVertices.begin(), cutGraphVertices.end(), vh.idx()) == cutGraphVertices.end()) {
            auto vh_pred = trimesh_.vertex_handle(vertexPredecessor[vh]);
            auto heh = trimesh_.find_halfedge(vh, vh_pred);
            auto oheh = trimesh_.opposite_halfedge_handle(heh);
            cutGraphVertices.push_back(vh.idx());
            complementHEdges.push_back(heh.idx());
            complementHEdges.push_back(oheh.idx());
            cutGraphWoBoundary.push_back(heh.idx());
            cutGraphWoBoundary.push_back(oheh.idx());
            vh = vh_pred;
            cutGraphHe[heh] = true;
            cutGraphHe[oheh] = true;
        } else {
            flag = false;
        }
    }
}

void DijkstraDistance::getDualGraph(const std::vector<int> &faces) {
    trimesh_.request_face_status();
    initDualGraphProp(faces);
    calculateDGDijkstra(faces);
    colorDualGraph();
    trimesh_.release_face_status();
}

void DijkstraDistance::initDualGraphProp(const std::vector<int> &faces) {
    double initValue = INT_MAX, zeroDist = 0.0;
    auto dualGraphDist = OpenMesh::FProp<double>(trimesh_, "dualGraphDist");
    auto dualGraphOrigin = OpenMesh::FProp<int>(trimesh_, "dualGraphOrigin");
    auto dualGraphPred = OpenMesh::FProp<int>(trimesh_, "dualGraphPred");
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");

    for (auto fh: trimesh_.faces()) {
        dualGraphDist[fh] = initValue;
        dualGraphOrigin[fh] = initValue;
        dualGraphPred[fh] = initValue;
        cutGraphFZone[fh] = zeroDist;
        trimesh_.status(fh).set_tagged(false);
    }
    auto fh = trimesh_.face_handle(faces[0]);
    dualGraphDist[fh] = zeroDist;
    dualGraphOrigin[fh] = faces[0];
    dualGraphPred[fh] = faces[0];
}

void DijkstraDistance::calculateDGDijkstra(const std::vector<int> &faces) {
    auto dualGraphDist = OpenMesh::FProp<double>(trimesh_, "dualGraphDist");
    auto dualGraphOrigin = OpenMesh::FProp<int>(trimesh_, "dualGraphOrigin");
    auto dualGraphPred = OpenMesh::FProp<int>(trimesh_, "dualGraphPred");
    auto faceSel = OpenMesh::FProp<bool>(trimesh_, "faceSel");
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
                && distance < dualGraphDist[*ff_it] && faceSel[*ff_it]) {
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

void DijkstraDistance::colorDualGraph() {
    auto dualGraphOrigin = OpenMesh::FProp<int>(trimesh_, "dualGraphOrigin");
    auto borderDualG = OpenMesh::EProp<int>(trimesh_, "borderDualG");
    for (auto fh: trimesh_.faces()) {
        for (TriMesh::FaceFaceIter ff_it = trimesh_.ff_iter(fh); ff_it.is_valid(); ++ff_it) {
            // border of 3d object
            if (fh.is_boundary() && (dualGraphOrigin[fh] != INT_MAX)) {
                for (TriMesh::FaceEdgeIter fhe_it = trimesh_.fe_iter(fh); fhe_it.is_valid(); ++fhe_it) {
                    if (fhe_it->is_boundary()) {
                        borderDualG[*fhe_it] = 1;
                    }
                }
            }
            // border of selection
            if (dualGraphOrigin[fh] != dualGraphOrigin[*ff_it]) {
                searchComEBetweenF(fh, *ff_it, 1);
                // inside dual graph
            } else if ((dualGraphOrigin[fh] == dualGraphOrigin[*ff_it]) && (dualGraphOrigin[fh] != INT_MAX)) {
                searchComEBetweenF(fh, *ff_it, 2);
                // outside dual graph
            } else {
                searchComEBetweenF(fh, *ff_it, 3);
            }
        }
    }
}

void DijkstraDistance::searchComEBetweenF(const OpenMesh::FaceHandle fh, const OpenMesh::SmartFaceHandle fh2,
                                          const int color) {
    auto borderDualG = OpenMesh::EProp<int>(trimesh_, "borderDualG");
    for (TriMesh::FaceEdgeIter fe_it = trimesh_.fe_iter(fh); fe_it.is_valid(); ++fe_it) {
        for (TriMesh::FaceEdgeIter fe2_it = trimesh_.fe_iter(fh2); fe2_it.is_valid(); ++fe2_it) {
            if (fe_it->idx() == fe2_it->idx()) {
                borderDualG[*fe_it] = color;
            }
        }
    }
}

std::vector<int>
DijkstraDistance::calculateDijkstra(const std::vector<int> HeConstraints, const double refDist,
                                    const bool includeBoundary) {
    auto distanceBaryCenter = OpenMesh::FProp<double>(trimesh_, "distanceBaryCenter");
    std::vector<int> includedFaces;
    std::vector<int> constraintFaces = transformHehToFaces(HeConstraints);
    dijkstraDistBaryCenter(includedFaces, refDist);
    if (includeBoundary) {
        includeBoundaryFaces(includedFaces, refDist);
    }
    return includedFaces;
}

void DijkstraDistance::dijkstraDistBaryCenter(std::vector<int> &includedFaces, const double refDist) {
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
        includedFaces.push_back(fh.idx());
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
    trimesh_.release_face_status();
}

std::vector<int> DijkstraDistance::transformHehToFaces(const std::vector<int> &constraintHeh) {
    trimesh_.request_face_status();
    auto origin_constraint = OpenMesh::FProp<int>(trimesh_, "origin_constraint");
    auto distanceBaryCenter = OpenMesh::FProp<double>(trimesh_, "distanceBaryCenter");
    auto predecessor_face = OpenMesh::FProp<int>(trimesh_, "predecessor_face");
    auto positionHessianMatrix = OpenMesh::FProp<int>(trimesh_, "positionHessianMatrix");
    auto periodJump = OpenMesh::HProp<int>(trimesh_, "periodJump");
    auto currentPJ = OpenMesh::FProp<int>(trimesh_, "currentPJ");
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
        currentPJ[fh] = minDistance;
        trimesh_.status(fh).set_tagged(false);
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
    trimesh_.release_face_status();
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
DijkstraDistance::getHeVectorOfSelection(const std::vector<int> &includedFaces) {
    std::vector<int> includedHe;
    for (int i: includedFaces) {
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
    for (auto eh: trimesh_.edges()) {
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

void DijkstraDistance::cleanMeshOfProps() {
//    auto solCoordSysUV = OpenMesh::VProp<std::vector<Point>>(trimesh_, "solCoordSysUV");
    if (OpenMesh::hasProperty<OpenMesh::VertexHandle, std::vector<Point>>(trimesh_, "solCoordSysUV")) {
        auto propH = OpenMesh::VProp<std::vector<Point>>(trimesh_, "solCoordSysUV").getRawProperty();
        trimesh_.remove_property(propH);
    }
    //    auto faceSel = OpenMesh::FProp<bool>(trimesh_, "faceSel");
    if (OpenMesh::hasProperty<OpenMesh::FaceHandle, bool>(trimesh_, "faceSel")) {
        auto propH = OpenMesh::FProp<bool>(trimesh_, "faceSel").getRawProperty();
        trimesh_.remove_property(propH);
    }
//    auto heColor = OpenMesh::HProp<int>(trimesh_, "heColor");
    if (OpenMesh::hasProperty<OpenMesh::HalfedgeHandle, int>(trimesh_, "heColor")) {
        auto propH = OpenMesh::HProp<int>(trimesh_, "heColor").getRawProperty();
        trimesh_.remove_property(propH);
    }
    //    auto vertexDist = OpenMesh::VProp<double>(trimesh_, "vertexDist");
    if (OpenMesh::hasProperty<OpenMesh::VertexHandle, double>(trimesh_, "vertexDist")) {
        auto propH = OpenMesh::VProp<double>(trimesh_, "vertexDist").getRawProperty();
        trimesh_.remove_property(propH);
    }
    //    auto vertexOrigin = OpenMesh::VProp<int>(trimesh_, "vertexOrigin");
    if (OpenMesh::hasProperty<OpenMesh::VertexHandle, int>(trimesh_, "vertexOrigin")) {
        auto propH = OpenMesh::VProp<int>(trimesh_, "vertexOrigin").getRawProperty();
        trimesh_.remove_property(propH);
    }
    //    auto vertexPredecessor = OpenMesh::VProp<int>(trimesh_, "vertexPredecessor");
    if (OpenMesh::hasProperty<OpenMesh::VertexHandle, int>(trimesh_, "vertexPredecessor")) {
        auto propH = OpenMesh::VProp<int>(trimesh_, "vertexPredecessor").getRawProperty();
        trimesh_.remove_property(propH);
    }
    //    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    if (OpenMesh::hasProperty<OpenMesh::VertexHandle, int>(trimesh_, "vertexAppearanceCG")) {
        auto propH = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG").getRawProperty();
        trimesh_.remove_property(propH);
    }
    //    auto cutGraphHe = OpenMesh::HProp<bool>(trimesh_, "cutGraphHe");
    if (OpenMesh::hasProperty<OpenMesh::HalfedgeHandle, bool>(trimesh_, "cutGraphHe")) {
        auto propH = OpenMesh::HProp<bool>(trimesh_, "cutGraphHe").getRawProperty();
        trimesh_.remove_property(propH);
    }
    //    auto dualGraphDist = OpenMesh::FProp<double>(trimesh_, "dualGraphDist");
    if (OpenMesh::hasProperty<OpenMesh::FaceHandle, double>(trimesh_, "dualGraphDist")) {
        auto propH = OpenMesh::FProp<double>(trimesh_, "dualGraphDist").getRawProperty();
        trimesh_.remove_property(propH);
    }
    //    auto dualGraphOrigin = OpenMesh::FProp<int>(trimesh_, "dualGraphOrigin");
    if (OpenMesh::hasProperty<OpenMesh::FaceHandle, int>(trimesh_, "dualGraphOrigin")) {
        auto propH = OpenMesh::FProp<int>(trimesh_, "dualGraphOrigin").getRawProperty();
        trimesh_.remove_property(propH);
    }
    //    auto dualGraphPred = OpenMesh::FProp<int>(trimesh_, "dualGraphPred");
    if (OpenMesh::hasProperty<OpenMesh::FaceHandle, int>(trimesh_, "dualGraphPred")) {
        auto propH = OpenMesh::FProp<int>(trimesh_, "dualGraphPred").getRawProperty();
        trimesh_.remove_property(propH);
    }
    //    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");
    if (OpenMesh::hasProperty<OpenMesh::FaceHandle, int>(trimesh_, "cutGraphFZone")) {
        auto propH = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone").getRawProperty();
        trimesh_.remove_property(propH);
    }
    //    auto borderDualG = OpenMesh::EProp<int>(trimesh_, "borderDualG");
    if (OpenMesh::hasProperty<OpenMesh::EdgeHandle, int>(trimesh_, "borderDualG")) {
        auto propH = OpenMesh::EProp<int>(trimesh_, "borderDualG").getRawProperty();
        trimesh_.remove_property(propH);
    }
    //    auto distanceBaryCenter = OpenMesh::FProp<double>(trimesh_, "distanceBaryCenter");
    if (OpenMesh::hasProperty<OpenMesh::FaceHandle, double>(trimesh_, "distanceBaryCenter")) {
        auto propH = OpenMesh::FProp<double>(trimesh_, "distanceBaryCenter").getRawProperty();
        trimesh_.remove_property(propH);
    }
//    auto origin_constraint = OpenMesh::FProp<int>(trimesh_, "origin_constraint");
    if (OpenMesh::hasProperty<OpenMesh::FaceHandle, int>(trimesh_, "origin_constraint")) {
        auto propH = OpenMesh::FProp<int>(trimesh_, "origin_constraint").getRawProperty();
        trimesh_.remove_property(propH);
    }
    //    auto predecessor_face = OpenMesh::FProp<int>(trimesh_, "predecessor_face");
    if (OpenMesh::hasProperty<OpenMesh::FaceHandle, int>(trimesh_, "predecessor_face")) {
        auto propH = OpenMesh::FProp<int>(trimesh_, "predecessor_face").getRawProperty();
        trimesh_.remove_property(propH);
    }
    //    auto periodJump = OpenMesh::HProp<int>(trimesh_, "periodJump");
    if (OpenMesh::hasProperty<OpenMesh::HalfedgeHandle, int>(trimesh_, "periodJump")) {
        auto propH = OpenMesh::HProp<int>(trimesh_, "periodJump").getRawProperty();
        trimesh_.remove_property(propH);
    }
//    auto positionHessianMatrix = OpenMesh::FProp<int>(trimesh_, "positionHessianMatrix");
    if (OpenMesh::hasProperty<OpenMesh::FaceHandle, int>(trimesh_, "positionHessianMatrix")) {
        auto propH = OpenMesh::FProp<int>(trimesh_, "positionHessianMatrix").getRawProperty();
        trimesh_.remove_property(propH);
    }
    //    auto currentPJ = OpenMesh::FProp<int>(trimesh_, "currentPJ");
    if (OpenMesh::hasProperty<OpenMesh::FaceHandle, int>(trimesh_, "currentPJ")) {
        auto propH = OpenMesh::FProp<int>(trimesh_, "currentPJ").getRawProperty();
        trimesh_.remove_property(propH);
    }
    //    auto referenceHeIdx = OpenMesh::FProp<int>(trimesh_, "referenceHeIdx");
    if (OpenMesh::hasProperty<OpenMesh::FaceHandle, int>(trimesh_, "referenceHeIdx")) {
        auto propH = OpenMesh::FProp<int>(trimesh_, "referenceHeIdx").getRawProperty();
        trimesh_.remove_property(propH);
    }
//    auto uVectorField = OpenMesh::FProp<Point>(trimesh_, "uVectorField");
    if (OpenMesh::hasProperty<OpenMesh::FaceHandle, Point>(trimesh_, "uVectorField")) {
        auto propH = OpenMesh::FProp<Point>(trimesh_, "uVectorField").getRawProperty();
        trimesh_.remove_property(propH);
    }
    //    auto uVectorFieldRotOne = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotOne");
    if (OpenMesh::hasProperty<OpenMesh::FaceHandle, Point>(trimesh_, "uVectorFieldRotOne")) {
        auto propH = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotOne").getRawProperty();
        trimesh_.remove_property(propH);
    }
    //    auto uVectorFieldRotTwo = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotTwo");
    if (OpenMesh::hasProperty<OpenMesh::FaceHandle, Point>(trimesh_, "uVectorFieldRotTwo")) {
        auto propH = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotTwo").getRawProperty();
        trimesh_.remove_property(propH);
    }
    //    auto constraint_angle = OpenMesh::FProp<double>(trimesh_, "constraint_angle");
    if (OpenMesh::hasProperty<OpenMesh::FaceHandle, double>(trimesh_, "constraint_angle")) {
        auto propH = OpenMesh::FProp<double>(trimesh_, "constraint_angle").getRawProperty();
        trimesh_.remove_property(propH);
    }
    //    auto crossFieldIdx = OpenMesh::VProp<double>(trimesh_, "crossFieldIdx");
    if (OpenMesh::hasProperty<OpenMesh::VertexHandle, double>(trimesh_, "crossFieldIdx")) {
        auto propH = OpenMesh::VProp<double>(trimesh_, "crossFieldIdx").getRawProperty();
        trimesh_.remove_property(propH);
    }
//    auto vertexPosUi = OpenMesh::VProp<int>(trimesh_, "vertexPosUi");
    if (OpenMesh::hasProperty<OpenMesh::VertexHandle, double>(trimesh_, "vertexPosUi")) {
        auto propH = OpenMesh::VProp<double>(trimesh_, "vertexPosUi").getRawProperty();
        trimesh_.remove_property(propH);
    }
//    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    if (OpenMesh::hasProperty<OpenMesh::VertexHandle, int>(trimesh_, "vertexPosVi")) {
        auto propH = OpenMesh::VProp<int>(trimesh_, "vertexPosVi").getRawProperty();
        trimesh_.remove_property(propH);
    }
//    auto uVectorFieldRotThree = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotThree");
    if (OpenMesh::hasProperty<OpenMesh::FaceHandle, Point>(trimesh_, "uVectorFieldRotThree")) {
        auto propH = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotThree").getRawProperty();
        trimesh_.remove_property(propH);
    }
//    auto uVectorFieldRotFour= OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotFour");
    if (OpenMesh::hasProperty<OpenMesh::FaceHandle, Point>(trimesh_, "uVectorFieldRotFour")) {
        auto propH = OpenMesh::FProp<Point>(trimesh_, "uVectorFieldRotFour").getRawProperty();
        trimesh_.remove_property(propH);
    }
    trimesh_.garbage_collection();
}
