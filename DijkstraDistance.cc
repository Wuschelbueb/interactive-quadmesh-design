#include "DijkstraDistance.hh"

void
DijkstraDistance::calcDijkstraWSingularities(std::vector<int> &complementHEdges, std::vector<int> &singularities,
                                             std::vector<int> &cutGraphWoBoundary) {
    trimesh_.request_vertex_status();
    std::vector<int> cutGraphVertices = createVerticesVector(complementHEdges, singularities);
    initProperties(cutGraphVertices, true);
    for (const int &i: singularities) {
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
    for (const int &i: complementHEdges) {
        auto heh = trimesh_.halfedge_handle(i);
        auto vht = trimesh_.to_vertex_handle(heh);
        auto vhf = trimesh_.from_vertex_handle(heh);
        auto findVhf = std::find(dualGraphVertices.begin(), dualGraphVertices.end(), vhf.idx());
        auto findVht = std::find(dualGraphVertices.begin(), dualGraphVertices.end(), vht.idx());
        if (findVhf == dualGraphVertices.end()) {
            dualGraphVertices.push_back(vhf.idx());
        }
        if (findVht == dualGraphVertices.end()) {
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
    int max = INT_MAX, zeroDist = 0.0, standardAppearance = 1;
    for (const auto &vh: trimesh_.vertices()) {
        vertexDist[vh] = max;
        trimesh_.status(vh).set_tagged(false);
        if (first_it) {
            vertexOrigin[vh] = max;
            vertexPredecessor[vh] = max;
            vertexAppearanceCG[vh] = standardAppearance;
        }
    }
    for (const auto &i: cutGraphVertices) {
        auto vh = trimesh_.vertex_handle(i);
        vertexDist[vh] = zeroDist;
        if (first_it) {
            vertexOrigin[vh] = vh.idx();
            vertexPredecessor[vh] = vh.idx();
        }
    }
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

int DijkstraDistance::vertexGetSmallestDist() {
    auto vertexDist = OpenMesh::VProp<double>(trimesh_, "vertexDist");
    double minDistance = DBL_MAX;
    int idx = INT_MAX;
    for (const auto &vh: trimesh_.vertices()) {
        if (!trimesh_.status(vh).tagged() && vertexDist[vh] < minDistance) {
            minDistance = vertexDist[vh];
            idx = vh.idx();
        }
    }
    return idx;
}

void DijkstraDistance::addPathToCutGraph(std::vector<int> &cutGraphVertices, std::vector<int> &complementHEdges,
                                         const int i, std::vector<int> &cutGraphWoBoundary) {
    auto vertexPredecessor = OpenMesh::VProp<int>(trimesh_, "vertexPredecessor");
    auto vertexOrigin = OpenMesh::VProp<int>(trimesh_, "vertexOrigin");
    auto cutGraphHe = OpenMesh::HProp<bool>(trimesh_, "cutGraphHe");;
    auto vh = trimesh_.vertex_handle(i);
    auto iter1 = std::find(cutGraphVertices.begin(), cutGraphVertices.end(), vh.idx());
    while (iter1 == cutGraphVertices.end()) {
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
        iter1 = std::find(cutGraphVertices.begin(), cutGraphVertices.end(), vh.idx());
    }
}

void DijkstraDistance::getDualGraph(const std::vector<int> &faces) {
    trimesh_.request_face_status();
    initDualGraphProp(faces);
    calculateDGDijkstra(faces);
    colorDualGraph(faces);
    trimesh_.release_face_status();
}

void DijkstraDistance::initDualGraphProp(const std::vector<int> &faces) {
    double initValue = INT_MAX, zeroDist = 0.0;
    auto dualGraphDist = OpenMesh::FProp<double>(trimesh_, "dualGraphDist");
    auto dualGraphOrigin = OpenMesh::FProp<int>(trimesh_, "dualGraphOrigin");
    auto dualGraphPred = OpenMesh::FProp<int>(trimesh_, "dualGraphPred");
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");

    for (const auto &fh: trimesh_.faces()) {
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
    for (const int &i: faces) {
        auto fh = trimesh_.face_handle(i);
        if (!trimesh_.status(fh).tagged() && dualGraphDist[fh] < minDistance) {
            minDistance = dualGraphDist[fh];
            idx = i;
        }
    }
    return idx;
}

void DijkstraDistance::colorDualGraph(const std::vector<int> &faces) {
    auto dualGraphOrigin = OpenMesh::FProp<int>(trimesh_, "dualGraphOrigin");
    auto faceSel = OpenMesh::FProp<bool>(trimesh_, "faceSel");
    auto borderDualG = OpenMesh::EProp<int>(trimesh_, "borderDualG");
    for (auto &i: faces) {
        auto fh = make_smart(trimesh_.face_handle(i), trimesh_);
        for (auto fhe_it: fh.halfedges()) {
            if (fhe_it.opp().is_boundary()) {
                borderDualG[fhe_it.edge()] = 1;
                continue;
            }
            if (!faceSel[fhe_it.opp().face()]) {
                borderDualG[fhe_it.edge()] = 1;
                continue;
            }
            borderDualG[fhe_it.edge()] = 2;
        }
    }
}

std::vector<int>
DijkstraDistance::calculateDijkstra(const std::vector<int> HeConstraints, const double refDist) {
    auto distanceBaryCenter = OpenMesh::FProp<double>(trimesh_, "distanceBaryCenter");
    std::vector<int> includedFaces;
    std::vector<int> constraintFaces = transformHehToFaces(HeConstraints);
    dijkstraDistBaryCenter(includedFaces, refDist);
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

std::vector<int>
DijkstraDistance::getAllHeFromFaces(const std::vector<int> &includedFaces) {
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

std::vector<int>
DijkstraDistance::getHeFromVertex(OpenMesh::VertexHandle selectedVertex, const std::vector<int> &originVertices) {
    trimesh_.request_halfedge_status();
    for (auto he: trimesh_.halfedges()) {
        trimesh_.status(he).set_tagged(false);
    }
    std::vector<int> constraints;
    for (auto i: originVertices) {
        auto vh = trimesh_.vertex_handle(i);
        for (auto voh_it = trimesh_.voh_iter(vh); voh_it.is_valid(); ++voh_it) {
            if (!voh_it->is_boundary() && !voh_it->tagged()) {
                trimesh_.status(*voh_it).set_tagged(true);
                constraints.push_back(voh_it->idx());
            }
        }
    }
    trimesh_.release_halfedge_status();
    return constraints;
}

void DijkstraDistance::cleanMeshOfProps() {
//    auto solCoordSysUV = OpenMesh::VProp<std::vector<OpenMesh::Vec2d>>(trimesh_, "solCoordSysUV");
    if (OpenMesh::hasProperty<OpenMesh::VertexHandle, std::vector<OpenMesh::Vec2d>>(trimesh_, "solCoordSysUV")) {
        auto propH = OpenMesh::VProp<std::vector<OpenMesh::Vec2d>>(trimesh_, "solCoordSysUV").getRawProperty();
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
    if (OpenMesh::hasProperty<OpenMesh::VertexHandle, int>(trimesh_, "vertexPosUi")) {
        auto propH = OpenMesh::VProp<int>(trimesh_, "vertexPosUi").getRawProperty();
        trimesh_.remove_property(propH);
    }
//    auto vertexPosVi = OpenMesh::VProp<int>(trimesh_, "vertexPosVi");
    if (OpenMesh::hasProperty<OpenMesh::VertexHandle, int>(trimesh_, "vertexPosVi")) {
        auto propH = OpenMesh::VProp<int>(trimesh_, "vertexPosVi").getRawProperty();
        trimesh_.remove_property(propH);
    }
    if (OpenMesh::hasProperty<OpenMesh::HalfedgeHandle, OpenMesh::Vec2d>(trimesh_, "quadTextr")) {
        auto propH = OpenMesh::HProp<OpenMesh::Vec2d>(trimesh_, "quadTextr").getRawProperty();
        trimesh_.remove_property(propH);
    }
    if (OpenMesh::hasProperty<OpenMesh::HalfedgeHandle, bool>(trimesh_, "boundaryHe")) {
        auto propH = OpenMesh::HProp<bool>(trimesh_, "boundaryHe").getRawProperty();
        trimesh_.remove_property(propH);
    }
    if (OpenMesh::hasProperty<OpenMesh::VertexHandle, int>(trimesh_, "vertexColor")) {
        auto propH = OpenMesh::VProp<int>(trimesh_, "vertexColor").getRawProperty();
        trimesh_.remove_property(propH);
    }
    trimesh_.request_halfedge_texcoords2D();
    for (auto he: trimesh_.halfedges()) {
        trimesh_.set_texcoord2D(he, {0.f, 0.f});
    }

    trimesh_.garbage_collection();
}

