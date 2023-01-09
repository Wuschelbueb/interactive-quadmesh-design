//
// Created by wuschelbueb on 04.01.23.
//

#include "CreateObjData.hh"

void CreateObjData::getStream() {
    trimesh_.release_halfedge_status();
    trimesh_.request_halfedge_status();
    trimesh_.request_vertex_status();
    trimesh_.request_face_status();
    getLambda();
    createData();
    trimesh_.release_vertex_status();
    trimesh_.release_face_status();
}

void CreateObjData::getLambda() {
    auto quadTextr = OpenMesh::HProp<OpenMesh::Vec2d>(trimesh_, "quadTextr");
    double maxX = 0, minX = 0, maxY = 0, minY = 0;
    for (auto he: trimesh_.halfedges()) {
        if (quadTextr[he][0] > maxX) {
            maxX = quadTextr[he][0];
        } else if (quadTextr[he][0] <= minX) {
            minX = quadTextr[he][0];
        }
        if (quadTextr[he][1] > maxY) {
            maxY = quadTextr[he][1];
        } else if (quadTextr[he][1] <= minY) {
            minY = quadTextr[he][1];
        }
    }
    // here min - max because openflipper saves coordinates as negatives
    double rangeX = maxX - minX;
    double rangeY = maxY - minY;
    if (rangeX > rangeY) {
        lambda = rangeX;
        oldMin = minX;
        oldMax = maxX;
        dataStream << "# went for X: " << lambda << " (" << rangeY << ") with minX: " << minX << " and maxX: " << maxX
                   << std::endl;
    } else {
        lambda = rangeY;
        oldMin = minY;
        oldMax = maxY;
        dataStream << "# went for Y: " << lambda << " (" << rangeX << ") with minY: " << minY << " and maxY: " << maxY
                   << std::endl;
    }
}

void CreateObjData::createData() {
    auto quadTextr = OpenMesh::HProp<OpenMesh::Vec2d>(trimesh_, "quadTextr");

    bool oneTime = true;
    std::map<int, int> connectVhToTexCoord;
    for (auto he: trimesh_.halfedges()) {
        if (he.to().idx() == centerVertex.idx() && oneTime) {
            dataStream << "# Selected Point: " << trimesh_.point(centerVertex)
                       << "\n# with Index: " << centerVertex.idx()
                       << "\n# and TexCoords: " << transformToTexCoord(quadTextr[he][0])
                       << " " << transformToTexCoord(quadTextr[he][1]) << std::endl;
            oneTime = false;
        }
    }
    resetVertexStatus();
    writeVertices(connectVhToTexCoord);
    resetVertexStatus();
    writeTexCoords(connectVhToTexCoord);
    resetVertexStatus();
    writeFaces(connectVhToTexCoord);
}

void CreateObjData::resetVertexStatus() {
    for (auto vh: trimesh_.vertices()) {
        trimesh_.status(vh).set_tagged(false);
    }
}

void CreateObjData::writeVertices(std::map<int, int> &connectVhToTexCoord) {
    int counter = 1;
    dataStream << "# object.obj\no object\n# Vertices\n";
    for (auto he: trimesh_.halfedges()) {
        if (he.to().tagged()) {
            continue;
        }
        dataStream << "# " << he.to().idx() << "\nv " << trimesh_.point(he.to()) << std::endl
                   << "vn " << trimesh_.normal(he.to())[0] << " " << trimesh_.normal(he.to())[1] << " "
                   << trimesh_.normal(he.to())[2] << std::endl;
        connectVhToTexCoord.insert({he.to().idx(), counter++});
        trimesh_.status(he.to()).set_tagged(true);
    }
}

void CreateObjData::writeTexCoords(std::map<int, int> &connectVhToTexCoord) {
    auto quadTextr = OpenMesh::HProp<OpenMesh::Vec2d>(trimesh_, "quadTextr");
    auto heColor = OpenMesh::HProp<int>(trimesh_, "heColor");
    int counter = 1;

    dataStream << "\n# Texture Coordinates\n";
    for (auto he: trimesh_.halfedges()) {
        if (!(!he.is_boundary() && !trimesh_.status(he.to()).tagged())) {
            continue;
        }
        auto heNb = connectVhToTexCoord[he.to().idx()];
        if (heColor[he] != 1) {
//            dataStream << "# " << heNb << "\nvt " << transformToTexCoord(quadTextr[he][0])
//                       << " " << transformToTexCoord(quadTextr[he][1]) << std::endl;
            dataStream << "# " << heNb << "\nvt " << quadTextr[he][0]
                       << " " << quadTextr[he][1] << std::endl;
        } else {
            dataStream << "vt 0 0" << std::endl;
        }
        trimesh_.status(he.to()).set_tagged(true);
    }
}

void CreateObjData::writeFaces(std::map<int, int> &connectVhToTexCoord) {
    auto quadTextr = OpenMesh::HProp<OpenMesh::Vec2d>(trimesh_, "quadTextr");
    auto heColor = OpenMesh::HProp<int>(trimesh_, "heColor");

    dataStream << "\n# Faces\n"; //Faces v/vt/vn
    for (auto he: trimesh_.halfedges()) {
        if (!(!trimesh_.is_boundary(he) && !trimesh_.status(he).tagged())) {
            continue;
        }
        auto heNb = connectVhToTexCoord[he.to().idx()];
        auto heNextNb = connectVhToTexCoord[he.next().to().idx()];
        auto hePrevNb = connectVhToTexCoord[he.prev().to().idx()];
        if (heColor[he] != 1) {
            dataStream << "\n# " << he.face().idx() << "\nf "
                       << heNb << "/" << heNb << "/" << heNb << " "
                       << heNextNb << "/" << heNextNb << "/" << heNextNb << " "
                       << hePrevNb << "/" << hePrevNb << "/" << hePrevNb;
        } else {
            dataStream << "\n# " << he.face().idx() << "\nf "
                       << heNb << "//" << heNb << " "
                       << heNextNb << "//" << heNextNb << " "
                       << hePrevNb << "//" << hePrevNb;
        }
        trimesh_.status(he).set_tagged(true);
        trimesh_.status(he.next()).set_tagged(true);
        trimesh_.status(he.prev()).set_tagged(true);
    }
    dataStream << "\n# End of File!";
}


double CreateObjData::transformToTexCoord(double value) {
    // NewValue = (((OldValue - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin
    value = (value - oldMin) / std::abs(lambda);
    return std::ceil(value * 100.0) / 100.0;
}