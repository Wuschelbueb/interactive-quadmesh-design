//
// Created by wuschelbueb on 04.01.23.
//

#include "CreateObjData.hh"

void CreateObjData::getStream() {
    trimesh_.release_halfedge_status();
    trimesh_.request_halfedge_status();
    getLambda();
    createData();
}

double CreateObjData::transformToTexCoord(double value) {
    return (value / (2 * lambda)) + 0.5;
}

void CreateObjData::getLambda() {
    auto quadTextr = OpenMesh::HProp<OpenMesh::Vec2d>(trimesh_, "quadTextr");
    double max = 0, min = 0;
    for (auto he: trimesh_.halfedges()) {
        if (quadTextr[he][0] > max) {
            max = quadTextr[he][0];
        } else if (quadTextr[he][0] <= min) {
            min = quadTextr[he][0];
        }
        if (quadTextr[he][1] > max) {
            max = quadTextr[he][1];
        } else if (quadTextr[he][1] <= min) {
            min = quadTextr[he][1];
        }
    }
    lambda = (std::abs(max) > std::abs(min)) ? max : min;
}

void CreateObjData::createData() {
    bool oneTime = true;
    auto quadTextr = OpenMesh::HProp<OpenMesh::Vec2d>(trimesh_, "quadTextr");
    auto heColor = OpenMesh::HProp<int>(trimesh_, "heColor");

    for (auto he: trimesh_.halfedges()) {
        if (he.to().idx() == centerVertex.idx() && oneTime) {
            objData << "# Point in 3D space: " << trimesh_.point(centerVertex)
                   << "\n# vt below is the centerVertex " << quadTextr[he]
                   << "\n# vt " << transformToTexCoord(quadTextr[he][0])
                   << " " << transformToTexCoord(quadTextr[he][1])
                   << "\n# lambda center: " << lambda / 2 << std::endl;
            oneTime = false;
        }
    }
    objData << "# object.obj\no object\n# Vertices\n";
    for (auto vh: trimesh_.vertices()) {
        objData<< "v " << trimesh_.point(vh) << std::endl;
    }
    objData<< "\n# Vertex Normals\n";
    for (auto fh: trimesh_.faces()) {
        objData<< "vn " << trimesh_.normal(fh) << std::endl;
    }
    objData<< "\n# Texture Coordinates\n";
    for (auto he: trimesh_.halfedges()) {
        if (!he.is_boundary() && heColor[he] != 1) {
            objData<< "vt " << transformToTexCoord(quadTextr[he][0])
                   << " " << transformToTexCoord(quadTextr[he][1]) << std::endl;
//            objData<< "vt " << quadTextr[he][0]
//                   << " " << quadTextr[he][1] << std::endl;
        } else if (!he.is_boundary() && heColor[he] == 1) {
            objData<< "vt 0 0" << std::endl;
        }
    }
    objData<< "\n# Faces\n"; //Faces v/vt/vn
    for (auto he: trimesh_.halfedges()) {
        if (!trimesh_.is_boundary(he) && !trimesh_.status(he).tagged() && heColor[he] != 1) {
            objData<< "f ";
            objData<< he.to().idx() + 1 << "/" << he.idx() + 1 << "/" << he.face().idx() + 1 << " ";
            objData<< he.next().to().idx() + 1 << "/" << he.next().idx() + 1 << "/" << he.face().idx() + 1 << " ";
            objData<< he.prev().to().idx() + 1 << "/" << he.prev().idx() + 1 << "/" << he.face().idx() + 1 << "\n";
        } else if (!trimesh_.is_boundary(he) && !trimesh_.status(he).tagged() && heColor[he] == 1) {
            objData<< "f ";
            objData<< he.to().idx() + 1 << "//" << he.face().idx() + 1 << " ";
            objData<< he.next().to().idx() + 1 << "//" << he.face().idx() + 1 << " ";
            objData<< he.prev().to().idx() + 1 << "//" << he.face().idx() + 1 << "\n";
        }
        trimesh_.status(he).set_tagged(true);
        trimesh_.status(he.next()).set_tagged(true);
        trimesh_.status(he.prev()).set_tagged(true);
    }
}
