//
// Created by wuschelbueb on 04.01.23.
//

#include "CreateObjData.hh"

void CreateObjData::getStream() {
    trimesh_.release_halfedge_status();
    trimesh_.request_halfedge_status();
    trimesh_.request_vertex_status();
    createData();
    trimesh_.release_vertex_status();
}

void CreateObjData::createData() {
    auto quadTextr = OpenMesh::HProp<OpenMesh::Vec2d>(trimesh_, "quadTextr");

    bool oneTime = true;
    std::map<int, int> mapVertexIdx;
    std::map<int, int> mapTexCoordToIdx;
    for (auto he: trimesh_.halfedges()) {
        if (he.to().idx() == centerVertex.idx() && oneTime) {
            dataStream << "# TexCoords: " << quadTextr[he][0] << " " << quadTextr[he][1] << std::endl;
            oneTime = false;
        }
    }
    resetVertexStatus();
    writeVertices(mapVertexIdx);
    resetVertexStatus();
    writeTexCoords(mapTexCoordToIdx);
    resetVertexStatus();
    writeFaces(mapVertexIdx, mapTexCoordToIdx);
    resetVertexStatus();
}

void CreateObjData::resetVertexStatus() {
    for (auto vh: trimesh_.vertices()) {
        trimesh_.status(vh).set_tagged(false);
    }
}

void CreateObjData::writeVertices(std::map<int, int> &mapVertexIdx) {
    int counter = 1;
    dataStream << "# object.obj\no object\n# Vertices\n";
    for (auto he: trimesh_.halfedges()) {
        if (he.is_boundary() || he.to().tagged()) {
            continue;
        }
        dataStream << "# " << he.to().idx() << "\nv " << trimesh_.point(he.to()) << std::endl
                   << "vn " << trimesh_.normal(he.to())[0] << " " << trimesh_.normal(he.to())[1] << " "
                   << trimesh_.normal(he.to())[2] << std::endl;
        mapVertexIdx.insert({he.to().idx(), counter++});

        trimesh_.status(he.to()).set_tagged(true);
    }
}

void CreateObjData::writeTexCoords(std::map<int, int> &mapTexCoordToIdx) {
    auto quadTextr = OpenMesh::HProp<OpenMesh::Vec2d>(trimesh_, "quadTextr");
    auto heColor = OpenMesh::HProp<int>(trimesh_, "heColor");
    int counter = 1;
    dataStream << "\n# Texture Coordinates\n";
    for (auto he: trimesh_.halfedges()) {
        if (he.is_boundary() || he.to().tagged() || heColor[he] == 1) {
            continue;
        }
        dataStream << "# " << he.to().idx() << "\nvt " << quadTextr[he][0]
                   << " " << quadTextr[he][1] << std::endl;
        mapTexCoordToIdx.insert({he.to().idx(), counter++});
        trimesh_.status(he.to()).set_tagged(true);
    }
}

void CreateObjData::writeFaces(const std::map<int, int> &mapVertexToIdx, const std::map<int, int> &mapTexCoordToIdx) {
    auto quadTextr = OpenMesh::HProp<OpenMesh::Vec2d>(trimesh_, "quadTextr");
    auto heColor = OpenMesh::HProp<int>(trimesh_, "heColor");

    dataStream << "\n# Faces (v/vt/n)\n"; //Faces v/vt/vn
    for (auto he: trimesh_.halfedges()) {
        if (he.is_boundary() || he.tagged()) {
            continue;
        }
        auto heNb = mapVertexToIdx.find(he.to().idx());
        auto heNextNb = mapVertexToIdx.find(he.next().to().idx());
        auto hePrevNb = mapVertexToIdx.find(he.prev().to().idx());
        auto heTex = mapTexCoordToIdx.find(he.to().idx());
        auto heNextTex = mapTexCoordToIdx.find(he.next().to().idx());
        auto hePrevTex = mapTexCoordToIdx.find(he.prev().to().idx());
        if (heColor[he] != 1) {
            dataStream << "\n# " << he.face().idx() << "\nf "
                       << heNb->second << "/" << heTex->second << "/" << heNb->second << " "
                       << heNextNb->second << "/" << heNextTex->second << "/" << heNextNb->second << " "
                       << hePrevNb->second << "/" << hePrevTex->second << "/" << hePrevNb->second;
        } else {
            dataStream << "\n# " << he.face().idx() << "\nf "
                       << heNb->second << "//" << heNb->second << " "
                       << heNextNb->second << "//" << heNextNb->second << " "
                       << hePrevNb->second << "//" << hePrevNb->second;
        }
        trimesh_.status(he).set_tagged(true);
        trimesh_.status(he.next()).set_tagged(true);
        trimesh_.status(he.prev()).set_tagged(true);
    }
    dataStream << "\n# placeholder";
    dataStream << "\n# End of File";
}