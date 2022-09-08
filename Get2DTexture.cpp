//
// Created by wuschelbueb on 27.06.22.
//

#include "Get2DTexture.h"

void Get2DTexture::initProperty() {
    auto quadTextr = OpenMesh::HProp<OpenMesh::Vec2d>(trimesh_, "quadTextr");
    for (auto he: trimesh_.halfedges()) {
        quadTextr[he] = {0, 0};
    }
}

void Get2DTexture::get2DTexture(OpenMesh::SmartHalfedgeHandle he) {
    auto quadTextr = OpenMesh::HProp<OpenMesh::Vec2d>(trimesh_, "quadTextr");
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    auto solCoordSysUV = OpenMesh::VProp<std::vector<OpenMesh::Vec2d>>(trimesh_, "solCoordSysUV");
    auto vh = he.to();
    int pos = 0;
    if (vertexAppearanceCG[vh] > 1) {
        pos = getPositionConstraintRow(vh, cutGraphFZone[he.face()]);
    }
    quadTextr[he] = solCoordSysUV[vh][pos];
}

int Get2DTexture::getPositionConstraintRow(OpenMesh::SmartVertexHandle &vh, int cutGraphZone) {
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    std::vector<int> sectors;
    for (TriMesh::VertexFaceIter vf_it = trimesh_.vf_begin(vh); vf_it.is_valid(); ++vf_it) {
        int zone = cutGraphFZone[*vf_it];
        if (zone != 0 && (std::find(sectors.begin(), sectors.end(), zone) == sectors.end())) {
            sectors.push_back(zone);
        }
    }
    try {
        if ((int) sectors.size() != vertexAppearanceCG[vh]) {
            throw 404;
        }
    } catch (int x) {
        std::cerr << "getPositionConstraintRow: vector doesn't have the same size as vertexAppearance property v: "
                  << sectors.size()
                  << " != p: " << vertexAppearanceCG[vh] << "\n";
    }
    std::sort(sectors.begin(), sectors.end());
    auto it1 = find(sectors.begin(), sectors.end(), cutGraphZone);
    int pos = it1 - sectors.begin();
    return pos;
}
