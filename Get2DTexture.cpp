//
// Created by wuschelbueb on 27.06.22.
//

#include "Get2DTexture.h"

void Get2DTexture::initProperty(OpenMesh::HPropHandleT<ACG::Vec2d> hp_texture) {
//    auto quad = OpenMesh::HProp<OpenMesh::Vec3f>(trimesh_, "quad");
    for (auto he: trimesh_.halfedges()) {
//        quad[he] = {0, 0, 0};
        trimesh_.property(hp_texture, he) = {0,0};
    }
}

void Get2DTexture::get2DTexture(TriMesh::FaceHalfedgeIter &fh_it, double &u, double &v) {
//    auto quad = OpenMesh::HProp<OpenMesh::Vec3f>(trimesh_, "quad");
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    auto solCoordSysUV = OpenMesh::VProp<std::vector<OpenMesh::Vec2d>>(trimesh_, "solCoordSysUV");
    auto vh = fh_it->to();
    OpenMesh::Vec2d uvCoord;
    if (vertexAppearanceCG[vh] > 1 && cutGraphFZone[fh_it->face()] != 0) {
        int pos = getPositionConstraintRow(vh, cutGraphFZone[fh_it->face()]);
        uvCoord = solCoordSysUV[vh][pos];
    } else {
        uvCoord = solCoordSysUV[vh][0];
    }
    u = uvCoord[0];
    v = uvCoord[1];
//    quad[*fh_it] = {u, v, 0};
//    auto test = OpenMesh::getProperty<OpenMesh::HalfedgeHandle, OpenMesh::Vec3f>(trimesh_, "quad");
//    std::cout << "halfedge " << fh_it->idx() << " with u: " << u << "\tv: " << v << std::endl
//              << "and point " << quad[*fh_it] << std::endl << "getName " << test.getName() << std::endl;
}

int Get2DTexture::getPositionConstraintRow(OpenMesh::SmartVertexHandle &vh, int cutGraphZone) {
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    std::vector<int> sectors;
    for (auto it = vh.outgoing_halfedges().begin(); it.is_valid(); ++it) {
        int zone = cutGraphFZone[it->face()];
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
