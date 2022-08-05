//
// Created by wuschelbueb on 27.06.22.
//

#include "Get2DTexture.h"

void Get2DTexture::initProperty(OpenMesh::HPropHandleT<ACG::Vec2d> hp_texture) {
//    auto quad = OpenMesh::HProp<OpenMesh::Vec3f>(trimesh_, "quad");
    for (auto he: trimesh_.halfedges()) {
//        quad[he] = {0, 0, 0};
        trimesh_.property(hp_texture, he) = {0, 0};
    }
}

void Get2DTexture::get2DTexture(TriMesh::FaceHalfedgeIter &fh_it, double &u, double &v) {
//    auto quad = OpenMesh::HProp<OpenMesh::Vec3f>(trimesh_, "quad");
    auto cutGraphFZone = OpenMesh::FProp<int>(trimesh_, "cutGraphFZone");
    auto vertexAppearanceCG = OpenMesh::VProp<int>(trimesh_, "vertexAppearanceCG");
    auto solCoordSysUV = OpenMesh::VProp<std::vector<OpenMesh::Vec2d>>(trimesh_, "solCoordSysUV");
    auto vh = fh_it->to();
//    std::cout << "vertex " << vh.idx() << " with app: " << vertexAppearanceCG[vh] << std::endl;
    OpenMesh::Vec2d uvCoord;
    if (vertexAppearanceCG[vh] > 1 && cutGraphFZone[fh_it->face()] != 0) {
        int pos = getPositionConstraintRow(vh, cutGraphFZone[fh_it->face()]);
//        std::cout << "position " << pos << std::endl;
        uvCoord = solCoordSysUV[vh][pos];
//        std::cout << "uvCoord (if) " << solCoordSysUV[vh][pos] << std::endl;
    } else {
        uvCoord = solCoordSysUV[vh][0];
//        std::cout << "uvCoord (else) " << solCoordSysUV[vh][0] << std::endl;
    }
    u = uvCoord[0];
    v = uvCoord[1];
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
