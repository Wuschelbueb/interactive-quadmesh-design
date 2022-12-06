//
// Created by wuschelbueb on 06.10.22.
//

#ifndef OPENFLIPPER_PATCHPREVIEW_HH
#define OPENFLIPPER_PATCHPREVIEW_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <OpenMesh/Core/Geometry/EigenVectorT.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include <MeshTools/MeshSelectionT.hh>
#include <QtOpenGL>
#include <functional>
#include <stdexcept>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <gmm/gmm.h>
#include <vector>
#include <float.h>
#include <cmath>

class PatchPreview {
public:
    /*
     * is the same as:
     * Crossfield(TriMesh &trimesh, std::vector<int> &heInRange) {
     *      trimesh_ = trimesh;
     *      heInRange_ = heInRange;
     * }
     */
    PatchPreview(TriMesh &trimesh)
            : trimesh_{trimesh} {
    }

    ~PatchPreview() {
    }

    void getCurvature();

private:

    TriMesh &trimesh_;

    void get2RingFaceNeighbours(std::vector<OpenMesh::FaceHandle> &twoRingFaces, OpenMesh::VertexHandle vh);

    void getParams();

    void calculateCurvature();

    double calculateArea(std::vector<OpenMesh::FaceHandle> twoRingfaces);

    std::vector<OpenMesh::HalfedgeHandle> getTwoRingHe(std::vector<OpenMesh::FaceHandle> twoRingFaces);

    OpenMesh::VertexHandle getSelectedVertex();
};


#endif //OPENFLIPPER_PATCHPREVIEW_HH
