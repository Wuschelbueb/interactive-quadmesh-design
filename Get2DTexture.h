//
// Created by wuschelbueb on 27.06.22.
//

#ifndef OPENFLIPPER_GET2DTEXTURE_H
#define OPENFLIPPER_GET2DTEXTURE_H

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <MeshTools/MeshSelectionT.hh>
#include <ACG/Geometry/Types/PlaneT.hh>
#include <ACG/Scenegraph/TransformNode.hh>
#include <ACG/Utils/ColorCoder.hh>
#include <CoMISo/Solver/ConstrainedSolver.hh>
#include <gmm/gmm.h>
#include <functional>
#include <iostream>
#include <vector>
#include <float.h>
#include <cmath>
#include "DijkstraDistance.hh"


class Get2DTexture {
public:
    using Point = OpenMesh::Vec3d;
//    typedef ACG::Vec2d  TexCoord2D;
    typedef gmm::wsvector<double> CVectorType;
public:

    /*
     * is the same as:>
     * Crossfield(TriMesh &trimesh, std::vector<int> &heInRange) {
     *      trimesh_ = trimesh;
     *      heInRange_ = heInRange;
     * }
     */
    Get2DTexture(TriMesh &trimesh)
            : trimesh_{trimesh} {
    }

    ~Get2DTexture() {
    }

    /**
     * initialize property to zero.\n
     * @param hp_texture property
     */
    void initProperty();

    /**
     * get vertex position in u,v system.
     * @param he halfedge handle
     * @param u u parameter of u,v system
     * @param v v parameter of u,v system
     */
    void get2DTexture(OpenMesh::SmartHalfedgeHandle he);


private:

    TriMesh &trimesh_;

    /**
     * get column position in constraint matrix.\n
     * is used when a vertex has multiple occurrences.\n
     * returns int which indicates which position in vector should be accessed
     * @param vh vertex handle
     * @param cutGraphZone face property
     * @return int
     */
    int getPositionConstraintRow(OpenMesh::SmartVertexHandle &vh, int cutGraphZone);

};


#endif //OPENFLIPPER_GET2DTEXTURE_H
