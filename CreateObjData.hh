//
// Created by wuschelbueb on 04.01.23.
//

#ifndef OPENFLIPPER_CREATEOBJDATA_HH
#define OPENFLIPPER_CREATEOBJDATA_HH

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

class CreateObjData {
public:
    CreateObjData(TriMesh &trimesh, OpenMesh::VertexHandle &selectedVertex)
            : trimesh_{trimesh}, centerVertex{selectedVertex} {
        getStream();
    }

    void getStream(std::string &data) {
        objData = dataStream.str();
        data = objData;
    }

private:
    TriMesh &trimesh_;
    OpenMesh::VertexHandle centerVertex;
    std::string objData;
    double lambda = 0;
    double oldMin, oldMax;
    std::stringstream dataStream;


    void getStream();

    void createData();

    double transformToTexCoord(double value);

    void getLambda();

    void resetVertexStatus();

    void writeTexCoords(std::map<int, int> &connectVhToTexCoord);

    void writeVertices(std::map<int, int> &connectVhToTexCoord);

    void writeFaces(std::map<int, int> &connectVhToTexCoord);
};


#endif //OPENFLIPPER_CREATEOBJDATA_HH
