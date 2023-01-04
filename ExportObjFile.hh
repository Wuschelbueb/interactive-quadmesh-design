//
// Created by wuschelbueb on 04.01.23.
//

#ifndef OPENFLIPPER_EXPORTOBJFILE_HH
#define OPENFLIPPER_EXPORTOBJFILE_HH

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

class ExportObjFile {
public:
    ExportObjFile(TriMesh &trimesh, OpenMesh::VertexHandle &selectedVertex)
        : trimesh_{trimesh}, centerVertex{selectedVertex} {
        getFile();
    }
private:
    TriMesh &trimesh_;
    OpenMesh::VertexHandle centerVertex;
    double lambda = 0;
    void getFile();
    void createFile();
    double transformToTexCoord(double value);
    void getLambda();

};


#endif //OPENFLIPPER_EXPORTOBJFILE_HH
