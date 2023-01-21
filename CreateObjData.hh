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


/**
 * class which creates an .obj file for the visualization.
 * creates, vertices, normals, texcoords, faces,
 * center of selection, minimum distance to nearest vertex of center
 * which then can be extracted by the openGL tool
 */
class CreateObjData {
public:
    CreateObjData(TriMesh &trimesh, OpenMesh::VertexHandle &selectedVertex)
            : trimesh_{trimesh}, centerVertex{selectedVertex} {
        getStream();
    }

    /**
     * is a public function which saves the dataStream as a string.
     * @param data pass by reference variable which save the dataStream as string
     */
    void getStream(std::string &data) {
        objData = dataStream.str();
        if (!objData.empty()) {
            objData.pop_back();
        }
        data += objData;
    }

    ~CreateObjData() {
    }

private:
    TriMesh &trimesh_;
    // vertex which is at the center of selection
    OpenMesh::VertexHandle centerVertex;
    // string of data
    std::string objData;
    // used to create data for .obj file
    std::stringstream dataStream;

    /**
     * gets called when object gets initialized.
     * requests all the necessary properties and
     * calls createData
     */
    void getStream();

    /**
     * creates the header comments of .obj file
     * calls the other functions as well
     */
    void createData();

    /**
     * resets the status of all vertices
     */
    void resetVertexStatus();

    /**
     * writes the vertex data to the stream
     * @param mapVertexIdx assigns nummeration for vertices
     * which later get used by writeFaces function
     */
    void writeVertices(std::map<int, int> &mapVertexIdx);


    /**
     * writes faces data to stream
     * uses both maps to:
     * 1. determine the order in which the faces have to be written
     * 2. which texCoords has to be assigned to which vertex
     * @param mapVertexToIdx contains the enummeration for vertices
     * @param mapTexCoordToIdx contains the enummeration for TexCoords
     */
    void writeFaces(const std::map<int, int> &mapVertexToIdx);
};


#endif //OPENFLIPPER_CREATEOBJDATA_HH
