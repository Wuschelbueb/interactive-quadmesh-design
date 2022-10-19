#include "MastersThesisPlugin.hh"
#include "DijkstraDistance.hh"
#include "Crossfield.hh"
#include "GlobalParametrization.hh"
#include "Get2DTexture.h"
#include "PatchPreview.hh"

void MastersThesisPlugin::initializePlugin() {
    tool_ = new MastersThesisToolbar();
    QSize size(300, 300);
    tool_->resize(size);

    connect(tool_->get_selection, SIGNAL(clicked()), this, SLOT(slot_get_boundary()));
    connect(tool_->getDualGraph, SIGNAL(clicked()), this, SLOT(slot_get_crossfield()));
    connect(tool_->getGlobalParam, SIGNAL(clicked()), this, SLOT(slot_get_global_param()));
    connect(tool_->get2DTexture, SIGNAL(clicked()), this, SLOT(slot_get_2d_texture()));

    emit addToolbox(tr("MastersThesis"), tool_);

}

void MastersThesisPlugin::pluginsInitialized() {
    emit addTexture(texture_name(), "quadTexture.png", 2);
    emit setTextureMode(texture_name(), "clamp=false,center=false,repeat=true,type=halfedgebased");
    emit switchTexture(texture_name());
}

void MastersThesisPlugin::slot_get_boundary() {
    const double refDist = tool_->dijkstra_distance->value();
    const bool inclBoundaryF = tool_->include_boundary_faces->isChecked();
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH);
         o_it != PluginFunctions::objectsEnd(); ++o_it) {
        // create mesh
        TriMeshObject *tri_obj = PluginFunctions::triMeshObject(*o_it);
        TriMesh *trimesh = tri_obj->mesh();
        // different creation of mesh: TriMesh *trimesh = PluginFunctions::triMesh(*o_it);

        if (trimesh) {
            DijkstraDistance dijkDistMesh{*trimesh};
//            PatchPreview patch{*trimesh};
//            patch.getCurvature();
            dijkDistMesh.cleanMeshOfProps();
            heConstraints = dijkDistMesh.getHeConstraints();
            std::vector<int> includedNodes = dijkDistMesh.calculateDijkstra(heConstraints, refDist, inclBoundaryF);
            includedHEdges = dijkDistMesh.getHeVectorOfSelection(includedNodes);
            dijkDistMesh.colorizeEdges(includedHEdges);
            // change layer of display
            // set draw mode
            PluginFunctions::triMeshObject(*o_it)->meshNode()->drawMode(
                    ACG::SceneGraph::DrawModes::SOLID_SMOOTH_SHADED | ACG::SceneGraph::DrawModes::EDGES_COLORED);
            emit updatedObject(tri_obj->id(), UPDATE_ALL);
        }
    }
}

void MastersThesisPlugin::slot_get_crossfield() {
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH);
         o_it != PluginFunctions::objectsEnd(); ++o_it) {
        // create mesh
        TriMeshObject *tri_obj = PluginFunctions::triMeshObject(*o_it);
        TriMesh *trimesh = tri_obj->mesh();
        if (trimesh) {
            Crossfield mesh{*trimesh, includedHEdges, heConstraints};
            mesh.getCrossfield();
            PluginFunctions::triMeshObject(*o_it)->meshNode()->drawMode(ACG::SceneGraph::DrawModes::EDGES_COLORED);
            emit updatedObject(o_it->id(), UPDATE_ALL);
        }
    }
}

void MastersThesisPlugin::slot_get_global_param() {
    double hValue = tool_->hValue->value();
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH);
         o_it != PluginFunctions::objectsEnd(); ++o_it) {
        // create mesh
        TriMeshObject *tri_obj = PluginFunctions::triMeshObject(*o_it);
        TriMesh *trimesh = tri_obj->mesh();
        if (trimesh) {
            GlobalParametrization mesh{*trimesh, hValue};
            mesh.getGlobalParam();
            PluginFunctions::triMeshObject(*o_it)->meshNode()->drawMode(ACG::SceneGraph::DrawModes::EDGES_COLORED);
            emit updatedObject(o_it->id(), UPDATE_ALL);
        }
    }
}

void MastersThesisPlugin::slot_get_2d_texture() {
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH);
         o_it != PluginFunctions::objectsEnd(); ++o_it) {
        // create mesh
        TriMeshObject *tri_obj = PluginFunctions::triMeshObject(*o_it);
        TriMesh *trimesh = tri_obj->mesh();
        if (trimesh) {
            auto quadTextr = OpenMesh::HProp<OpenMesh::Vec2d>(*trimesh, "quadTextr");
            const char *propertyName = quadTextr.getName().c_str();

            Get2DTexture mesh{*trimesh};
            mesh.initProperty();
            for (auto he: trimesh->halfedges()) {
                if (!he.is_boundary()) {
                    mesh.get2DTexture(he);
                }
            }
            //todo doesn't work, fix
            emit switchTexture(propertyName, o_it->id());
            emit updatedTextures(propertyName, o_it->id());
            PluginFunctions::triMeshObject(*o_it)->meshNode()->drawMode(
                    ACG::SceneGraph::DrawModes::SOLID_2DTEXTURED_FACE_SHADED);
            emit updatedObject(o_it->id(), UPDATE_ALL);
        }
    }
}

#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(tetctplugin, TetCTPlugin);
#endif