#include "MastersThesisPlugin.hh"
#include "DijkstraDistance.hh"
#include "Crossfield.hh"
#include "GlobalParametrization.hh"

typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;

void MastersThesisPlugin::initializePlugin() {
    tool_ = new MastersThesisToolbar();
    QSize size(300, 300);
    tool_->resize(size);

    connect(tool_->get_selection, SIGNAL(clicked()), this, SLOT(slot_get_boundary()));
    connect(tool_->getDualGraph, SIGNAL(clicked()), this, SLOT(slot_get_dualGraph()));

    emit addToolbox(tr("MastersThesis"), tool_);

}

void MastersThesisPlugin::pluginsInitialized() {}

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
            heConstraints = dijkDistMesh.getHeConstraints();
            std::vector<int> includedNodes = dijkDistMesh.calculateDijkstra(heConstraints, refDist, inclBoundaryF);
            includedHEdges = dijkDistMesh.getHeVectorOfSelection(includedNodes);
            dijkDistMesh.colorizeEdges(includedHEdges);
            // change layer of display
            PluginFunctions::triMeshObject(*o_it)->meshNode()->drawMode(ACG::SceneGraph::DrawModes::EDGES_COLORED);
            emit updatedObject(o_it->id(), UPDATE_ALL);
        }
    }

}

void MastersThesisPlugin::slot_get_dualGraph() {
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
            GlobalParametrization mesh2{*trimesh};
            mesh2.getGlobalParam();
        }
    }
}

#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(tetctplugin, TetCTPlugin);
#endif