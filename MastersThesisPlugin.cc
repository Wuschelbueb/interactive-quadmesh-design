#include "MastersThesisPlugin.hh"
#include "DijkstraDistance.hh"
#include "Crossfield.hh"

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
            DijkstraDistance mesh{*trimesh};
            constrained_vertices = mesh.calculateDijkstra(refDist);
            if (inclBoundaryF)
                mesh.includeBoundaryFaces(constrained_vertices, refDist);
            constrained_HEdges = mesh.getHEinRange(constrained_vertices, refDist, inclBoundaryF);
            mesh.colorizeArea(constrained_HEdges);
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
            Crossfield mesh{*trimesh, constrained_HEdges};
            mesh.getCrossfield();
            emit updatedObject(o_it->id(), UPDATE_ALL);
        }
    }
}

#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(tetctplugin, TetCTPlugin);
#endif