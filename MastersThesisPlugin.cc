#include "MastersThesisPlugin.hh"
#include "DijkstraDistance.hh"
#include "Crossfield.hh"
#include "GlobalParametrization.hh"
#include "Get2DTexture.h"

void MastersThesisPlugin::initializePlugin() {
    tool_ = new MastersThesisToolbar();
    QSize size(300, 300);
    tool_->resize(size);

    connect(tool_->get_selection, SIGNAL(clicked()), this, SLOT(slot_get_boundary()));
    connect(tool_->getDualGraph, SIGNAL(clicked()), this, SLOT(slot_get_dualGraph()));
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
            dijkDistMesh.cleanMeshOfProps();
            emit updateView();
            emit updatedObject(o_it->id(), UPDATE_ALL);
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
        }
    }
}

void MastersThesisPlugin::slot_get_global_param() {
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH);
         o_it != PluginFunctions::objectsEnd(); ++o_it) {
        // create mesh
        TriMeshObject *tri_obj = PluginFunctions::triMeshObject(*o_it);
        TriMesh *trimesh = tri_obj->mesh();
        if (trimesh) {
            GlobalParametrization mesh{*trimesh};
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
            Get2DTexture mesh{*trimesh};
            std::cout << texture_name() << std::endl;

            // copy scale texture coordinates to mesh property and set view mode
//          double scale = tool_->uvScaleSpinBox->value();
            OpenMesh::HPropHandleT<ACG::Vec2d> hp_texture;
            if (!trimesh->get_property_handle(hp_texture, texture_name())) {
                trimesh->add_property(hp_texture, texture_name());
            }
            mesh.initProperty(hp_texture);
//            mesh.initProperty();
//          double u_off = tool_->u_offset_sb->value();
//          double v_off = tool_->v_offset_sb->value();
            for (TriMesh::FaceIter f_it = trimesh->faces_begin(); f_it != trimesh->faces_end(); ++f_it) {
                for (TriMesh::FaceHalfedgeIter fh_it = trimesh->fh_iter(*f_it); fh_it.is_valid(); ++fh_it) {
                    double u, v;
//                    std::cout << "new iteration with he: " << fh_it->idx() << "\n";
                    mesh.get2DTexture(fh_it, u, v);
//                    std::cout << "halfedge " << fh_it->idx() << "\tu " << u << "\tv " << v << std::endl;
                    trimesh->property(hp_texture, *fh_it) = {u, v};
                }
            }
//            auto test = OpenMesh::getProperty<OpenMesh::HalfedgeHandle, OpenMesh::Vec3f>(trimesh, "quad");
            emit switchTexture(texture_name(), o_it->id());
            emit updatedTextures(texture_name(), o_it->id());

            PluginFunctions::setDrawMode(ACG::SceneGraph::DrawModes::SOLID_2DTEXTURED_FACE_SHADED);
        }
    }
}

#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(tetctplugin, TetCTPlugin);
#endif