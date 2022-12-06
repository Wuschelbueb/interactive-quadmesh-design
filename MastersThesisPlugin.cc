#include "MastersThesisPlugin.hh"
#include "DijkstraDistance.hh"
#include "Crossfield.hh"
#include "GlobalParametrization.hh"
#include "Get2DTexture.h"
#include "PatchPreview.hh"
#include "myShaderNode.hh"


void MastersThesisPlugin::initializePlugin() {
    tool_ = new MastersThesisToolbar();
    QSize size(300, 300);
    tool_->resize(size);

    connect(tool_->get_selection, SIGNAL(clicked()), this, SLOT(slot_get_boundary()));
    connect(tool_->getDualGraph, SIGNAL(clicked()), this, SLOT(slot_get_crossfield()));
    connect(tool_->getGlobalParam, SIGNAL(clicked()), this, SLOT(slot_get_global_param()));
    connect(tool_->get2DTexture, SIGNAL(clicked()), this, SLOT(slot_get_2d_texture()));
    connect(tool_->testButton, SIGNAL(clicked()), this, SLOT(slot_select_point()));

    emit addToolbox(tr("MastersThesis"), tool_);

}

void MastersThesisPlugin::pluginsInitialized() {
    emit addTexture(texture_name(), "quadTexture.png", 2);
    emit setTextureMode(texture_name(), "clamp=false,center=false,repeat=true,type=halfedgebased");
    emit switchTexture(texture_name());
    emit addPickMode("MastersThesisPlugin");
}

void MastersThesisPlugin::slot_select_point() {
    if (tool_->testButton->isChecked()) {
        // Picking mode of our plugin shall be activated
        // set OpenFlipper's action mode to picking
        PluginFunctions::actionMode(Viewer::PickingMode);
        // Now activate our picking mode
        PluginFunctions::pickMode("MastersThesisPlugin");
    } else {
        // Picking mode shall be deactivated
        PluginFunctions::actionMode(Viewer::ExamineMode);
    }
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
            // calculates the curvature of a patch as a preview
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
//            PluginFunctions::triMeshObject(*o_it)->meshNode()->drawMode(
//                    ACG::SceneGraph::DrawModes::NONE);
//            PluginFunctions::triMeshObject(*o_it)->myShaderNode()->drawMode() //doesn't work
            emit updatedObject(tri_obj->id(), UPDATE_ALL);
        }
    }

}

void MastersThesisPlugin::slotPickModeChanged(const std::string &_mode) {
    // Set button checked if pick mode is our
    // plugin's pick mode
    tool_->testButton->setChecked(_mode == "MastersThesisPlugin");
}

void MastersThesisPlugin::slotMouseEvent(QMouseEvent *_event) {
    //figure out where identification mode is
    //todo use identification mode to extract (selected) point!
    //todo create a vector between selected vertex and extracted point
    if (PluginFunctions::pickMode() == "MastersThesisPlugin" &&
        PluginFunctions::actionMode() == Viewer::PickingMode) {
        // If double click has been performed
        if (_event->type() == QEvent::MouseButtonDblClick) {
            size_t node_idx, target_idx;
            OpenMesh::Vec3d hitPoint;
            // Get picked object's identifier
            if (PluginFunctions::scenegraphPick(ACG::SceneGraph::PICK_ANYTHING, _event->pos(), node_idx,
                                                target_idx, &hitPoint)) {
                BaseObjectData *object;
                // Get picked object
                if (PluginFunctions::getPickedObject(node_idx, object)) {
                    // Show small dialog window;
                    QDialog *dlg = new QDialog;
                    QGridLayout *grid = new QGridLayout;
                    QLabel *label = new QLabel;
                    QString str = QString("Point: [");
                    str += QString::number(hitPoint[0]);
                    str += QString(", ");
                    str += QString::number(hitPoint[1]);
                    str += QString(", ");
                    str += QString::number(hitPoint[2]);
                    str += QString("]");
                    clickedPoint = hitPoint;
                    label->setText(str);
                    grid->addWidget(label, 0, 0);
                    dlg->setLayout(grid);
                    dlg->show();
                    // Set last selected object
                    activeObject_ = node_idx;
                } else {
                    emit log(LOGINFO, "Picking failed");
                }
            }
            return;
        }
        // Continue traversing scene graph
        ACG::SceneGraph::MouseEventAction action(_event, PluginFunctions::viewerProperties().glState());
        PluginFunctions::traverse(action);
    }
}

void MastersThesisPlugin::slot_get_crossfield() {
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH);
         o_it != PluginFunctions::objectsEnd(); ++o_it) {
        // create mesh
        TriMeshObject *tri_obj = PluginFunctions::triMeshObject(*o_it);
        TriMesh *trimesh = tri_obj->mesh();
        if (trimesh) {
            std::vector<int> vertices = MeshSelection::getVertexSelection(trimesh);
            OpenMesh::VertexHandle vh = trimesh->vertex_handle(vertices[0]);
            refVector = clickedPoint - trimesh->point(vh);
            Crossfield mesh{*trimesh, includedHEdges, heConstraints, refVector};
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
Q_EXPORT_PLUGIN2(mastersthesisplugin, MastersThesisPlugin);
#endif