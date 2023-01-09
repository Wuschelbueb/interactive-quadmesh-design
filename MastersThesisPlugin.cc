#include "MastersThesisPlugin.hh"
#include "DijkstraDistance.hh"
#include "Crossfield.hh"
#include "GlobalParametrization.hh"
#include "Get2DTexture.hh"
#include "PatchPreview.hh"
#include "CreateObjData.hh"
#include "ACG/Scenegraph/LineNode.hh"


void MastersThesisPlugin::initializePlugin() {
    tool_ = new MastersThesisToolbar();
    QSize size(300, 300);
    tool_->resize(size);

    connect(tool_->calculationButton, SIGNAL(clicked()), this, SLOT(slot_calculate_quad_mesh()));
    connect(tool_->selectionButton, SIGNAL(clicked()), this, SLOT(slot_select_point()));
    connect(tool_->previewButton, SIGNAL(clicked()), this, SLOT(slot_get_preview_dijkstra()));
    connect(tool_->saveObjectFileButton, SIGNAL(clicked()), this, SLOT(slot_save_object_file()));

    emit addToolbox(tr("MastersThesis"), tool_);

}

void MastersThesisPlugin::pluginsInitialized() {
    emit addTexture("quadTextr", "quadTexture.png", 2);
    emit setTextureMode("quadTextr", "clamp=false,center=false,repeat=true,type=halfedgebased");
    emit addPickMode("Vertex Selection");

}

void MastersThesisPlugin::slot_select_point() {
    if (tool_->selectionButton->isChecked()) {
        // Picking mode of our plugin shall be activated
        // set OpenFlipper's action mode to picking
        PluginFunctions::actionMode(Viewer::PickingMode);
        // Now activate our picking mode
        PluginFunctions::pickMode("Vertex Selection");
    } else {
        // Picking mode shall be deactivated
        PluginFunctions::actionMode(Viewer::ExamineMode);
    }
}

void MastersThesisPlugin::slotPickModeChanged(const std::string &_mode) {
    // change visualization of 3D object
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH);
         o_it != PluginFunctions::objectsEnd(); ++o_it) {
        TriMeshObject *tri_obj = PluginFunctions::triMeshObject(*o_it);
        tri_obj->meshNode()->drawMode(ACG::SceneGraph::DrawModes::WIREFRAME
                                      | ACG::SceneGraph::DrawModes::SOLID_SMOOTH_SHADED
                                      | ACG::SceneGraph::DrawModes::POINTS_COLORED);
        emit updatedObject(tri_obj->id(), UPDATE_COLOR);
    }
    // Set button checked if pick mode is our
    // plugin's pick mode
//    tool_->selectionButton->setChecked(_mode == "Vertex Selection");
}

void MastersThesisPlugin::slotMouseEvent(QMouseEvent *_event) {
    // selects vertex which is at center of patch
    if (PluginFunctions::pickMode() == ("Vertex Selection") &&
        PluginFunctions::actionMode() == Viewer::PickingMode &&
        tool_->selectionButton->isChecked()) {
        // handle mouse events
        if (_event->button() == Qt::LeftButton) {
            size_t node_idx, target_idx;
            ACG::Vec3d hit_point;
            // pick vertices
            if (PluginFunctions::scenegraphPick(ACG::SceneGraph::PICK_VERTEX, _event->pos(),
                                                node_idx, target_idx, &hit_point)) {
                BaseObjectData *obj;
                if (PluginFunctions::getPickedObject(node_idx, obj)) {
                    // is picked object a triangle mesh?
                    TriMeshObject *tri_obj = PluginFunctions::triMeshObject(obj);

                    if (tri_obj) {
                        auto targetedVh = tri_obj->mesh()->vertex_handle(target_idx);
                        selectedVertex = targetedVh;
                        selectedVertexAsPoint = tri_obj->mesh()->point(targetedVh);
                        if (targetedVh == TriMesh::InvalidVertexHandle) {
                            return;
                        }
                        // size of selection
                        tri_obj->materialNode()->set_point_size(18);
                        // updates vertex selection and deselects old one
                        for (auto vh: tri_obj->mesh()->vertices()) {
                            if ((int) target_idx == vh.idx()) {
                                tri_obj->mesh()->status(vh).set_selected(true);
                                tri_obj->mesh()->set_color(vh, ACG::Vec4f(0, 1, 0, 1));
                            } else {
                                tri_obj->mesh()->status(vh).set_selected(false);
                                tri_obj->mesh()->set_color(vh, ACG::Vec4f(1, 1, 1, 0));
                            }
                        }
                        // visualization
                        tri_obj->meshNode()->drawMode(ACG::SceneGraph::DrawModes::WIREFRAME
                                                      | ACG::SceneGraph::DrawModes::SOLID_SMOOTH_SHADED |
                                                      ACG::SceneGraph::DrawModes::POINTS_COLORED);

                        tri_obj->materialNode()->enable_alpha_test(0.8);
                        // updates color of vertices
                        emit updatedObject(tri_obj->id(), UPDATE_COLOR);
                        return;
                    }
                }
            }
        }
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
                    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH);
                         o_it != PluginFunctions::objectsEnd(); ++o_it) {
                        auto tri_obj = PluginFunctions::triMeshObject(*o_it);
                        auto trimesh = tri_obj->mesh();
                        ACG::SceneGraph::LineNode *lineNode;
                        //create line node    connect(tool_->calculationButton, SIGNAL(clicked()), this, SLOT(slot_calculate_quad_mesh()));

                        if (!tri_obj->getAdditionalNode(lineNode, name(), "Cross field direction")) {
                            lineNode = new ACG::SceneGraph::LineNode(ACG::SceneGraph::LineNode::LineSegmentsMode,
                                                                     tri_obj->manipulatorNode(),
                                                                     "Cross field direction");
                            tri_obj->addAdditionalNode(lineNode, name(), "Cross field direction");

                            //creates the line
                            lineNode->clear_points();
                            lineNode->set_color(OpenMesh::Vec4f(1.0f, 0.0f, 0.0f, 1.0f));
                            lineNode->set_line_width(3);
                            lineNode->add_line(selectedVertexAsPoint, hitPoint);
                            lineNode->alwaysOnTop() = true;
                            clickedPoint = hitPoint;
                        } else {
                            //creates the line
                            lineNode->clear_points();
                            lineNode->set_color(OpenMesh::Vec4f(1.0f, 0.0f, 0.0f, 1.0f));
                            lineNode->set_line_width(3);
                            lineNode->add_line(selectedVertexAsPoint, hitPoint);
                            lineNode->alwaysOnTop() = true;
                            clickedPoint = hitPoint;
                        }
                    }
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
    emit updateView();
}

void MastersThesisPlugin::slotUpdateTexture(QString _textureName, int _identifier) {
    if (_textureName != "quadTextr") {
        return;
    }

    BaseObjectData *object;
    if (!PluginFunctions::getObject(_identifier, object)) {
        return;
    }

    if (object->dataType(DATA_TRIANGLE_MESH)) {
        TriMesh *mesh = PluginFunctions::triMesh(object);
        if (_textureName == "quadTextr") {
            emit updatedTextures("quadTextr", _identifier);
        }
    }
    if (object->dataType(DATA_POLY_MESH)) {
        TriMesh *mesh = PluginFunctions::triMesh(object);
        if (_textureName == "quadTextr") {
            emit updatedTextures("quadTextr", _identifier);
        }
    }
}

void MastersThesisPlugin::slot_get_preview_dijkstra() {
    const double refDist = tool_->dijkstra_distance->value();
//    const bool inclBoundaryF = tool_->include_boundary_faces->isChecked();
    const bool inclBoundaryF = false;
    PluginFunctions::actionMode(Viewer::ExamineMode);
    tool_->selectionButton->setChecked(false);
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
            originHalfedges = dijkDistMesh.getHeFromVertex(selectedVertex);
            std::vector<int> includedFaces = dijkDistMesh.calculateDijkstra(originHalfedges, refDist, inclBoundaryF);
            includedHalfedges = dijkDistMesh.getAllHeFromFaces(includedFaces);
            dijkDistMesh.colorizeEdges(includedHalfedges);
            PluginFunctions::triMeshObject(*o_it)->meshNode()->drawMode(
                    ACG::SceneGraph::DrawModes::SOLID_SMOOTH_SHADED | ACG::SceneGraph::DrawModes::EDGES_COLORED |
                    ACG::SceneGraph::DrawModes::POINTS_COLORED);
            emit updatedObject(tri_obj->id(), UPDATE_ALL);
        }
    }
}

void MastersThesisPlugin::slot_calculate_quad_mesh() {
    double hValue = tool_->hValue->value();
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH);
         o_it != PluginFunctions::objectsEnd(); ++o_it) {
        // create mesh
        TriMeshObject *tri_obj = PluginFunctions::triMeshObject(*o_it);
        TriMesh *trimesh = tri_obj->mesh();
        // different creation of mesh: TriMesh *trimesh = PluginFunctions::triMesh(*o_it);
        PluginFunctions::actionMode(Viewer::ExamineMode);
        auto quadTextr = OpenMesh::HProp<OpenMesh::Vec2d>(*trimesh, "quadTextr");
        auto heColor = OpenMesh::HProp<int>(*trimesh, "heColor");
        refVector = clickedPoint - selectedVertexAsPoint;

        if (trimesh) {
            Crossfield crossFieldMesh{*trimesh, includedHalfedges, originHalfedges, refVector};
            crossFieldMesh.getCrossfield();
            ACG::SceneGraph::LineNode *lineNode;
            if (tri_obj->getAdditionalNode(lineNode, name(), "Cross field direction")) {
                lineNode = new ACG::SceneGraph::LineNode(ACG::SceneGraph::LineNode::LineSegmentsMode,
                                                         tri_obj->manipulatorNode(),
                                                         "Cross field direction");
                tri_obj->removeAdditionalNode(lineNode, name(), "Cross field direction");
            }

            GlobalParametrization gpMesh{*trimesh, hValue};
            gpMesh.getGlobalParam();

            Get2DTexture twoDimMesh{*trimesh};
            twoDimMesh.initProperty();
            for (auto he: trimesh->halfedges()) {
                if (!he.is_boundary() && heColor[he] != 1) {
                    twoDimMesh.setQuadTexHeProperty(he);
                }
            }
            CreateObjData objDataCreation{*trimesh, selectedVertex};
            objDataCreation.getStream(objData);

            // switch to the right texture
            emit switchTexture("quadTextr", o_it->id());
            // select texture from menu
            tri_obj->setObjectDrawMode(ACG::SceneGraph::DrawModes::SOLID_2DTEXTURED_FACE_SHADED,
                                       PluginFunctions::ALL_VIEWERS);
        }
    }
}


void MastersThesisPlugin::slot_save_object_file() {
    // Save selection button has been clicked

    QString filename = QFileDialog::getSaveFileName(0, tr("Save As"), "object.obj", tr(".obj ( *.obj )"));

    if (filename.isEmpty()) {
        return;
    } else {
        QFile file(filename);
        if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
            QMessageBox::information(0, tr("Unable to open file"),
                                     file.errorString());
            return;
        }
        QTextStream out(&file);
        out << QString::fromStdString(objData);
        file.close();
    }
}

#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(mastersthesisplugin, MastersThesisPlugin
);
#endif