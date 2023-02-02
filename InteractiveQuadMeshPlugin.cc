#include "InteractiveQuadMeshPlugin.hh"
#include "DijkstraDistance.hh"
#include "Crossfield.hh"
#include "GlobalParametrization.hh"
#include "Get2DTexture.hh"
#include "PatchPreview.hh"
#include "CreateObjData.hh"
#include "ACG/Scenegraph/LineNode.hh"


void InteractiveQuadMeshPlugin::initializePlugin() {
    tool_ = new InteractiveQuadMeshToolbar();
    QSize size(300, 300);
    tool_->resize(size);

    connect(tool_->calculationButton, SIGNAL(clicked()), this, SLOT(slot_calculate_quad_mesh()));
    connect(tool_->selectionButton, SIGNAL(clicked()), this, SLOT(slot_select_point()));
    connect(tool_->previewButton, SIGNAL(clicked()), this, SLOT(slot_get_preview_dijkstra()));
    connect(tool_->saveObjectFileButton, SIGNAL(clicked()), this, SLOT(slot_save_object_file()));

    emit addToolbox(tr("Interactive Quad Mesh Design"), tool_);

}

void InteractiveQuadMeshPlugin::pluginsInitialized() {
    emit addTexture("quadTextr", "quadTexture.png", 2);
    emit setTextureMode("quadTextr", "clamp=false,center=false,repeat=true,type=halfedgebased");
    emit addPickMode("Vertex Selection");

}

void InteractiveQuadMeshPlugin::slot_select_point() {
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

void InteractiveQuadMeshPlugin::slotPickModeChanged(const std::string &_mode) {
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

void InteractiveQuadMeshPlugin::slotMouseEvent(QMouseEvent *_event) {
    // selects vertex which is at center of patch
    if (PluginFunctions::pickMode() == ("Vertex Selection") &&
        PluginFunctions::actionMode() == Viewer::PickingMode &&
        tool_->selectionButton->isChecked()) {
        // handle mouse events

        // If double click has been performed
        if (_event->button() == Qt::LeftButton && (_event->modifiers() & Qt::ControlModifier)) {
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
                        ACG::SceneGraph::LineNode *lineNode;
                        //create line node    connect(tool_->calculationButton, SIGNAL(clicked()), this, SLOT(slot_calculate_quad_mesh()));

                        if (!tri_obj->getAdditionalNode(lineNode, name(), "Cross field direction")) {
                            lineNode = new ACG::SceneGraph::LineNode(ACG::SceneGraph::LineNode::LineSegmentsMode,
                                                                     tri_obj->manipulatorNode(),
                                                                     "Cross field direction");
                            tri_obj->addAdditionalNode(lineNode, name(), "Cross field direction");

                            //creates the line
                            lineNode->clear_points();
                            lineNode->set_color(OpenMesh::Vec4f(0.0f, 0.0f, 1.0f, 1.0f));
                            lineNode->set_line_width(3);
                            lineNode->add_line(selectedOriginPoint, hitPoint);
                            lineNode->alwaysOnTop() = true;
                            clickedPointU = hitPoint;
                        } else {
                            //creates the line
                            lineNode->clear_points();
                            lineNode->set_color(OpenMesh::Vec4f(0.0f, 0.0f, 1.0f, 1.0f));
                            lineNode->set_line_width(3);
                            lineNode->add_line(selectedOriginPoint, hitPoint);
                            lineNode->alwaysOnTop() = true;
                            clickedPointU = hitPoint;
                        }
                    }
                } else {
                    emit log(LOGINFO, "Picking failed");
                }
            }
            return;
        }
        if (_event->button() == Qt::LeftButton && (_event->modifiers() & Qt::ShiftModifier)) {
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
                        ACG::SceneGraph::LineNode *lineNode;
                        //create line node    connect(tool_->calculationButton, SIGNAL(clicked()), this, SLOT(slot_calculate_quad_mesh()));

                        if (!tri_obj->getAdditionalNode(lineNode, name(), "Second Direction")) {
                            lineNode = new ACG::SceneGraph::LineNode(ACG::SceneGraph::LineNode::LineSegmentsMode,
                                                                     tri_obj->manipulatorNode(),
                                                                     "Second Direction");
                            tri_obj->addAdditionalNode(lineNode, name(), "Second Direction");

                            //creates the line
                            lineNode->clear_points();
                            lineNode->set_color(OpenMesh::Vec4f(1.0f, 0.0f, 0.0f, 1.0f));
                            lineNode->set_line_width(3);
                            lineNode->add_line(selectedOriginPoint, hitPoint);
                            lineNode->alwaysOnTop() = true;
                            clickedPointV = hitPoint;
                        } else {
                            //creates the line
                            lineNode->clear_points();
                            lineNode->set_color(OpenMesh::Vec4f(1.0f, 0.0f, 0.0f, 1.0f));
                            lineNode->set_line_width(3);
                            lineNode->add_line(selectedOriginPoint, hitPoint);
                            lineNode->alwaysOnTop() = true;
                            clickedPointV = hitPoint;
                        }
                    }
                } else {
                    emit log(LOGINFO, "Picking failed");
                }
            }
            return;
        }

        if (_event->button() == Qt::LeftButton) {
            size_t node_idx, target_idx;
            ACG::Vec3d hit_point;
            // pick vertices
            if (PluginFunctions::scenegraphPick(ACG::SceneGraph::PICK_ANYTHING,
                                                _event->pos(), node_idx, target_idx, &hit_point)) {
                BaseObjectData *object(0);
                PluginFunctions::getPickedObject(node_idx, object);
                if (!object) return;
                TriMeshObject *tri_obj = PluginFunctions::triMeshObject(object);

                if (tri_obj) {
                    if (PluginFunctions::scenegraphPick(ACG::SceneGraph::PICK_FACE, _event->pos(), node_idx, target_idx,
                                                        &hit_point)) {
                        auto targetedFh = make_smart(tri_obj->mesh()->face_handle(target_idx), tri_obj->mesh());
                        float distance = INT_MAX;
                        selectedOriginPoint = hit_point;
                        originVertices.clear();
                        for (auto he: targetedFh.halfedges()) {
                            originVertices.push_back(he.to().idx());
                            ACG::Vec3d distVector = selectedOriginPoint - tri_obj->mesh()->point(he.to());
                            if (distVector.norm() < distance) {
                                distance = distVector.norm();
                                selectedVertex = he.to();
                            }
                        }

//                        std::cout << "selectedOrigin " << selectedOriginPoint << std::endl;
                        // visualization
                        tri_obj->meshNode()->drawMode(ACG::SceneGraph::DrawModes::WIREFRAME
                                                      | ACG::SceneGraph::DrawModes::SOLID_SMOOTH_SHADED);

                        tri_obj->materialNode()->enable_alpha_test(0.8);
                        // updates color of vertices
                        emit updatedObject(tri_obj->id(), UPDATE_COLOR);
                        return;
                    }
                }
            }

        }
        // Continue traversing scene graph
        ACG::SceneGraph::MouseEventAction action(_event, PluginFunctions::viewerProperties().glState());
        PluginFunctions::traverse(action);
    }
    emit updateView();
}

void InteractiveQuadMeshPlugin::slotUpdateTexture(QString _textureName, int _identifier) {
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

void InteractiveQuadMeshPlugin::slot_get_preview_dijkstra() {
    elementM = tool_->sizeM->value();
    elementN = tool_->sizeN->value();
    refVectorU = clickedPointU - selectedOriginPoint;
    refVectorV = clickedPointV - selectedOriginPoint;
    double angle = acos((refVectorU | refVectorV) / refVectorU.norm() * refVectorV.norm());
    double restAngle;
    // step one calculate length of quads
    if (angle > M_PI / 2) {
        restAngle = angle - M_PI / 2;
    } else {
        restAngle = M_PI / 2 - angle;
    }
    refVectorVLengthAt90Deg = cos(restAngle) * refVectorV.norm();
    double distanceU = elementM * refVectorU.norm();
    double distanceV = elementN * refVectorVLengthAt90Deg;

    if (distanceU > distanceV) {
        refDist = distanceU * 1.5;
    } else {
        refDist = distanceV * 1.5;
    }
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
            originHalfedges = dijkDistMesh.getHeFromVertex(selectedVertex, originVertices);
            std::vector<int> includedFaces = dijkDistMesh.calculateDijkstra(originHalfedges, refDist);
            includedHalfedges = dijkDistMesh.getAllHeFromFaces(includedFaces);
            dijkDistMesh.colorizeEdges(includedHalfedges);
            PluginFunctions::triMeshObject(*o_it)->meshNode()->drawMode(
                    ACG::SceneGraph::DrawModes::SOLID_SMOOTH_SHADED | ACG::SceneGraph::DrawModes::EDGES_COLORED |
                    ACG::SceneGraph::DrawModes::POINTS_COLORED);
            emit updatedObject(tri_obj->id(), UPDATE_ALL);
        }
    }
}

void InteractiveQuadMeshPlugin::slot_calculate_quad_mesh() {
//    const double hValue = tool_->hValue->value();

    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH);
         o_it != PluginFunctions::objectsEnd(); ++o_it) {
        // create mesh
        TriMeshObject *tri_obj = PluginFunctions::triMeshObject(*o_it);
        TriMesh *trimesh = tri_obj->mesh();
        // different creation of mesh: TriMesh *trimesh = PluginFunctions::triMesh(*o_it);
        PluginFunctions::actionMode(Viewer::ExamineMode);
        auto quadTextr = OpenMesh::HProp<OpenMesh::Vec2d>(*trimesh, "quadTextr");
        auto heColor = OpenMesh::HProp<int>(*trimesh, "heColor");

        if (trimesh) {
            objData.clear();
            Crossfield crossFieldMesh{*trimesh, includedHalfedges, originHalfedges, refVectorU};
            crossFieldMesh.getCrossfield();
            std::string initData = setNbOfQuads(elementM, elementN);
            objData += initData;

            GlobalParametrization gpMesh{*trimesh, refVectorU.norm(), refVectorVLengthAt90Deg};
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


void InteractiveQuadMeshPlugin::slot_save_object_file() {
    // Save selection button has been clicked
    QString filename = QFileDialog::getSaveFileName(0, tr("Save As"), "object.obj", tr(".obj ( *.obj )"));

    if (filename.isEmpty()) {
        return;
    } else {
        QFile file(filename);
        if (!file.open(QIODevice::WriteOnly)) {
            QMessageBox::information(0, tr("Unable to open file"),
                                     file.errorString());
            return;
        }
        QTextStream out(&file);
        out << QString::fromStdString(objData);

        file.close();
        file.disconnect();
    }
}

std::string
InteractiveQuadMeshPlugin::setNbOfQuads(const int elementsM, const int elementsN) {
    std::stringstream tempStream;
    tempStream << "# nbQuadsU: " << elementsM * 2 << "\n# nbQuadsV: " << elementsN * 2 << std::endl;
    return tempStream.str();
}

#if QT_VERSION < 0x050000
Q_EXPORT_PLUGIN2(mastersthesisplugin, MastersThesisPlugin
);
#endif