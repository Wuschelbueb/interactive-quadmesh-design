#ifndef OPENFLIPPER_MASTERSTHESISPLUGIN_H
#define OPENFLIPPER_MASTERSTHESISPLUGIN_H

#include <QObject>
#include <QMessageBox>
#include <OpenFlipper/common/Types.hh>
#include <OpenFlipper/BasePlugin/MouseInterface.hh>
#include <OpenFlipper/BasePlugin/PickingInterface.hh>
#include <OpenFlipper/BasePlugin/BaseInterface.hh>
#include <OpenFlipper/BasePlugin/ToolboxInterface.hh>
#include <OpenFlipper/BasePlugin/LoggingInterface.hh>
#include <OpenFlipper/BasePlugin/LoadSaveInterface.hh>
#include <OpenFlipper/BasePlugin/TextureInterface.hh>
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>
#include "InteractiveQuadMeshToolbar.hh"

class InteractiveQuadMeshPlugin
        : public QObject,
          BaseInterface,
          MouseInterface,
          ToolboxInterface,
          TextureInterface,
          LoadSaveInterface,
          PickingInterface,
          LoggingInterface {
Q_OBJECT
    Q_INTERFACES(BaseInterface)
    Q_INTERFACES(TextureInterface)
    Q_INTERFACES(ToolboxInterface)
    Q_INTERFACES(LoggingInterface)
    Q_INTERFACES(LoadSaveInterface)
    Q_INTERFACES(MouseInterface)
    Q_INTERFACES(PickingInterface)


#if QT_VERSION >= 0x050000
    Q_PLUGIN_METADATA(IID "org.OpenFlipper.Plugins.Plugin-MastersThesis")
#endif

signals:

    //BaseInterface
    void updateView();

    //PickingInterface
    void addPickMode(const std::string &_mode);

    void addHiddenPickMode(const std::string &_mode);

    //LoggingInterface
    void log(Logtype _type, QString _message);

    void log(QString _message);

    // ToolboxInterface
    void addToolbox(QString _name, QWidget *_widget);

    // LoadSaveInterface
    void addEmptyObject(DataType _type, int &_id);

    // BackupInterface
    void createBackup(int _objectid, QString _name, UpdateType _type = UPDATE_ALL);

    void updatedObject(int _id, const UpdateType &_type);

    // Texture Interface
    void addTexture(QString _textureName, QString _filename, uint dimension);

    void updatedTextures(QString, int);

    void setTextureMode(QString _textureName, QString _mode);

    void switchTexture(QString _textureName);

    void switchTexture(QString _textureName, int _id);


private slots:

    // BaseInterface
    void initializePlugin();

    void pluginsInitialized();

    // MouseInterface
    void slotMouseEvent(QMouseEvent *_event);

    // Update texture
    void slotUpdateTexture(QString _textureName, int _identifier);

    //PickingInterface
    void slotPickModeChanged(const std::string &_mode);


public :

    ~InteractiveQuadMeshPlugin() {}

    InteractiveQuadMeshPlugin() :
            activeObject_(-1) {}

    QString name() { return QString("Simple plugin"); };

    QString description() { return QString("Does actually nothing but works!"); };

public slots:

    /**
     * get first selection with the help of dijkstra algorithm.\n
     */
    void slot_calculate_quad_mesh();

    /**
     * get first selection with the help of dijkstra algorithm.\n
     */
    void slot_get_preview_dijkstra();

    /**
     * select a point for crossfield
     */
    void slot_select_point();

    /**
     * open dialog to save obj file
     */
    void slot_save_object_file();

    QString version() { return QString("1.0"); };


private:
    // initialize toolbar plugin
    InteractiveQuadMeshToolbar *tool_;
    //store outgoing he of selected Vertex
    std::vector<int> originHalfedges;
    // vector based on included faces
    std::vector<int> includedHalfedges;
    // point which gives direction of crossfield
    ACG::Vec3d clickedPointU;
    ACG::Vec3d clickedPointV;
    // selected starting vertex
    OpenMesh::VertexHandle selectedVertex;
    // point of selectedVertex
    ACG::Vec3d selectedVertexAsPoint;
    // vector between selectedVertex and clickedPoint
    ACG::Vec3d refVectorU;
    ACG::Vec3d refVectorV;
    double refVectorVLengthAt90Deg = 1.0;
    double refDist = 0;
    // stringstream containing .obj data
    std::string objData;

    int elementM = 0;
    int elementN = 0;

    /**
     * helper function which gets the number of quads which can be displayed
     * if the number is higher there is a chance that artifacts pop up.
     */
    std::string setNbOfQuads(const int elementsM, const int elementsN);

    // Last picked object
    int activeObject_;

};

#endif //OPENFLIPPER_MASTERSTHESISPLUGIN_H