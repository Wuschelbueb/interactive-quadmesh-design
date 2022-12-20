#ifndef OPENFLIPPER_MASTERSTHESISPLUGIN_H
#define OPENFLIPPER_MASTERSTHESISPLUGIN_H

#include <QObject>
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
#include "myShaderNode.hh"
#include "MastersThesisToolbar.hh"

class MastersThesisPlugin
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

    // Texture interface
    // tell OpenFlipper about the texture we want to use
    void addTexture(QString _textureName, QString _filename, uint dimension);

    // tell OpenFlipper that our texture (coordinates) have changed
    void updatedTextures(QString, int);
    void updatedTexture(QString);

    // tell OpenFlipper which texture settings we want to use
    void setTextureMode(QString _textureName, QString _mode);

    // tell OpenFlipper to use the texture with name _textureName
    void switchTexture(QString _textureName);

    /** \brief emit this signal if you want to switch the texture of a specific object
     * This signal can be called from any thread.\n
     */
    void switchTexture(QString _textureName, int _id);


private slots:

    // BaseInterface
    void initializePlugin();
    void pluginsInitialized();

    //MouseInterface
    void slotMouseEvent(QMouseEvent *_event);

    //PickingInterface
    void slotPickModeChanged(const std::string &_mode);

    const char *texture_name() const { return "quadTextr"; }


public :

    ~MastersThesisPlugin() {}

    MastersThesisPlugin() :
            activeObject_(-1) {}

    QString name() { return QString("Simple plugin"); };

    QString description() { return QString("Does actually nothing but works!"); };

public slots:

    /**
     * get first selection with the help of dijkstra algorithm.\n
     */
    void slot_get_boundary();

    /**
     * get crossfield for each triangle.\n
     */
    void slot_get_crossfield();

    /**
     * select a point for crossfield
     */
     void slot_select_point();

    /**
     * create global parametrization for selection.\n
     */
    void slot_get_global_param();

    /**
     * put a quad mesh texture over tri-mesh.\n
     */
    void slot_get_2d_texture();

    QString version() { return QString("1.0"); };


private:
    MastersThesisToolbar *tool_;
    myShaderNode* myShaderNode_;

    //store selected vertices
    std::vector<int> originHalfedges;
    std::vector<int> includedHalfedges;
    ACG::Vec3d clickedPoint;
    OpenMesh::VertexHandle selectedVertex;
    ACG::Vec3d selectedVertexAsPoint;
    ACG::Vec3d refVector;

    // Last picked object
    int activeObject_;

};

#endif //OPENFLIPPER_MASTERSTHESISPLUGIN_H