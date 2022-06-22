#ifndef OPENFLIPPER_MASTERSTHESISPLUGIN_H
#define OPENFLIPPER_MASTERSTHESISPLUGIN_H

#include <QObject>
#include <OpenFlipper/common/Types.hh>
#include <OpenFlipper/BasePlugin/BaseInterface.hh>
#include <OpenFlipper/BasePlugin/ToolboxInterface.hh>
#include <OpenFlipper/BasePlugin/LoggingInterface.hh>
#include <OpenFlipper/BasePlugin/LoadSaveInterface.hh>
#include <OpenFlipper/BasePlugin/TextureInterface.hh>
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include "MastersThesisToolbar.hh"

class MastersThesisPlugin : public QObject, BaseInterface, ToolboxInterface, LoggingInterface, LoadSaveInterface {
Q_OBJECT
    Q_INTERFACES(BaseInterface)
    Q_INTERFACES(ToolboxInterface)
    Q_INTERFACES(LoggingInterface)
    Q_INTERFACES(LoadSaveInterface)
//    Q_INTERFACES(TextureInterface)


#if QT_VERSION >= 0x050000
    Q_PLUGIN_METADATA(IID "org.OpenFlipper.Plugins.Plugin-MastersThesis")
#endif

signals:

    void updateView();

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
//    void addTexture(QString _textureName, QString _filename, uint dimension);
//
//    // tell OpenFlipper that our texture (coordinates) have changed
//    void updatedTextures(QString, int);
//
//    void updatedTexture(QString);
//
//    // tell OpenFlipper which texture settings we want to use
//    void setTextureMode(QString _textureName, QString _mode);
//
//    // tell OpenFlipper to use the texture with name _textureName
//    void switchTexture(QString _textureName);

    /** \brief emit this signal if you want to switch the texture of a specific object
     * This signal can be called from any thread.\n
     */
//    void switchTexture(QString _textureName, int _id);


private slots:

    // initialization functions
    void initializePlugin();

    void pluginsInitialized();


public :

    ~MastersThesisPlugin() {}

    QString name() { return QString("Simple plugin"); };

    QString description() { return QString("Does actually nothing but works!"); };

public slots:

    void slot_get_boundary();

    void slot_get_dualGraph();

    void slot_get_global_param();

private:
    MastersThesisToolbar *tool_;

    //store selected vertices
    std::vector<int> heConstraints;
    std::vector<int> includedHEdges;


};

#endif //OPENFLIPPER_MASTERSTHESISPLUGIN_H
