#ifndef OPENFLIPPER_MASTERSTHESISPLUGIN_H
#define OPENFLIPPER_MASTERSTHESISPLUGIN_H

#include <QObject>
#include <OpenFlipper/common/Types.hh>
#include <OpenFlipper/BasePlugin/BaseInterface.hh>
#include <OpenFlipper/BasePlugin/ToolboxInterface.hh>
#include <OpenFlipper/BasePlugin/LoggingInterface.hh>
#include <OpenFlipper/BasePlugin/LoadSaveInterface.hh>
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include "MastersThesisToolbar.hh"

class MastersThesisPlugin : public QObject, BaseInterface, ToolboxInterface, LoggingInterface, LoadSaveInterface {
Q_OBJECT
    Q_INTERFACES(BaseInterface)
    Q_INTERFACES(ToolboxInterface)
    Q_INTERFACES(LoggingInterface)
    Q_INTERFACES(LoadSaveInterface)


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

private:
    MastersThesisToolbar *tool_;

    //store selected vertices
    std::vector<int> heConstraints;
    std::vector<int> includedHEdges;


};

#endif //OPENFLIPPER_MASTERSTHESISPLUGIN_H
