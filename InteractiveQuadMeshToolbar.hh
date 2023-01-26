#pragma once

#if QT_VERSION >= 0x050000
#include <QtWidgets>
#else
#include <QtGui>
#endif

#include "ui_InteractiveQuadMeshToolbarBase.h"

class InteractiveQuadMeshToolbar : public QWidget, public Ui::InteractiveQuadMeshToolbarBase
{
    Q_OBJECT

public:
    InteractiveQuadMeshToolbar(QWidget * parent = 0);
};
