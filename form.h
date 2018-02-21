#ifndef FORM_H
#define FORM_H

#include <QComboBox>
#include <QPushButton>
#include <QSlider>
#include <QSpinBox>
#include <QTabWidget>
#include <QTimer>
#include <QWidget>

#include <QtCharts/QtCharts>
QT_CHARTS_USE_NAMESPACE

constexpr double kRangeX = 10.0;
constexpr double kRangeT = 1.0;
constexpr int kNxMin = 32;
constexpr int kNxMax = 256;
constexpr int kNtMin = 1;
constexpr int kNtMax = 1000;

class Form : public QWidget
{
    Q_OBJECT

public:
    Form(QWidget *parent = 0);
    ~Form();

    enum InitialProfile {Gauss, SuperGauss, Rectangle, Delta};
    Q_ENUM(InitialProfile)

    enum MethodType {Explicit, Implicit, CrankNicolson};
    Q_ENUM(MethodType)

private slots:

private:
    QLabel *labelInitial;
    QComboBox *comboBoxInitial;
    QLabel *labelSizeX_1, *labelSizeX_2, *labelSizeT_1, *labelSizeT_2, *labelNX_1, *labelNX_2, *labelNT_1, *labelNT_2;
    QLabel *labelSizeX, *labelSizeT;
    QSlider *sliderNX, *sliderNT;
    QSpinBox *spinBoxNX, *spinBoxNT;
    QLabel *labelStepX_1, *labelStepX_2, *labelStepX;
    QLabel *labelStepT_1, *labelStepT_2, *labelStepT;
    QLabel *labelCFL_1, *labelCFL_2, *labelCFL;
    QPushButton *pushButtonSolve;
    QTabWidget *tabWidgetMethods;
    QWidget *widgetExplicit, *widgetImplicit, *widgetCrankNicolson;
};

#endif // FORM_H
