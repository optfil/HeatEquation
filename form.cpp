#include "form.h"

Form::Form(QWidget *parent)
    : QWidget(parent)
{
    labelInitial = new QLabel(tr("Temperature profile"));
    comboBoxInitial = new QComboBox();
    comboBoxInitial->addItem(tr("Gauss"), QVariant(Gauss));
    comboBoxInitial->addItem(tr("SuperGauss"), QVariant(SuperGauss));
    comboBoxInitial->addItem(tr("Rectangle"), QVariant(Rectangle));
    comboBoxInitial->addItem(tr("Delta"), QVariant(Delta));

    labelSizeX_1 = new QLabel(tr("Grid size"));
    labelSizeX_2 = new QLabel(tr(" L = "));
    labelSizeX_2->setAlignment(Qt::AlignRight);
    labelSizeT_1 = new QLabel(tr("Integration time"));
    labelSizeT_2 = new QLabel(tr(" T = "));
    labelSizeT_2->setAlignment(Qt::AlignRight);
    labelNX_1 = new QLabel(tr("Number of spatial points"));
    labelNX_2 = new QLabel(tr("NX = "));
    labelNX_2->setAlignment(Qt::AlignRight);
    labelNT_1 = new QLabel(tr("Number of temporal points"));
    labelNT_2 = new QLabel(tr("NT = "));
    labelNT_2->setAlignment(Qt::AlignRight);

    labelSizeX = new QLabel(QString::number(kRangeX, 'f', 1));
    labelSizeT = new QLabel(QString::number(kRangeT, 'f', 1));

    sliderNX = new QSlider(Qt::Horizontal);
    sliderNX->setRange(1, 4);
    sliderNX->setSingleStep(1);
    sliderNX->setPageStep(1);
    sliderNX->setTickInterval(1);
    sliderNX->setTickPosition(QSlider::TicksBelow);
    sliderNX->setValue(1);

    sliderNT = new QSlider(Qt::Horizontal);
    sliderNT->setRange(kNtMin, kNtMax);
    sliderNT->setTickInterval(10);
    sliderNT->setTickPosition(QSlider::TicksBelow);
    sliderNT->setValue(kNtMin);

    spinBoxNX = new QSpinBox();
    spinBoxNX->setMinimum(0);
    spinBoxNX->setMaximum(kNxMax);
    spinBoxNX->setValue(kNxMin);
    spinBoxNX->setSingleStep(kNxMin);

    spinBoxNT = new QSpinBox();
    spinBoxNT->setMinimum(kNtMin);
    spinBoxNT->setMaximum(kNtMax);
    spinBoxNT->setValue(kNtMin);
    spinBoxNT->setSingleStep(1);

    labelStepX_1 = new QLabel(tr("Spatial step"));
    labelStepX_2 = new QLabel(tr("dx = "));
    labelStepX_2->setAlignment(Qt::AlignRight);
    labelStepX = new QLabel();

    labelStepT_1 = new QLabel(tr("Time step"));
    labelStepT_2 = new QLabel(tr("dt = "));
    labelStepT_2->setAlignment(Qt::AlignRight);
    labelStepT = new QLabel();

    labelCFL_1 = new QLabel(tr("Alpha"));
    labelCFL_2 = new QLabel(tr("α = "));
    labelCFL_2->setAlignment(Qt::AlignRight);
    labelCFL = new QLabel();

    pushButtonSolve = new QPushButton(tr("Start"));

    widgetExplicit = new QWidget();
    widgetImplicit = new QWidget();
    widgetCrankNicolson = new QWidget();

    tabWidgetMethods = new QTabWidget();
    tabWidgetMethods->addTab(widgetExplicit, tr("Explicit"));
    tabWidgetMethods->addTab(widgetImplicit, tr("Implicit"));
    tabWidgetMethods->addTab(widgetCrankNicolson, tr("Crank-Nicolson"));
}

Form::~Form()
{

}
