#include "form.h"

Form::Form(QWidget *parent)
    : QWidget(parent), param(nullptr)
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

    labelAlpha_1 = new QLabel(tr("Alpha"));
    labelAlpha_2 = new QLabel(tr("α = "));
    labelAlpha_2->setAlignment(Qt::AlignRight);
    labelAlpha = new QLabel();

    pushButtonSolve = new QPushButton(tr("Start"));

    widgetExplicit = new QWidget();
    widgetImplicit = new QWidget();
    widgetCrankNicolson = new QWidget();

    tabWidgetMethods = new QTabWidget();
    tabWidgetMethods->addTab(widgetExplicit, tr("Explicit"));
    tabWidgetMethods->addTab(widgetImplicit, tr("Implicit"));
    tabWidgetMethods->addTab(widgetCrankNicolson, tr("Crank-Nicolson"));

    QGridLayout *layoutNxNt = new QGridLayout();
    layoutNxNt->addWidget(labelInitial, 0, 0, 1, 1);
    layoutNxNt->addWidget(comboBoxInitial, 0, 1, 1, 3);
    layoutNxNt->addWidget(labelSizeX_1, 1, 0, 1, 1);
    layoutNxNt->addWidget(labelSizeX_2, 1, 1, 1, 1);
    layoutNxNt->addWidget(labelSizeX, 1, 2, 1, 1);
    layoutNxNt->addWidget(labelSizeT_1, 2, 0, 1, 1);
    layoutNxNt->addWidget(labelSizeT_2, 2, 1, 1, 1);
    layoutNxNt->addWidget(labelSizeT, 2, 2, 1, 1);
    layoutNxNt->addWidget(labelNX_1, 3, 0, 1, 1, Qt::AlignBaseline);
    layoutNxNt->addWidget(labelNX_2, 3, 1, 1, 1, Qt::AlignBaseline);
    layoutNxNt->addWidget(spinBoxNX, 3, 2, 1, 1);
    layoutNxNt->addWidget(sliderNX, 3, 3, 1, 1);
    layoutNxNt->addWidget(labelNT_1, 4, 0, 1, 1, Qt::AlignBaseline);
    layoutNxNt->addWidget(labelNT_2, 4, 1, 1, 1, Qt::AlignBaseline);
    layoutNxNt->addWidget(spinBoxNT, 4, 2, 1, 1);
    layoutNxNt->addWidget(sliderNT, 4, 3, 1, 1);
    layoutNxNt->addWidget(labelStepX_1, 5, 0, 1, 1);
    layoutNxNt->addWidget(labelStepX_2, 5, 1, 1, 1);
    layoutNxNt->addWidget(labelStepX, 5, 2, 1, 1);
    layoutNxNt->addWidget(labelStepT_1, 6, 0, 1, 1);
    layoutNxNt->addWidget(labelStepT_2, 6, 1, 1, 1);
    layoutNxNt->addWidget(labelStepT, 6, 2, 1, 1);
    layoutNxNt->addWidget(labelAlpha_1, 7, 0, 1, 1);
    layoutNxNt->addWidget(labelAlpha_2, 7, 1, 1, 1);
    layoutNxNt->addWidget(labelAlpha, 7, 2, 1, 1);
    layoutNxNt->addWidget(pushButtonSolve, 6, 3, 2, 1);

    QVBoxLayout *layoutParam = new QVBoxLayout;
    //layoutParam->addWidget(chartView);
    layoutParam->addLayout(layoutNxNt);

    QHBoxLayout *layoutMain = new QHBoxLayout();
    layoutMain->addLayout(layoutParam);
    layoutMain->addWidget(tabWidgetMethods);

    setLayout(layoutMain);

    connect(comboBoxInitial, SIGNAL(currentIndexChanged(int)), this, SLOT(selectionChanged()));
    connect(sliderNX, SIGNAL(valueChanged(int)), this, SLOT(update_nx_from_slider(int)));
    connect(sliderNT, SIGNAL(valueChanged(int)), this, SLOT(update_nt(int)));
    connect(spinBoxNX, SIGNAL(valueChanged(int)), this, SLOT(update_nx(int)));
    connect(spinBoxNT, SIGNAL(valueChanged(int)), this, SLOT(update_nt(int)));
}

Form::~Form()
{
    delete param;
}

void Form::update_nx_from_slider(int log_n)
{
    int new_nx = static_cast<int>(std::round(std::pow(2.0, n+4)));

    spinBoxNX->blockSignals(true);
    sliderNX->blockSignals(true);

    spinBoxNX->setValue(new_nx);
    spinBoxNX->setSingleStep(new_nx);

    spinBoxNX->blockSignals(false);
    sliderNX->blockSignals(false);

    initiateState();
    updateSpectrum();
}

void Form::update_nx(int n)
{
    int old_nx = param->get_nx();
    if (n == 0)
    {
        n = std::max(old_nx/2, kNxMin);
    }

    int new_nx_log = static_cast<int>(std::round(std::log2(static_cast<double>(n))));
    int new_nx = static_cast<int>(std::round(std::pow(2.0, new_nx_log)));

    spinBoxNX->blockSignals(true);
    sliderNX->blockSignals(true);

    sliderNX->setValue(new_nx_log-3);
    spinBoxNX->setValue(new_nx);
    spinBoxNX->setSingleStep(new_nx);

    spinBoxNX->blockSignals(false);
    sliderNX->blockSignals(false);

    initiateState();
    updateSpectrum();
}

void Form::update_nt(int n)
{
    spinBoxNT->blockSignals(true);
    sliderNT->blockSignals(true);

    sliderNT->setValue(n);
    spinBoxNT->setValue(n);

    spinBoxNT->blockSignals(false);
    sliderNT->blockSignals(false);

    initiateState();
}

void Form::selectionChanged()
{
    initiateState();
    updateSpectrum();
}

void Form::updateLabels()
{
    labelStepX->setText(QString::number(param->get_dx(), 'f', 3));
    labelStepT->setText(QString::number(param->get_dt(), 'f', 3));
    labelAlpha->setText(QString::number(param->get_alpha(), 'f', 3));
}

void Form::initiateState()
{}

void Form::updateSpectrum()
{}
