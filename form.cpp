#include "form.h"

#include <complex>

static double initial(double x, Form::InitialProfile profile, double ampl = 1.0)
{
    switch (profile)
    {
    case Form::Gauss:
        return ampl * std::exp(-std::pow(x / (0.1*kRangeX), 2.0));
    case Form::SuperGauss:
        return ampl * std::exp(-std::pow(x / (0.1*kRangeX), 8.0));
    case Form::Rectangle:
        return ampl * ((std::abs(x) < 0.1*kRangeX) ? 1.0 : 0.0);
    case Form::Delta:
        return ampl * ((std::abs(x) < 1e-10*kRangeX) ? 1.0 : 0.0);
    default:
        return 0;
    }
}

static std::pair<double, double> dispersion_diffusion(double q_N, double alpha, Form::MethodType type)
{
    std::complex<double> lambda;
    double kappa = 2.0*M_PI*q_N;
    switch (type)
    {
    case Form::Explicit:
        lambda = 1.0 - 2.0 * alpha * (1.0 - std::cos(kappa));
        break;
    case Form::Implicit:
        lambda = 1.0 / (1.0 + 2.0 * alpha * (1.0 - std::cos(kappa)));
        break;
    case Form::CrankNicolson:
        lambda = (1.0 - alpha * (1.0 - std::cos(kappa))) / (1.0 + alpha * (1.0 - std::cos(kappa)));
        break;
    default:
        lambda = 1.0;
        break;
    }

    lambda = std::log(lambda);

    return std::make_pair(std::imag(lambda), -std::real(lambda));
}

static void setGrid(QValueAxis* ax)
{
    ax->setGridLineVisible(true);
    QPen pen = ax->gridLinePen();
    pen.setWidth(2);
    pen.setColor(Qt::gray);
    ax->setGridLinePen(pen);
}

Form::Form(QWidget *parent)
    : QWidget(parent), param_(nullptr), t_cur_(0.0)
{
    timer = new QTimer();
    timer->setInterval(10);

    seriesInitial = new QLineSeries();
    seriesInitial->setColor(Qt::blue);
    seriesInitial->setPen(QPen(seriesInitial->pen().brush(), 3));
    seriesExplicitIdealDispersion = new QLineSeries();
    seriesExplicitIdealDispersion->setColor(Qt::blue);
    seriesExplicitIdealDispersion->setPen(QPen(seriesExplicitIdealDispersion->pen().brush(), 3));
    seriesImplicitIdealDispersion = new QLineSeries();
    seriesImplicitIdealDispersion->setColor(Qt::blue);
    seriesImplicitIdealDispersion->setPen(QPen(seriesImplicitIdealDispersion->pen().brush(), 3));
    seriesCrankNicolsonIdealDispersion = new QLineSeries();
    seriesCrankNicolsonIdealDispersion->setColor(Qt::blue);
    seriesCrankNicolsonIdealDispersion->setPen(QPen(seriesCrankNicolsonIdealDispersion->pen().brush(), 3));
    seriesExplicitIdealDissipation = new QLineSeries();
    seriesExplicitIdealDissipation->setColor(Qt::blue);
    seriesExplicitIdealDissipation->setPen(QPen(seriesExplicitIdealDissipation->pen().brush(), 3));
    seriesImplicitIdealDissipation = new QLineSeries();
    seriesImplicitIdealDissipation->setColor(Qt::blue);
    seriesImplicitIdealDissipation->setPen(QPen(seriesImplicitIdealDissipation->pen().brush(), 3));
    seriesCrankNicolsonIdealDissipation = new QLineSeries();
    seriesCrankNicolsonIdealDissipation->setColor(Qt::blue);
    seriesCrankNicolsonIdealDissipation->setPen(QPen(seriesCrankNicolsonIdealDissipation->pen().brush(), 3));
    seriesExplicitDispersion = new QLineSeries();
    seriesExplicitDispersion->setColor(Qt::red);
    seriesExplicitDispersion->setPen(QPen(seriesExplicitDispersion->pen().brush(), 3));
    seriesImplicitDispersion = new QLineSeries();
    seriesImplicitDispersion->setColor(Qt::red);
    seriesImplicitDispersion->setPen(QPen(seriesImplicitDispersion->pen().brush(), 3));
    seriesCrankNicolsonDispersion = new QLineSeries();
    seriesCrankNicolsonDispersion->setColor(Qt::red);
    seriesCrankNicolsonDispersion->setPen(QPen(seriesCrankNicolsonDispersion->pen().brush(), 3));
    seriesExplicitDissipation = new QLineSeries();
    seriesExplicitDissipation->setColor(Qt::red);
    seriesExplicitDissipation->setPen(QPen(seriesExplicitDissipation->pen().brush(), 3));
    seriesImplicitDissipation = new QLineSeries();
    seriesImplicitDissipation->setColor(Qt::red);
    seriesImplicitDissipation->setPen(QPen(seriesImplicitDissipation->pen().brush(), 3));
    seriesCrankNicolsonDissipation = new QLineSeries();
    seriesCrankNicolsonDissipation->setColor(Qt::red);
    seriesCrankNicolsonDissipation->setPen(QPen(seriesCrankNicolsonDissipation->pen().brush(), 3));

    seriesExplicitIdealDispersion->append(QList<QPointF>() << QPointF(0.0, 0.0) << QPointF(0.5, 0.0));
    seriesImplicitIdealDispersion->append(QList<QPointF>() << QPointF(0.0, 0.0) << QPointF(0.5, 0.0));
    seriesCrankNicolsonIdealDispersion->append(QList<QPointF>() << QPointF(0.0, 0.0) << QPointF(0.5, 0.0));

    QChart *chartInitial = new QChart();
    chartInitial->addSeries(seriesInitial);
    chartInitial->setTitle(tr("Initial profile"));
    chartInitial->legend()->hide();

    QValueAxis *axisXInitial = new QValueAxis;
    axisXInitial->setLineVisible(false);
    setGrid(axisXInitial);
    axisXInitial->setLabelsVisible(false);
    axisXInitial->setRange(-kRangeX/2, kRangeX/2);
    chartInitial->addAxis(axisXInitial, Qt::AlignBottom);
    seriesInitial->attachAxis(axisXInitial);
    QValueAxis *axisYInitial = new QValueAxis;
    axisYInitial->setLineVisible(false);
    setGrid(axisYInitial);
    axisYInitial->setLabelsVisible(false);
    axisYInitial->setRange(0.0, 1.003);
    chartInitial->addAxis(axisYInitial, Qt::AlignLeft);
    seriesInitial->attachAxis(axisYInitial);

    chartView = new QChartView(chartInitial);
    chartView->setRenderHint(QPainter::Antialiasing);

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

    QChart *explicitDispersionChart = new QChart();
    //explicitDispersionChart->addSeries(spectrumExplicitDispersion);
    explicitDispersionChart->addSeries(seriesExplicitIdealDispersion);
    explicitDispersionChart->addSeries(seriesExplicitDispersion);
    explicitDispersionChart->setTitle(tr("Dispersion"));
    explicitDispersionChart->legend()->hide();

    QValueAxis *axisXExplicitDispersion = new QValueAxis;
    axisXExplicitDispersion->setLineVisible(false);
    setGrid(axisXExplicitDispersion);
    axisXExplicitDispersion->setTitleText("ϰ / ϰ_N");
    axisXExplicitDispersion->setTitleFont(QFont("Times New Roman", 14));
    axisXExplicitDispersion->setTickCount(3);
    axisXExplicitDispersion->setRange(0.0, 0.5);
    explicitDispersionChart->addAxis(axisXExplicitDispersion, Qt::AlignBottom);
    seriesExplicitIdealDispersion->attachAxis(axisXExplicitDispersion);
    seriesExplicitDispersion->attachAxis(axisXExplicitDispersion);
    //QValueAxis *axisSpectrumXExplicitDispersion = new QValueAxis;
    //axisSpectrumXExplicitDispersion->setLineVisible(false);
    //axisSpectrumXExplicitDispersion->setLabelsVisible(false);
    //explicitDispersionChart->addAxis(axisSpectrumXExplicitDispersion, Qt::AlignBottom);
    //spectrumExplicitDispersion->attachAxis(axisSpectrumXExplicitDispersion);
    QValueAxis *axisYExplicitDispersion = new QValueAxis;
    axisYExplicitDispersion->setLineVisible(false);
    setGrid(axisYExplicitDispersion);
    axisYExplicitDispersion->setTitleText("Ω⋅Δt / π");
    axisYExplicitDispersion->setTitleFont(QFont("Times New Roman", 14));
    axisYExplicitDispersion->setTickCount(5);
    axisYExplicitDispersion->setRange(-0.5, 1.5);
    explicitDispersionChart->addAxis(axisYExplicitDispersion, Qt::AlignLeft);
    seriesExplicitIdealDispersion->attachAxis(axisYExplicitDispersion);
    seriesExplicitDispersion->attachAxis(axisYExplicitDispersion);
    //spectrumExplicitDispersion->attachAxis(axisYExplicitDispersion);

    explicitDispersion = new QChartView();
    explicitDispersion->setRenderHint(QPainter::Antialiasing);
    explicitDispersion->setChart(explicitDispersionChart);

    QChart *explicitDissipationChart = new QChart();
    //explicitDissipationChart->addSeries(spectrumExplicitDissipation);
    explicitDissipationChart->addSeries(seriesExplicitIdealDissipation);
    explicitDissipationChart->addSeries(seriesExplicitDissipation);
    explicitDissipationChart->setTitle(tr("Dissipation"));
    explicitDissipationChart->legend()->hide();

    QValueAxis *axisXExplicitDissipation = new QValueAxis;
    axisXExplicitDissipation->setLineVisible(false);
    setGrid(axisXExplicitDissipation);
    axisXExplicitDissipation->setTitleText("ϰ / ϰ_N");
    axisXExplicitDissipation->setTitleFont(QFont("Times New Roman", 14));
    axisXExplicitDissipation->setTickCount(3);
    axisXExplicitDissipation->setRange(0.0, 0.5);
    explicitDissipationChart->addAxis(axisXExplicitDissipation, Qt::AlignBottom);
    seriesExplicitIdealDissipation->attachAxis(axisXExplicitDissipation);
    seriesExplicitDissipation->attachAxis(axisXExplicitDissipation);
    //QValueAxis *axisSpectrumXExplicitDissipation = new QValueAxis;
    //axisSpectrumXExplicitDissipation->setLineVisible(false);
    //axisSpectrumXExplicitDissipation->setLabelsVisible(false);
    //explicitDissipationChart->addAxis(axisSpectrumXExplicitDissipation, Qt::AlignBottom);
    //spectrumExplicitDissipation->attachAxis(axisSpectrumXExplicitDissipation);
    QValueAxis *axisYExplicitDissipation = new QValueAxis;
    axisYExplicitDissipation->setLineVisible(false);
    setGrid(axisYExplicitDissipation);
    axisYExplicitDissipation->setTitleText("γ⋅Δt / (4π²α)");
    axisYExplicitDissipation->setTitleFont(QFont("Times New Roman", 14));
    axisYExplicitDissipation->setTickCount(5);
    axisYExplicitDissipation->setRange(-0.1, 0.3);
    explicitDissipationChart->addAxis(axisYExplicitDissipation, Qt::AlignLeft);
    seriesExplicitIdealDissipation->attachAxis(axisYExplicitDissipation);
    seriesExplicitDissipation->attachAxis(axisYExplicitDissipation);
    //spectrumExplicitDissipation->attachAxis(axisYExplicitDissipation);

    explicitDissipation = new QChartView();
    explicitDissipation->setRenderHint(QPainter::Antialiasing);
    explicitDissipation->setChart(explicitDissipationChart);

    QChart *explicitSolutionChart = new QChart();
    explicitSolutionChart->setTitle(tr("Solution"));
    explicitSolutionChart->legend()->hide();
    QValueAxis *axisXExplicitSolution = new QValueAxis;
    axisXExplicitSolution->setLineVisible(false);
    setGrid(axisXExplicitSolution);
    axisXExplicitSolution->setLabelsVisible(false);
    axisXExplicitSolution->setRange(-0.5*kRangeX, 0.5*kRangeX);
    explicitSolutionChart->addAxis(axisXExplicitSolution, Qt::AlignBottom);
    QValueAxis *axisYExplicitSolution = new QValueAxis;
    axisYExplicitSolution->setLineVisible(false);
    setGrid(axisYExplicitSolution);
    axisYExplicitSolution->setLabelsVisible(false);
    axisYExplicitSolution->setRange(-0.2, 1.2);
    explicitSolutionChart->addAxis(axisYExplicitSolution, Qt::AlignLeft);

    explicitSolution = new QChartView();
    explicitSolution->setRenderHint(QPainter::Antialiasing);
    explicitSolution->setChart(explicitSolutionChart);

    QVBoxLayout *explicitLeft = new QVBoxLayout();
    explicitLeft->addWidget(explicitDispersion);
    explicitLeft->addWidget(explicitDissipation);
    QHBoxLayout *explicitMain = new QHBoxLayout();
    explicitMain->addLayout(explicitLeft);
    explicitMain->addWidget(explicitSolution);
    widgetExplicit->setLayout(explicitMain);

    widgetImplicit = new QWidget();

    QChart *implicitDispersionChart = new QChart();
    //implicitDispersionChart->addSeries(spectrumImplicitDispersion);
    implicitDispersionChart->addSeries(seriesImplicitIdealDispersion);
    implicitDispersionChart->addSeries(seriesImplicitDispersion);
    implicitDispersionChart->setTitle(tr("Dispersion"));
    implicitDispersionChart->legend()->hide();

    QValueAxis *axisXImplicitDispersion = new QValueAxis;
    axisXImplicitDispersion->setLineVisible(false);
    setGrid(axisXImplicitDispersion);
    axisXImplicitDispersion->setTitleText("ϰ / ϰ_N");
    axisXImplicitDispersion->setTitleFont(QFont("Times New Roman", 14));
    axisXImplicitDispersion->setTickCount(3);
    axisXImplicitDispersion->setRange(0.0, 0.5);
    implicitDispersionChart->addAxis(axisXImplicitDispersion, Qt::AlignBottom);
    seriesImplicitIdealDispersion->attachAxis(axisXImplicitDispersion);
    seriesImplicitDispersion->attachAxis(axisXImplicitDispersion);
    //QValueAxis *axisSpectrumXImplicitDispersion = new QValueAxis;
    //axisSpectrumXImplicitDispersion->setLineVisible(false);
    //axisSpectrumXImplicitDispersion->setLabelsVisible(false);
    //implicitDispersionChart->addAxis(axisSpectrumXImplicitDispersion, Qt::AlignBottom);
    //spectrumImplicitDispersion->attachAxis(axisSpectrumXImplicitDispersion);
    QValueAxis *axisYImplicitDispersion = new QValueAxis;
    axisYImplicitDispersion->setLineVisible(false);
    setGrid(axisYImplicitDispersion);
    axisYImplicitDispersion->setTitleText("Ω⋅Δt / π");
    axisYImplicitDispersion->setTitleFont(QFont("Times New Roman", 14));
    axisYImplicitDispersion->setTickCount(5);
    axisYImplicitDispersion->setRange(-0.5, 1.5);
    implicitDispersionChart->addAxis(axisYImplicitDispersion, Qt::AlignLeft);
    seriesImplicitIdealDispersion->attachAxis(axisYImplicitDispersion);
    seriesImplicitDispersion->attachAxis(axisYImplicitDispersion);
    //spectrumImplicitDispersion->attachAxis(axisYImplicitDispersion);

    implicitDispersion = new QChartView();
    implicitDispersion->setRenderHint(QPainter::Antialiasing);
    implicitDispersion->setChart(implicitDispersionChart);

    QChart *implicitDissipationChart = new QChart();
    //implicitDissipationChart->addSeries(spectrumImplicitDissipation);
    implicitDissipationChart->addSeries(seriesImplicitIdealDissipation);
    implicitDissipationChart->addSeries(seriesImplicitDissipation);
    implicitDissipationChart->setTitle(tr("Dissipation"));
    implicitDissipationChart->legend()->hide();

    QValueAxis *axisXImplicitDissipation = new QValueAxis;
    axisXImplicitDissipation->setLineVisible(false);
    setGrid(axisXImplicitDissipation);
    axisXImplicitDissipation->setTitleText("ϰ / ϰ_N");
    axisXImplicitDissipation->setTitleFont(QFont("Times New Roman", 14));
    axisXImplicitDissipation->setTickCount(3);
    axisXImplicitDissipation->setRange(0.0, 0.5);
    implicitDissipationChart->addAxis(axisXImplicitDissipation, Qt::AlignBottom);
    seriesImplicitIdealDissipation->attachAxis(axisXImplicitDissipation);
    seriesImplicitDissipation->attachAxis(axisXImplicitDissipation);
    //QValueAxis *axisSpectrumXImplicitDissipation = new QValueAxis;
    //axisSpectrumXImplicitDissipation->setLineVisible(false);
    //axisSpectrumXImplicitDissipation->setLabelsVisible(false);
    //implicitDissipationChart->addAxis(axisSpectrumXImplicitDissipation, Qt::AlignBottom);
    //spectrumImplicitDissipation->attachAxis(axisSpectrumXImplicitDissipation);
    QValueAxis *axisYImplicitDissipation = new QValueAxis;
    axisYImplicitDissipation->setLineVisible(false);
    setGrid(axisYImplicitDissipation);
    axisYImplicitDissipation->setTitleText("γ⋅Δt");
    axisYImplicitDissipation->setTitleFont(QFont("Times New Roman", 14));
    axisYImplicitDissipation->setTickCount(5);
    axisYImplicitDissipation->setRange(-0.1, 0.3);
    implicitDissipationChart->addAxis(axisYImplicitDissipation, Qt::AlignLeft);
    seriesImplicitIdealDissipation->attachAxis(axisYImplicitDissipation);
    seriesImplicitDissipation->attachAxis(axisYImplicitDissipation);
    //spectrumImplicitDissipation->attachAxis(axisYImplicitDissipation);

    implicitDissipation = new QChartView();
    implicitDissipation->setRenderHint(QPainter::Antialiasing);
    implicitDissipation->setChart(implicitDissipationChart);

    QChart *implicitSolutionChart = new QChart();
    implicitSolutionChart->setTitle(tr("Solution"));
    implicitSolutionChart->legend()->hide();
    QValueAxis *axisXImplicitSolution = new QValueAxis;
    axisXImplicitSolution->setLineVisible(false);
    setGrid(axisXImplicitSolution);
    axisXImplicitSolution->setLabelsVisible(false);
    axisXImplicitSolution->setRange(-0.5*kRangeX, 0.5*kRangeX);
    implicitSolutionChart->addAxis(axisXImplicitSolution, Qt::AlignBottom);
    QValueAxis *axisYImplicitSolution = new QValueAxis;
    axisYImplicitSolution->setLineVisible(false);
    setGrid(axisYImplicitSolution);
    axisYImplicitSolution->setLabelsVisible(false);
    axisYImplicitSolution->setRange(-0.2, 1.2);
    implicitSolutionChart->addAxis(axisYImplicitSolution, Qt::AlignLeft);

    implicitSolution = new QChartView();
    implicitSolution->setRenderHint(QPainter::Antialiasing);
    implicitSolution->setChart(implicitSolutionChart);

    QVBoxLayout *implicitLeft = new QVBoxLayout();
    implicitLeft->addWidget(implicitDispersion);
    implicitLeft->addWidget(implicitDissipation);
    QHBoxLayout *implicitMain = new QHBoxLayout();
    implicitMain->addLayout(implicitLeft);
    implicitMain->addWidget(implicitSolution);
    widgetImplicit->setLayout(implicitMain);

    widgetCrankNicolson = new QWidget();

    QChart *crankNicolsonDispersionChart = new QChart();
    //crankNicolsonDispersionChart->addSeries(spectrumCrankNicolsonDispersion);
    crankNicolsonDispersionChart->addSeries(seriesCrankNicolsonIdealDispersion);
    crankNicolsonDispersionChart->addSeries(seriesCrankNicolsonDispersion);
    crankNicolsonDispersionChart->setTitle(tr("Dispersion"));
    crankNicolsonDispersionChart->legend()->hide();

    QValueAxis *axisXCrankNicolsonDispersion = new QValueAxis;
    axisXCrankNicolsonDispersion->setLineVisible(false);
    setGrid(axisXCrankNicolsonDispersion);
    axisXCrankNicolsonDispersion->setTitleText("ϰ / ϰ_N");
    axisXCrankNicolsonDispersion->setTitleFont(QFont("Times New Roman", 14));
    axisXCrankNicolsonDispersion->setTickCount(3);
    axisXCrankNicolsonDispersion->setRange(0.0, 0.5);
    crankNicolsonDispersionChart->addAxis(axisXCrankNicolsonDispersion, Qt::AlignBottom);
    seriesCrankNicolsonIdealDispersion->attachAxis(axisXCrankNicolsonDispersion);
    seriesCrankNicolsonDispersion->attachAxis(axisXCrankNicolsonDispersion);
    //QValueAxis *axisSpectrumXCrankNicolsonDispersion = new QValueAxis;
    //axisSpectrumXCrankNicolsonDispersion->setLineVisible(false);
    //axisSpectrumXCrankNicolsonDispersion->setLabelsVisible(false);
    //crankNicolsonDispersionChart->addAxis(axisSpectrumXCrankNicolsonDispersion, Qt::AlignBottom);
    //spectrumCrankNicolsonDispersion->attachAxis(axisSpectrumXCrankNicolsonDispersion);
    QValueAxis *axisYCrankNicolsonDispersion = new QValueAxis;
    axisYCrankNicolsonDispersion->setLineVisible(false);
    setGrid(axisYCrankNicolsonDispersion);
    axisYCrankNicolsonDispersion->setTitleText("Ω⋅Δt / π");
    axisYCrankNicolsonDispersion->setTitleFont(QFont("Times New Roman", 14));
    axisYCrankNicolsonDispersion->setTickCount(5);
    axisYCrankNicolsonDispersion->setRange(-0.5, 1.5);
    crankNicolsonDispersionChart->addAxis(axisYCrankNicolsonDispersion, Qt::AlignLeft);
    seriesCrankNicolsonIdealDispersion->attachAxis(axisYCrankNicolsonDispersion);
    seriesCrankNicolsonDispersion->attachAxis(axisYCrankNicolsonDispersion);
    //spectrumCrankNicolsonDispersion->attachAxis(axisYCrankNicolsonDispersion);

    crankNicolsonDispersion = new QChartView();
    crankNicolsonDispersion->setRenderHint(QPainter::Antialiasing);
    crankNicolsonDispersion->setChart(crankNicolsonDispersionChart);

    QChart *crankNicolsonDissipationChart = new QChart();
    //crankNicolsonDissipationChart->addSeries(spectrumCrankNicolsonDissipation);
    crankNicolsonDissipationChart->addSeries(seriesCrankNicolsonIdealDissipation);
    crankNicolsonDissipationChart->addSeries(seriesCrankNicolsonDissipation);
    crankNicolsonDissipationChart->setTitle(tr("Dissipation"));
    crankNicolsonDissipationChart->legend()->hide();

    QValueAxis *axisXCrankNicolsonDissipation = new QValueAxis;
    axisXCrankNicolsonDissipation->setLineVisible(false);
    setGrid(axisXCrankNicolsonDissipation);
    axisXCrankNicolsonDissipation->setTitleText("ϰ / ϰ_N");
    axisXCrankNicolsonDissipation->setTitleFont(QFont("Times New Roman", 14));
    axisXCrankNicolsonDissipation->setTickCount(3);
    axisXCrankNicolsonDissipation->setRange(0.0, 0.5);
    crankNicolsonDissipationChart->addAxis(axisXCrankNicolsonDissipation, Qt::AlignBottom);
    seriesCrankNicolsonIdealDissipation->attachAxis(axisXCrankNicolsonDissipation);
    seriesCrankNicolsonDissipation->attachAxis(axisXCrankNicolsonDissipation);
    //QValueAxis *axisSpectrumXCrankNicolsonDissipation = new QValueAxis;
    //axisSpectrumXCrankNicolsonDissipation->setLineVisible(false);
    //axisSpectrumXCrankNicolsonDissipation->setLabelsVisible(false);
    //crankNicolsonDissipationChart->addAxis(axisSpectrumXCrankNicolsonDissipation, Qt::AlignBottom);
    //spectrumCrankNicolsonDissipation->attachAxis(axisSpectrumXCrankNicolsonDissipation);
    QValueAxis *axisYCrankNicolsonDissipation = new QValueAxis;
    axisYCrankNicolsonDissipation->setLineVisible(false);
    setGrid(axisYCrankNicolsonDissipation);
    axisYCrankNicolsonDissipation->setTitleText("γ⋅Δt");
    axisYCrankNicolsonDissipation->setTitleFont(QFont("Times New Roman", 14));
    axisYCrankNicolsonDissipation->setTickCount(5);
    axisYCrankNicolsonDissipation->setRange(-0.1, 0.3);
    crankNicolsonDissipationChart->addAxis(axisYCrankNicolsonDissipation, Qt::AlignLeft);
    seriesCrankNicolsonIdealDissipation->attachAxis(axisYCrankNicolsonDissipation);
    seriesCrankNicolsonDissipation->attachAxis(axisYCrankNicolsonDissipation);
    //spectrumCrankNicolsonDissipation->attachAxis(axisYCrankNicolsonDissipation);

    crankNicolsonDissipation = new QChartView();
    crankNicolsonDissipation->setRenderHint(QPainter::Antialiasing);
    crankNicolsonDissipation->setChart(crankNicolsonDissipationChart);

    QChart *crankNicolsonSolutionChart = new QChart();
    crankNicolsonSolutionChart->setTitle(tr("Solution"));
    crankNicolsonSolutionChart->legend()->hide();
    QValueAxis *axisXCrankNicolsonSolution = new QValueAxis;
    axisXCrankNicolsonSolution->setLineVisible(false);
    setGrid(axisXCrankNicolsonSolution);
    axisXCrankNicolsonSolution->setLabelsVisible(false);
    axisXCrankNicolsonSolution->setRange(-0.5*kRangeX, 0.5*kRangeX);
    crankNicolsonSolutionChart->addAxis(axisXCrankNicolsonSolution, Qt::AlignBottom);
    QValueAxis *axisYCrankNicolsonSolution = new QValueAxis;
    axisYCrankNicolsonSolution->setLineVisible(false);
    setGrid(axisYCrankNicolsonSolution);
    axisYCrankNicolsonSolution->setLabelsVisible(false);
    axisYCrankNicolsonSolution->setRange(-0.2, 1.2);
    crankNicolsonSolutionChart->addAxis(axisYCrankNicolsonSolution, Qt::AlignLeft);

    crankNicolsonSolution = new QChartView();
    crankNicolsonSolution->setRenderHint(QPainter::Antialiasing);
    crankNicolsonSolution->setChart(crankNicolsonSolutionChart);

    QVBoxLayout *crankNicolsonLeft = new QVBoxLayout();
    crankNicolsonLeft->addWidget(crankNicolsonDispersion);
    crankNicolsonLeft->addWidget(crankNicolsonDissipation);
    QHBoxLayout *crankNicolsonMain = new QHBoxLayout();
    crankNicolsonMain->addLayout(crankNicolsonLeft);
    crankNicolsonMain->addWidget(crankNicolsonSolution);
    widgetCrankNicolson->setLayout(crankNicolsonMain);

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
    layoutParam->addWidget(chartView);
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
    connect(tabWidgetMethods, SIGNAL(currentChanged(int)), this, SLOT(updateDispersionDiffusion()));
    connect(pushButtonSolve, SIGNAL(clicked(bool)), this, SLOT(Solve()));
    connect(timer, SIGNAL(timeout()), this, SLOT(Tick()));

    initiateState();
    updateSpectrum();
}

Form::~Form()
{
    delete param_;
}

void Form::update_nx_from_slider(int n)
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
    int old_nx = param_->get_nx();
    if (n == 0)
    {
        n = std::max(old_nx/2, kNxMin);
    }

    int new_nx_log = static_cast<int>(std::round(std::log2(static_cast<double>(n))));
    int new_nx = static_cast<int>(std::round(std::pow(2.0, new_nx_log)));

    spinBoxNX->blockSignals(true);
    sliderNX->blockSignals(true);

    sliderNX->setValue(new_nx_log-4);
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
    labelStepX->setText(QString::number(param_->get_dx(), 'f', 3));
    labelStepT->setText(QString::number(param_->get_dt(), 'f', 3));
    labelAlpha->setText(QString::number(param_->get_alpha(), 'f', 3));
}

void Form::initiateState()
{
    delete param_;
    param_ = new Parameters(spinBoxNX->value()+1, spinBoxNT->value(), kRangeX, kRangeT);

    InitialProfile profile = comboBoxInitial->currentData().value<InitialProfile>();
    state_.resize(param_->get_nx());
    tmp_state_.resize(state_.size());
    tdma_u_.resize(state_.size()-1);
    tdma_v_.resize(state_.size()-1);

    double ampl = (profile == Delta) ? kRangeX*0.1/param_->get_dx() : 1.0;

    for (decltype(state_.size()) i = 0; i < state_.size(); ++i)
        state_[i] = initial((double(i) - state_.size()/2) * param_->get_dx(), profile, ampl);

    QList<QPointF> init_data;
    for (decltype(state_.size()) i = 0; i < state_.size(); ++i)
        init_data.append(QPointF((double(i) - state_.size()/2) * param_->get_dx(), state_[i]));
    seriesInitial->clear();
    seriesInitial->append(init_data);

    updateLabels();
    updateDispersionDiffusion();
    cleanSolution();
}

void Form::updateDispersionDiffusion()
{
    method_ = static_cast<MethodType>(tabWidgetMethods->currentIndex());

    QList<QPointF> ideal_diff_data, explicit_disp_data, explicit_diff_data, implicit_disp_data, implicit_diff_data, crank_nicolson_disp_data, crank_nicolson_diff_data;
    std::pair<double, double> coeffs;
    for (int i = 0; i < param_->get_nx()/2+1; ++i)
    {
        double xi = static_cast<double>(i) / (param_->get_nx()-1);
        ideal_diff_data.append(QPointF(xi, xi*xi));

        coeffs = dispersion_diffusion(xi, param_->get_alpha(), Form::Explicit);
        explicit_disp_data.append(QPointF(xi, coeffs.first));
        explicit_diff_data.append(QPointF(xi, coeffs.second));

        coeffs = dispersion_diffusion(xi, param_->get_alpha(), Form::Implicit);
        implicit_disp_data.append(QPointF(xi, coeffs.first));
        implicit_diff_data.append(QPointF(xi, coeffs.second));

        coeffs = dispersion_diffusion(xi, param_->get_alpha(), Form::CrankNicolson);
        crank_nicolson_disp_data.append(QPointF(xi, coeffs.first));
        crank_nicolson_diff_data.append(QPointF(xi, coeffs.second));
    }

    for (int i = 0; i < param_->get_nx()/2+1; ++i)
    {
        explicit_disp_data[i].setY(explicit_disp_data[i].y() / M_PI);
        explicit_diff_data[i].setY(explicit_diff_data[i].y() / (4.0*M_PI*M_PI*param_->get_alpha()));
        implicit_disp_data[i].setY(implicit_disp_data[i].y() / M_PI);
        implicit_diff_data[i].setY(implicit_diff_data[i].y() / (4.0*M_PI*M_PI*param_->get_alpha()));
        crank_nicolson_disp_data[i].setY(crank_nicolson_disp_data[i].y() / M_PI);
        crank_nicolson_diff_data[i].setY(crank_nicolson_diff_data[i].y() / (4.0*M_PI*M_PI*param_->get_alpha()));
    }

    seriesExplicitIdealDissipation->clear();
    seriesExplicitIdealDissipation->append(ideal_diff_data);
    seriesExplicitDispersion->clear();
    seriesExplicitDispersion->append(explicit_disp_data);
    seriesExplicitDissipation->clear();
    seriesExplicitDissipation->append(explicit_diff_data);

    seriesImplicitIdealDissipation->clear();
    seriesImplicitIdealDissipation->append(ideal_diff_data);
    seriesImplicitDispersion->clear();
    seriesImplicitDispersion->append(implicit_disp_data);
    seriesImplicitDissipation->clear();
    seriesImplicitDissipation->append(implicit_diff_data);

    seriesCrankNicolsonIdealDissipation->clear();
    seriesCrankNicolsonIdealDissipation->append(ideal_diff_data);
    seriesCrankNicolsonDispersion->clear();
    seriesCrankNicolsonDispersion->append(crank_nicolson_disp_data);
    seriesCrankNicolsonDissipation->clear();
    seriesCrankNicolsonDissipation->append(crank_nicolson_diff_data);
}

void Form::updateSpectrum()
{}

void Form::cleanSolution()
{
    explicitSolution->chart()->removeAllSeries();
    implicitSolution->chart()->removeAllSeries();
    crankNicolsonSolution->chart()->removeAllSeries();
}

void Form::Solve()
{
    pushButtonSolve->setEnabled(false);
    tabWidgetMethods->setEnabled(false);
    comboBoxInitial->setEnabled(false);
    spinBoxNX->setEnabled(false);
    spinBoxNT->setEnabled(false);
    sliderNX->setEnabled(false);
    sliderNT->setEnabled(false);

    initiateState();
    updateSpectrum();

    showState();

    t_cur_ = 0.0;
    timer->start();
}

void Form::Tick()
{
    static int t_index = 1;

    if (t_cur_ < kRangeT + 1e-3*param_->get_dt())
    {
        t_cur_ += param_->get_dt();

        switch (method_)
        {
        case Explicit:
            tmp_state_.front() = state_.front();
            tmp_state_.back() = state_.back();
            for (decltype(state_.size()) i = 1; i < state_.size()-1; ++i)
                tmp_state_[i] = state_[i] + param_->get_alpha() * (state_[i+1] - 2.0*state_[i] + state_[i-1]);
            break;
        case Implicit:
            tdma_u_[0] = 0.0;
            tdma_v_[0] = state_[0];
            double inv_denominator;
            for (decltype(state_.size()) i = 1; i < state_.size()-1; ++i)
            {
                inv_denominator = 1.0 / (param_->get_alpha() * tdma_u_[i-1] - (2.0*param_->get_alpha()+1));
                tdma_u_[i] = -param_->get_alpha() * inv_denominator;
                tdma_v_[i] = (-state_[i] - param_->get_alpha() * tdma_v_[i-1]) * inv_denominator;
            }
            tmp_state_[state_.size()-1] = state_[state_.size()-1];
            for (decltype(state_.size()) i = state_.size()-2; i > 0; --i)
                tmp_state_[i] = tdma_u_[i] * tmp_state_[i+1] + tdma_v_[i];
            tmp_state_[0] = tdma_u_[0] * tmp_state_[1] + tdma_v_[0];

            break;
        case CrankNicolson:/*
            for (decltype(state_.size()) i = 1; i < state_.size()-1; ++i)
                tmp_state_[i] = (1.0 - param->get_alpha()*param->get_alpha()) * state_[i] - 0.5*param->get_alpha() * (state_[i+1] - state_[i-1]) + 0.5*param->get_alpha()*param->get_alpha() * (state_[i+1] + state_[i-1]);*/
            break;
        }

        state_ = tmp_state_;

        if (t_cur_ > kRangeT / 5.0 * t_index)
        {
            ++t_index;
            showState();
        }

        if (*std::max_element(&state_[0], &state_[static_cast<int>(state_.size()*0.4)]) > 3.0 || *std::min_element(&state_[0], &state_[static_cast<int>(state_.size()*0.4)]) < -3.0)
        {
            t_index = 1;
            showState();
            finishCalculation();
        }
    }
    else
    {
        t_index = 1;

        finishCalculation();
    }
}

void Form::finishCalculation()
{
    timer->stop();
    pushButtonSolve->setEnabled(true);
    tabWidgetMethods->setEnabled(true);
    comboBoxInitial->setEnabled(true);
    spinBoxNX->setEnabled(true);
    spinBoxNT->setEnabled(true);
    sliderNX->setEnabled(true);
    sliderNT->setEnabled(true);
}

void Form::showState()
{
    QChart *chart = nullptr;
    switch(method_)
    {
    case Explicit:
        chart = explicitSolution->chart();
        break;
    case Implicit:
        chart = implicitSolution->chart();
        break;
    case CrankNicolson:
        chart = crankNicolsonSolution->chart();
        break;
    }

    for (auto& series: chart->series())
        series->setOpacity(0.5);

    QLineSeries *series = new QLineSeries();
    chart->addSeries(series);
    series->attachAxis(chart->axisX());
    series->attachAxis(chart->axisY());

    QList<QPointF> data;
    for (decltype(state_.size()) i = 0; i < state_.size(); ++i)
        data << QPointF(((double)i - state_.size()/2) * param_->get_dx(), state_[i]);
    series->append(data);
}
