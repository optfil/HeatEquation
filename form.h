#ifndef FORM_H
#define FORM_H

#include <vector>

#include <QComboBox>
#include <QPushButton>
#include <QSlider>
#include <QSpinBox>
#include <QTabWidget>
#include <QTimer>
#include <QWidget>

#include <QtCharts/QtCharts>
QT_CHARTS_USE_NAMESPACE

#include "parameters.h"

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
    void update_nx_from_slider(int log_n);
    void update_nx(int n);
    void update_nt(int n);
    void selectionChanged();
    void updateLabels();
    void updateSpectrum();
    void initiateState();
    void updateDispersionDiffusion();
    void Solve();
    void Tick();

private:
    QChartView *chartView;
    QLabel *labelInitial;
    QComboBox *comboBoxInitial;
    QLabel *labelSizeX_1, *labelSizeX_2, *labelSizeT_1, *labelSizeT_2, *labelNX_1, *labelNX_2, *labelNT_1, *labelNT_2;
    QLabel *labelSizeX, *labelSizeT;
    QSlider *sliderNX, *sliderNT;
    QSpinBox *spinBoxNX, *spinBoxNT;
    QLabel *labelStepX_1, *labelStepX_2, *labelStepX;
    QLabel *labelStepT_1, *labelStepT_2, *labelStepT;
    QLabel *labelAlpha_1, *labelAlpha_2, *labelAlpha;
    QPushButton *pushButtonSolve;
    QTabWidget *tabWidgetMethods;
    QWidget *widgetExplicit, *widgetImplicit, *widgetCrankNicolson;
    QChartView *explicitDispersion, *explicitDissipation, *explicitSolution;
    QChartView *implicitDispersion, *implicitDissipation, *implicitSolution;
    QChartView *crankNicolsonDispersion, *crankNicolsonDissipation, *crankNicolsonSolution;
    QLineSeries *seriesInitial;
    QLineSeries *seriesExplicitIdealDispersion, *seriesImplicitIdealDispersion, *seriesCrankNicolsonIdealDispersion;
    QLineSeries *seriesExplicitIdealDissipation, *seriesImplicitIdealDissipation, *seriesCrankNicolsonIdealDissipation;
    QLineSeries *seriesExplicitDispersion, *seriesImplicitDispersion, *seriesCrankNicolsonDispersion;
    QLineSeries *seriesExplicitDissipation, *seriesImplicitDissipation, *seriesCrankNicolsonDissipation;

    QTimer *timer;

    Parameters *param_;
    MethodType method_;
    std::vector<double> state_, tmp_state_, tdma_u_, tdma_v_;
    double t_cur_;

    void showState();
    void finishCalculation();
    void cleanSolution();
};

#endif // FORM_H
