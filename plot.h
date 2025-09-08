#ifndef plotData_H
#define plotData_H

#include <QMainWindow>
#include <QWidget>
#include <QVector>
#include "qcustomplot.h"
#include "io.h"

QT_BEGIN_NAMESPACE
class QAction;
class QGroupBox;
class QLabel;
class QMenu;
class QMenuBar;
class QPushButton;
class QFileDialog;
class QListWidget;
class QRadioButton;
class QString;
class QCPGraph;
QT_END_NAMESPACE

struct plotCurve
{
    dVector xData; // original (non-concatenated x values); this should be removed from here and put into graphWindow
    dVector yData; // data points
    dVector yFit;   // fit to data
    dVector yError; // error bars (0=none, 1=bars, 2=envelope)
    bVector ignore; // if true, ignore the point when auto-scaling

    bool visible=true;    // allows curves to be turned off (e.g., single run mode)
    bool enabled=true;    // allows curves to be always invisible (but potentially available for export)
    bool histogram=false;
    bool exportCurve=true; // export curve to file? (e.g., ROIs)
    double scaleFactorX=1;
    double scaleFactorY=1;    // this is a property of a single curve
    int pointSizeForData=8;

    QString fileName;      // file source for this curve
    QString legend;
    QColor colorForData=Qt::black;
    QColor colorForFit=Qt::red;
    int lineThicknessForData=1;
    int lineThicknessForFit=2;
    bool dashedForData=false;
    bool dashedForFit=false;
    QCPScatterStyle::ScatterShape pointStyle;
    int errorBars=0;      // 0=none, 1=bars, 2=envelop
};

class plotData : public QMainWindow
{
    Q_OBJECT
private:
    QCustomPlot *_qcplot;              // the QCustomPlot plotting class object
    QWidget *_parentPage;
    QStatusBar *_statusBar=nullptr;

    QVector<QVector<plotCurve>> _listOfCurves;      // [nFiles][nCurvesPerFile]
    QVector<QVector<QCPGraph *>> _listOfDataGraphs; // [nFiles][nCurvesPerFile]

    // Plotting mode and style
    int _plotID;                       // ID different plot surfacers (only the time plot gets some connections)
    QCPRange _autoScaleXRange={0.,0.};
    bool _pressToMoveTracer=true;
    QCPItemTracer *_positionTracer;
    // keep track of position
    int _currentFileInit=0;
    int _currentFile=0;
    int _currentPoint=0;
    // to enable a feature of jumping back to the last point (anywhere), keep the track of the last position
    int _lastFile=0;
    int _lastPoint=0;

    // these are handled in class plot
    bool _autoScale=true;
    bool _normalizeRunMode=false;
    bool _concatenateRunMode=true;
    bool _singleRunMode=true; // if true, show only 1 run and potentially also the fit

    bool _cursorOnPlot=false;
    double _xAxis2Ratio;    // ratio of yAxis2 max to yAxis1 max
    double _yAxis2Ratio;    // ratio of yAxis2 max to yAxis1 max
    int _legendPosition=2;  // 1 = top left, 2 = top right, 3 = bottom left, 4 = bottom right

    int whichConcatenatedFile(double x);
    int whichTracerTimePoint(int iFile, double xPosition);
    void interpretMousePosition(QMouseEvent *event);
    int getTotalDataPoints();
    void concatenateRuns();
    void normalizeCurves();
    void reScaleAxes(QCPRange *xRange, QCPRange *yRange);
    void writeTemporaryStream(bool newDirectory, QString regionName, QTextStream &tempStream);

public:
    plotData(int plotID);
//    virtual ~plotData() {};

    void plotDataAndFit(bool newData);
    void plotDataCurve(plotCurve curve, plotCurve errorBarTargetCurve);
    void plotFitCurve(plotCurve curve);
    void plotIgnoreCurve(plotCurve curve);
    inline void setTracerDisplay(bool state) {_positionTracer->setVisible(state);}
    void writeGraphAllCurves(bool newDirectory, bool newFile, QString fileName, QString regionName);
    void writeGraphJustData(bool newDirectory, bool newGraph, QString fileName, QString regionName);

    inline QCustomPlot *getPlotSurface() const {return _qcplot;}
    inline bool getSingleRunMode() const {return _singleRunMode;}
    inline double getXPosition(int iFile) {return _listOfCurves[iFile][0].xData[_currentPoint];}
    inline double getXPosition(int iFile, int iTime) {return _listOfCurves[iFile][0].xData[iTime];}
    inline double getXPosition(int iFile, int iCurve, int iTime) {return _listOfCurves[iFile][iCurve].xData[iTime];}
    void setCurrentTimeAndPoint(int iFile, double xPosition);
    void setCurrentTimeAndPoint(int iFile, int iTime);
    inline int getNumberFiles() {return _listOfCurves.size();}
    inline int getNumberCurves(int iFile) {return _listOfCurves[iFile].size();}
    inline bool isLegendOn() {return _qcplot->legend->visible();}

    // setting up a series of curves
    void init(int nFiles);
    void conclude(bool singleRun);

    void addCurve(int iFile, QString legend, dVector xData, dVector yData);
    void addCurve(int iFile, QString legend, dVector xData, dVector yData, dVector yFit);
    void addCurve(int iFile, QString legend, dVector xData, dVector yData, bVector ignore);
    void addCurve(int iFile, QString legend, dVector xData, dVector yData, dVector yFit, bVector ignore);
    void addCurve(int iFile, QString legend, dVector xData, dVector yData, dVector yError, dVector yFit);
    void addCurve(int iFile, QString legend, dVector xData, dVector yData, dVector yError, dVector yFit, bVector ignore);

    inline void setDataPointSize(int pointSizeForData) {_listOfCurves[_currentFileInit].last().pointSizeForData=pointSizeForData;
        if (_listOfCurves[_currentFileInit].last().pointStyle == QCPScatterStyle::ssNone)
            _listOfCurves[_currentFileInit].last().pointStyle = QCPScatterStyle::ssDisc;}
    inline void setDataLineThickness(int thickness) {_listOfCurves[_currentFileInit].last().lineThicknessForData=thickness;}
    inline void setFitLineThickness(int thickness) {_listOfCurves[_currentFileInit].last().lineThicknessForFit=thickness;}
    inline void setDataLineDashed(bool state) {_listOfCurves[_currentFileInit].last().dashedForData=state;}
    inline void setFitLineDashed(bool state) {_listOfCurves[_currentFileInit].last().dashedForFit=state;}
    inline void setPointStyle(QCPScatterStyle::ScatterShape pointStyle) {_listOfCurves[_currentFileInit].last().pointStyle=pointStyle;}
    inline void setErrorBars(int errorBars) {_listOfCurves[_currentFileInit].last().errorBars=errorBars;}
    inline void setColorForData(QColor color) {_listOfCurves[_currentFileInit].last().colorForData=color;}
    inline void setColorForFit(QColor color) {_listOfCurves[_currentFileInit].last().colorForFit=color;}
    inline void setVisible(bool visible) {_listOfCurves[_currentFileInit].last().visible=visible;}
    inline void setEnabled(bool enabled) {_listOfCurves[_currentFileInit].last().enabled=enabled;}
    inline void setHistogram(bool histogram) {_listOfCurves[_currentFileInit].last().histogram=histogram;}
    inline void setExport(bool exportCurve) {_listOfCurves[_currentFileInit].last().exportCurve=exportCurve;}
    inline void setLegend(QString legend) {_listOfCurves[_currentFileInit].last().legend=legend;}
    // set scaleFactor for single curve; set _xAxis2Ratio; if necessary, this can be overridden by a later call of setYAxisRatio
    inline void setScaleFactorX(double scaleFactor) {_listOfCurves[_currentFileInit].last().scaleFactorX=scaleFactor; _xAxis2Ratio=1./scaleFactor;}
    // set scaleFactor for single curve; set _yAxis2Ratio; if necessary, this can be overridden by a later call of setYAxisRatio
    inline void setScaleFactorY(double scaleFactor) {_listOfCurves[_currentFileInit].last().scaleFactorY=scaleFactor; _yAxis2Ratio=1./scaleFactor;}
    void setScaleFactorYRelative(int iCurveRef);
    inline void setXAxisRatio(double ratio) {_xAxis2Ratio = ratio;}
    void setLegendPostition(int iPos);

    inline void setQCStatusBar(QStatusBar *statusBar ) {_statusBar = statusBar;}
    inline void setVisibleXAxis(bool visible) {_qcplot->xAxis->setVisible(visible);}
    inline void setLabelXAxis(QString label) {_qcplot->xAxis->setLabel(label);}
    inline void setLabelYAxis(QString label) {_qcplot->yAxis->setLabel(label);}
    inline void setLabelXAxis2(QString label) {_qcplot->xAxis2->setLabel(label);}
    inline void setLabelYAxis2(QString label) {_qcplot->yAxis2->setLabel(label);}
    inline void setMainPage(QWidget *main) {_parentPage = main;}
    inline void setAutoScaleXRange(QCPRange range) {_autoScaleXRange = range;}
    inline void setLegendOn(bool state) {_qcplot->legend->setVisible(state);}
    inline void setAxisRanges(QCPRange xRange, QCPRange yRange)
    {
        _qcplot->xAxis->setRange(xRange);
        _qcplot->yAxis->setRange(yRange);
        _autoScale = false;
    }
    inline void setAutoScale( bool value ) {_autoScale = value;}
    inline void setSingleRunMode( bool mode ) {_singleRunMode = mode;}
    inline void setConcatenatedRunMode(bool state) {_concatenateRunMode = state; concatenateRuns();}
    inline void setNormalizeRunMode(bool state) {_normalizeRunMode = state;}

    // getters
    inline int getLastCurveIndex() {return _listOfCurves[_currentFileInit].count()-1;}
    inline bool containsFit(plotCurve curve)    {return curve.yFit.size()   == curve.xData.size();}
    inline bool containsError(plotCurve curve)  {return curve.yError.size() == curve.xData.size();}
    inline bool containsIgnore(plotCurve curve) {return curve.ignore.size() == curve.xData.size() && curve.ignore.contains(true);}
    inline QString getLabelXAxis() {return _qcplot->xAxis->label();}
    double getMaxYAbs(int iGraph);
    inline void setXRange(QCPRange range) {_qcplot->xAxis->setRange(range); _qcplot->replot();}
    inline void setYRange(QCPRange range) {_qcplot->yAxis->setRange(range); _qcplot->replot();}
    inline bool isAutoScale() {return _autoScale;}
    inline double getXAxisRatio() {return _xAxis2Ratio;}
    inline double getYAxisRatio() {return _yAxis2Ratio;}
    inline int getNumberPoints(int iFile) {return _listOfCurves[iFile][0].xData.size();}
    inline int getNumberPoints(int iFile, int iCurve) {return _listOfCurves[iFile][iCurve].xData.size();}
    inline double getXData(int iFile, int iCurve, int iPoint) {return _listOfCurves[iFile][iCurve].xData[iPoint];}
    inline double getYData(int iFile, int iCurve, int iPoint) {return _listOfCurves[iFile][iCurve].yData[iPoint];}
    inline dVector getYData(int iFile, int iCurve) {return _listOfCurves[iFile][iCurve].yData;}
    inline int getCurrentFile()  {return _currentFile;}
    inline int getCurrentPoint() {return _currentPoint;}
private slots:
    void axisLabelDoubleClick(QCPAxis *axis, QCPAxis::SelectablePart part);
    void legendDoubleClick(QCPLegend *legend, QCPAbstractLegendItem *item);
    void contextMenuRequest(QPoint pos);
    void makeSelectedGraphDotted();
    void makeSelectedGraphBiggerPoints();
    void makeSelectedGraphSmallerPoints();
    void removeSelectedGraph();
    void writeSelectedGraph();
    void moveLegend();
    void plotMousePress(QMouseEvent *event);
    void plotMouseMove(QMouseEvent *event);
    void writeOneGraph();
    void changeSelectMethod();
    void popUpChangeXRange();
    void popUpChangeYRange();

    void lastTimePoint();
    void setXAxis2Range(QCPRange Range);
    void setYAxis2Range(QCPRange Range);

    void keyboardPlus();
    void keyboardMinus();

public slots:
    void setXZoom();
    void setYZoom();
    void autoScale(bool state);
    void setSelectPoints(bool pointNotCurve);
    inline void showLegend() {_qcplot->legend->setVisible(true); _qcplot->replot();}
    inline void hideLegend() {_qcplot->legend->setVisible(false); _qcplot->replot();}
signals:
    void autoScaleRanges(bool state);
    void changedPointFromGraph(int iDataCurve, int iTime);
};

#endif // plotData_H
