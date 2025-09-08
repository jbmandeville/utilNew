#include <QtWidgets>
#include <QFrame>
#include <QDebug>
#include <QVector>
#include <QObject>
#include <QFileDialog>
#include <QString>
#include <QPen>
#include <QListWidget>
#include <QListWidgetItem>
#include <QSize>

#include "plot.h"

/*
plotData::~plotData()
{
}
*/

plotData::plotData(int plotID)
{
    FUNC_ENTER;
//    _parentPage = containingPage;

    _plotID = plotID;

    _qcplot = new QCustomPlot();
    _qcplot->setObjectName(QStringLiteral("plot"));
    _qcplot->xAxis->setLabel("time");
    _qcplot->yAxis->setLabel("y");
    _qcplot->xAxis2->setLabel("");
    _qcplot->yAxis2->setLabel("");
    _qcplot->xAxis->setLabelFont(QFont("Arial", 20));
    _qcplot->yAxis->setLabelFont(QFont("Helvetica", 20));
    _qcplot->xAxis->setTickLabelFont(QFont("Arial", 16));
    _qcplot->yAxis->setTickLabelFont(QFont("Helvetica",16));
    _qcplot->xAxis2->setVisible(true);
    _qcplot->xAxis2->setTickLabels(false);
    _qcplot->yAxis2->setVisible(true);
    _qcplot->yAxis2->setTickLabels(false);
//    _qcplot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables | QCP::iSelectAxes);
//    _qcplot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    _qcplot->setAutoAddPlottableToLegend(true);
    _qcplot->legend->setVisible(false);
    _qcplot->legend->setFont(QFont("Arial", 20));

    // add the phase tracer (red circle) which sticks to the graph data (and gets updated in bracketDataSlot by timer event):
    _positionTracer = new QCPItemTracer(_qcplot);
    _positionTracer->setStyle(QCPItemTracer::tsCircle);
    _positionTracer->setPen(QPen(Qt::red));
    _positionTracer->setBrush(Qt::cyan);
    _positionTracer->setSize(16);
//    _positionTracer->setVisible(false);
    _positionTracer->setGraphKey(0);
    _positionTracer->setVisible(true);
    _positionTracer->position->setAxes(_qcplot->xAxis,_qcplot->yAxis);

    _xAxis2Ratio = _yAxis2Ratio = 0.;

    setSelectPoints(true);

    connect(_qcplot, SIGNAL(axisClick(QCPAxis*,QCPAxis::SelectablePart,QMouseEvent*)), this,
            SLOT(axisLabelDoubleClick(QCPAxis*,QCPAxis::SelectablePart)));
//    connect(_qcplot, SIGNAL(axisDoubleClick(QCPAxis*,QCPAxis::SelectablePart,QMouseEvent*)), this,
//            SLOT(axisLabelDoubleClick(QCPAxis*,QCPAxis::SelectablePart)));
    connect(_qcplot, SIGNAL(legendDoubleClick(QCPLegend*,QCPAbstractLegendItem*,QMouseEvent*)), this,
            SLOT(legendDoubleClick(QCPLegend*,QCPAbstractLegendItem*)));
    connect(_qcplot, SIGNAL(mousePress(QMouseEvent*)),this,SLOT(plotMousePress(QMouseEvent*)));
    connect(_qcplot, SIGNAL(mouseMove(QMouseEvent*)), this,SLOT(plotMouseMove(QMouseEvent*)));

    connect(_qcplot->xAxis, SIGNAL(rangeChanged(QCPRange)), this, SLOT(setXAxis2Range(QCPRange)));
    connect(_qcplot->yAxis, SIGNAL(rangeChanged(QCPRange)), this, SLOT(setYAxis2Range(QCPRange)));

    // setup policy and connect slot for context menu popup:
    _qcplot->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(_qcplot, SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(contextMenuRequest(QPoint)));

    if ( _plotID == 0 )
    { // for the time-series plot only
        QAction *plusAction = new QAction(this);
        plusAction->setShortcut(Qt::Key_Equal);
        connect(plusAction, SIGNAL(triggered()), this, SLOT(keyboardPlus()));
        QAction *minusAction = new QAction(this);
        minusAction->setShortcut(Qt::Key_Minus);
        connect(minusAction, SIGNAL(triggered()), this, SLOT(keyboardMinus()));

        _qcplot->addAction(plusAction);
        _qcplot->addAction(minusAction);
    }

    FUNC_EXIT;
}

void plotData::init(int nFiles)
{
    _listOfCurves.resize(nFiles);
    for (int jFile=0; jFile<nFiles; jFile++)
        _listOfCurves[jFile].clear();  // curves will be appended
    setLegendOn(false);
}

void plotData::keyboardMinus()
{
    FUNC_ENTER << _currentFile << _currentPoint;
    int iPoint = _currentPoint - 1;
    int iFile  = _currentFile;
    if ( iPoint < 0 )
    { // go to the previous file
        iFile--;
        if ( iFile < 0 ) iFile  = getNumberFiles() - 1;
        iPoint = getNumberPoints(iFile,0) - 1;
    }
    setCurrentTimeAndPoint(iFile, iPoint);
    emit changedPointFromGraph(_currentFile,_currentPoint);
    FUNC_EXIT << _currentFile << _currentPoint;
}

void plotData::keyboardPlus()
{
    FUNC_ENTER << _currentFile << _currentPoint;
    int iPoint = _currentPoint + 1;
    int iFile  = _currentFile;
    if ( iPoint >= getNumberPoints(iFile) )
    { // go to the next file
        iFile++;
        if ( iFile >= getNumberFiles() ) iFile  = 0;
        iPoint=0;
    }
    setCurrentTimeAndPoint(iFile, iPoint);
    emit changedPointFromGraph(_currentFile,_currentPoint);
    FUNC_EXIT << _currentFile << _currentPoint;
}

///////////////////////////////////////////////////////////////////
///////////////////////// slots ///////////////////////////////////
///////////////////////////////////////////////////////////////////
void plotData::setXZoom()
{
    _qcplot->setCursor(Qt::OpenHandCursor);
    _qcplot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    _qcplot->axisRect(0)->setRangeDrag(Qt::Horizontal);
    _qcplot->axisRect(0)->setRangeZoom(Qt::Horizontal);
    autoScale(false);
}
void plotData::setYZoom()
{
    _qcplot->setCursor(Qt::OpenHandCursor);
    _qcplot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    _qcplot->axisRect(0)->setRangeDrag(Qt::Vertical);
    _qcplot->axisRect(0)->setRangeZoom(Qt::Vertical);
    autoScale(false);
}
void plotData::setSelectPoints(bool pointNotCurve)
{
    if ( pointNotCurve )
    {
        _pressToMoveTracer = true;
        _qcplot->setCursor(Qt::CrossCursor);
        if ( _autoScale )
            _qcplot->setInteractions(QCP::iNone);
        else
            _qcplot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    }
    else
    {
        _pressToMoveTracer = false;
        _qcplot->setCursor(Qt::ArrowCursor);
        if ( _autoScale )
            _qcplot->setInteractions(QCP::iSelectPlottables);
        else
            _qcplot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    }
}
void plotData::lastTimePoint()
{
    FUNC_ENTER << _lastFile << _lastPoint << _currentFile << _currentPoint;
    setCurrentTimeAndPoint(_lastFile, _lastPoint);
    emit changedPointFromGraph(_currentFile,_currentPoint);
}

void plotData::setXAxis2Range( QCPRange Range1 )
{
    QCPRange Range2;
    Range2.lower = Range1.lower * _xAxis2Ratio;
    Range2.upper = Range1.upper * _xAxis2Ratio;
    FUNC_INFO << "Range2" << _xAxis2Ratio << Range1 << Range2;
    _qcplot->xAxis2->setRange(Range2);
}
void plotData::setYAxis2Range( QCPRange Range1 )
{
    QCPRange Range2;
    Range2.lower = Range1.lower * _yAxis2Ratio;
    Range2.upper = Range1.upper * _yAxis2Ratio;
    _qcplot->yAxis2->setRange(Range2);
}
void plotData::axisLabelDoubleClick(QCPAxis *axis, QCPAxis::SelectablePart part)
{
  // Set an axis label by double clicking on it
  if (part == QCPAxis::spAxisLabel) // only react when the actual axis label is clicked, not tick label or axis backbone
  {
    bool ok;
    QString newLabel = QInputDialog::getText(this, "QCustomPlot example", "New axis label:", QLineEdit::Normal, axis->label(), &ok);
    if (ok)
    {
      axis->setLabel(newLabel);
      _qcplot->replot();
    }
  }
}
void plotData::legendDoubleClick(QCPLegend *legend, QCPAbstractLegendItem *item)
{
  // Rename a graph by double clicking on its legend item
  Q_UNUSED(legend)
  if (item) // only react if item was clicked (user could have clicked on border padding of legend where there is no item, then item is 0)
  {
    QCPPlottableLegendItem *plItem = qobject_cast<QCPPlottableLegendItem*>(item);
    bool ok;
    QString newName = QInputDialog::getText(this, "QCustomPlot example", "New graph name:", QLineEdit::Normal, plItem->plottable()->name(), &ok);
    if (ok)
    {
      plItem->plottable()->setName(newName);
      _qcplot->replot();
    }
  }
}
void plotData::contextMenuRequest(QPoint pos)
{
    QMenu *menu = new QMenu(this);
    menu->setAttribute(Qt::WA_DeleteOnClose);

    if (_qcplot->legend->selectTest(pos, false) >= 0 ) // context menu on legend requested
        //  if (_qcplot->legend->selectTest(pos, false) >= 0 && _qcplot->legend->visible() ) // context menu on legend requested
    {
        if  ( _qcplot->legend->visible() )
            menu->addAction("Hide legend", this, SLOT(hideLegend()));
        else
            menu->addAction("Show legend", this, SLOT(showLegend()));

        menu->addAction("Move to top left", this, SLOT(moveLegend()))->setData(static_cast<int>((Qt::AlignTop|Qt::AlignLeft)));
        menu->addAction("Move to top center", this, SLOT(moveLegend()))->setData(static_cast<int>((Qt::AlignTop|Qt::AlignHCenter)));
        menu->addAction("Move to top right", this, SLOT(moveLegend()))->setData(static_cast<int>((Qt::AlignTop|Qt::AlignRight)));
        menu->addAction("Move to bottom right", this, SLOT(moveLegend()))->setData(static_cast<int>((Qt::AlignBottom|Qt::AlignRight)));
        menu->addAction("Move to bottom left", this, SLOT(moveLegend()))->setData(static_cast<int>((Qt::AlignBottom|Qt::AlignLeft)));
    }
    else if (_qcplot->selectedGraphs().size() > 0)
    {
        menu->addAction("Write selected graph", this, SLOT(writeSelectedGraph()));
        menu->addAction("Remove selected graph", this, SLOT(removeSelectedGraph()));
        menu->addAction("LineStyle dotted", this, SLOT(makeSelectedGraphDotted()));
        if ( _qcplot->selectedGraphs().at(0)->scatterStyle().size() != 0. )
        {
            menu->addAction("Bigger points", this, SLOT(makeSelectedGraphBiggerPoints()));
            menu->addAction("Smaller points", this, SLOT(makeSelectedGraphSmallerPoints()));
        }
    }
    else
    {
        menu->addAction("Change x range", this, SLOT(popUpChangeXRange()));
        menu->addAction("Change y range", this, SLOT(popUpChangeYRange()));
        if ( _pressToMoveTracer )
            menu->addAction("Select one curve", this, SLOT(changeSelectMethod()));
        else
            menu->addAction("Select one point", this, SLOT(changeSelectMethod()));
        menu->addAction("Write all curves to table file", this, SLOT(writeOneGraph()));
    }
    /*
  if (_qcplot->selectedGraphs().size() > 0)
  {
      menu->addAction("Remove selected graph", this, SLOT(removeSelectedGraph()));
      menu->addAction("LineStyle dotted", this, SLOT(makeSelectedGraphDotted()));
      menu->addAction("Bigger points", this, SLOT(makeSelectedGraphBiggerPoints()));
      menu->addAction("Smaller points", this, SLOT(makeSelectedGraphSmallerPoints()));
  }
*/
    menu->popup(_qcplot->mapToGlobal(pos));
}

void plotData::popUpChangeXRange()
{
    QCPRange xRange = _qcplot->xAxis->range();
    bool ok;
    QString lower, upper, range;
    range = lower.setNum(xRange.lower) + " " + upper.setNum(xRange.upper);
    QString rangeString = QInputDialog::getText(this, "X Range:", "min,max (use comma/space)", QLineEdit::Normal, range, &ok);
    if ( ok )
    {
        QRegularExpression whiteSpaceComma("[,\\s]");// match a comma or a space
        QStringList valueList = rangeString.split(whiteSpaceComma,Qt::SkipEmptyParts); // clazy:exclude=use-static-qregularexpression
        if ( valueList.size() != 2 ) return;
        lower = valueList[0];  upper = valueList[1];
        bool ok1, ok2;
        double lowerValue = lower.toDouble(&ok1);
        double upperValue = upper.toDouble(&ok2);
        if ( ok1 && ok2 )
        {
            xRange.lower = qMin(lowerValue,upperValue);
            xRange.upper = qMax(lowerValue,upperValue);
            _qcplot->xAxis->setRange(xRange);
            _qcplot->replot();
        }
        _autoScale = false;
        emit autoScaleRanges(_autoScale);
    }
}
void plotData::popUpChangeYRange()
{
    QCPRange yRange = _qcplot->yAxis->range();
    bool ok;
    QString lower, upper, range;
    range = lower.setNum(yRange.lower) + " " + upper.setNum(yRange.upper);
    QString rangeString = QInputDialog::getText(this, "Y Range:", "min,max (use comma/space)", QLineEdit::Normal, range, &ok);
    if ( ok )
    {
        QRegularExpression whiteSpaceComma("[,\\s]");// match a comma or a space
        QStringList valueList = rangeString.split(whiteSpaceComma,Qt::SkipEmptyParts); // clazy:exclude=use-static-qregularexpression
        if ( valueList.size() != 2 ) return;
        lower = valueList[0];  upper = valueList[1];
        bool ok1, ok2;
        double lowerValue = lower.toDouble(&ok1);
        double upperValue = upper.toDouble(&ok2);
        if ( ok1 && ok2 )
        {
            yRange.lower = qMin(lowerValue,upperValue);
            yRange.upper = qMax(lowerValue,upperValue);
            _qcplot->yAxis->setRange(yRange);
            _qcplot->replot();
        }
        _autoScale = false;
        emit autoScaleRanges(_autoScale);
    }
}

void plotData::changeSelectMethod()
{
    setSelectPoints(!_pressToMoveTracer);
}

void plotData::makeSelectedGraphDotted()
{
  if (_qcplot->selectedGraphs().size() > 0)
    {
        QCPScatterStyle currentScatter = _qcplot->selectedGraphs().at(0)->scatterStyle();
        QPen currentPen = _qcplot->selectedGraphs().at(0)->pen();
        //      QBrush currentBrush = _qcplot->selectedGraphs().at(0)->brush();
        currentPen.setStyle(Qt::DotLine);
        //      currentPen.setColor(Qt::green);
        //      currentScatter.setBrush(Qt::green);
        //      currentScatter.setPen(currentPen);
        //      currentScatter.setSize(5);
        //      currentScatter.setShape(QCPScatterStyle::ssCircle);
        _qcplot->selectedGraphs().at(0)->setScatterStyle(currentScatter);
        _qcplot->selectedGraphs().at(0)->setPen(currentPen);
        _qcplot->replot();
    }
}
void plotData::makeSelectedGraphBiggerPoints()
{
  if (_qcplot->selectedGraphs().size() > 0)
  {
      QCPScatterStyle currentScatter = _qcplot->selectedGraphs().at(0)->scatterStyle();
      currentScatter.setSize(1.5*currentScatter.size());
      _qcplot->selectedGraphs().at(0)->setScatterStyle(currentScatter);
      _qcplot->replot();
  }
}
void plotData::makeSelectedGraphSmallerPoints()
{
  if (_qcplot->selectedGraphs().size() > 0)
  {
        QCPScatterStyle currentScatter = _qcplot->selectedGraphs().at(0)->scatterStyle();
      currentScatter.setSize(0.75*currentScatter.size());
      _qcplot->selectedGraphs().at(0)->setScatterStyle(currentScatter);
      _qcplot->replot();
  }
}
void plotData::removeSelectedGraph()
{
  if (_qcplot->selectedGraphs().size() > 0)
  {
    _qcplot->removeGraph(_qcplot->selectedGraphs().at(0));
    _qcplot->replot();
  }
}

void plotData::writeOneGraph()
{
    QString fileName = "/graph.dat";
    QFileDialog fileDialog;
    QString fullFileName;
    fullFileName = fileDialog.getSaveFileName(this,
                                              "Name of file",
                                              QDir::currentPath()+fileName,
                                              tr("Text files (*.txt *.dat *.roi)"));
    if ( fullFileName.isEmpty() ) return;
    writeGraphAllCurves(true,true,fullFileName,"");
}

void plotData::writeTemporaryStream(bool newDirectory, QString regionName, QTextStream &tempStream)
{
    FUNC_ENTER;
    if ( newDirectory ) tempStream << "directory " << QDir::currentPath() << "\n";
    // Write the headers; assume all files have the same curves to be exported
    FUNC_INFO << "nfiles" << getNumberFiles();
    for (int jFile=0; jFile<getNumberFiles(); jFile++)
    {
        tempStream << "new File\n";
        // Write headers
        tempStream << "index time ";
        for (int jCurve=0; jCurve<getNumberCurves(jFile); jCurve++)
        {
            plotCurve currentCurve = _listOfCurves[jFile][jCurve];
            FUNC_INFO << "file" << jFile << "curve" << jCurve << currentCurve.legend;
            FUNC_INFO << "yData" << currentCurve.yData;
            if ( currentCurve.exportCurve )
            {
                if ( currentCurve.legend == "data" && !regionName.isEmpty() )
                    tempStream << regionName << " ";
                else
                    tempStream << currentCurve.legend << " ";
                if ( currentCurve.errorBars != 0 )
                    tempStream << currentCurve.legend << "_err ";
                if ( currentCurve.scaleFactorY != 1. )
                    tempStream << currentCurve.legend << "_scaled ";
            }
        } // jCurve
        tempStream << "\n";
        // loop over time points and write all curves for this file in columns
        for (int jt=0; jt<_listOfCurves[jFile][0].yData.size(); jt++)
        {
            tempStream << jt << " " << _listOfCurves[jFile][0].xData[jt] << " ";
            for (int jCurve=0; jCurve<_listOfCurves[jFile].size(); jCurve++)
            {
                plotCurve currentCurve = _listOfCurves[jFile][jCurve];
                if ( currentCurve.exportCurve ) // if correct file and export
                {
//                    if ( jt==0) FUNC_INFO << "output" << currentCurve.legend << currentCurve.yData[jt];
                    if ( currentCurve.scaleFactorY != 1. )
                        tempStream << currentCurve.yData[jt] / currentCurve.scaleFactorY << " ";
                    else
                        tempStream << currentCurve.yData[jt] << " ";
                    if ( currentCurve.errorBars != 0 )
                        tempStream << currentCurve.yError[jt] << " ";
                } // write curve
            } // jCurve
            tempStream << "\n";
        } // jt
        tempStream << "\n";
    } // jFile
    tempStream.flush();//flush the stream into the file
    FUNC_EXIT;
}


void plotData::writeGraphAllCurves(bool newDirectory, bool newFile, QString fileName, QString regionName)
{
    FUNC_ENTER << newDirectory << newFile << fileName << regionName;
    QTemporaryFile tempFile;
    if (!tempFile.open()) return;
    QTextStream tempStream(&tempFile);
    writeTemporaryStream(newDirectory,  regionName, tempStream);
    tempStream.seek(0);

    QFile file(fileName);
    QTextStream outStream(&file);

    if ( newFile )
    { // open the output file as write only and do a copy/paste
        if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return;
    }
    else
    { // not new file: append
        if (!file.open(QIODevice::Append | QIODevice::Text)) return;
        if ( !newDirectory ) outStream << "new ROI\n";
    }
    outStream << tempFile.readAll(); // just copy/paste
    file.close();
    tempFile.close();
}

void plotData::writeGraphJustData(bool newDirectory, bool newFile, QString fileName, QString regionName)
{
    FUNC_ENTER << newDirectory << newFile << fileName << regionName;
    QTemporaryFile tempFile;
    if (!tempFile.open()) return;
    QTextStream tempStream(&tempFile);
    writeTemporaryStream(newDirectory,  regionName, tempStream);
    tempStream.seek(0);

    if ( newFile ) // straightforward: just take 3 first words
    { // open the output file as write only and take the 1st 3 words (index time region)
        QFile file(fileName);
        QTextStream outStream(&file);
        if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return;
         // just take the 1st 3 words of the temporary file (index time region)
        QString line = tempStream.readLine();
        while (!tempStream.atEnd() )
        {
            QRegularExpression whiteSpaceComma("[,\\s]");// match a comma or a space
            QStringList stringList = line.split(whiteSpaceComma,Qt::SkipEmptyParts); // clazy:exclude=use-static-qregularexpression
            FUNC_INFO << "lineOrig" << line << "stringList" << stringList;
            if ( stringList.count() < 3 ) // blank line or "new file"
                outStream << line << "\n";
            else // index time region
                outStream << stringList.at(0) << " " << stringList.at(1) << " " << stringList.at(2) << "\n";
            line = tempStream.readLine();
        }
        file.close();
    }
    else // take all words and append the data column
    { // not new file: append
        // open the output file 1st as read and then as write
        QFile file(fileName);
        QTextStream outStream(&file);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return;
        QTemporaryFile tempFileOriginal;
        QTextStream tempStreamOriginal(&tempFileOriginal);
        if (!tempFileOriginal.open()) return;
        tempStreamOriginal << outStream.readAll(); // first make a copy of the existing data.
        tempStreamOriginal.seek(0);
        file.close();

        if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return;
        outStream.setDevice(&file);  outStream.seek(0);
        if ( !newDirectory )
        { // only the 1st/original has the directory info
            QString lineOriginal = tempStreamOriginal.readLine();  // to match prior write
            outStream << lineOriginal << "\n";
        }
        while (!tempStreamOriginal.atEnd())
        {
            QString lineOriginal = tempStreamOriginal.readLine();
            QString lineNew      = tempStream.readLine();
            FUNC_INFO << "original" << lineOriginal;
            FUNC_INFO << "new line" << lineNew;
            QRegularExpression whiteSpaceComma("[,\\s]");// match a comma or a space
            QStringList stringListOrig = lineOriginal.split(whiteSpaceComma,Qt::SkipEmptyParts); // clazy:exclude=use-static-qregularexpression
            QStringList stringListNew  = lineNew.split(whiteSpaceComma,Qt::SkipEmptyParts); // clazy:exclude=use-static-qregularexpression
            if ( stringListNew.count() < 3 ) // blank line or "new file"
                outStream << lineNew << "\n";
            else
            {
                // index time region, but use only region
                //                    FUNC_INFO << "new line:" << lineOriginal << " " << stringList.at(2) << "\n";
                if ( !stringListOrig.at(0).compare("directory") )
                    outStream << lineOriginal << "\n";
                else
                    outStream << lineOriginal << " " << stringListNew.at(2) << "\n";
            }
            //                FUNC_INFO << "end it?" << tempStreamOriginal.atEnd();
        }
        file.close();
        tempFileOriginal.close();

    }
    tempFile.close();
}

void plotData::writeSelectedGraph()
{
    if (_qcplot->selectedGraphs().size() == 1)
    {
        QString fileName;
        QFileDialog fileDialog;
        fileName = fileDialog.getSaveFileName(this,
                                              "Name of file",
                                              QDir::currentPath(),
                                              "(*.dat)");
        if ( fileName.isEmpty() ) return;
        QFile file(fileName);
        QTextStream out(&file);
        if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
            return;
        out << "x y\n";
        int nTime = _qcplot->selectedGraphs().at(0)->dataCount();
        for (int jt=0; jt<nTime; jt++)
            out << _qcplot->selectedGraphs().at(0)->dataMainKey(jt)   << " "
                << _qcplot->selectedGraphs().at(0)->dataMainValue(jt) << "\n";
        file.close();
    }
}
void plotData::setLegendPostition(int iPos)
{
    FUNC_ENTER << iPos;
    _legendPosition = iPos;
    Qt::Alignment align;
    if ( _legendPosition == 1 )
        align =  Qt::AlignLeft | Qt::AlignTop;
    else if ( _legendPosition == 2 )
        align =  Qt::AlignRight | Qt::AlignTop;
    else if ( _legendPosition == 3 )
        align =  Qt::AlignLeft | Qt::AlignBottom;
    else if ( _legendPosition == 4 )
        align =  Qt::AlignRight | Qt::AlignBottom;
    else
    {
        _legendPosition = 2;
        align =  Qt::AlignRight | Qt::AlignTop;
    }
    _qcplot->axisRect()->insetLayout()->setInsetAlignment(0, align);
    _qcplot->replot();
}
void plotData::moveLegend()
{
    FUNC_ENTER;
    if (QAction* contextAction = qobject_cast<QAction*>(sender())) // make sure this slot is really called by a context menu action, so it carries the data we need
    {
        bool ok;
        int dataInt = contextAction->data().toInt(&ok);
        if (ok)
        {
            _qcplot->axisRect()->insetLayout()->setInsetAlignment(0, static_cast<Qt::Alignment>(dataInt));
            _qcplot->replot();
            if ( (dataInt & Qt::AlignLeft) && (dataInt & Qt::AlignTop) )
                _legendPosition = 1;
            else if ( (dataInt & Qt::AlignRight) && (dataInt & Qt::AlignTop) )
                _legendPosition = 2;
            else if ( (dataInt & Qt::AlignLeft) && (dataInt & Qt::AlignBottom) )
                _legendPosition = 3;
            else if ( (dataInt & Qt::AlignRight) && (dataInt & Qt::AlignBottom) )
                _legendPosition = 4;
            else
                _legendPosition = 2;
        }
    }
}

void plotData::setCurrentTimeAndPoint(int iFile, int iPoint)
{ // xxx needs work
    FUNC_ENTER << iFile << iPoint;
    if ( _plotID == 0 ) FUNC_ENTER << iFile << iPoint << _plotID;
    if ( iFile < 0  || iFile >= getNumberFiles() ) return;
    if ( iPoint < 0 || iPoint >= getNumberPoints(iFile,0) ) return;

    if ( _plotID == 0 ) FUNC_INFO << "before set" << _lastFile << _lastPoint << _currentFile << _currentPoint;
    bool notSamePoint = (iFile != _currentFile) || (iPoint != _currentPoint);
    if ( notSamePoint )
    {
        _lastFile = _currentFile;
        _lastPoint = _currentPoint;
    }

    // set the file/point indices
    _currentFile  = iFile;
    _currentPoint = iPoint;
    if ( _plotID == 0 ) FUNC_INFO << "set values" << _currentFile << _currentPoint;

    double x = _listOfCurves[_currentFile][0].xData.at(_currentPoint);
    double y = _listOfCurves[_currentFile][0].yData.at(_currentPoint);
    QString message;
    message = QString("(scan,pt) = (%1 , %2)    (x,y) = (%3 , %4)").arg(0+1).arg(_currentPoint+1).arg(x).arg(y);
    if ( _plotID == 0 ) FUNC_INFO << "position" << message;
    if ( _plotID == 0 ) FUNC_INFO << "exit" << _lastFile << _lastPoint << _currentFile << _currentPoint << "\n";
}

void plotData::setCurrentTimeAndPoint(int iFile, double xPosition)
{ // xPosition is actual x value on axis, with or without concatenation
    FUNC_ENTER << "***" << iFile << xPosition << "_listOfDataGraphs.size" << _listOfDataGraphs.size() << "plotID" << _plotID;
    FUNC_INFO << "size2" << _listOfDataGraphs[iFile].size();
    int iPoint = whichTracerTimePoint(iFile, xPosition);
    setCurrentTimeAndPoint(iFile, iPoint);
}

void plotData::plotMousePress(QMouseEvent *event)
{
    double x = _qcplot->xAxis->pixelToCoord(event->pos().x());
    // double y = _qcplot->yAxis->pixelToCoord(event->pos().y());
    if ( event->buttons() == Qt::LeftButton && _pressToMoveTracer )
    {
        if ( !_concatenateRunMode )
        {
            // then don't identify the curve; keep the same one
            int iFile = _currentFile;
            int iPoint = whichTracerTimePoint(iFile,x);
            FUNC_INFO << "setCurrentTimeAndPoint" << iFile << iPoint;
            setCurrentTimeAndPoint(iFile, iPoint);  // store last time point info
            emit changedPointFromGraph(_currentFile,_currentPoint);
        }
        else
        {
            // Need to indentify which curve and which time point
            int iFile = whichConcatenatedFile(x);
            int iPoint = whichTracerTimePoint(iFile,x);
            FUNC_INFO << "setCurrentTimeAndPoint" << iFile << iPoint;
            setCurrentTimeAndPoint(iFile, iPoint);   // store last time point info
            emit changedPointFromGraph(_currentFile,_currentPoint);
        }
        _qcplot->replot();
    }
    FUNC_EXIT;
}

int plotData::getTotalDataPoints()
{
    int nPointsTotal = 0;
    for (int jFile=0; jFile<getNumberFiles(); jFile++)
        nPointsTotal += getNumberPoints(jFile,0);
    return nPointsTotal;
}

int plotData::whichConcatenatedFile( double x )
{
    double firstX, lastX;
    for (int jFile=0; jFile<getNumberFiles(); jFile++)
    {
        plotCurve dataCurve = _listOfCurves[jFile][0];
        firstX = dataCurve.xData[0];
        lastX  = dataCurve.xData[getNumberPoints(jFile,0)-1];
        if ( firstX == lastX )
        {
            firstX -= 0.5;  lastX += 0.5;
        }
//        FUNC_INFO << "x, firstX, lastX" << x << firstX << lastX << _listOfCurves[jCurve].xData.size();
//        if ( _listOfCurves[jCurve].xData.size() == 1 || (x >= firstX && x<= lastX) )
        if ( x >= firstX && x<= lastX ) return jFile;
    }
    return 0;
}

int plotData::whichTracerTimePoint(int iFile, double xPosition)
{
    double distanceMin = 1.e10;
    int iPoint=0;
    for (int jPoint=0; jPoint<getNumberPoints(iFile,0); jPoint++)
    {
        double distance = qAbs(xPosition - _listOfCurves[iFile][0].xData[jPoint]);
        if ( distance < distanceMin )
        {
            distanceMin = distance;
            iPoint = jPoint;
        }
    }
    return iPoint;
}

void plotData::plotMouseMove(QMouseEvent *event)
{
    interpretMousePosition(event);
    // bool leftButtonMove = (event->type() == QEvent::MouseMove) && (event->buttons() == Qt::LeftButton);

    if ( _cursorOnPlot )
    {
        double x = _qcplot->xAxis->pixelToCoord(event->pos().x());
        double y = _qcplot->yAxis->pixelToCoord(event->pos().y());
//        double y2 = _qcplot->yAxis2->pixelToCoord(event->pos().y());
        QString message = QString("%1 , %2").arg(x).arg(y);
        _qcplot->setToolTip(message);
        if ( _statusBar )
            _statusBar->showMessage(message);
    }
    else
    {
        QCPItemPosition *position = _positionTracer->position;
        QString message;
        message = QString("(scan,pt) = (%1 , %2)    (x,y) = (%3 , %4)").arg(_currentFile+1).arg(_currentPoint+1)
                .arg(position->coords().x()).arg(position->coords().y());
        if ( _statusBar ) _statusBar->showMessage(message);
    }
}

void plotData::interpretMousePosition(QMouseEvent *event)
{
    qreal x = event->localPos().x();
    qreal y = event->localPos().y();
    _cursorOnPlot = false;

    if ( x >= _qcplot->axisRect()->left() && x <=_qcplot->axisRect()->rect().right() &&
         y >= _qcplot->axisRect()->top()  && y <=_qcplot->axisRect()->rect().bottom() )
        _cursorOnPlot = true;
}

double plotData::getMaxYAbs(int iGraph)
{
    plotCurve *curve = &(_listOfCurves[0][iGraph]);
    double dMax = 0.;
    for ( int jPoint=0; jPoint<curve->yData.size(); jPoint++ )
    {
        if ( qAbs(curve->yData[jPoint]) > dMax ) dMax = qAbs(curve->yData[jPoint]);
    }
    return dMax;
}

void plotData::conclude(bool singleRun)
{
    FUNC_ENTER;

    if ( _currentFile < 0  || _currentFile  >= getNumberFiles() )              _currentFile = 0;
    if ( _currentPoint < 0 || _currentPoint >= getNumberPoints(_currentFile) ) _currentPoint = 0;
    bool somePoints = getTotalDataPoints() > 200 && getTotalDataPoints() < 800;
    bool manyPoints = getTotalDataPoints() > 800;
    for (int jFile=0; jFile<getNumberFiles(); jFile++)
    {
        bool visible = true;
        if ( singleRun && jFile != _currentFile ) visible = false;
        if ( singleRun )
        {
            somePoints = getNumberPoints(jFile,0) > 200 && getNumberPoints(jFile,0) < 800;
            manyPoints = getNumberPoints(jFile,0) > 800;
        }
        for (int jCurve=0; jCurve<getNumberCurves(jFile); jCurve++ )
        {
            plotCurve *dataCurve = &(_listOfCurves[jFile][jCurve]);
            if ( dataCurve->pointSizeForData > 0 && somePoints )
                dataCurve->pointSizeForData = qMin(4,dataCurve->pointSizeForData);
            else if ( dataCurve->pointSizeForData > 0 && manyPoints )
                dataCurve->pointSizeForData = qMin(2,dataCurve->pointSizeForData);
            dataCurve->visible = visible;
            if ( dataCurve->scaleFactorX != 1. )
            {
                FUNC_INFO << "scale factor" << jFile << dataCurve->scaleFactorY << _yAxis2Ratio;
                for (int jt=0; jt<getNumberPoints(jFile,jCurve); jt++)
                    dataCurve->xData[jt] *= dataCurve->scaleFactorX;
            }
            FUNC_INFO << "file, curve" << jFile << jCurve << "scaleFactorY" << dataCurve->scaleFactorY;
            if ( dataCurve->scaleFactorY != 1. )
            {
                FUNC_INFO << "rescale";
                for (int jt=0; jt<getNumberPoints(jFile,jCurve); jt++)
                    dataCurve->yData[jt] *= dataCurve->scaleFactorY;
                for (int jt=0; jt<dataCurve->yError.size(); jt++)
                    dataCurve->yError[jt] *= dataCurve->scaleFactorY;
            }
        }
    }

    if ( _normalizeRunMode )   normalizeCurves();
    if ( _concatenateRunMode ) concatenateRuns();
    FUNC_EXIT;
}

void plotData::concatenateRuns()
{
    FUNC_ENTER << _concatenateRunMode;

    for (int jFile=0; jFile<getNumberFiles(); jFile++)
    {
        if ( jFile > 0 && _concatenateRunMode ) // for the 1st file, use the original xData (no offset due to concatenation)
        {
            plotCurve *dataCurve     = &(_listOfCurves[jFile][0]);
            plotCurve *lastDataCurve = &(_listOfCurves[jFile-1][0]);
            double offset;
            int nTime  = dataCurve->xData.size();
            if ( nTime == 1 )
                offset = 0.;
            else
            {
                double delta = lastDataCurve->xData[nTime-1] - lastDataCurve->xData[nTime-2];
                offset = lastDataCurve->xData.last() + delta;
            }
            for (int jCurve=0; jCurve<getNumberCurves(jFile); jCurve++)
            {
                for ( int jt=0; jt<getNumberPoints(jFile,jCurve); jt++ )
                    _listOfCurves[jFile][jCurve].xData[jt] += offset;
            }
        } //jFile >0
    } // jFile
    FUNC_EXIT;
}

void plotData::normalizeCurves()
{
    double yMax = 0.;
    for (int jFile=0; jFile<getNumberFiles(); jFile++)
    {
        plotCurve *dataCurve = &(_listOfCurves[jFile][0]);
        for (int jTime=0; jTime<dataCurve->yData.size(); jTime++)
        {
            double yData = dataCurve->yData.at(jTime);
            if ( yData > yMax) yMax = yData;
        }
    }

    if ( yMax != 0. )
    {
        for (int jFile=0; jFile<getNumberFiles(); jFile++)
        {
            for (int jCurve=0; jCurve<getNumberCurves(jFile); jCurve++)
            {
                for (int jTime=0; jTime<getNumberPoints(jFile,jCurve); jTime++)
                {
                    _listOfCurves[jFile][jCurve].yData[jTime] /= yMax;
                    _listOfCurves[jFile][jCurve].yError[jTime] /= yMax;
                }
            }
        } // jCurve
    } // jFile
}

void plotData::addCurve(int iFile, QString legend, dVector xData, dVector yData)
{
    dVector yError, yFit;
    bVector ignore;
    addCurve(iFile, legend, xData, yData, yError, yFit, ignore);
}

void plotData::addCurve(int iFile, QString legend, dVector xData, dVector yData, dVector yFit)
{
    dVector yError;
    bVector ignore;
    addCurve(iFile, legend, xData, yData, yError, yFit, ignore);
}

void plotData::addCurve(int iFile, QString legend, dVector xData, dVector yData, bVector ignore)
{
    dVector yError, yFit;
    addCurve(iFile, legend, xData, yData, yError, yFit, ignore);
}

void plotData::addCurve(int iFile, QString legend, dVector xData, dVector yData, dVector yFit, bVector ignore)
{
    dVector yError;
    addCurve(iFile, legend, xData, yData, yError, yFit, ignore);
}

void plotData::addCurve(int iFile, QString legend, dVector xData, dVector yData, dVector yError, dVector yFit)
{
    bVector ignore;
    addCurve(iFile, legend, xData, yData, yError, yFit, ignore);
}

void plotData::addCurve(int iFile, QString legend, dVector xData, dVector yData, dVector yError, dVector yFit, bVector ignore)
{ // xData and yData dimensions must match.  yError, yFit, and ignore dimensions can either match or be zero
    FUNC_ENTER << iFile << getNumberFiles();
    if ( iFile < 0 || iFile >= getNumberFiles() )
    {
        QString errorString = QString("Error: file index %1 is out of range (0 to %2) when adding curve").arg(iFile).arg(getNumberFiles()-1);
        qInfo() << errorString;
        qFatal("Programming terminated: fix code.");
    }
    else if ( yData.size() != xData.size() )
    {
        QString errorString = QString("Use vectors with consistent x and y sizes in plotData::addCurve(): %1, %2").
                              arg(xData.size()).arg(yData.size());
        qInfo() << errorString;
        qFatal("Programming terminated: fix code.");
    }
    else if ( yError.size() != xData.size() && yError.size() != 0 )
    {
        QString errorString = QString("Use vectors with consistent x and error sizes in plotData::addCurve(): %1, %2").
                              arg(xData.size()).arg(yError.size());
        qInfo() << errorString;
        qFatal("Programming terminated: fix code.");
    }
    else if ( yFit.size() != xData.size() && yFit.size() != 0 )
    {
        QString errorString = QString("Use vectors with consistent x and fit sizes in plotData::addCurve(): %1, %2").
                              arg(xData.size()).arg(yFit.size());
        qInfo() << errorString;
        qFatal("Programming terminated: fix code.");
    }
    else if ( ignore.size() != xData.size() && ignore.size() != 0 )
    {
        QString errorString = QString("Use vectors with consistent x and ignore sizes in plotData::addCurve(): %1, %2").
                              arg(xData.size()).arg(ignore.size());
        qInfo() << errorString;
        qFatal("Programming terminated: fix code.");
    }
    else
        _currentFileInit = iFile;
    plotCurve newCurve;
    newCurve.legend = legend;
    newCurve.visible = newCurve.enabled   = newCurve.exportCurve = true;
    newCurve.dashedForData = newCurve.dashedForFit = newCurve.histogram = newCurve.errorBars   = false;
    newCurve.scaleFactorX = newCurve.scaleFactorY = 1.;
    newCurve.pointStyle   = QCPScatterStyle::ssCircle;
    // default values for data
    if ( _listOfCurves[_currentFileInit].size() != 0 )
        newCurve.colorForData = Qt::gray;
    _listOfCurves[_currentFileInit].append(newCurve);
    _listOfCurves[_currentFileInit].last().xData = xData;
    _listOfCurves[_currentFileInit].last().yData = yData;
    _listOfCurves[_currentFileInit].last().yError = yError;
    _listOfCurves[_currentFileInit].last().yFit = yFit;
    _listOfCurves[_currentFileInit].last().ignore = ignore;

    FUNC_EXIT;
}

void plotData::setScaleFactorYRelative(int iCurveRef)
// set the scale factor for working curve [_currentFile].last() to the scale factor [_currentFile][iCurveRef]
{
    FUNC_ENTER;
    plotCurve currentCurve = _listOfCurves[_currentFile].last();
    plotCurve refCurve     = _listOfCurves[_currentFile][iCurveRef];

    double maxRef = 0.;
    for (int j=0; j< refCurve.yData.size(); j++)
        maxRef = qMax(maxRef,refCurve.yData[j]);
    double maxThis = 0.;
    for (int j=0; j< currentCurve.yData.size(); j++)
        maxThis = qMax(maxThis,currentCurve.yData[j]);

    double scaleFactorY = 1.;
    if ( maxThis != 0. ) scaleFactorY = maxRef / maxThis;
    setScaleFactorY(scaleFactorY);
}

void plotData::plotDataAndFit(bool newData)
{
    FUNC_ENTER;

    // clear the graph (the following can return the number cleared)
    _qcplot->clearGraphs();
    QCPRange Range1 = _qcplot->yAxis->range();
    _qcplot->yAxis2->setRange(Range1);
    _qcplot->yAxis2->setTickLabels(false);

    _listOfDataGraphs.resize(getNumberFiles());
    // plot curves in reverse order so that the first curves (e.g., jCurve=0) are not occluded by later curves
    for ( int jFile=0; jFile<getNumberFiles(); jFile++ )
    {
        _listOfDataGraphs[jFile].resize(getNumberCurves(jFile));
        for (int jCurve=getNumberCurves(jFile)-1; jCurve>=0; jCurve--)
        {
            plotCurve currentCurve = _listOfCurves[jFile][jCurve];
            plotCurve errorBarTargetCurve;
            if ( currentCurve.errorBars == 2 && jCurve+1 < _listOfCurves[jFile].size() )
            {
                FUNC_INFO << "curve and target" << jCurve << jCurve+1;
                errorBarTargetCurve = _listOfCurves[jFile][jCurve+1];
            }
            FUNC_INFO << "curve" << currentCurve.legend << jFile << jCurve;
            plotDataCurve(currentCurve, errorBarTargetCurve);
            _listOfDataGraphs[jFile][jCurve] = _qcplot->graph();

            if ( containsFit(currentCurve) )    plotFitCurve(currentCurve);
            if ( containsIgnore(currentCurve) ) plotIgnoreCurve(currentCurve);
        } //jcurve
    } // jfile

    // set the tracer position
    FUNC_INFO << _currentFile << _currentPoint;
    // concurrently, update the tracer to match
    QCPGraph *graph = _listOfDataGraphs[_currentFile][0];
    plotCurve curve = _listOfCurves[_currentFile][0];
    FUNC_INFO << "set Position to" << curve.xData[_currentPoint];
    _positionTracer->setGraph(graph);
    _positionTracer->setGraphKey(curve.xData[_currentPoint]);
    _positionTracer->updatePosition();

    FUNC_INFO << "autoScale and plot";
    if ( _autoScale && newData )
        autoScale(true);
    else
        _qcplot->replot();

    FUNC_EXIT;
}

void plotData::plotDataCurve(plotCurve curve, plotCurve errorBarTargetCurve)
{
    FUNC_ENTER;
    if ( curve.visible && curve.enabled )
    {
        // Fill color for error envelopes
        QColor fillColor;
        int hue; int saturation; int value;
        curve.colorForData.getHsv(&hue, &saturation, &value);
        fillColor.setHsv(hue,saturation/4,value);

        QPen myPen;
        // Pen for lines & points
        if ( curve.errorBars == 2 )
            myPen.setColor(fillColor);
        else
            myPen.setColor(curve.colorForData);
        myPen.setWidth(curve.lineThicknessForData);
        if ( curve.dashedForData )
            myPen.setStyle(Qt::DotLine);
        else
            myPen.setStyle(Qt::SolidLine);

        // Points
        QCPScatterStyle myScatter;
        myScatter.setShape(curve.pointStyle);
        myScatter.setPen(myPen);
        myScatter.setBrush(curve.colorForData);
        myScatter.setSize(curve.pointSizeForData);

        _qcplot->addGraph();
        // set the curve name and color
        _qcplot->graph()->setName(curve.legend);
        if (curve.legend.compare("none") == 0  )
        {
            int nLegendItems = _qcplot->legend->itemCount();
            _qcplot->legend->removeItem(nLegendItems-1);
        }
        // Points
        _qcplot->graph()->setPen(myPen);
        _qcplot->graph()->setScatterStyle(myScatter);
        // lines
        if ( curve.lineThicknessForData == 0 )
            _qcplot->graph()->setLineStyle(QCPGraph::lsNone);
        if ( curve.histogram )
            _qcplot->graph()->setLineStyle(QCPGraph::lsStepCenter);
        // Set the data
        FUNC_INFO << "set data" << curve.legend << curve.xData.size() << curve.yData.size();
//        FUNC_INFO << "set xData" << curve.xData;
//        FUNC_INFO << "set yData" << curve.yData;
         _qcplot->graph()->setData(curve.xData,curve.yData);
        if ( _qcplot->graph() == _listOfDataGraphs[_currentFile][0] )
             _positionTracer->setGraph(_qcplot->graph());  // set this when graphs are defined
        //                if ( curve.errorBars == 2 && jCurve != 0 )
        if ( curve.errorBars == 1 )
        {
            FUNC_INFO << "errorBars == 1";
            // error bars
            QCPErrorBars *errorBars = new QCPErrorBars(_qcplot->xAxis, _qcplot->yAxis);
            errorBars->removeFromLegend();
            errorBars->setAntialiased(false);
            // errorBars->setDataPlottable(_qcplot->graph());
            int iGraph = _qcplot->graphCount() - 1; // iGraph-1 is  current graph
            errorBars->setDataPlottable(_qcplot->graph(iGraph));
            errorBars->setPen(myPen);
            // Set the data
            errorBars->setData(curve.yError);
        }
        else if ( curve.errorBars == 2 )
        {
            FUNC_INFO << "errorBars == 2";
            if ( errorBarTargetCurve.errorBars == 2 ) // jCurve+1 due to reverse order plotting; otherwise jCurve-1
            {
                int iGraph = _qcplot->graphCount() - 2; // iGraph-1 is  current graph
                _qcplot->graph()->setBrush(QBrush(QBrush(fillColor)));
                _qcplot->graph()->setChannelFillGraph(_qcplot->graph(iGraph));
            }
        } // error bars
        if ( curve.scaleFactorX != 1. )
        {
            FUNC_INFO << "scaleFactorX";
            _qcplot->xAxis2->setVisible(true);
            _qcplot->xAxis2->setTickLabels(true);
            _qcplot->xAxis2->setTickLabelColor(curve.colorForData);
            _qcplot->xAxis2->setTickLabelFont(QFont("Helvetica",16));
            QCPRange Range1 = _qcplot->xAxis->range();
            FUNC_INFO << "Range1" << Range1;
            setXAxis2Range(Range1);
            //                _qcplot->xAxis2->setLabel(curve.legend);
            _qcplot->xAxis2->setLabelFont(QFont("Helvetica", 20));
            _qcplot->xAxis2->setLabelColor((curve.colorForData));
        }
        if ( curve.scaleFactorY != 1. )
        {
            FUNC_INFO << "scaleFactorY" << curve.scaleFactorY;
            _qcplot->yAxis2->setVisible(true);
            _qcplot->yAxis2->setTickLabels(true);
            _qcplot->yAxis2->setTickLabelColor(curve.colorForData);
            _qcplot->yAxis2->setTickLabelFont(QFont("Helvetica",16));
            QCPRange range1 = _qcplot->yAxis->range();
            setYAxis2Range(range1);
            FUNC_INFO << "range" << range1.lower << range1.upper;
            //                _qcplot->yAxis2->setLabel(curve.legend);
            _qcplot->yAxis2->setLabelFont(QFont("Helvetica", 20));
            _qcplot->yAxis2->setLabelColor((curve.colorForData));
        }
    } // visible
    FUNC_EXIT;
}

void plotData::plotFitCurve(plotCurve curve)
{
    FUNC_ENTER;
    if ( curve.visible && curve.enabled )
    {
        _qcplot->addGraph();
        // set the curve name and color
        _qcplot->graph()->setName("fit");
        // Points
        _qcplot->graph()->setScatterStyle(QCPScatterStyle::ssNone);

        // lines
        QPen myPen;
        // Pen for lines & points
        if ( curve.colorForFit == Qt::red ) FUNC_INFO << "red for" << curve.legend;
        myPen.setColor(curve.colorForFit);
        myPen.setWidth(curve.lineThicknessForFit);
        if ( curve.dashedForFit )
            myPen.setStyle(Qt::DotLine);
        else
            myPen.setStyle(Qt::SolidLine);
        _qcplot->graph()->setPen(myPen);
        if ( curve.lineThicknessForFit == 0 ) _qcplot->graph()->setLineStyle(QCPGraph::lsNone);
        FUNC_INFO << "thickness" << curve.lineThicknessForFit;

        // Set the data
        _qcplot->graph()->setData(curve.xData,curve.yFit);
    } // visible
    FUNC_EXIT;
}

void plotData::plotIgnoreCurve(plotCurve curve)
{
    FUNC_ENTER;
    if ( curve.visible && curve.enabled )
    {
        _qcplot->addGraph();
        // set the curve name and color
        _qcplot->graph()->setName("ignore");

        dVector xData, yData;
        for (int j=0; j<curve.xData.size(); j++)
        {
            if ( curve.ignore.at(j) )
            {
                xData.append(curve.xData.at(j));
                yData.append(curve.yData.at(j));
            }
        }

        // lines
        QPen myPen;
        // Pen for lines & points
        if ( curve.colorForFit == Qt::red ) FUNC_INFO << "red for ignored" << curve.legend;
        myPen.setColor(curve.colorForFit);
        myPen.setWidth(curve.lineThicknessForData);
        // Points
        QCPScatterStyle myScatter;
        myScatter.setShape(curve.pointStyle);
        myScatter.setPen(myPen);
        myScatter.setBrush(curve.colorForFit);
        myScatter.setSize(2*curve.pointSizeForData);
        myScatter.setShape(QCPScatterStyle::ssStar);
        _qcplot->graph()->setLineStyle(QCPGraph::lsNone);
        _qcplot->graph()->setScatterStyle(myScatter);
        // Set the data
        _qcplot->graph()->setData(xData,yData);
    } // visible
    FUNC_EXIT;
}

void plotData::autoScale(bool state)
{
    FUNC_INFO << "plotData::autoScale" << state;
    _autoScale = state;
    if ( _autoScale )
    {
        QCPRange xRange, yRange;
        reScaleAxes(&xRange, &yRange);
        FUNC_INFO << "plotData::autoScale yRange" << yRange.lower << yRange.upper;
        _qcplot->yAxis->setRange(yRange);
        setXAxis2Range(xRange);
        setYAxis2Range(yRange);
        if ( _autoScaleXRange.lower == _autoScaleXRange.upper )
            _qcplot->xAxis->setRange(xRange);
        else
            _qcplot->xAxis->setRange(_autoScaleXRange);
        _qcplot->replot();
        setSelectPoints(_pressToMoveTracer);
    }
}

void plotData::reScaleAxes(QCPRange *xRange, QCPRange *yRange)
{ // don't use the QCustomPlot function rescaleAxes (_qcplot->rescaleAxes()) so we can flag ignored points in the re-scaling
    xRange->lower=1.e10;  xRange->upper=-1.e10;
    yRange->lower=1.e10;  yRange->upper=-1.e10;
    FUNC_ENTER << "Data::reScaleAxes enter";

    for ( int jFile=0; jFile<getNumberFiles(); jFile++ )
    {
        for (int jCurve=0; jCurve<getNumberCurves(jFile); jCurve++ )
        {
            plotCurve currentCurve = _listOfCurves[jFile][jCurve];
            if ( currentCurve.visible && currentCurve.enabled )
            {
                int nPoints = currentCurve.yData.size();
                for ( int jPoint=0; jPoint<nPoints; jPoint++ )
                {
                    bool errorBars = currentCurve.errorBars == 1; // for errorBars=2 (envelop, data are passed as two line curves)
                    bool includePoint = currentCurve.ignore.size() == 0;
                    if ( currentCurve.ignore.size() > 0 ) includePoint = !currentCurve.ignore[jPoint];
                    if ( includePoint )
                    {   // Do not ignore this point
                        //                    FUNC_INFO << "Data::reScaleAxes 4";
                        double yMin, yMax;
                        if ( errorBars )
                        {
                            yMin = currentCurve.yData[jPoint] - currentCurve.yError[jPoint];
                            yMax = currentCurve.yData[jPoint] + currentCurve.yError[jPoint];
                        }
                        else
                            yMin = yMax = currentCurve.yData[jPoint];
                        //                    FUNC_INFO << "graph" << jFile << "point" << jPoint << "yMin, yMax" << yMin << yMax;
                        if ( yMin < yRange->lower ) yRange->lower = yMin;
                        if ( yMax > yRange->upper ) yRange->upper = yMax;
                    }
                    // Do not ignore the x value of "ignored" points.
                    if ( currentCurve.xData[jPoint] < xRange->lower ) xRange->lower = currentCurve.xData[jPoint];
                    if ( currentCurve.xData[jPoint] > xRange->upper ) xRange->upper = currentCurve.xData[jPoint];
                }
            }
        } // jCurve
    } // jFile
    if ( xRange->lower ==  1.e10 ) xRange->lower = -10.;
    if ( xRange->upper == -1.e10 ) xRange->upper =  10.;
    // expland the x and y ranges by 10% total
    double extra = (xRange->upper - xRange->lower) * 0.05;
    xRange->lower -= extra;
    xRange->upper += extra;

    if ( yRange->lower ==  1.e10 ) yRange->lower = -10.;
    if ( yRange->upper == -1.e10 ) yRange->upper =  10.;
    extra = (yRange->upper - yRange->lower) * 0.05;
    yRange->lower -= extra;
    yRange->upper += extra;
    FUNC_EXIT << "Data::reScaleAxes exit";
}
