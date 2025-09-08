#ifndef GENERALGLM_H
#define GENERALGLM_H

#include "io.h"

enum fMRIEventShapes // used both in fMRI and timePage classes
{
    Shape_none,
    Shape_separator1,

    Shape_Gamma,
    Shape_Sigmoid,
    Shape_separator2,

    Shape_Square,
    Shape_RampUp,
    Shape_RampDown,
    Shape_FIR,
    Shape_FR,
    Shape_separator3,

    Shape_ASL,
    Shape_Sine,
    Shape_Cosine,
    Shape_separator4,

    Shape_Table,
    Shape_Nuisance,
    Shape_separator5,

    Shape_Baseline0,
    Shape_Baseline1,
    Shape_Baseline2,
    Shape_Baseline3,
    Shape_Baseline4,
    Shape_Baseline5,
    Shape_Baseline6,
    Shape_Baseline7,
    Shape_Baseline8,
    Shape_Baseline9,
    Shape_Baseline10,
    Shape_separator6,

    Shape_Interaction
};
enum PETWeightingModels
{
    Weights_Uniform,
    Weights_11C_Noiseless,
    Weights_11C,
    Weights_18F_Noiseless,
    Weights_18F
};
enum PETEventTypes
{
    Type_R1,
    Type_k2,
    Type_k2a,
    Type_plasma,
    Type_challenge
};
class GeneralGLM
{
#define MaxConditions 100

private:
    int _nTime=0;
    int _currentCondition=0;
    iVector _nContrasts;     // [nConditions]; = 1 for T test or >= 1 for F test
    QVector<dMatrix> _fCovar;          // [nConditions][nContrasts][nContrasts]
    QVector<dMatrix> _fCovarInv;       // [nConditions][nContrasts][nContrasts]

    bool _initialized=false; // has the GLM been initialized successfully in constructor?
    bool _prepared=false;    // has the pseudo-inverse has been created?

    bool _basisFunctionsChanged=false;  // use this in a set/get fashion to determine when basis functions have changed
    dMatrix _XTWXm1_cc;      // used for variance for OLS
    int _nIncluded;          // sum of non-zero weights
    dVector _weight_t;
    dVector _data_t;
    dVector _fit_t;
    dVector _fitErr_t;
    dVector _beta_c;         // coefficient attached to a basis function
    dVector _sem_c;          // SEM of coefficient
    double _dof=0.;          // DOF for fit
    double _dofEffective=0.; // potentially override _dof by using non-zero _dofEffective (e.g, smoothing to boost dof alla Worsley 2002)
    double _average=0.;      // average across time course

protected:
    int _nCoeff=0;
    dMatrix _X_tc;           // Design matrix
    dMatrix _XTWXm1XTW_ct;   // used for Beta = (Xt*W*X)^-1 * Xt * W for either OLS or WLS (weights=1)
    double _sigma2=0.;       // SOS for fit; for 2nd-order GLM, sigma2 is incorporated into weights (=1/sigma2), so set _sigma2=1;
                             // otherwise, determine _sigma2 post-hoc from data-fit, so that 1/weights = _sigma2 * normalized weights
    double _resStDev=0.;     // residual standard deviation: sqrt(res)

public:
    // conditions
    QStringList _conditionList; // nConditions = _conditionList.size()
    iVector _conditionShape;    // [nConditions]: matches _basisShape but with condition index (not coeff index)

    // vector of length _nCoeff (= # events)
    cVector _basisID;           // [_nCoeff]; points to eventID for the respective basis function
    iVector _basisShape;        // [_nCoeff]; R1, k2, k2a, challenge, or shape for fMRI
    QVector<dMatrix> _contrastMatrix;  // [nConditions][_nContrasts][_nCoeff]
    double _sinusoidPeriod=10.;

    double _FIRresolution=0.;
    int _FIRlength=0;  // set
    iVector _FIRContrast;

    iVector _nuisanceLag;           // [_nCoeff]
    dVector _nuisanceContrast;      // [_nCoeff]

    // conditions
    double _fStatistic;      // F statitics for condition
    double _cnr;             // CNR for condition
    double _pRaw;            // raw p-value for condition
    double _effectSize;      // magnitude for condition
    double _effectStDev;       // SEM for condition
    iMatrix _indexCoeffInCondition; // [nConditions][# coeff in condition]; keep track of which coefficients are a condition

    // functions
    void init(int nTime, int nCoeff);  // call this first
    void setWeights(dVector weight_t);
    void setVariances(dVector variance_t);
    void setOLS();
    void replaceBasisFunction( int iBasis, dVector X_t);   // define/replace existing basis function
    void addBasisFunction( dVector X_t);
    void addOrInsertBasisFunction( int iBasis, dVector X_t);
    void fitWLS(dVector data) { fitWLS(data, true);}  // a convenience function for most situations (error not known a-priori)
    void fitWLS(dVector data, bool computeSigma2);
    void calculatePseudoInverse();
    void removeCoefficient(int iCoeff, dVector &data);
    void removeCoefficient(int iCoeff, dVector &data, dVector &fit);
    void orthogonalize(int iCoeff, dVector &data);
    void orthogonalizeExistingBasisFunction(int iCoeff, int iCoeffMod);
    void calculateCovarianceMatrix();
    void calculateFStat (int lCondition);
    void defineConditions(QString conditionString);
    void decodeConditions(bool defineMatrices );
    void evaluateCurrentCondition();

    // setters
    void setNumberContrasts(int iCondition, int nContrasts );
    inline void setPrepared(bool state) {_prepared = state;}
    inline void setBasisFunctionsChanged(bool state) {_basisFunctionsChanged=state;}  // allows recoloring of basis functions
    inline void setCurrentCondition(int lCondition) {_currentCondition = qMax(0,qMin(lCondition,_conditionList.size()-1));}
    inline void setDOFEffective(double dof) {_dofEffective=dof;}
    inline void setBasisID(int iBasis, QChar ID) {_basisID[iBasis] = ID;}

    // getters
    inline bool getBasisFunctionsChanged() {return _basisFunctionsChanged;}
    QString getConditionString();
    inline QChar getBasisID(int iBasis) {return _basisID[iBasis];}
    inline int getNumberTimePoints()    {return _nTime;}
    inline double getSigma2()           {return _sigma2;}
    inline double getResidualSD()       {return _resStDev;}
    inline double getAIC()              {double var = _sigma2*_dof/_weight_t.size();
                                         return 2*_nCoeff + _nTime * qLn(var);}
    inline double getAverage()          {return _average;}
    inline double getBeta(int iCoeff)   {return _beta_c[iCoeff];}
    inline double getSEM(int iCoeff)    {return _sem_c[iCoeff];}
    inline double getVar(int iCoeff)    {return _sem_c[iCoeff] * _sem_c[iCoeff];}
    inline double getCovar(int iC1, int iC2) {return _sigma2 * _XTWXm1_cc[iC1][iC2];}
    inline double getDOF()              {return _dof;}
    inline dVector getDataAll() {return _data_t;}
    inline dVector getFitAll() {return _fit_t;}
    inline double getData(int iTime)    {return _data_t[iTime];}
    inline int getFitSize() {return _fit_t.size();}
    inline double getFit(int iTime)     {return _fit_t[iTime];}
    inline double getFitErr(int iTime)  {return _fitErr_t[iTime];}
    inline double getRegressor(int iCoeff, int iTime) {return _beta_c[iCoeff] * _X_tc[iTime][iCoeff];}
    inline bool isReady()               {return _prepared;}
    inline int getNumberCoefficients()  {return _nCoeff;}
    inline int getCurrentCondition()    {if ( _currentCondition< 0 || _currentCondition>=_conditionList.size() ) _currentCondition=0; return _currentCondition;}
    inline double getEffectSize()       {return _effectSize;}  // call "setCurrentCondition" and "evaluateCurrentCondition" prior to this
    inline double getEffectStDev()      {return _effectStDev;} // ...
    inline double getpRaw()             {return _pRaw;}        // ...
    inline double getCNR()              {return _cnr;}         // ...
    inline double getBasisPoint(int iBasis, int iTime) {return _X_tc[iTime][iBasis];}
    inline int getNumberConditions()                {return _conditionList.size();}
    inline int getNumberContrasts(int iCondition)   {return _nContrasts[iCondition];}
    inline QString getConditionName(int lCondition) {return _conditionList.at(lCondition);}
    inline int getCurrentConditionShape() {return _conditionShape[_currentCondition];}
    inline iVector getFIRContrast() {return _FIRContrast;}
    inline int getNuisancelag(int iCoeff)  {return _nuisanceLag[iCoeff];}
    inline bool isValidID(QChar ID) {return getEventIndex(ID) >= 0;}

    QChar getEventChar(int indexEvent);
    int getEventIndex(QChar eventID);
    int getNumberCoefficientsInCondition(int iCondition);
    int getOnlyCoefficientInCondition(int iCondition);
    int getFirstCoefficientInCondition(int iCondition);
    int getEventCoefficient(QChar eventID);
    int getEventCoefficient(QChar eventID, int iOccurrence);
    iVector getEventCoefficients(QChar eventID);
    double getWeight(int iTime);
    inline dVector getWeights() {return _weight_t;}
};

class SimplePolynomial : public GeneralGLM
{ // basis functions = {1,x,x^2};  use this for cases to fit a simple polynomial without normalizating to the range x=(-1,1)
private:
    int _nPoly;  // usually this is the same as _nCoeff in generalGLM, but it could be different if non-poly basis functions get tacked on
    dVector _xInputVector;
    double baselineBasisFunction(int iPoly, double x);

public:
    void define (int nCoeff, int nTime);
    void define (int nCoeff, dVector xVector);
};


////////////////////////////////// Polynomial GLM //////////////////////////////////
class PolynomialGLM : public GeneralGLM
{ // basis functions up to order 6; input vectors get normalized to range x=(-1,1), so retrieve fits from this class
private:
    int _nPoly;  // usually this is the same as _nCoeff in generalGLM, but it could be different if non-poly basis functions get tacked on
    dVector _xInputVector;
    int _centralPoint=-1;  // for value >= 0, this is LOESS, and the central point is the point of interest
//    inline double baselineBasisFunction(int iPoly, int x) {return utilMath::polynomialLegendre(iPoly, static_cast<double>(x));}
    double getBasisDerivative(int iPoly, double x);
    inline double getNormalizedTime(int iTime) {return getNormalizedTime(_xInputVector[iTime]);}
    double getNormalizedTime(double time);
    // most general
    void define(int nPoly, int nCoeff, dVector xVector, int iCenter);
public:
    // polynomial-only entries
    inline void define (int nCoeff, int nTime)          {define(nCoeff,nCoeff,nTime);}
    inline void define (int nCoeff, dVector xVector)    {define(nCoeff,nCoeff,xVector);}
    // allow defining functions that provide for additional coefficients to be added after polynomial terms
    void define(int nPoly, int nCoeff, int nTime);
    void define(int nPoly, int nCoeff, dVector xVector);   // for non-evenly spaced points
    // provide for asymmetric LOESS
    inline void define(int nCoeff, dVector xVector, int iCenter) {define(nCoeff,nCoeff,xVector,iCenter);}
    double getDerivative(int iTime);
    double getDerivative(double time);
    double getFitInterpolation(double time); // provids a mechanism for getting non-indexed (continuous) fit values
    inline double getFitAtCentralPoint() {return getFit(_centralPoint);}
};

////////////////////////////////// LOESS Polynomial vector //////////////////////////////////
class LOESS
{
private:
    QVector<PolynomialGLM> _polyLOESS;  // vector of [nTime] polynomials, usually of order 3
    dVector _xPoints;   // x points used in to produce fit
    dVector _yPoints;   // y points used in to produce fit

public:
    void defineAndFit(dVector yVector, double smoothingScale, bool asymmetric, bool linkWidthToTime);
    void defineAndFit(dVector yVector, double smoothingScale, bool asymmetric, bool linkWidthToTime, int onset);
    void defineAndFit(dVector timeVector, dVector yVector, double smoothingScale, bool asymmetric, bool linkWidthToTime);
    void defineAndFit(dVector timeVector, dVector yVector, double smoothingScale, bool asymmetric, bool linkWidthToTime, int onset);
    void defineAndFitPolar(dVector thetaVector, dVector radiusVector, double smoothingScale);
    dVector fit(dVector data);
    double getFitAtCentralPoint(int iTime) {return _polyLOESS[iTime].getFitAtCentralPoint();}
    double getFirstDerivative(int iTime);
    double getSecondDerivative(int iTime);
    double getFitInterpolation(double xBetween);
    inline int getNumberPoints() {return _xPoints.size();}
    inline double getX(int iX) {return _xPoints[iX];}
    inline double operator() (int iX) {return getFitAtCentralPoint(iX);}
};

////////////////////////////////// Polynomial GLM //////////////////////////////////
class GaussianPlusPolynomial : public GeneralGLM
{
    // Refernece:
    // A New Algorithm for Fitting a Gaussian Function Riding on the Polynomial Background
    // Roonizi, IEEE SIGNAL PROCESSING LETTERS, VOL. 20, NO. 11, NOVEMBER 2013
private:
    // data (input)
    dVector _histogram;
    dVector _fit;
    double _xOffset;
    double _xStep;
    // fit: polynomial plus gaussian
    // f(x) = alpha_k * x^k + mag * exp{ (x-mu)^2 / 2/sigma2 }
    dVector _poly;
    double _mag;
    double _mu;
    double _sigma2;

    double baselineBasisFunction(int iPoly, double x);
    dVector getPhi (int order);

public:
        void fit(int nPoly, double xOffset, double xStep, dVector histogram);

        // getters
        double getGaussianMu();
        double getGaussianSigma2();
        double getPolynomialBeta(int iCoeff);
        inline double getGaussianMag() {return _mag;};
        inline dVector getFit() {return _fit;};
};

#endif // GENERALGLM_H
