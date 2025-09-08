#include <QDebug>
#include <QRegularExpression>
#include "generalglm.h"

QString GeneralGLM::getConditionString()
{
    QString string="";
    for (int jString=0; jString<_conditionList.size(); jString++)
    {
        string.append(_conditionList.at(jString));
        string.append(" ");
    }
    return string;
}

double GeneralGLM::getWeight(int iTime)
{
    if ( _nTime == 0 || iTime >= _nTime )
        return 1.;
    else
        return _weight_t[iTime];
};


void GeneralGLM::calculateCovarianceMatrix()
{ // call this after conditions have been defined
    for (int jCondition=0; jCondition<_conditionList.size(); jCondition++)
    {
        // Define the F convariance matrix for each condition.
        _fCovar[jCondition].resize(_nContrasts[jCondition]);
        _fCovarInv[jCondition].resize(_nContrasts[jCondition]);
        for (int jc1=0; jc1<_nContrasts[jCondition]; jc1++)
        {
            _fCovar[jCondition][jc1].fill(0.,_nContrasts[jCondition]);
            _fCovarInv[jCondition][jc1].fill(0.,_nContrasts[jCondition]);
            for (int jc2=0; jc2<_nContrasts[jCondition]; jc2++)
            {
                _fCovar[jCondition][jc1][jc2] = 0.;
                for (int jCoeff=0; jCoeff<_nCoeff; jCoeff++)
                {
                    for (int jCoeff1=0; jCoeff1<_nCoeff; jCoeff1++)
                        _fCovar[jCondition][jc1][jc2]        += _contrastMatrix[jCondition][jc1][jCoeff]
                                * _XTWXm1_cc[jCoeff][jCoeff1] * _contrastMatrix[jCondition][jc2][jCoeff1];
                }
            }
        }
        // Make a copy for the inverse
        bool allZero=true;
        for (int jc1=0; jc1<_nContrasts[jCondition]; jc1++)
        {
            for (int jc2=0; jc2<_nContrasts[jCondition]; jc2++)
            {
                _fCovarInv[jCondition][jc1][jc2] = _fCovar[jCondition][jc1][jc2];
                allZero &= _fCovar[jCondition][jc1][jc2] == 0.;
            }
        }
        if ( allZero ) return; // don't bother trying to invert a zero matrix
        // Now invert the matrix.
        if ( ! utilMath::invertSquareMatrix(_fCovarInv[jCondition]) ) return;
    }
}

int GeneralGLM::getNumberCoefficientsInCondition(int iCondition)
{ // find the number of coefficients in the condition AND assign a vector matching them to the basis coefficient list
    int nCoefficientsInCondition = 0;
    // _nContrasts[iCondition] should be 1 contrast for a T test, but more for an F test
    for ( int jContrast=0; jContrast<_nContrasts[iCondition]; jContrast++ )
    {
        for ( int jCoeff=0; jCoeff<_nCoeff; jCoeff++ )
        {
            if ( _contrastMatrix[iCondition][jContrast][jCoeff] != 0. )
                nCoefficientsInCondition++;
        }
    }
    return nCoefficientsInCondition;
}

int GeneralGLM::getOnlyCoefficientInCondition(int iCondition)
{
    int iCoeff = -1;
    if ( getNumberCoefficientsInCondition(iCondition) == 1 )
    {
        // should be 1 contrast for a T test, but more for an F test
        for ( int jContrast=0; jContrast<_nContrasts[iCondition]; jContrast++ )
        {
            for ( int jCoeff=0; jCoeff<_nCoeff; jCoeff++ )
            {
                if ( _contrastMatrix[iCondition][jContrast][jCoeff] != 0. )
                    iCoeff = jCoeff;
            }
        }
    }
    return iCoeff;
}

int GeneralGLM::getFirstCoefficientInCondition(int iCondition)
{
    int iCoeff = -1;
    // should be 1 contrast for a T test, but more for an F test
    for ( int jContrast=0; jContrast<_nContrasts[iCondition] || iCoeff >=0; jContrast++ )
    {
        for ( int jCoeff=0; jCoeff<_nCoeff; jCoeff++ )
        {
            if ( _contrastMatrix[iCondition][jContrast][jCoeff] != 0. )
            {
                iCoeff = jCoeff;
                break;
            }
        }
    }
    return iCoeff;
}

void GeneralGLM::setNumberContrasts(int iCondition, int nContrasts )
{
//    FUNC_ENTER << iCondition << nContrasts;
    // This function needs to be called if either the condition changes or the # ncoefficients changes (call through defineConditions)
    if ( iCondition >= _conditionList.size() )
    {
        decodeConditions(false);  // this excludes bad conditions in the string and reset the string to valid conditions
        decodeConditions(true);
    }
    _nContrasts[iCondition] = nContrasts;
    _contrastMatrix[iCondition].resize(_nContrasts[iCondition]);
    for (int jContrast=0; jContrast<_nContrasts[iCondition]; jContrast++)
        _contrastMatrix[iCondition][jContrast].fill(0.,_nCoeff);
//    FUNC_EXIT;
}

void GeneralGLM::evaluateCurrentCondition()
{
    FUNC_ENTER;
    int currentCondition = getCurrentCondition();
    FUNC_INFO << currentCondition << getNumberConditions();
    if ( getNumberConditions() != 0 && currentCondition >= 0 )
        calculateFStat(currentCondition);
}

void GeneralGLM::calculateFStat (int lCondition)
{
    FUNC_ENTER;
    if ( _sigma2 == 0. )
    {
        _cnr = _effectSize = _effectStDev = 0.;
        _pRaw = 1.;
        return;
    }

    // Make sure the condition is within range.
    if ( lCondition < 0 || lCondition >= _conditionList.size() )
    {
        decodeConditions(false);  // this excludes bad conditions in the string and reset the string to valid conditions
        decodeConditions(true);
    }

    dVector CBeta;

    CBeta.fill(0.,_nContrasts[lCondition]);
    int lContrastMax = 0;  double CBetaMax  = 0.;
    // Calculate the product CBeta
    for (int jc1=0; jc1<_nContrasts[lCondition]; jc1++)
    {
        for (int jCoeff=0; jCoeff<_nCoeff; jCoeff++)
            CBeta[jc1] += _contrastMatrix[lCondition][jc1][jCoeff] * _beta_c[jCoeff];
        if ( qAbs(CBeta[jc1]) > qAbs(CBetaMax) )
        {
            CBetaMax     = CBeta[jc1];
            lContrastMax = jc1;
        }
    }
    // Determine the effect size and standard deviation.
    _effectSize  = CBetaMax;
    _effectStDev = qSqrt(_fCovar[lCondition][lContrastMax][lContrastMax] * _sigma2);

    // Compute the numerator for F
    double numerator = 0.;
    for (int jc1=0; jc1<_nContrasts[lCondition]; jc1++)
    {
        for (int jc2=0; jc2<_nContrasts[lCondition]; jc2++)
            numerator += CBeta[jc1] * _fCovarInv[lCondition][jc1][jc2] * CBeta[jc2];
    }

    // Degrees of freedom
    double dof1 = _dof;
    double dof2 = _nContrasts[lCondition];
    if ( dof2 < 0. ) qInfo() << "dof2 =" << dof2 << lCondition;
    // potentially override _dof by using an effective DOF due to regularization/smoothing
    if ( _dofEffective != 0. ) dof1 = _dofEffective;
    if ( dof1 == 0 )
    {   // This could happen rarely if no points are included (e.g., less included points than coefficients)
        // Failure to account for this chance will lead to a fatal error in "ibeta()" below (dof1=0)
        _fStatistic = _cnr = 0.;
        _pRaw = 1.;
        return;
    }
    ///////////////////////////////////////////////
    // The covariance matrix, _fCovarInv, is formed from (_XTWXm1_cc)^-1
    // As such, _sigma2 can be included 2 ways:
    // 1) as long as the weights are normalized, _sigma2 can be determined post-hoc from data vs. fit
    // 2) in some cases (e.g., 2nd order GLM), weights are known explicitly (1/sigma2), so set _sigma2=1 in formula below
    // Note that in any case,  "weights/_sigma2" should = 1/sigma2_actual, so weight normalization is NOT arbitrary
    ///////////////////////////////////////////////
    _fStatistic = numerator / _sigma2 / dof2;
    _cnr        = qSqrt(_fStatistic);
    if ( CBetaMax < 0. ) _cnr *= -1.;

    double x = dof1/(dof1+dof2*_fStatistic);
    if (x < 0. || x > 1.0) x = 1.;
    double beta = utilMath::ibeta(0.5*dof1,0.5*dof2,x); // cephes
    if ( dof2 == 1. )
        _pRaw = beta;
    else
    {
        _pRaw = 2.*beta;
        if ( _pRaw == 2. ) _pRaw = 1.;
        if ( _pRaw > 1.0 ) _pRaw = 2. - _pRaw;
    }

    double tooSmall = exp(-100.);
    if ( _pRaw < tooSmall ) _pRaw = tooSmall;

    FUNC_EXIT;
    return;
}

void GeneralGLM::removeCoefficient(int iCoeff, dVector &data)
{
    removeCoefficient(iCoeff, data, _fit_t);
    _data_t = data;
}

void GeneralGLM::removeCoefficient(int iCoeff, dVector &data, dVector &fit)
{
    if ( _nTime != data.size() )
    {
        qWarning() << "Error: the # time points (" << data.size() << ") does not match the GLM (" << _nTime << ").";
        exit(1);
    }
    for (int jt=0; jt<_nTime; jt++)
    {
        // remove coefficients
        fit[jt] -= _X_tc[jt][iCoeff] * _beta_c[iCoeff];
        data[jt]-= _X_tc[jt][iCoeff] * _beta_c[iCoeff];
    }
}

void GeneralGLM::orthogonalize(int iCoeff, dVector &data)
{ // orthogonalize data wrt basis[iCoeff]
//    return;
    double dotProduct = 0.;  double mag2 = 0.;
    for (int jt=0; jt<_nTime; jt++)
    {
        dotProduct += data[jt] * _X_tc[jt][iCoeff];
        mag2       += _X_tc[jt][iCoeff] * _X_tc[jt][iCoeff];
    }
    // Subtract the portion of the input vector that is parallel to the given coefficient vector
    for (int jt=0; jt<_nTime; jt++)
        data[jt] -= dotProduct / mag2 * _X_tc[jt][iCoeff];
}

void GeneralGLM::orthogonalizeExistingBasisFunction(int iCoeff, int iCoeffMod)
{ // orthogonalize basis[iCoeffMod] wrt coefficient iCoeff
//    return;
    if ( iCoeff == iCoeffMod ) return;  // do nothing if they are the same coefficient
    double dotProduct = 0.;  double mag2 = 0.;
    for (int jt=0; jt<_nTime; jt++)
    {
        dotProduct += _X_tc[jt][iCoeffMod] * _X_tc[jt][iCoeff];
        mag2       += _X_tc[jt][iCoeff]    * _X_tc[jt][iCoeff];
    }
    // Subtract the portion of the input vector that is parallel to the given coefficient vector
    for (int jt=0; jt<_nTime; jt++)
        _X_tc[jt][iCoeffMod] -= dotProduct / mag2 * _X_tc[jt][iCoeff];
}

void GeneralGLM::setOLS()
{
    if ( _nTime == 0 )
        qFatal("Error: GeneralGLM not initialized (nTime=0) in function setOLS");
    _weight_t.fill(1.,_nTime);
}


void GeneralGLM::init(int nTime , int nCoeff)
{ // this should be called first to initialize the GLM
    if ( nTime != _nTime )
    {
        _nTime  = _nIncluded = nTime;
        _data_t.fill(0.,_nTime);
        _fit_t.fill(0.,_nTime);
        FUNC_INFO << "resize1 _fit_t to" << _fit_t.size();
        _fitErr_t.fill(0.,_nTime);
        _X_tc.resize(_nTime);
    }

    _contrastMatrix.resize(MaxConditions);
    _fCovar.resize(MaxConditions);
    _fCovarInv.resize(MaxConditions);

    if ( nCoeff != _nCoeff )
    {
        _nCoeff = nCoeff;
        _XTWXm1_cc.resize(_nCoeff);
        _XTWXm1XTW_ct.resize(_nCoeff);
        for (int jc=0; jc<_nCoeff; jc++)
        {
            _XTWXm1_cc[jc].fill(0.,_nCoeff);
            _XTWXm1XTW_ct[jc].fill(0.,_nTime);
        }
    }

    _dof = _nIncluded - _nCoeff;   // this should be changed if time points are ignored
    for ( int jt=0; jt<_nTime; jt++ )
        _X_tc[jt].resize(_nCoeff);

    _initialized = true;
}

void GeneralGLM::setVariances( dVector variance_t)
{
    int nTime = variance_t.size();
    dVector weight_t; weight_t.resize(nTime);
    for (int jt=0; jt<nTime; jt++)
    {
        if ( variance_t[jt] > 0. )
            weight_t[jt] = 1. / variance_t[jt];
        else
            weight_t[jt] = 1.;
    }
    /*
    // Normalize weights
    double average=0.;
    for (int jt=0; jt<nTime; jt++)
        average += weight_t[jt];
    average /= static_cast<double>(nTime);
    for (int jt=0; jt<nTime; jt++)
        weight_t[jt] /= average;
        */
    setWeights(weight_t);
}

void GeneralGLM::setWeights( dVector weight_t)
{
//    FUNC_ENTER << "reset weights" << weight_t;
    if ( weight_t.size() != _nTime )
    {
        qWarning() << "Error: weight size " << weight_t.size() << " not = time size " << _nTime;
        exit(1);
    }
    _nIncluded = _nTime;
    _weight_t.resize(_nTime);
    for ( int jt=0; jt<_nTime; jt++ )
    {
        _weight_t[jt] = weight_t[jt];
        if ( _weight_t[jt] == 0. ) _nIncluded--;
    }
    _dof = _nIncluded - _nCoeff;
    _prepared = false;
//    FUNC_EXIT << "_weight_t" << _weight_t;
}

void GeneralGLM::replaceBasisFunction( int iBasis, dVector X_t)
{
    if ( X_t.size() != _nTime )
    {
        qWarning() << "Error: time dimension " << X_t.size() << " not equal to GLM time dimension " << _nTime;
        exit(1);
    }
    else if ( iBasis < 0 || iBasis >= _nCoeff )
    {
        qWarning() << "Error: basis function index " << iBasis << " out of range in replaceBasisFunction";
        exit(1);
    }
    for ( int jt=0; jt<_nTime; jt++ )
    {
        _X_tc[jt][iBasis] = X_t[jt];
//        if ( _weight_t[jt] == 0. ) _X_tc[jt][iBasis] = 0.;
    }
    _basisFunctionsChanged=true;
    _prepared = false;
}

void GeneralGLM::addBasisFunction( dVector X_t)
{
    if ( X_t.size() != _nTime )
    {
        qWarning() << "Error: time dimension " << X_t.size() << " not equal to GLM time dimension " << _nTime;
        exit(1);
    }
    else
    {
        int nTime = X_t.size();
        init(nTime,_nCoeff+1);
        for (int jt=0; jt<X_t.size(); jt++)
        {
            _X_tc[jt][_nCoeff-1] = X_t[jt];
//            if ( _weight_t[jt] == 0. ) _X_tc[jt][_nCoeff-1] = 0.;
        }
    }
    _basisFunctionsChanged=true;
    _prepared = false;
}

void GeneralGLM::addOrInsertBasisFunction( int iBasis, dVector X_t)
{
    if ( X_t.size() != _nTime )
    {
        qWarning() << "Error: time dimension " << X_t.size() << " not equal to GLM time dimension " << _nTime;
        exit(1);
    }
    else if ( iBasis < 0 || iBasis > _nCoeff )
    {
        qWarning() << "Error: basis function index " << iBasis << " out of range in addBasisFunction";
        exit(1);
    }
    else if ( iBasis == _nCoeff )
    { // add basis function
        int nTime = X_t.size();
        init(nTime,_nCoeff+1);
    }
    for (int jt=0; jt<X_t.size(); jt++)
        _X_tc[jt][iBasis] = X_t[jt];
    _basisFunctionsChanged=true;
    _prepared = false;
    FUNC_EXIT;
}
void GeneralGLM::fitWLS(dVector data, bool computeSigma2)
{ // computeSigma2: if true, find sigma2 from residue; else assume sigma2 is KNOWN and used in weighting matrix
    // If any sigma values are zero, ignore the voxel.
    _data_t = data;
//    FUNC_INFO << "data" << _data_t;
    bool allZeroes=true;
    for (int jt=0; jt<_nTime; jt++)
        allZeroes &= (data[jt] == 0.);
    if ( allZeroes )
    {
        FUNC_INFO << "ncoeff" << _nCoeff << "ntime" << _nTime;
        // Set everything to zero
        _average = _sigma2 = 0.;
        _beta_c.fill(0.,_nCoeff);   _sem_c.fill(0.,_nCoeff);
        _fit_t.fill(0.,_nTime);     _fitErr_t.fill(0.,_nTime);
        FUNC_INFO << "resize2 _fit_t to" << _fit_t.size();
        return;
    }

    if ( !_prepared )
        calculatePseudoInverse();
    FUNC_INFO << 2;
    // Calculate Beta_c = (X_transpose * W * X)^-1 * X_transpose * W * Y
    for (int jc=0; jc<_nCoeff; jc++)
    {
        _beta_c[jc] = 0.;
        // Calculate Beta, summing over time
        for (int jt=0; jt<_nTime; jt++)
            _beta_c[jc] += _XTWXm1XTW_ct[jc][jt] * data[jt];
    }
    FUNC_INFO << 3;
    // Compute the fit: Y_t = X_tc * beta_c
    FUNC_INFO << "nTime" << _nTime;
    for (int jt=0; jt<_nTime; jt++)
    {
        _fit_t[jt] = _fitErr_t[jt] = 0.;
        for (int jc=0; jc<_nCoeff; jc++)
            _fit_t[jt] += _X_tc[jt][jc] * _beta_c[jc];
    }

    FUNC_INFO << 4;
    if ( computeSigma2 )
    {
        // Compute sigma
        _sigma2 = _average = 0;
        for (int jt=0; jt<_nTime; jt++)
        {
            if ( _weight_t[jt] != 0. )
            {
                _sigma2 += SQR(data[jt]-_fit_t[jt]);
                _average += data[jt];
            }
        }
        _resStDev = qSqrt(_sigma2 / static_cast<double>(_nTime-1));
        _sigma2  /= _dof;
        _average /= _dof;
        if ( _nTime == _nCoeff )
        {
            _sigma2 = 0.;
            _sem_c.fill(0.,_nCoeff);
        }
        else
        {
            // Compute the SEM for each coefficient
            for (int jc=0; jc<_nCoeff; jc++)
                _sem_c[jc] = sqrt(_sigma2*_XTWXm1_cc[jc][jc]);
        }
    }
    else
    {
        // Set _sigma2=1. This variable is only used in post-hoc determination of variance.
        // In 2nd-order GLM, _sigma2 is known/set explicitly in the weighting matrix.
        _sigma2 = 1.;
        if ( _nTime == _nCoeff )
            _sem_c.fill(0.,_nCoeff);
        else
        {
            // Compute the SEM for each coefficient
            for (int jc=0; jc<_nCoeff; jc++)
                _sem_c[jc] = sqrt(_XTWXm1_cc[jc][jc]);
        }
    }
    FUNC_INFO << 5;
    // Compute the fit error
    _fitErr_t.fill(0.,_nTime);
    FUNC_INFO << 6;
    for (int jt=0; jt<_nTime; jt++)
    {
        for (int jc=0; jc<_nCoeff; jc++)
            _fitErr_t[jt] += _X_tc[jt][jc] * _sem_c[jc];
    }
    FUNC_EXIT;
}

void GeneralGLM::defineConditions(QString conditionString )
{ // must be called AFTER GLM_define_basis_functions()
//    FUNC_ENTER;
    //  DefineContrastMatrix:
    //    - Generally will be true,
    //    - except when one wants to interpret the conditionString to create a list of files without knowing about basis functions
    //  Set up the conditional tests based upon "condition_string".
    //  This string can contain more than 1 condition.
    //  Example: "1 2 3 4 1234 12-34 12-33 1-3,2-4
    //  yields 7 conditions, including 7 T tests:
    //    1-4) event 1 versus baseline, 2 versus baseline, ...
    //    5) event 1+2+3+4 versus baseline
    //    6) event 1+2-3-4 versus baseline
    //    7) event 1+2-3-3 versus baseline (e.g., [1 1 -2 ..] in contrast vector)
    //  ... and 1 F test
    //    1) events 1-3 or 2-4

    _conditionList = conditionString.split(QRegularExpression("[\\s]"),Qt::SkipEmptyParts); // clazy:exclude=use-static-qregularexpression
    decodeConditions(false);  // this excludes bad conditions in the string and reset the string to valid conditions
    decodeConditions(true);
//    FUNC_EXIT;
}

void GeneralGLM::decodeConditions( bool defineMatrices )
{ // generally call this function 2ce: 1) to exclude bad conditions with argument false, then define conditions with argument true
//    FUNC_ENTER << _conditionList;
    if ( _conditionList.isEmpty() )
        return;
    else if ( defineMatrices )
    {
        _nContrasts.fill(1,_conditionList.size());
        _conditionShape.resize(_conditionList.size());
        _indexCoeffInCondition.resize(_conditionList.size());
    }

    QStringList validConditionList;
    int nConditions = _conditionList.size();
    for (int jCondition=0; jCondition<nConditions; jCondition++)
    {
        QString condition = _conditionList.at(jCondition);
        // allocate contrast matrix and set to 0
        if ( defineMatrices )
        {
            _indexCoeffInCondition[jCondition].resize(0);
            setNumberContrasts(jCondition,condition.count(',') + 1); // # of contrasts = # commas + 1 in the condition string
        }
        // Go through the condition string and decode it
        bool allGoodCoefficientsInCondition = true;
        int iContrast=0;
        int sign=1;
        bool anyBPType = false;
        for ( int jChar=0; jChar<condition.size(); jChar++)
        {
            QChar eventID = condition.at(jChar);
            if ( eventID == ',' )
            {
                iContrast++;
                sign = 1;
            }
            else if ( eventID == '-')
                sign = -1;
            else
            {
                iVector iCoeffVector = getEventCoefficients(eventID);
                FUNC_INFO << "eventID" << eventID << "iCoeffVector" << iCoeffVector;
                for (int jCoeffInVector=0; jCoeffInVector<iCoeffVector.size(); jCoeffInVector++)
                {
                    int iCoeff = iCoeffVector[jCoeffInVector];
                    allGoodCoefficientsInCondition &= (iCoeff >= 0);
                    if ( iCoeff >= 0 && defineMatrices )
                    {
                        _indexCoeffInCondition[jCondition].append(iCoeff);
                        _conditionShape[jCondition] = _basisShape[iCoeff];
                        anyBPType = anyBPType || (_basisShape[iCoeff] == Type_k2a);
                        if ( _conditionShape[jCondition] == Shape_FIR || _conditionShape[jCondition] == Shape_FR )
                        {
                            if ( _FIRContrast[jCoeffInVector] == 1 )
                                _contrastMatrix[jCondition][iContrast][iCoeff] += sign;
                        }
                        else if ( _conditionShape[jCondition] == Shape_Nuisance )
                        {
                            if ( _nuisanceContrast[iCoeff] != 0. )
                                _contrastMatrix[jCondition][iContrast][iCoeff] += sign * _nuisanceContrast[iCoeff];
                        }
                        else
                            _contrastMatrix[jCondition][iContrast][iCoeff] += sign;
                    }
                } // jCoeffInVector
                // BP takes precedence over dBP (challenge), so that BP+challenge yield BP, not dBP and occupancy.
                if ( anyBPType ) _conditionShape[jCondition] = Type_k2a;
            } // apply sign
        } // jChar
        if ( allGoodCoefficientsInCondition && ! validConditionList.contains(condition) ) validConditionList.append(condition);
    } // jCondition
    // reset the condition list
    _conditionList = validConditionList;

    // Go back through the contrast matrix and make the sum of all positives = 1 within a contrast (and same for the negatives)
    // Thus, the contrast ab-c = "(a+b)-c") becomes "(a+b)/2 - c"
    if ( defineMatrices )
    {
        for ( int jCondition=0; jCondition<nConditions; jCondition++)
        {
            if ( _conditionShape[jCondition] != Shape_Nuisance )
            {
                for ( int jContrast=0; jContrast<_nContrasts[jCondition]; jContrast++ )
                {
//                    double sumPos=0;  double sumNeg=0;
                    iVector shapesPosUnique, shapesNegUnique, shapesPosAll, shapesNegAll;
                    for ( int jCoeff=0; jCoeff<_nCoeff; jCoeff++)
                    {
                        double cValue = _contrastMatrix[jCondition][jContrast][jCoeff];
                        if ( cValue > 0. ) shapesPosAll.append(_basisShape[jCoeff]);
                        if ( cValue < 0. ) shapesNegAll.append(_basisShape[jCoeff]);
                        if ( cValue > 0. && !shapesPosUnique.contains(_basisShape[jCoeff]) ) shapesPosUnique.append(_basisShape[jCoeff]);
                        if ( cValue < 0. && !shapesNegUnique.contains(_basisShape[jCoeff]) ) shapesNegUnique.append(_basisShape[jCoeff]);
                    }
                    for ( int jType =0; jType<shapesPosUnique.count(); jType++)
                    {
                        int iType = shapesPosUnique.at(jType);
                        int nSum  = shapesPosAll.count(iType);
                        for ( int jCoeff=0; jCoeff<_nCoeff; jCoeff++)
                        {
                            double cValue = _contrastMatrix[jCondition][jContrast][jCoeff];
                            if ( cValue > 0. && _basisShape[jCoeff] == iType )
                                _contrastMatrix[jCondition][jContrast][jCoeff] /= static_cast<double>(nSum);
                        }
                    }
                    for ( int jType =0; jType<shapesNegUnique.count(); jType++)
                    {
                        int iType = shapesNegUnique.at(jType);
                        int nSum  = shapesNegAll.count(iType);
                        for ( int jCoeff=0; jCoeff<_nCoeff; jCoeff++)
                        {
                            double cValue = _contrastMatrix[jCondition][jContrast][jCoeff];
                            if ( cValue < 0. && _basisShape[jCoeff] == iType )
                                _contrastMatrix[jCondition][jContrast][jCoeff] /= static_cast<double>(nSum);
                        }
                    }
//                    qInfo() << "contrastMatrix[" << jCondition << "] =" << _contrastMatrix[jCondition][jContrast];
                } // jContrast
            } // ! Shape_Nuisance
        } // jContrast
    } // jCondition

    // The conditions have been defined, so calculate the covariance matrix.
    if ( defineMatrices ) calculateCovarianceMatrix();
//    FUNC_EXIT << _conditionList;
}

QChar GeneralGLM::getEventChar(int indexEvent)
{
    // indexEvent   char
    //  0-8         1-9
    //  9-34        a-z
    // 35-60        A-Z

    int iAscii;
    if ( indexEvent < 0)
        iAscii = 32;               // should only occur for erroneous entry
    else if ( indexEvent < 9 )
        iAscii = indexEvent + 49;  // event index 0 = '1' = ascii 49
    else if ( indexEvent < 35)
        iAscii = indexEvent + 88;  // event index 9 = 'a' = ascii 97
    else
        iAscii = indexEvent + 30;  // event index 35 = 'A' = ascii 65
    QChar qAscii = static_cast<QChar>(iAscii);
    return qAscii;
}

int GeneralGLM::getEventIndex(QChar eventID)
{
    int iAscii = eventID.toLatin1();
    int indexEvent;
    if ( iAscii >= 49 && iAscii <= 57 )
        indexEvent = iAscii - 49;  // '1' is 49 and goes into index 0, 9->8
    else if ( iAscii >= 97 && iAscii <= 122 )
        indexEvent = iAscii - 88;  // 'a' is 97 and goes into index 9
    else if ( iAscii >= 65 && iAscii <= 90 )
        indexEvent = iAscii - 30;  // 'A' is 65 and goes into index 35
    else
        indexEvent=-1; // invalid ID
    return indexEvent;
}

int GeneralGLM::getEventCoefficient(QChar eventID)
{
    for (int jCoeff=0; jCoeff<_basisID.size(); jCoeff++)
    {
        if ( _basisID[jCoeff] == eventID )
            return(jCoeff);
    }
    return -1;
}
int GeneralGLM::getEventCoefficient(QChar eventID, int iOccurrence)
{
    iVector eventCoefficients = getEventCoefficients(eventID);
    if ( iOccurrence >= eventCoefficients.size() )
        return eventCoefficients[0];
    else
        return eventCoefficients[iOccurrence];
}
iVector GeneralGLM::getEventCoefficients(QChar eventID )
{
    int nCoeffPerEvent=0;
    for (int jCoeff=0; jCoeff<_basisID.size(); jCoeff++)
    {
        if ( _basisID[jCoeff] == eventID ) nCoeffPerEvent++;
        FUNC_INFO << "eventID" << eventID << "basisID" << _basisID[jCoeff];
    }
    iVector eventCoefficients;
    if ( nCoeffPerEvent == 0 )
    {
        eventCoefficients.resize(1);
        eventCoefficients[0]=-1;
    }
    else
    {
        eventCoefficients.resize(nCoeffPerEvent);
        nCoeffPerEvent=0;
        for (int jCoeff=0; jCoeff<_basisID.size(); jCoeff++)
        {
            if ( _basisID[jCoeff] == eventID )
            {
                eventCoefficients[nCoeffPerEvent] = jCoeff;
                nCoeffPerEvent++;
            }
        }
    }
    return eventCoefficients;
}

void GeneralGLM::calculatePseudoInverse()
{
    FUNC_ENTER;
    // Zero all data
    _beta_c.fill(0.,_nCoeff);
    _sem_c.fill(0.,_nCoeff);
    _XTWXm1_cc.resize(_nCoeff);
    _XTWXm1XTW_ct.resize(_nCoeff);
    for (int jc=0; jc<_nCoeff; jc++)
    {
        _XTWXm1_cc[jc].fill(0.,_nCoeff);
        _XTWXm1XTW_ct[jc].fill(0.,_nTime);
    }
    if ( _weight_t.size() == 0 ) _weight_t.fill(1.,_nTime);
    FUNC_INFO << "here" << 1;
    FUNC_INFO << "nCoeff" << _nCoeff << _nTime;

    // Calculate (X_transpose * W * X) by summing over the time index.
    for (int jc=0; jc<_nCoeff; jc++)
    {
        for (int jc1=0; jc1<_nCoeff; jc1++)
        {
            FUNC_INFO << "_XTWXm1_cc size" << _XTWXm1_cc.size() << "weight size" << _weight_t.size();
            for (int jt=0; jt<_nTime; jt++)
                _XTWXm1_cc[jc][jc1] += _X_tc[jt][jc] * _weight_t[jt] * _X_tc[jt][jc1];
        }
    }
    FUNC_INFO << "_XTWXm1_cc" << _XTWXm1_cc;

    // Invert to get (X_transpose * W * X)^-1
    if ( ! utilMath::invertSquareMatrix(_XTWXm1_cc) ) return;

    FUNC_INFO << "here" << 2;
    // Calculate (X_transpose * W * X)^-1 * X_transpose
    for (int jc=0; jc<_nCoeff; jc++)
    {
        for (int jt=0; jt<_nTime; jt++)
        {
            for (int jc1=0; jc1<_nCoeff; jc1++)
                _XTWXm1XTW_ct[jc][jt] += _XTWXm1_cc[jc][jc1] * _X_tc[jt][jc1]* _weight_t[jt];
        }
    }
    _prepared = true;
    FUNC_EXIT;
}

void SimplePolynomial::define(int nCoeff, int nTime)
{
    FUNC_ENTER << nCoeff << nTime;
    dVector xVector;  xVector.resize(nTime);
    double duration = qMax(1.,static_cast<double>(nTime-1));
    // convert to a vector on the range -1,1
    for ( int jt=0; jt<nTime; jt++ )
        xVector[jt] = 2.*static_cast<double>(jt)/duration - 1.;  // range=(-1,1) with center point at x=0 or x=1/(n-1)
    define(nCoeff, xVector);
}
void SimplePolynomial::define(int nCoeff, dVector xVector)
{ // This defines a simple polynomial GLM for fitting a single scan
    FUNC_ENTER << nCoeff << xVector;
    int nTime = xVector.size();
    if ( nTime < nCoeff )
      {
        qInfo() << "Error: the # of points (" << nTime << ") must be >= # coefficients" << nCoeff;
        return;
      }
    else if ( nTime == nCoeff )
        nCoeff = nTime;
    nCoeff = qMin(nCoeff,4);
    // The following sets a flag; failure above will leave the flag unset
    FUNC_INFO << "init" << nTime << nCoeff;
    init(nTime,nCoeff);

    // Create the basis functions.
    for (int jPoly=0; jPoly<nCoeff; jPoly++)
    {
        dVector X_t;
        for ( int jt=0; jt<nTime; jt++ )
            X_t.append( baselineBasisFunction(jPoly,xVector[jt]) );
        addOrInsertBasisFunction(jPoly,X_t);
    }
    // set uniform weights as default
    dVector weights;  weights.fill(1.,nTime);
    setWeights(weights);
    calculatePseudoInverse();
    FUNC_INFO << "get fit size0" << getFitSize();
}
double SimplePolynomial::baselineBasisFunction(int iPoly, double x)
{
    double value;
    if ( iPoly == 0 )
        value = 1.;
    else if ( iPoly == 1 )
        value = x;
    else
        value = qPow(x,static_cast<double>(iPoly));
    return value;
}

void PolynomialGLM::define(int nPoly, int nCoeff, int nTime)
{
//    FUNC_ENTER << nPoly << nCoeff << nTime;
    dVector xVector;  xVector.resize(nTime);
    double duration = qMax(1.,static_cast<double>(nTime-1));
    // convert to a vector on the range -1,1
    for ( int jt=0; jt<nTime; jt++ )
        xVector[jt] = 2.*static_cast<double>(jt)/duration - 1.;  // range=(-1,1) with center point at x=0 or x=1/(n-1)
    define(nPoly, nCoeff, xVector);
}

void PolynomialGLM::define(int nPoly, int nCoeff, dVector xVector)
{
    int iCenter = -1;
    define(nPoly, nCoeff, xVector, iCenter);
}

void PolynomialGLM::define(int nPoly, int nCoeff, dVector xVector, int iCenter)
{ // This defines a simple polynomial GLM for fitting a single scan
    _nPoly = nPoly;  // usually this is the same as _nCoeff in generalGLM, but it could differ if additional basis functions are added
    _xInputVector = xVector;
    _centralPoint = iCenter;
    int nTime = _xInputVector.size();
    if ( nTime < nCoeff )
      {
        qInfo() << "Error: the # of points (" << nTime << ") must be >= # coefficients" << nCoeff;
        return;
      }
    else if ( nTime < nCoeff )
        nCoeff = nTime;
    nCoeff = qMin(nCoeff,6);
    // The following sets a flag; failure above will leave the flag unset
    init(nTime,nCoeff);

    // Create the basis functions.
    dVector X_t;
    X_t.resize(nTime);
    for (int jPoly=0; jPoly<nCoeff; jPoly++)
    {
        for ( int jt=0; jt<nTime; jt++ )
            X_t[jt] = utilMath::polynomialLegendre(jPoly,getNormalizedTime(jt));
        addOrInsertBasisFunction(jPoly,X_t);
    }
    // set uniform weights as default
    dVector weights;  weights.fill(1.,nTime);
    setWeights(weights);
    calculatePseudoInverse();
}
double PolynomialGLM::getNormalizedTime(double time)
{ // convert time (via the index into _xInputVector) to the range (-1,1)
    double minusOneToOne = 0.;
    if ( _centralPoint < 0 )
    {
        double duration = _xInputVector[_xInputVector.size()-1] - _xInputVector[0];  // duration = max - min
        if ( duration != 0. )
        {
            double scale = 2./duration;  double offset = 2.*_xInputVector[0]/duration - 1.;
            minusOneToOne = scale * time + offset;
        }
    }
    else
    { // LOESS; generally #pts should be odd, so there is a central point
        double x0 = _xInputVector[_centralPoint]; double xLow = _xInputVector[0]; double xHigh = _xInputVector[_xInputVector.size()-1];
        double maxHalfWidth = qMax(qAbs(x0-xLow),qAbs(xHigh-x0));
        if ( maxHalfWidth == 0.)
            minusOneToOne = 0.;
        else
            minusOneToOne = (time-x0)/maxHalfWidth;
    }
    return minusOneToOne;
}
double PolynomialGLM::getBasisDerivative(int iPoly, double x)
{ // return Legendre polynomial of x with order iPoly
    double value;
    if ( iPoly == 0 )
        value = 0.;
    else if ( iPoly == 1 )
        value = 1.;
    else if ( iPoly == 2 )
        value = 3*x;
    else if ( iPoly == 3 )
        value = 1.5 * (5.*x*x - 1.);
    else if ( iPoly == 4 )
        value = 2.5 * (7.*x*x*x - 3.*x);
    else // if ( iPoly == 5 )
        value = 1.875 * (21.*x*x*x*x - 14.*x*x + 1.);
    return value;
}
double PolynomialGLM::getFitInterpolation(double time)
{ // if _nCoeff > _nPoly, the extra terms will need to be computed outside of this class based upon the extra terms parametric definition
    double x = getNormalizedTime(time);
//    FUNC_ENTER << _nPoly << _nCoeff << time << x;
    double fitValue=0.;
    for (int jPoly=0; jPoly<_nPoly; jPoly++)
        fitValue += utilMath::polynomialLegendre(jPoly, x) * getBeta(jPoly);
//    FUNC_EXIT << fitValue;
    return fitValue;
}

double PolynomialGLM::getDerivative(int iTime)
{ // dBasis/dt = dBasis/dt' * dt'/dt
    double derivative = 0.;
    for (int jPoly=0; jPoly<_nCoeff; jPoly++)
    {
        double polyDeriv = getBasisDerivative(jPoly,getNormalizedTime(iTime));
        derivative += polyDeriv;
    }
    double duration = _xInputVector[_xInputVector.size()-1] - _xInputVector[0];  // duration = max - min
    return derivative * 2./duration;  // dBasis/dt = dBasis/dt' * dt'/dt
}

double PolynomialGLM::getDerivative(double time)
{ // dBasis/dt = dBasis/dt' * dt'/dt
    double derivative = 0.;
    for (int jPoly=0; jPoly<_nCoeff; jPoly++)
    {
        double polyDeriv = getBasisDerivative(jPoly,getNormalizedTime(time));
        derivative += polyDeriv;
    }
    double duration = _xInputVector[_xInputVector.size()-1] - _xInputVector[0];  // duration = max - min
    return derivative * 2./duration;  // dBasis/dt = dBasis/dt' * dt'/dt
}

void LOESS::defineAndFit(dVector yVector, double smoothingScale, bool asymmetric, bool linkWidthToTime, int onset)
{
    dVector xVector;
    for (int j=0; j<yVector.size(); j++)
        xVector.append(j+0.5);
    defineAndFit(xVector, yVector, smoothingScale, asymmetric, linkWidthToTime, onset);
}

void LOESS::defineAndFit(dVector yVector, double smoothingScale, bool asymmetric, bool linkWidthToTime)
{
    dVector xVector;
    for (int j=0; j<yVector.size(); j++)
        xVector.append(j+0.5);
    defineAndFit(xVector, yVector, smoothingScale, asymmetric, linkWidthToTime);
}

void LOESS::defineAndFit(dVector timeVector, dVector yVector, double smoothingScale, bool asymmetric, bool linkWidthToTime)
{
    int onset = 0;
    defineAndFit(timeVector, yVector, smoothingScale, asymmetric, linkWidthToTime, onset);
}

void LOESS::defineAndFit(dVector timeVector, dVector yVector, double smoothingScale, bool asymmetric, bool linkWidthToTime, int onset)
{
    // timeVector, yVector = x and y
    // 0 < smoothingScale <= 0.5, where 0.5 means use half the data on either side
    // asymmetric: allow asymmetric smoothing windows to account for boundaries; otherwise, force smaller windows near boundaries
    // linkWidthToTime: if true, the smoothing scale increases with time
//    FUNC_ENTER << timeVector.size() << yVector.size() << smoothingScale << asymmetric << linkWidthToTime;
    _xPoints = timeVector;
    _yPoints = yVector;
    double smoothingScaleNonzero = qMax(0.01, smoothingScale);
    double smoothingHalfWidth = smoothingScaleNonzero;  // smoothingWidth is the HALF Width
    // Define local quadratic smoothing for each this run: (one function for each time point)
    int nTime = timeVector.size();
    _polyLOESS.resize(nTime);
    for ( int jt=0; jt<nTime; jt++)
    {
//        FUNC_INFO << "jt" << jt;
        double time0 = timeVector[jt];
        if ( linkWidthToTime )
            smoothingHalfWidth = qMax(0.01,(time0-onset) * smoothingScaleNonzero);    // link smoothing width to time, e.g., smoothingScale=0.5
        else
            smoothingHalfWidth = (timeVector.last() - timeVector.first()) * smoothingScaleNonzero;
//        FUNC_INFO << "timeVector" << timeVector;
//        FUNC_INFO << "calc" << timeVector.last() << timeVector.first() << smoothingScaleNonzero;
//        FUNC_INFO << "smoothingHalfWidth" << smoothingHalfWidth;
        int iLowSide=0;  double width=0.;
        while ( jt - iLowSide >= 0 && width < smoothingHalfWidth )
        {
            double time = timeVector[jt-iLowSide];
            width = qAbs(time-time0);
            iLowSide++;
        }
        iLowSide--;
        int iHighSide=0;  width = 0.;
        while ( jt + iHighSide < nTime && width < smoothingHalfWidth )
        {
            double time = timeVector[jt+iHighSide];
            width = qAbs(time-time0);
            iHighSide++;
        }
        iHighSide--;
        // Set weights
        dVector weight;  weight.clear();
        int iHalfWidth = qMax(1,qMin(iLowSide,iHighSide));
        dVector localTime, localYData;
        int iLocal=0;
        int iCenter=0;
        if ( !asymmetric ) iLowSide = iHighSide = iHalfWidth;
        for ( int jRel=-iLowSide; jRel<=iHighSide; jRel++)
        {
            int iTime = jt + jRel;
            if ( iTime >= 0 && iTime < nTime )
            {
                double time = timeVector[iTime];
                double distance = qAbs(time-time0) / smoothingHalfWidth;
                if ( distance < 0.99 )
                {
                    localTime.append(time);
                    localYData.append(yVector[iTime]);
                    weight.append(qPow(1-qPow(distance,3),3));  // use traditional tri-cube weight function
                    if ( jRel == 0 ) iCenter = iLocal;
                    iLocal++;
                }
            }
        }
//        FUNC_INFO << localTime;
//        FUNC_INFO << localYData;
//        FUNC_INFO << weight;
        //             define       (nCoeff,          nTime)
        _polyLOESS[jt].define(qMin(weight.size(),3),localTime,iCenter);  // qMax: don't define 3 terms with only 2 points
        _polyLOESS[jt].setWeights(weight);
        _polyLOESS[jt].fitWLS(localYData,false);
//        FUNC_EXIT << "poly[" << jt << "] = " << _polyLOESS[jt].getNumberTimePoints();
    }
}

void LOESS::defineAndFitPolar(dVector thetaVector, dVector radiusVector, double smoothingScale)
{ // assume range (0, 2*PI)
//    FUNC_ENTER << thetaVector;
    _xPoints = thetaVector;
    _yPoints = radiusVector;
    // Define local quadratic smoothing for each this run: (one function for each theta point)
    int nTheta = thetaVector.size();
    _polyLOESS.resize(nTheta);

    dVector expandedThetaVector  = thetaVector;
    dVector expandedRadiusVector = radiusVector;
    // append high side using low values
    int addHighSide=0;  double theta0 = 2*PI;  double width=0.;
    double smoothingHalfWidth = smoothingScale * PI; // 0.5 means ± 90 degrees, or half the data
    while ( width < smoothingHalfWidth )
    {
        double theta  = thetaVector[addHighSide] + 2.*PI;
        double radius = radiusVector[addHighSide];
        width = qAbs(theta-theta0);
        expandedThetaVector.append(theta);
        expandedRadiusVector.append(radius);
        addHighSide++;
    }
    // append low side using high values
    int addLowSide=0;  theta0 = 0;  width=0.;
    while ( width < smoothingHalfWidth )
    {
        double theta  = thetaVector[nTheta-1-addLowSide] - 2.*PI;
        double radius = radiusVector[nTheta-1-addLowSide];
        width = qAbs(theta-theta0);
        expandedThetaVector.prepend(theta);
        expandedRadiusVector.prepend(radius);
        addLowSide++;
    }

    for ( int jTheta=0; jTheta<nTheta; jTheta++)
    {
        int jThetaShift = jTheta + addLowSide;
        double theta0 = expandedThetaVector[jThetaShift];
        int iLowSide=0;  double theta = theta0;  double width=0.;
        while ( jThetaShift - iLowSide >= 0 && width < smoothingHalfWidth )
        {
            theta = expandedThetaVector[jThetaShift-iLowSide];
            width = qAbs(jThetaShift-theta0);
            iLowSide++;
        }
        iLowSide--;     iLowSide = qMax(1,iLowSide);
        int iHighSide=0;  theta = theta0;  width = 0.;
        while ( jThetaShift + iHighSide < expandedThetaVector.size() && width < smoothingHalfWidth )
        {
            theta = expandedThetaVector[jThetaShift+iHighSide];
            width = qAbs(theta-theta0);
            iHighSide++;
        }
        iHighSide--;    iHighSide = qMax(1,iHighSide);
        // Set weights
        dVector weight;  weight.clear();
        dVector localtheta, localYData;
        int iLocal=0;
        int iCenter=0;
        for ( int jRel=-iLowSide; jRel<=iHighSide; jRel++)
        {
            int iTheta = jThetaShift + jRel;
            if ( iTheta >= 0 && iTheta < expandedThetaVector.size() )
            {
                theta = expandedThetaVector[iTheta];
                double distance = qAbs(theta-theta0) / smoothingHalfWidth;
                if ( distance < 0.99 )
                {
                    localtheta.append(theta);
                    localYData.append(expandedRadiusVector[iTheta]);
                    weight.append(qPow(1-qPow(distance,3),3));  // use traditional tri-cube weight function
                    if ( jRel == 0 ) iCenter = iLocal;
                    iLocal++;
                }
            }
        }
        //             define       (nCoeff,          nTheta)
        _polyLOESS[jTheta].define(qMin(weight.size(),3),localtheta,iCenter);  // qMax: don't define 3 terms with only 2 points
        _polyLOESS[jTheta].setWeights(weight);
        _polyLOESS[jTheta].fitWLS(localYData,false);
//        FUNC_EXIT << "poly[" << jTheta << "] = " << _polyLOESS[jTheta].getNumberTimePoints();
    }
}

double LOESS::getFitInterpolation(double xBetween)
{
    double minDistance=1.e10;  int iClosest = -1;
    for ( int jX=0; jX<_polyLOESS.size(); jX++)
    {
        double x = _xPoints[jX];
        double distance = qAbs(x-xBetween);
        if ( distance < minDistance)
        {
            minDistance = distance;
            iClosest = jX;
        }
    }
    return _polyLOESS[iClosest].getFitInterpolation(xBetween);
}

dVector LOESS::fit(dVector data)
{
    dVector fit = data;
    for (int jt=0; jt<data.size(); jt++)
    {
        int nLocal  = _polyLOESS[jt].getNumberTimePoints();
        int iHalfWidth = nLocal/2;
        dVector localData;  localData.resize(nLocal);
        int jLocal=0;
        for ( int jRel=-iHalfWidth; jRel<=iHalfWidth; jRel++, jLocal++ )
        {
            int iTime = jt + jRel;
            if ( iTime >= 0 && iTime < data.size() )
                localData[jLocal] = data[iTime];
            else
                localData[jLocal] = 0.;  // should also have weight=0.
        }
        _polyLOESS[jt].fitWLS(localData,true);
        fit[jt] = getFitAtCentralPoint(jt);
    }
    return fit;
}

double LOESS::getFirstDerivative(int iTime)
{ // presume define() and fit() preceeded this
//    FUNC_ENTER << iTime;
    double deriv = 0.;
    if ( iTime > 0 && iTime < _polyLOESS.size()-1 )
        deriv = (getFitAtCentralPoint(iTime+1) - getFitAtCentralPoint(iTime-1)) / (_xPoints[iTime+1] - _xPoints[iTime-1]);
    else if ( iTime == 0 )
        deriv = (getFitAtCentralPoint(1) - getFitAtCentralPoint(0)) / (_xPoints[1] - _xPoints[0]);
    else // iTime == _polyLOESS.size()-1
    {
        int it = _polyLOESS.size()-1;
        deriv= (getFitAtCentralPoint(it) - getFitAtCentralPoint(it-1)) / (_xPoints[it] - _xPoints[it-1]);
    }
    return deriv;
}

double LOESS::getSecondDerivative(int iTime)
{ // presume define() and fit() preceeded this
//    FUNC_ENTER << iTime;
    double deriv = 0.;
    if ( iTime > 0 && iTime < _polyLOESS.size()-1 )
    {
        double x1 = _xPoints[iTime-1];  double y1 = getFitAtCentralPoint(iTime-1);
        double x2 = _xPoints[iTime];    double y2 = getFitAtCentralPoint(iTime);
        double x3 = _xPoints[iTime+1];  double y3 = getFitAtCentralPoint(iTime+1);
        double x21 = x2 - x1;           double x31 = x3 - x1;         double x32 = x3 - x2;
        deriv = 2.*y1/x21/x31 - 2.*y2/x32/x21 + 2.*y3/x32/x31;
    }
    else if ( iTime == 0 )
        deriv = getSecondDerivative(1);
    else
        deriv = getSecondDerivative(_polyLOESS.size()-2);
    return deriv;
}

void GaussianPlusPolynomial::fit (int nPoly, double xOffset, double xStep, dVector histogram)
{
    // The function to fit is
    //     f(x) = alpha_k * x^k + mag * exp{ (x-mu)^2 / 2/sigma2 }
    // BUT we are fitting a modified function that linearize the equation:
    //
    //     y_i = X_ij * beta_j
    //           with (i,j) both have same dimensions (coefficients),
    //           X_ij = sum_over_bins(phi_i * phi_j)
    //           y_i  = sum_over_bins(phi_i * f=data)
    //
    //       phi_i(x) = x^i  for 0<i<=N+1
    //     phi_N+2(x) = sum_over_bins_0_x{u*f(u)*du), where f(u) is taken from data
    //     phi_N+3(x) = sum_over_bins_0_x{f(u)*du)
    int nBins = histogram.size();
    int nCoeff = nPoly + 4;  // polynomial order + 1 + A, mu, sigma2
    // The following sets a flag; failure above will leave the flag unset
    init(nCoeff,nCoeff);

    // save the histogram and polynomial order
    _histogram = histogram;
    _poly.resize(nPoly);
    _xOffset = xOffset;
    _xStep   = xStep;

    // make the design matrix X_tc = X_ij
    for (int jCoeff=0; jCoeff<nCoeff; jCoeff++)
    {
        dVector phi1 = getPhi(jCoeff);
        dVector column; column.fill(0.,nCoeff);
        for (int jCoeff2=0; jCoeff2<nCoeff; jCoeff2++)
        {
            dVector phi2 = getPhi(jCoeff2);
            for (int jBin=0; jBin<nBins; jBin++)
                column[jCoeff2] += phi1[jBin] * phi2[jBin];
        }
        addOrInsertBasisFunction(jCoeff,column);
    }

    // set uniform weights as default
    dVector weights;  weights.fill(1.,nCoeff);
    setWeights(weights);
    calculatePseudoInverse();

    // create the modified y vector and then fit it.
    dVector yVector;
    for (int jCoeff=0; jCoeff<nCoeff; jCoeff++)
    {
        dVector phi = getPhi(jCoeff);
        double sum = 0.;
        for (int jBin=0; jBin<nBins; jBin++)
            sum += _histogram[jBin] * phi[jBin];
        yVector.append(sum);
    }
   fitWLS(yVector);

    // Now beta_c are available in the modified framework; use getters to convert them
    _mu = getGaussianMu();
    _sigma2 = getGaussianSigma2();
    for (int jPoly=0; jPoly<nPoly; jPoly++)
        _poly[jPoly] = getPolynomialBeta(jPoly);
    // Calculate mag by regressing out polynomial and then finding best mag
    dVector gaussianData = _histogram;
    for (int jBin=0; jBin<nBins; jBin++)
    {
        double x = xOffset + xStep * jBin;
        for (int jPoly=0; jPoly<nPoly; jPoly++)
            gaussianData[jBin] -= _poly[jPoly] * utilMath::polynomialLegendre(jPoly,x);
    }
    double numerator=0.;  double denominator=0.;
    for (int jBin=0; jBin<nBins; jBin++)
    {
        double x = xOffset + xStep * jBin;
        double gaussShape = qExp( -SQR(x-_mu)/2./_sigma2 );
        numerator   += gaussianData[jBin] * gaussShape;
        denominator += gaussShape*gaussShape;
    }
    _mag = numerator / denominator;
    // assign the fit
    _fit.resize(nBins);
    for (int jBin=0; jBin<nBins; jBin++)
    {
        double x = xOffset + xStep * jBin;
        _fit[jBin] = _mag * qExp( -SQR(x-_mu)/2./_sigma2 );
        FUNC_INFO << "x, diff, arg, fit" << x << x-_mu << SQR(x-_mu)/2./_sigma2 << _fit[jBin];
        for (int jPoly=0; jPoly<nPoly; jPoly++)
            _fit[jBin] += _poly[jPoly] * baselineBasisFunction(jPoly,x);
    }
    FUNC_EXIT << "fit" << _fit;
}
dVector GaussianPlusPolynomial::getPhi (int order)
{
    int nPoly = _poly.size();
    int nBins = _histogram.size();
    dVector phi;
    if ( order <= nPoly+1 )
    {
        for ( int jBin=0; jBin<nBins; jBin++ )
        {
            double x = _xOffset + _xStep * jBin;
            phi.append( baselineBasisFunction(order,x) );  // x^k
        }
    }
    else if ( order == nPoly+2)
    {
        for ( int jBin=0; jBin<nBins; jBin++ )
        {
            double sum = 0.;
            for (int jSum=0; jSum<jBin; jSum++)
            {
                double x = _xOffset + _xStep * jSum;
                sum += x * _histogram[jSum] * _xStep;  // integral_0_x{u*f(u)*du)
            }
            phi.append(sum);
        }
    }
    else // if ( order == nPoly+3)
    {
        for ( int jBin=0; jBin<nBins; jBin++ )
        {
            double sum = 0.;
            for (int jSum=0; jSum<jBin; jSum++)
                sum += _histogram[jSum] * _xStep;  // integral_0_x{f(u)*du)
            phi.append(sum);
        }
    }
    return phi;
}
double GaussianPlusPolynomial::baselineBasisFunction(int iPoly, double x)
{
    double value;
    if ( iPoly == 0 )
        value = 1.;
    else if ( iPoly == 1 )
        value = x;
    else
        value = qPow(x,static_cast<double>(iPoly));
    return value;
}
double GaussianPlusPolynomial::getGaussianMu()
{
    int nPoly = _poly.size();
    double numerator   =   getBeta(nPoly+3);
    double denomenator = - getBeta(nPoly+2);
    return numerator / denomenator;
}
double GaussianPlusPolynomial::getGaussianSigma2()
{
    int nPoly = _poly.size();
    FUNC_INFO << nPoly+2;
    double numerator   = 1.;
    double denomenator = - getBeta(nPoly+2);
    return numerator / denomenator;
}
double GaussianPlusPolynomial::getPolynomialBeta(int iCoeff)
{
    int nPoly = _poly.size();
    if ( iCoeff < 0 || iCoeff >= nPoly )
        qFatal("Error: iCoeff out of range in GaussianPlusPolynomial::getPolynomialBeta");
    else
    {
        double sigma2 = getGaussianSigma2();
        double mu     = getGaussianMu();
        if ( iCoeff == nPoly-1 || iCoeff == nPoly-2 )
        {
            double alpha_nm1 = static_cast<double>(nPoly+1)*getBeta(nPoly+1)*sigma2;
            if ( iCoeff == nPoly-1 )
                return alpha_nm1;
            else
            {
                double alpha_nm2 = mu * alpha_nm1 + static_cast<double>(nPoly)*getBeta(nPoly)*sigma2;
                return alpha_nm2;
            }
        }
        else
        {
            int k = iCoeff + 2;
            double alpha_km2 = static_cast<double>(k) * sigma2 * ( getBeta(k) - getPolynomialBeta(k) )
                    + mu * getPolynomialBeta(k-1);
            return alpha_km2;
        }
    }
}
