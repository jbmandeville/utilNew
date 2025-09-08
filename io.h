#ifndef IO_H
#define IO_H

#include <QtGlobal>
#include <QComboBox>
#include <QString>
#include <QVector>
#include <QProcess>
#include <QtMath>
#include <QElapsedTimer>
#include <QDebug>
#include <QDir>
#ifdef USE_FFTW
   #include <fftw3.h>
#endif

#if QT_VERSION < QT_VERSION_CHECK(5, 14, 0)
namespace Qt
{
    static auto endl = ::endl;
    static auto SkipEmptyParts = QString::SkipEmptyParts;
    void enableComboBoxItem(QComboBox *comboBox, int itemNumber, bool enable);
    void visibleComboBoxItem(QComboBox *comboBox, int itemNumber, bool enable);
    int locateStringInComboxBox(QComboBox *combo, QString text);
}
#endif

#define FUNC_ENTER qDebug() << Q_FUNC_INFO << "enter"
#define FUNC_EXIT  qDebug() << Q_FUNC_INFO << "exit"
#define FUNC_INFO  qDebug() << Q_FUNC_INFO
#define QINFO      qInfo()  << Q_FUNC_INFO

#define PI 3.141592654
#define PI_180 (PI/180.)
#define RAND_0_1 (rand()/2147483647.)

#define NCOLORS 15 // do not include undefined value
enum overlayColor : short
{   //                 R   G   B
    Color_Cyan,    //  0  255 255
    Color_Teal,    //  0  128 128
    Color_Purple,  // 128 128  0
    Color_Magenta, // 255  0  255
    Color_Pink,    // 255 192 203
    Color_Red,     // 255  0   0
    Color_Lime,    //  0  255  0
    Color_Green,   //  0  128  0
    Color_Sky,     //  0  191 255
    Color_Blue,    //  0   0  255
    Color_Yellow,  // 255 255  0
    Color_Orange,  // 255 165  0
    Color_Tan,     // 210 180 140
    Color_Turquoise,//  0 206 209
    Color_Brown,   // 139 69  19
    Color_Undefined // this must be last, so that the enum values match the combobox
};

enum fileStorageTypes
{
    BFLOAT,
    BSHORT,
    BLONG,
    UINT8,
    USHORT,
    BDOUBLE,
    NIFTI_NII,
    NIFTI_NII_GZ,
    NIFTI_IMG,
    TypeUnknown
};

enum fileCategory
{ // this category mixes data types and parameter interpretation (!?)
    category_none,
    category_timeSeries,
    category_anatomy,
    category_atlas,
    category_summary,
    category_preProcessingSource,
    category_color,
    category_pValue,
    category_TValue,
    category_signalValue,
    category_BPValue,
    category_Occupancy,
    category_T2,
    category_R2,
    category_ADC,
    category_vector,
    category_RGB,
    category_complex,
    category_tensor,
    category_chi2,
    fileCategory_transient // files generated within the program but not saved
};

struct iPoint2D {union {int x; int lower;}; union {int y; int upper;}; };
struct iPoint3D {int x, y, z;};
struct iPoint4D {int x, y, z, t;};
struct bPoint3D {bool x, y, z;};
struct Point2D  {float x, y;};
struct Point3D  {float x, y, z;};
struct Point4D  {float x, y, z, t;};
struct dPoint2D {union {double x; double lower; double value;}; union {double y; double upper; double error;}; };
struct dPoint3D {double x, y, z;};
struct dPoint4D {double x, y, z, t;};
struct Mat33 {double m[3][3];};
struct Mat44 {double m[4][4];};
struct fComplex {float real, imag;};
struct dComplex {double real, imag;};
struct int8{unsigned int data : 8;};

typedef QVector<short>              shortVector;
typedef QVector<int>                iVector;
typedef QVector<bool>               bVector;
typedef QVector<QChar>              cVector;
typedef QVector<QString>            sVector;
typedef QVector<double>             dVector;
typedef QVector<float >             fVector;
typedef QVector<iPoint2D>          i2Vector;
typedef QVector<iPoint3D>          i3Vector;
typedef QVector<iPoint4D>          i4Vector;
typedef QVector<dPoint2D>          d2Vector;
typedef QVector<dPoint3D>          d3Vector;
typedef QVector<fComplex>          fcVector;
typedef QVector<dComplex>          dcVector;
typedef QVector<Point3D>            f3Vector;
typedef QVector<QVector<QString>>   sMatrix;
typedef QVector<QVector<bool>>      bMatrix;
typedef QVector<QVector<QChar>>     cMatrix;
typedef QVector<QVector<int>>       iMatrix;
typedef QVector<QVector<float>>     fMatrix;
typedef QVector<QVector<double>>    dMatrix;
typedef QVector<QVector<fComplex>> fcMatrix;
typedef QVector<QVector<Point3D>>  f3Matrix;
typedef QVector<QVector<dPoint2D>> d2Matrix;
typedef QVector<QVector<dPoint3D>> d3Matrix;
typedef QVector<QVector<iPoint2D>> i2Matrix;
typedef QVector<QVector<QVector<int>>>    iMatrix3;
typedef QVector<QVector<QVector<double>>> dMatrix3;
struct fVoxel      {int x, y, z, index3d=0; double binaryValue; overlayColor color=Color_Undefined;};  // can hold a voxel location, value, and optionally a color for display
struct fVoxelLight {int x, y, z; int t=0; double binaryValue;};
typedef QVector<fVoxel>      VoxelSet;
typedef QVector<fVoxelLight> VoxelSetLight;
struct kernelVoxel {int x, y, z, t; double value;};
typedef QVector<kernelVoxel> kernelSet;
struct OVLFile
{  // typical overlay-list format: [name] [path] [optional color]
    iPoint3D dimensions={0,0,0};     // keeps image dimensions (not always available)
    VoxelSet voxelList;              // vector of (x,y,z,value)
    QString name;                    // e.g. putamen
    QString path;                    // e.g. $TemplateDir/putamen.ovl
    overlayColor color = Color_Cyan; // Attach a color to an overlay file
    int IDRegion;                    // indicates region from atlas, with ID as given (e.g region 21); set to 0 for non-atlas regions
    bool saved=true;
    bool loaded=false;
};
struct OVLLight  // lightweight version of above
{  // typical overlay-list format: [name] [path] [optional color]
    VoxelSetLight voxelList;         // vector of (x,y,z,value); sometimes (x,y,z,t,value)
    int IDRegion;                    // indicates region from atlas, with ID as given (e.g region 21); set to 0 for non-atlas regions
                                     // IDRegion also could hold a 3d/4d index into voxel data
};
typedef QVector<OVLFile>  ovlVector;
typedef QVector<OVLLight> ovlVectorLight;
typedef struct {QVector<double> x, y, fit; QString name;} GraphVector;

struct basisFunction
{
    int IDRegion;          // index in basis function vector
    iPoint3D center;       // for some types of symmetric basis functions
    iPoint3D neighborLow;  // for basis functions, one can keep track of neighbors (low  side)
    iPoint3D neighborHigh; // for basis functions, one can keep track of neighbors (high side)
    VoxelSet voxelList;    // vector of fVoxels = (x,y,z,index3d,value,color)
    double dSum=0.;        // sum of overlay (not always computed)
    double scalarBeta=1.;  // a scaling factor (beta) applied to basis function
    double vectorBeta[3]={1.,1.,1.};
};
typedef QVector<basisFunction> bfVector;
struct parameterMaps
{
    QString parameterTypeName;  // condition-type name (e.g., "P" or "BP")
    QString mapTip;             // e.g., "T2 is the transverse relaxation time for MRI signal"
    // note that
    // 1) the file-name for output maps will always be "parameterTypeName"-"conditionName" (e.g., P-1) or just "conditionName" if "parameterTypeName" is empty
    // 2) the parameter name (e.g, shown in GUI) is "conditionName"; this could be a GLM condition or just "T2"
    fileCategory category;   // e.g. category_pValue
    dVector scalarValues;    // usually 1 or 2 (2 = mean,error) but could be more (e.g., signal, stdev, snr)
    dPoint3D vectorValue={0.,0.,0.};   // e.g. FA
    dPoint2D colorBarRange;
};
struct conditionMaps
{   // dimension this to [nConditions]; this holds condition values for a single overlay
    // For each overlay, all of the following are the same
    QString conditionName;  // e.g., "1" or "a-b" for GLM conditions or just something like "T2"
    QString conditionTip;   // e.g., "T2 is the transverse relaxation time for MRI signal"
    int DOF;
    QVector<parameterMaps> parameterMap;  // e.g., (P,T,S)
};

struct imageHeader
{ // all info to interpret a file should be stored in the header in order to facilitate voxel-reading in a different thread
    QString fileName;       // original file name (if read from file), including directory path
    iPoint4D dim={0,0,0,0}; // # voxels in x, y, z, t
    dPoint4D resolution={0.,0.,0.,0.};    // resolutions in x, y, z
    iPoint4D points={0,0,0,0};        // (nx, nx*ny, nx*ny*nz, nx*ny*nz*nt)
    dMatrix ijk_to_xyz;     // 4x4 transformation from data indices (i,j,k) to space (x,y,z)
    dMatrix xyz_to_ijk;     // 4x4 inverse transformation
    int byte_order;         // byte order 0 = Mac, 1 = PC
    bool offDiagonal=false ;// if true, the transformation matrix does not have parallel voxels & coordinates
    int fileType=NIFTI_NII; // e.g., bshort, bfloat, blong, nifti_nii, nifti_img
    int dataType=0;         // 0 = float, 1 = real-imaginary, 2 = magnitude-phase, 3 = vector, 4 = rbg
    int storageType=BFLOAT; // storage type: bshort, bfloat, blong
    double DOF=0.;          // degrees of freedom for statistics; default value 0
    int voxelOffset;        // offset to the first voxel data in the file
    float scaleSlope;       // scaling factor for data (data --> slope*data + offset); ignore if 0; default 0
    float scaleOffset;      // offset for data; default 0
    double noiseLevel=0.;   // st Dev of noise; determined by testing edges of image
};

struct volumeData
{
    fVector  f1;            // 3-dimensional volume data vector [i3d]
    fcVector fc;            // 3-dimensional complex vector
    f3Vector f3;            // 3-dimensional 3-vector (could also be RGB)
    // (0=x,1=y,2=z,3=phi_x,4=phi_y,5=phi_z,6=r_displace,7=COM.x,8=COM.y,9=COM.z)
    dVector mcParameters={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
};

// Define a structure for a stack of images
struct imageData
{
    imageHeader hdr;
    QVector<volumeData> timePoints;  // allocated as a time vector per file
    fVector averageVolume;
    shortVector color;                  // 3-dimensional color (for bitmaps/overlays)
    bool empty=false;                   // if true, all data is zero

//    iPoint4D dim={0,0,0,0}; // # voxels in x, y, z, w=t
//    Point4D resolution;    // resolutions in x, y, z, t
//    Point4D origin;        // origins in x, y, z, t
//    Point3D direction;     // directions in x, y, z
//    iPoint4D points;       // (nx, nx*ny, nx*ny*nz, nx*ny*nz*nt)
//    Mat44 ijk_to_xyz;      // transformation from data indices (i,j,k) to space (x,y,z)
//    Mat44 xyz_to_ijk;      // inverse transformation
};

struct ROI_data
{
    QString fileName;// Vector of file names for each run in the ROI
    QString name;    // some identifier
    iPoint3D voxel;  // voxel for 1-point ROI
    VoxelSet voxelList; // vector of (x,y,z,value)
    double nspace=0; // # points in space (non-integer due to potential non-binary voxel weights)
    double dt=0;     // time step
    double dof;
    dVector xSignal; // x values, including potentially a known time step read from files
    dVector ySignal; // original signal time series
};

namespace utilIO
{
    template <class T> T* allocate_vector(int n)
    {
        T* a = new T[n];
        return a;
    }
    template <class T> void delete_vector(T *vec)
    {
        delete vec;
    }

    int machine_byte_order();
    void swap(short *ptr);
    void swap(unsigned short *ptr);
    void swap(float *ptr);
    void swap(int  *ptr);
    void swap_16( void *ptr );
    void swap_32( void *ptr );

    inline int index1d(iPoint3D dim, fVoxel voxel)   {return voxel.z*dim.y*dim.x + voxel.y*dim.x + voxel.x;}
    inline int index1d(iPoint3D dim, iPoint3D voxel) {return voxel.z*dim.y*dim.x + voxel.y*dim.x + voxel.x;}
    inline int index1d(iPoint3D dim, int iX, int iY, int iZ) {return iZ*dim.y*dim.x + iY*dim.x + iX;}
    inline int index1d(iPoint4D dim, iPoint4D voxel)
    {
        return voxel.t*dim.z*dim.y*dim.x + voxel.z*dim.y*dim.x + voxel.y*dim.x + voxel.x;
    }
    inline int index1d(iPoint4D dim, int iX, int iY, int iZ, int iT)
    {
        return iT*dim.z*dim.y*dim.x + iZ*dim.y*dim.x + iY*dim.x + iX;
    }
    inline iPoint3D index3d(const iPoint3D dim, const int index3d)
    {
        iPoint3D voxel;
        voxel.z = index3d/(dim.x*dim.y);
        voxel.y = (index3d-voxel.z*dim.x*dim.y)/dim.x;
        voxel.x = index3d - voxel.z * dim.x*dim.y - voxel.y*dim.x;
        return voxel;
    };
    inline bool voxelInBounds(iPoint3D dim, iPoint3D voxel) {return
                voxel.x>=0    && voxel.y>=0    && voxel.z>=0 &&
                voxel.x<dim.x && voxel.y<dim.y && voxel.z<dim.z;}
    inline bool voxelInBounds(iPoint3D dim, fVoxel voxel) {return
                voxel.x>=0    && voxel.y>=0    && voxel.z>=0 &&
                voxel.x<dim.x && voxel.y<dim.y && voxel.z<dim.z;}

    int readOverlayFile(VoxelSet &voxelList, QString fileName, overlayColor color, iPoint3D dimSpace);

    void delayMS( int millisecondsToWait );
    void elapsedTime(const QString text, QElapsedTimer *timer);
    QString readTimeTableFile(QString fileName, QStringList &columnNames, dMatrix &table);
    int reflectBoundary(int index1d, int pad, int dim);

}

namespace tk
{
// band matrix solver
class band_matrix
{
private:
    dMatrix m_upper;  // upper band
    dMatrix m_lower;  // lower band
public:
    // constructors/destructors
    band_matrix() {};
    inline band_matrix(int dim, int n_u, int n_l) {resize(dim, n_u, n_l);}
    ~band_matrix() {};                       // destructor
    void resize(int dim, int n_u, int n_l);  // init with dim,n_u,n_l
    inline int dim() const {if ( m_upper.size() > 0) return m_upper[0].size(); else return 0;}
    int num_upper() const
    {
        return m_upper.size()-1;
    }
    int num_lower() const
    {
        return m_lower.size()-1;
    }
    // access operator
    double & operator () (int i, int j);            // write
    double   operator () (int i, int j) const;      // read
    // we can store an additional diogonal (in m_lower)
    // second diag (used in LU decomposition), saved in m_lower
    inline double saved_diag(int i) const {Q_ASSERT( (i>=0) && (i<dim()) ); return m_lower[0][i];}
    inline double &saved_diag(int i) {Q_ASSERT( (i>=0) && (i<dim()) ); return m_lower[0][i];}

    void lu_decompose();
    dVector r_solve(const dVector& b) const;
    dVector l_solve(const dVector& b) const;
    dVector lu_solve(const dVector& b,bool is_lu_decomposed=false);
};

// spline interpolation
class spline
{
public:
    enum bd_type {
        first_deriv = 1,
        second_deriv = 2
    };

private:
    dVector m_x,m_y;            // x,y coordinates of points
    // interpolation parameters
    // f(x) = a*(x-x_i)^3 + b*(x-x_i)^2 + c*(x-x_i) + y_i
    dVector m_a,m_b,m_c;        // spline coefficients
    double  m_b0, m_c0;                     // for left extrapol
    bd_type m_left, m_right;
    double  m_left_value, m_right_value;
    bool    m_force_linear_extrapolation;

public:
    // set default boundary condition to be zero curvature at both ends
    spline(): m_left(second_deriv), m_right(second_deriv),
        m_left_value(0.0), m_right_value(0.0),
        m_force_linear_extrapolation(false)
    {
        ;
    }

    // optional, but if called it has to come be before set_points()
    void set_boundary(bd_type left, double left_value,
                      bd_type right, double right_value,
                      bool force_linear_extrapolation=false);
    void set_points(const dVector& x,
                    const dVector& y, bool cubic_spline=true);
    double operator() (double x) const;
};

}

namespace utilMath
{
#define SQR(a) ((a) == 0.0 ? 0.0 : (a)*(a))
#define NMAT 3

    dComplex complexMultiply(dComplex dc1, dComplex dc2, bool star);
    double computePhase(dComplex dc);
    double computeAmplitude(dComplex dc);
    double dotProduct(dVector v1, dVector v2);
    double dotProduct(dPoint3D v1, dPoint3D v2);
    double vectorDistance(dPoint3D v1, dPoint3D v2);
#ifdef USE_FFTW
    void FFTW1D_Volume( iPoint3D dim, iPoint3D pad, dcVector &volume, int iDimFFT, bool forward );
    void FFTW1DVolumeInThread( iPoint3D dim, dcVector &volume, int iDimFFT, bool forward, iPoint3D pad,
                               dComplex *dcDataX, fftw_complex *fftwDataX, fftw_plan planBackwardX, fftw_plan planForwardX,
                               dComplex *dcDataY, fftw_complex *fftwDataY, fftw_plan planBackwardY, fftw_plan planForwardY,
                               dComplex *dcDataZ, fftw_complex *fftwDataZ, fftw_plan planBackwardZ, fftw_plan planForwardZ);
#endif
    Mat44 invertMat44(Mat44 matrix44 );
    bool invertSquareMatrix(dMatrix &dSquareMatrix);
    void LUDecomposeAndSolve(dVector YVector, dMatrix X, dVector &beta);
    void LUDecompose(dMatrix XMat, dMatrix &lower, dMatrix &upper);
    void LUDecomposeNew(dMatrix XMat, dMatrix &lower, dMatrix &upper);
    void LUSolve(dVector YVec, dMatrix lower, dMatrix upper, dVector &beta);
    void CholeskyDecomposeAndSolve(dVector YVec, dMatrix XMat, dVector &beta);
    void CholeskyDecompose(dMatrix matrix, dMatrix &lower, dMatrix &lowerT);

    double deltaKernel(double x);
    double LanczosKernel(double x, double width);
    double GaussKernel(double x, double fwhm);
    double inverse_Gauss(double x, double fwhm);
    void swapX( dComplex *dcData, int lDim, int lAdd );
    bool ParabolicInterpolation(dVector xParabola, dVector yParabola, double &xMax, double &yMax);
    void fuzzyBinning(double value, double min, double max, int nBins, iPoint2D &iBin, dPoint2D &weightBin );
    double polynomialLegendre(int iPoly, double x);

    dMatrix alignTransformForward(bool consistentCoordinates,
                                  dMatrix affine,         dMatrix sourceParallelXYZtoIJK, dMatrix sourceOriginalXYZtoIJK,
                                  dMatrix targetIJKtoXYZ, dMatrix resliceMatrix,          dMatrix &transformDistort);
    dMatrix alignTransformInverse(bool consistentCoordinates,
                                  dMatrix affineInverse,  dMatrix sourceParallelIJKtoXYZ, dMatrix sourceOriginalIJKtoXYZ,
                                  dMatrix targetXYZtoIJK, dMatrix resliceMatrixInverse,   dMatrix &transformDistort);
    void matrixTransform(dMatrix matrix, dVector &input);
    void matrixTransform(dMatrix matrix, dPoint4D &v4d);
    void matrixTransform(dMatrix matrix, fVoxel &v4d);
    void matrixTransform(dMatrix matrix, dPoint3D &v3d);
    dPoint3D matrixTransform(dMatrix matrix, iPoint3D i3d);
    dMatrix matrixProduct( dMatrix X1, dMatrix X2 );
    void matrixProductTranspose1( dMatrix X1, dMatrix X2, dMatrix &X1X2 );
    void matrixProductTranspose2( dMatrix X1, dMatrix X2, dMatrix &X1X2 );
    void matrixPseudoInverse( dMatrix X_12, dMatrix &piX_21 );
    double traceOfSquareMatrix(dMatrix mat);

    dPoint3D normalizedVector(dPoint3D vec);
    dPoint3D vectorStep(dPoint3D rOrig, dPoint3D unitVector, float step);
    basisFunction createBSplineShape(iPoint3D grid);
    void defineBasis1DBSpline3Dim(iPoint3D grid, dMatrix &basis1d);
    double basisValue1dBSpline(double xAbs, double xRel);
    int OtsuThreshold(dVector histogram);

    void eigen_decomposition(double A[NMAT][NMAT], double V[NMAT][NMAT], double d[NMAT]);
    inline double hypot2(double x, double y) {return qSqrt(x*x+y*y);}
    void tred2(double V[NMAT][NMAT], double d[NMAT], double e[NMAT]);
    void tql2(double V[NMAT][NMAT], double d[NMAT], double e[NMAT]);

    double ibeta(double aa, double bb, double xx);
    double incbcf(double a, double b, double x);
    double incbd(double a, double b, double x);
    double pseries(double a, double b, double x);

    void topDownMergeSort(dVector &array);
    void topDownSplitMerge(int iBegin, int iEnd, dVector &array, dVector &arrayWorking);
    void topDownMerge(int iBegin, int iMiddle, int iEnd, dVector &array, dVector &arrayWorking);
    inline double medianValue(dVector array) {if (array.isEmpty()) return 0.; else {topDownMergeSort(array); return array[array.size()/2];}}
    void copyArray(int iBegin, int iEnd, dVector &array, dVector &arrayWorking);

    void topDownMergeSort(dVector &array, iVector &indexArray);
    void topDownSplitMerge(int iBegin, int iEnd, dVector &array, dVector &arrayWorking, iVector &indexArray, iVector &indexWorking);
    void topDownMerge(int iBegin, int iMiddle, int iEnd, dVector &array, dVector &arrayWorking, iVector &indexArray, iVector &indexWorking);
    void copyArray(int iBegin, int iEnd, dVector &array, dVector &arrayWorking, iVector &indexArray, iVector &indexWorking);

    VoxelSet createVoxelList(int iThread, int nThreads, iPoint3D dim);
    iMatrix setupSliceProcessing(int nZ, int &nThreads);

    inline double randomFlat(double low, double high) {return RAND_0_1*(high-low)+low;}
    double randomGaussian(double fwhm, double cutoff);

}
namespace utilString
{
    QString createMapName(QString conditionName, QString parameterTypeName);
    QString createMapFileName(QString conditionName, QString parameterTypeName, bool nifti);
    void errorMessage(QString errorText);
    int decodeSelectionList(QString inputList, bool usePlusOneSyntax, iVector &includeVolume);
    int decodeSelectionList(QString inputList, bool usePlusOneSyntax, bVector &includeVolume);
    QString recodeSelectionList(iVector includeVolume, bool usePlusOneSyntax);
    QString recodeSelectionList(bVector includeVolume, bool usePlusOneSyntax);
    void replaceEnvironmentVariables(QStringList &inputStringList);
    QString replaceEnvironmentVariables(QString inputString);
    QString insertEnvVariables(QString inputString, QStringList templateDirectories);
    QString getFileNameExtension(QString fileName);
    QString getFileNameWithoutExtension(QString fileName);
    QString getFileNameWithoutDirectory(QString fileName);
    QString getFileBaseName(QString fileName);
    QString getDirectoryName(QString fileName);
    bool fileHasExtension(QString fileName);
    int extension_type(const QString extension );
    QString attachExtension(QString fileName, int fileType);
    QString ensureExtension(QString fileName);
}
#endif // IO_H
