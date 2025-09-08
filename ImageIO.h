#ifndef IMAGEIO
#define IMAGEIO

#include <QMainWindow>
#include <QProgressBar>
#include <QMap>
#include <QString>
#include <QVector>
#include <QVector2D>
#include <QVector3D>
#include <QVector4D>
#include <QFile>
#include <QFileInfo>
#include <QDir>
#include <QWidget>
#include <QRunnable>
#include <QDebug>
#include "io.h"
#ifdef USE_FFTW
    #include "fftw3.h"
#endif

using namespace utilIO;

/*
inline QDataStream &operator>>(QDataStream &in, QVector<short> &vector)
{
   in.readRawData(reinterpret_cast<char*>(&vector), sizeof(vector));
   return in;
}
*/

enum downSampleTypes
{
    downSample1,
    downSample2,
    downSample4
};
enum CostFunctionChoice
{
    CostFunction_MI,   // mutual information
    CostFunction_WLS,  // weighted least squares
    CostFunction_RV,   // residual variance
    CostFunction_NMI,  // normalized mutual information
    CostFunction_NC    // normalized correlation
};
enum vectorScalars
{
    vectorScalar_Magnitude,
    vectorScalar_x,
    vectorScalar_y,
    vectorScalar_z
};

#define NCOSTBINS 100 // 1% per bin

class bitMap;

class ImageIO
{
    struct NiftiDataType
    {
        int code, size, swapsize;
        QString name;
    };
private:
    QString _imageFileName; // File name exactly as passed
    QFileInfo _fileInfo;
    int _fileID=-1;        // ID files by a number (usually location in series)
    bVector _loaded;       // [nThreads]; load/reslice/map file in separate threads using initLoaded; setLoaded; isLoaded;
    fileCategory _category=category_none; // categorize datasets for presentation grouping, thresholds, etc.
    int _swap;             // True if byte order on disk is different to CPU
    int _iTime=0;          // points to current volume index
    dPoint3D _xdOrigin={0.,0.,0.};
    iPoint3D _xdDirection={1,1,1};
    vectorScalars _vectorScalar=vectorScalar_Magnitude;
    QElapsedTimer _timer;

    int read_xdisplay_header();
    int decode_xdisplay_header_keyword( const QString line);
    int readNiftiHeader( bool verbose );

    dMatrix nifti_quatern_to_mat44( float qb, float qc, float qd,
                                    float qx, float qy, float qz,
                                    float dx, float dy, float dz, float qfac );
    void nifti_mat44_to_quatern( dMatrix R, float *qb, float *qc, float *qd,
                                 float *qx, float *qy, float *qz,
                                 float *dx, float *dy, float *dz, float *qfac );
    Mat33 nifti_mat33_polar( Mat33 A );
    Mat33 nifti_mat33_inverse( Mat33 R );
    float nifti_mat33_determ( Mat33 R );
    float nifti_mat33_rownorm( Mat33 A );
    float nifti_mat33_colnorm( Mat33 A );

    void swap_16( void *ptr );
    void swap_32( void *ptr );
    dMatrix general_affine_to_mat44( float *srow_x, float *srow_y, float *srow_z );

    static void SwapBytes(size_t n, int siz, void *ar);
    static void SwapNiftiHeader(struct nifti_1_header *h);

    void dilateOrErode(bool dilate, int iTime);

    double spatialPoly1D(int iPower, int iPos, iPoint2D range);

protected:  // can be accessed from derived classes
    void CheckVoxelIndices(fVoxel voxel, const QString callingFunction);
    void CheckVoxelIndices(iPoint3D voxel, const QString callingFunction);
    void CheckVoxelIndices(int iX, int iY, int iZ, const QString callingFunction);
    void CheckVoxelIndices( int i3d, const QString callingFunction );

public:
    imageData _data; // make this public until I figure out how to get it for swapping
    // ImageIO( const ImageIO &imageStack);
    // Public I/O functions
    int readFileHeader(QString fileName, bool verbose );
    int readImageData();
    int readImageData(QProgressBar *progressBar);
    void resliceImageData(ImageIO *templateImage);
    void resliceImageData(ImageIO *templateImage, QProgressBar *progressBar);
    dPoint3D convertVoxelToSpace(iPoint3D iVoxel);
    dPoint3D convertVoxelToSpace( fVoxel voxel);
    dPoint3D convertVoxelToSpace(dPoint3D rVoxel);
    iPoint3D convertSpaceToVoxel(dPoint3D space);
    dPoint3D convertSpaceToVoxelFraction(dPoint3D space);
    double interpolateValueFromSpaceCoord(dPoint3D space, dPoint3D fwhm, iPoint3D iWidth, iPoint3D wrap,
                                          iPoint3D kernelType, int iTime, bool &inVolume);
    double interpolateValueFromSpaceCoord(dPoint3D space, dPoint3D fwhm, iPoint3D iWidth, iPoint3D wrap,
                                          iPoint3D kernelType, int iTime);
    double interpolateValueFromSpaceCoordBilinear(dPoint3D space, int iTime);
    double interpolateValueFromrIndex(dPoint3D rIndex, dPoint3D fwhm, iPoint3D iWidth,
                                      iPoint3D wrap, iPoint3D kernelType, int iTime);
    double interpolateValueFromrIndex(dPoint4D rIndex, dPoint3D fwhm, iPoint3D iWidth, iPoint3D wrap,
                                      iPoint3D kernelType, int iTime, bool &inVolume);
    double interpolateValueFromrIndex(dPoint4D rIndex, dPoint3D fwhm, iPoint3D iWidth, iPoint3D wrap,
                                      iPoint3D kernelType, int iTime);
    double interpolateValueFromrIndex(dPoint3D rIndex, dPoint3D fwhm, iPoint3D iWidth, iPoint3D wrap,
                                      iPoint3D kernelType, int iTime, bool &inVolume);
    double interpolateValueFromrIndexBilinear(dPoint3D rIndex, int iTime);
    void flattenByGLM(iPoint3D dimPoly, bitMap mask);
#ifdef USE_FFTW
    void smoothVolumeInThread(int iTime, dPoint3D fwhm, iPoint3D pad,
                              dComplex *dcDataX, fftw_complex *fftwDataX, fftw_plan planBackwardX, fftw_plan planForwardX,
                              dComplex *dcDataY, fftw_complex *fftwDataY, fftw_plan planBackwardY, fftw_plan planForwardY,
                              dComplex *dcDataZ, fftw_complex *fftwDataZ, fftw_plan planBackwardZ, fftw_plan planForwardZ);
    void smoothVolume(int iTime, double fwhm);
    void smoothVolume(int iTime, dPoint3D fwhm);
    void smoothVolumes(double fwhm);
    void smoothVolumes(dPoint3D fwhm);
#endif
    void smoothVolumeNoFFT(int iTime, double fwhm);
    void smoothVolumeNoFFT(int iTime, dPoint3D fwhm);
    void smoothVolumesNoFFT(double fwhm);
    void smoothVolumesNoFFT(dPoint3D fwhm);
    void dotProductWithKernelStride1(ImageIO *inputImage, int iTime, iPoint4D kernelSize, kernelSet kernel);
    void poolStride2(ImageIO *inputImage, int iTime, int iStride, kernelSet kernel, bool maxPool);
    void logTransform(bool forward, bitMap mask);
    void logTransform(bool forward);
    inline void dilate(int iTime) {dilateOrErode(true,iTime);}
    inline void erode(int iTime)  {dilateOrErode(false,iTime);}
    inline void close(int iTime)  {dilate(iTime); erode(iTime);}
    inline void open(int iTime)   {erode(iTime); dilate(iTime);}
    double calculateCostFunction(CostFunctionChoice costFunction,
                                 int iVolume, ImageIO *templateVolume, int iVolumeTemplate,
                                 bool includeZeroes);
    double calculateCostFunction(CostFunctionChoice costFunction,
                                 int iVolume, ImageIO *templateVolume, int iVolumeTemplate,
                                 bitMap *mask);
    double calculateCostFunction(CostFunctionChoice costFunction, int iVolume,
                                 ImageIO *templateVolume, int iVolumeTemplate);

    double getMinDistanceToBoundaryWeight(iPoint3D voxel);
    double calculateAverage(int iTime);
    double calculateAverage(int iTime, bitMap *mask);
    double calculateAverage(int iTime, VoxelSet voxelList);
    double calculateAverage(int iTime, OVLFile overlay);
    double calculateAverage(int iTime, ImageIO *templateFile, int iTimeTemplate, bool includeZeroes);
    dPoint2D getThreshold(int iTime);
    dPoint2D getThreshold();
    void resample( iPoint3D dimensionNew, dPoint3D resolutionNew );
    void upSampleOrDownSample(iPoint3D dimensionNew);
    void edgeFilter(int iTime);
    void edgeFilter3D(int iTime);
    void edgeFilter2D(int iTime);
    void cropUsingOverlay(bitMap mask);

    void NLMeans(int iTime, bool dim3D, int patchVoxels, int searchVoxels);

    int writeNiftiFile(QString dirName, bool useCompression);
    int writeNiftiData(QString fileName, bool useCompression);
    int writeNiftiHeader(QString fileName, bool useCompression);

    void setParallelCoordinates(dPoint3D resolution, dPoint3D origin, iPoint3D direction);
    void setParallelCoordinates(dPoint4D resolution, dPoint3D origin, iPoint3D direction);
    void setIjkToXyz(dMatrix ijkToxyz);

    // Public setters
    inline void clearFileName()        {setFileName("none");}
    inline void setFileID(int id) {_fileID = id;}
    inline void setFileName(QString fileName) {_imageFileName = fileName; _fileInfo.setFile(fileName);}
    void setFileName(QString rootName, int fileType);
    void setFileRootName(QString root) {setFileName(root + ".nii");}
    void resetDirName(QString newDirName);
    void resetFileName(ImageIO *templateImage, QString newName);
    inline void setCategory(fileCategory category) {_category=category;}
    inline void setVectorScalarMode(vectorScalars scalarMode) {_vectorScalar = scalarMode;}
    inline void defineNewMapSeries(QString mapName, fileCategory category, ImageIO *templateImage, int nTime, double DOF)
    { _category=category; copyTemplateDimensions(templateImage, nTime); setFileName(mapName); setDOF(DOF);}
    inline void setMotionCorrectionParameters(int iTime, dVector parameters) {_data.timePoints[iTime].mcParameters = parameters;}

    void setResolution(dPoint4D resolution);
    void setResolution(dPoint3D resolution);
    void setTimeResolution(double dt) {_data.hdr.resolution.t = dt;}
    inline void setStorageFormat(fileStorageTypes iType) {_data.hdr.storageType = iType;}
    inline void setNoiseLevel(double value) {_data.hdr.noiseLevel = value;}

    void initLoaded (int nThreads);
    void setLoaded();
    void setCurrentVolume(int iTime);
    void setStackValue(const int i3d, const int iT, const double value);
    void setStackValue(const int iX, const int iY, const int iZ, const int iT, const double value);
    void setStackValue(const iPoint3D voxel, const int iT, const double value);
    void setStackValue(const fVoxel voxel, const int iT);
    void setStackValue(const iPoint4D voxel, const double value);

    void setStackValue(const int i3d, const int iT, const Point3D value);
    void setStackValue(const iPoint3D voxel, const int iT, const Point3D value);

    void setParallelCenteredCoordinates();
    void setParallelCenteredCoordinates(iPoint4D dimension, dPoint4D resolution);
    void setParallelCenteredCoordinates(iPoint4D dimension, dPoint3D resolution);
    void setParallelCenteredCoordinates(iPoint3D dimension, dPoint3D resolution);
    void setParallelCOMCoordinates(dPoint3D COMIndex);

    void copyTemplate(ImageIO *templateImage);
    void copyTemplate(ImageIO *templateImage, QString name);
    void copyTemplateVolume(ImageIO *templateImage, int iVolume, QString name);
    void copyTemplateVolume(ImageIO *templateImage, int iVolume);
    void copyTemplateDimensions(ImageIO *templateImage, int nTime);
    void setDimensions(const iPoint4D dim);
    void setDimensions(const iPoint3D dim, const int nt);
    void setDimensions(const iPoint3D dim);
    void setDimensions(const int nx, const int nY, const int nZ, const int nT);
    void setOrigin(dPoint3D origin);

    void appendImageFileToTimeData(ImageIO imageFile);
    void deleteFirstVolume();

    inline void setImageHeader(imageHeader *hdr) {_data.hdr = *hdr;}
    inline void setDOF(double dof) {_data.hdr.DOF = dof;}

    // Public Getters
    bool voxelInBounds(iPoint3D voxel);
    bool voxelInBounds(iPoint4D voxel);
    bool voxelInBounds(fVoxel voxel);
    bool voxelInBounds(fVoxelLight voxel);
    inline bool voxelInBounds(int iVoxel) {return iVoxel >= 0 && iVoxel < nx();}

    inline bool isAllocated() {return _data.timePoints.size() != 0;}
    QString getFileHeaderName();
    inline int getFileType() {return utilString::extension_type(getExtension());}
    inline int getFileID() {return _fileID;}
    inline vectorScalars getVectorScalarMode() {return _vectorScalar;}
    inline int getCurrentVolumeIndex() {return _iTime;}
    inline imageHeader *getImageDataHeader() {return &_data.hdr;}
    inline double getDOF() {return _data.hdr.DOF;}
    bool isLoaded();
    bool noneLoaded();
    inline dVector getMotionCorrectionParameters(int iTime) {return _data.timePoints[iTime].mcParameters;}
    inline iPoint3D getImageDimensions() const
    {
        iPoint3D dim;
        dim.x = _data.hdr.dim.x;  dim.y = _data.hdr.dim.y;  dim.z = _data.hdr.dim.z;
        return dim;
    }
    inline dPoint3D getImageResolution() const
    {
        dPoint3D res;
        res.x = _data.hdr.resolution.x;  res.y = _data.hdr.resolution.y;  res.z = _data.hdr.resolution.z;
        return res;
    }
    inline iPoint4D getDimensions() const { return _data.hdr.dim;}
    inline dMatrix getIjkToXyz() const {return _data.hdr.ijk_to_xyz;}
    inline dMatrix getXyzToIjk() const {return _data.hdr.xyz_to_ijk;}

    inline int nx() const { return _data.hdr.dim.x; }
    inline int ny() const { return _data.hdr.dim.y; }
    inline int nz() const { return _data.hdr.dim.z; }
    inline int nt() const { return _data.timePoints.size(); }
    inline double dx() const { return _data.hdr.resolution.x; }
    inline double dy() const { return _data.hdr.resolution.y; }
    inline double dz() const { return _data.hdr.resolution.z; }
    inline double dt() const { return _data.hdr.resolution.t; }

    inline int index3dFromVoxel(const iPoint3D voxel)
    {
        return voxel.z * _data.hdr.points.y + voxel.y * _data.hdr.points.x + voxel.x;
    }
    inline int index3dFromVoxel(const fVoxel voxel)
    {
        return voxel.z * _data.hdr.points.y + voxel.y * _data.hdr.points.x + voxel.x;
    }
    inline int index3dFromVoxel(const fVoxelLight voxel)
    {
        return voxel.z * _data.hdr.points.y + voxel.y * _data.hdr.points.x + voxel.x;
    }
    inline int index3dFromVoxel(const int iX, const int iY, const int iZ)
    {
        return iZ * _data.hdr.points.y + iY * _data.hdr.points.x + iX;
    }

    inline dPoint4D getResolution() const {return _data.hdr.resolution;}
    dPoint3D getOrigin();
    dPoint3D getFOV();
    iPoint3D getDirection();
    iPoint3D determineDirectionMap(dMatrix ijk_to_xyz);
    inline bool hasOffDiagonalTransformation() {return _data.hdr.offDiagonal;}

    inline int voxelsPerSlice() const  { return _data.hdr.points.y; }
    inline int voxelsPerVolume() const { return _data.hdr.points.z; }
    inline int voxelsTotal() const     { return _data.hdr.points.t; }
    inline fVector getAverageVolume()  {return  _data.averageVolume;}

    inline int datatype() const { return _data.hdr.dataType; }
    inline fileCategory getCategory() {return _category;}
    inline bool is3Vector() {return _category == category_vector || _category == category_RGB;}
    inline double getNoiseLevel() {return _data.hdr.noiseLevel;}

    inline QString getReferencedName() {return _imageFileName;}               // exactly as it was passed
    inline QString getAbsoluteName()   {return _fileInfo.absoluteFilePath();} // absolute file name with full path
    inline QString getDirectory()      {return _fileInfo.absolutePath();}     // absolute path without file name
    inline QString getFileName()       {FUNC_ENTER << _fileInfo.fileName(); return _fileInfo.fileName();}         // full name with extension but not path
    QString getBaseName();                                                    // file name without path or extension
    inline QString getExtension()      {return _fileInfo.suffix();}           // extension (e.g, nii, gz, bfloat, ...)
    QString getRelativeName(QDir dir)  {return dir.relativeFilePath(getAbsoluteName());} // name relative to some directory

    dPoint3D coordinateFromVoxel(iPoint3D voxel);
    dPoint3D coordinateFromVoxel(int ix, int iY, int iZ);
    dcVector convertToComplexVolume(int iTime);

    // 4D data getters
    double getStackValue(int iX, int iY, int iZ, int it);
    double getStackValue(int i2d, int iZ, int it);
    double getStackValue(iPoint3D, int it);
    double getStackValue(fVoxel voxel, int it);
    double getStackValue(fVoxelLight voxelLight, int it);
    double getStackValue(iPoint4D voxel4D);
    double getStackValue(int i3d, int iTime);
    Point3D getStackVectorValue(iPoint3D voxel, int iTime);
    Point3D getStackVectorValue(fVoxel voxel, int iTime);

    // 3D data getters for volumes with a single time point
    double getVolumeValue(int iX, int iY, int iZ);
    double getVolumeValue(int i3d);
    double getVolumeValue(iPoint3D voxel);
    double getVolumeValue(fVoxel voxel);
    Point3D getVolumeVectorValue(int i3d);
    Point3D getVolumeVectorValue(iPoint3D voxel);

    dPoint3D computeCOM(int iVolume);
    dPoint3D computeCOMIndex(int iVolume);
    dPoint3D computeCOMIndex(int iVolume, bool inverse);
    dPoint3D computeCOMStdev(int iVolume, bool inverse, dPoint3D COM);
    d2Vector createProjectionInXYPlane(dPoint2D COM, int iZ, int angleDeg, double radiusStep, d2Vector &coordinates);
    dPoint2D getMinAndMax(int iTime);
    void fillHistogram(int iTime, double lowerLimit, double upperLimit, dVector &histogram);
    void fillHistogram(int iTime, double lowerLimit, double upperLimit, bitMap overlay, dVector &histogram);
    double getAverageSignalAfterThreshold(int iTime, double peakFraction);
    void multiplyConstant(double multiplier);
    void multiplyConstant(int iTime, double multiplier);
    void addImageData(ImageIO *addFile, int iTime, int iTimeAdd);
};

class bitMap : public ImageIO
{
private:
    VoxelSet _voxelSet;             // (this is a list of voxels, not a volume)
    bool _colorInVoxelSet[NCOLORS]; // keep track of which colors are used to speed display of color overlay maps
    int _overlayInputLogic=0;  // 0=OR, 1=AND, 2=NOT

    inline bool voxelSurrounded(fVoxel voxel) {iPoint3D voxel3d={voxel.x,voxel.y,voxel.z}; return voxelSurrounded(voxel3d);}
    bool voxelSurrounded(iPoint3D voxel);

public:
    void defineBitMap(iPoint3D dim);
    void defineBitMap(iPoint3D dim, dPoint3D res);
    void defineBitMap(iPoint3D dim, dMatrix ijkToXyz, dMatrix xyzToIjk);
    void destroyBitMap();
    void clearBitMap();
    void setAllBits();

    inline void setOverlayLogicToOR()  {_overlayInputLogic=0;}
    inline void setOverlayLogicToAND() {_overlayInputLogic=1;}
    inline void setOverlayLogicToNOT() {_overlayInputLogic=2;}
    int readOverlayFile(QString fileName, overlayColor color);
    int readBitmapFile(QString fileName, overlayColor color );
    int writeOverlayFile(QString fileName);
    int writeBitmapFile(QString fileName, bool useCompression);
    void setOrAppendBitMap( bitMap overlay, bool clear);
    void setOverlayData(OVLFile ovl, bool clear);
    void setRectangularOverlay( iPoint3D center, iPoint3D width, bool clear);
    void mirrorOverlay();
    void moveOverlay(iPoint3D move);
    void growOverlayField(const iPoint3D seedVoxel, ImageIO *colorFile, double threshold);
    void clusterOverlayField(const iPoint3D seedVoxel);
    iPoint3D getNearestNonzeroVoxel(const iPoint3D seedVoxel);
    void setBitMapVoxelIfTouchingConnectedNeighbor(const fVoxel seedVoxel);
    void dilateOverlayUsingThreshold(ImageIO *file, double thresholdRatio,
                                     bool aboveThreshold, bool twoDimensional, bool fillSurroundedVoxels);
    void fillSurroundedVoxels();
    void medianFilterOverlay(bool twoDimensional, int iHalfWidth);
    void dilateOverlay(bool twoDimensional);
    void erodeOverlay(bool medianFilter, bool twoDimensional);
    inline void dilateOverlay() {dilateOverlay(false);}
    inline void erodeOverlay()  {erodeOverlay(false,false);}
    inline void openOverlay()  {FUNC_ENTER; erodeOverlay(false,false); dilateOverlay(false);}
    inline void closeOverlay() {dilateOverlay(false); erodeOverlay(false,false);}
    inline void openOverlay(bool twoDimensional)  {erodeOverlay(false,twoDimensional); dilateOverlay(twoDimensional);}
    inline void closeOverlay(bool twoDimensional) {dilateOverlay(twoDimensional); erodeOverlay(false,twoDimensional);}

    void outerShell(bitMap source);
    bitMap findPerimeter();
    bitMap constraintField(double cutOffVoxels);

    void selectOverlayFromThresholdGrayScale(ImageIO *signalFile, bool above, double value, bool binary);
    void selectOverlayFromThresholdColorScale(ImageIO *signalFile, bool above, double value);
    void applyThresholds(bool reset, ImageIO *sourceFile, int iTime, double lowThresh, double highThresh);

    // Setters
    void setOverlayDimensions(const iPoint3D voxel);
    void setOverlayDimensions(const int nx, const int ny, const int nz);
    void setBitMapValue(const int i3d, const float value, bool updateVoxelSet);
    void setBitMapValue(const int iX, const int iY, const int iZ, const float value, overlayColor color, bool updateVoxelSet);
    void setBitMapValue(const iPoint3D voxel, const float value, overlayColor color, bool updateVoxelSet);
    void setBitMapValue(const iPoint3D voxel, const float value, bool updateVoxelSet);
    void setBitMapValue(fVoxel ovl, bool updateVoxelSet);
    void setBitMapValue(fVoxelLight voxelLight, bool updateVoxelSet);
    void setBitMapValueBroadBrush(fVoxel voxel, bool updateVoxelSet);
    void setBitMapValueBroadBrush(const iPoint3D voxel3d, const float value, overlayColor color, bool updateVoxelSet);
    void updateBitMapVoxelSet(overlayColor inputColor, iPoint2D rangeX, iPoint2D rangeY, iPoint2D rangeZ);
    void updateBitMapVoxelSet(overlayColor inputColor);
    inline void updateBitMapVoxelSet() {updateBitMapVoxelSet(Color_Undefined);}
    inline void updateBitMapVoxelSet(iPoint2D rangeX, iPoint2D rangeY, iPoint2D rangeZ)
    { updateBitMapVoxelSet(Color_Undefined, rangeX, rangeY, rangeZ);}
    inline void setColorUsedInOverlay(int iColor, bool state) {_colorInVoxelSet[iColor] = state;}

    // Getters
    inline iPoint3D getDimensions() const {return getImageDimensions();}
    inline dPoint3D getResolution() const {return getImageResolution();}
    inline bool someBitsAreSet() {return (_voxelSet.size() != 0);}
    inline int getBitsSetInBitMap() {return _voxelSet.size();}
    inline VoxelSet getVoxelSet() {return _voxelSet;}
    inline fVoxel getBitMapVoxel(int iVoxel) {return _voxelSet[iVoxel];}
    inline iPoint3D getBitMapVoxel3D(int iVoxel) {fVoxel voxel=_voxelSet[iVoxel]; return {voxel.x,voxel.y,voxel.z};}
    double getNumberVoxelsInOverlay();
    double getVolumeInOverlay();
    overlayColor getBitMapColor(int iX, int iY, int iZ);
    overlayColor getBitMapColor(iPoint3D voxel);
    overlayColor getBitMapColor(fVoxel voxel);
    inline bool colorIsUsedInOverlay(int iColor) {return _colorInVoxelSet[iColor];}
    double getBitMapValue(int i3d);
    double getBitMapValue(int iX, int iY, int iZ);
    double getBitMapValue(iPoint3D voxel);
    double getBitMapValue(fVoxel voxel);
    inline bool isVoxelNonZero(int iVoxel)          {return getBitMapValue(iVoxel)!= 0.f;}
    inline bool isVoxelNonZero(iPoint3D voxel)      {return getBitMapValue(voxel) != 0.f;}
    inline bool isVoxelZero(int iVoxel)             {return getBitMapValue(iVoxel)== 0.f;}
    inline bool isVoxelZero(iPoint3D voxel)         {return getBitMapValue(voxel) == 0.f;}
    inline bool isVoxelAbovePointP5(iPoint3D voxel) {return getBitMapValue(voxel) > 0.5;}
    inline bool isVoxelAbovePointP5(int iVoxel)     {return getBitMapValue(iVoxel) > 0.5;}
    inline bool isVoxelUnity(iPoint3D voxel)        {return getBitMapValue(voxel)  == 1.;}
    inline bool isVoxelUnity(int iVoxel)            {return getBitMapValue(iVoxel) == 1.;}
    double getMedianSignal(ImageIO *file, bool aboveThreshold, double thresholdRatio);
    iPoint3D computeCOMIndex();
    d3Vector computeCOMs();
};

class voxelAverager : public QObject, public QRunnable
{
    Q_OBJECT   // this can create a problem of "missing vtable"; fix by "run qmake" from Qt Creator build menu
public:
    voxelAverager(bVector includeVolume, ImageIO *sourceFile);
    void run(); // if public, it can be run directly in the same thread to test thread safety
private:
    ImageIO _avFile;
//    int _nThreads;
    bVector _includeVolume;
    int _nAveraged;
    ImageIO _sourceFile;
signals:
    void progressVoxelAverager(int iProgress);
    void finishedVoxelAverager(ImageIO avFile);
public slots:
};

class overlayReader : public QObject, public QRunnable
{
    Q_OBJECT
public:
    OVLFile _file;
    overlayReader(OVLFile overlayFile, iPoint3D dimensions, bool currentSpaceShown);
private:
    void run();
    iPoint3D _dimensions={0,0,0};
    bool _currentSpaceShown;
signals:
    void successfullyReadOverlay(OVLFile imageFile);
    void errorReadingOverlay(QString errorText);
public slots:
};

class voxelReader : public QObject, public QRunnable
{
    Q_OBJECT
public:
    ImageIO _imageFile;
    voxelReader(ImageIO imageFile);
    //    int readImageData(QString imageFileName, imageHeader header);
private:
    void run();
    int machine_byte_order();
    void bracketGoodSignal(imageHeader *hdr);
signals:
    void finishedReadingVoxels(ImageIO imageFile);
    void readProgress(int);
public slots:
    //    void readImageData();
};

class voxelSlicer : public QObject , public QRunnable
{
    Q_OBJECT   // this can create a problem of "missing vtable"; fix by "run qmake" from Qt Creator build menu
public:
    voxelSlicer(int firstSlice, int lastSlice, ImageIO imageOriginal, ImageIO imageReslice);
    void run();
private:
    int _iFileType;          // 0 = anatomy, 1 = atlas, 2 = color
    int _firstSlice;         // 1st slice to map
    int _lastSlice;          // last slice to map
    ImageIO  _imageOriginal; // original data that will be resliced
    ImageIO _imageReslice;   // the image data to modify; header must be correct already at initialization
signals:
    void finishedSlice(int whichSlice);
    void finished(int firstSlice, int lastSlice, ImageIO imageFile);
public slots:
};

#endif /* _NIFTI_IO_HEADER_ */
