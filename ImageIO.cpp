#include <QDebug>
#include <QBuffer>
#include <QtMath>
#include <QMutex>
#include <QtDebug>
#include <QMessageBox>
#include <QRegularExpression>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "ImageIO.h"
#include "nifti1.h"                  /*** NIFTI-1 header specification ***/
#include "znzlib.h"

//using namespace utilIO;

voxelAverager::voxelAverager(bVector includeVolume, ImageIO *sourceFile)
{
    FUNC_ENTER << includeVolume;
//    _nThreads = nThreads;
    _includeVolume = includeVolume;
    _sourceFile = *sourceFile;
    _avFile.copyTemplateDimensions(sourceFile,1);
    if ( includeVolume.size() != _sourceFile.nt() )
    {
        qInfo() << "Programming error 1 in voxelAverager::voxelAverager" << includeVolume.size() << _sourceFile.nt();
        exit(1);
    }
    _nAveraged = 0;
    for ( int jt=0; jt<_sourceFile.nt(); jt++)
        if ( _includeVolume[jt] ) _nAveraged++;
    if ( _nAveraged == 0 )
    {
        qInfo() << "Programming error 2 in voxelAverager::voxelAverager";
        exit(1);
    }
    FUNC_EXIT << _nAveraged;
}
void voxelAverager::run()
{
    FUNC_ENTER;
    int iDiv = _avFile.voxelsTotal()/10 + 1;  // so the total number of outputs will be about 10
    FUNC_INFO << "voxels" << _avFile.voxelsTotal() << "iDiv" << iDiv;
    for ( int jVoxel=0; jVoxel<_avFile.voxelsPerVolume(); ++jVoxel)
    {
        _avFile.setStackValue(jVoxel,0,0.);
        double av=0.;
        for (int jt=0; jt<_sourceFile.nt(); jt++)
        {
            if ( _includeVolume[jt] )
                av += static_cast<float>(_sourceFile.getStackValue(jVoxel,jt));
        }
        av /= static_cast<float>(_nAveraged);
        _avFile.setStackValue(jVoxel,0,av);
        if ( jVoxel%iDiv == 0 ) emit progressVoxelAverager(0);
    } // loop over voxels
    emit progressVoxelAverager(1);
    emit finishedVoxelAverager(_avFile);
    FUNC_EXIT;
}
/*
voxelAverager::voxelAverager(int nThreads, VoxelSet voxelList, bVector includeVolume, ImageIO *sourceFile)
{
    FUNC_ENTER << nThreads;
    _nThreads = nThreads;
    _voxelList = voxelList;  // list of voxels to map for this thread
    _includeVolume = includeVolume;
    _sourceFile = *sourceFile;
    _nTime = sourceFile->nt();
    if ( includeVolume.size() != sourceFile->nt() )
    {
        qInfo() << "Programming error 1 in voxelAverager::voxelAverager" << includeVolume.size() << _nTime;
        exit(1);
    }
    _nAveraged = 0;
    for ( int jt=0; jt<_nTime; jt++)
        if ( _includeVolume[jt] ) _nAveraged++;
    if ( _nAveraged == 0 )
    {
        qInfo() << "Programming error 2 in voxelAverager::voxelAverager";
        exit(1);
    }
    FUNC_EXIT << _nAveraged;
}
void voxelAverager::run()
{
    FUNC_ENTER << _voxelList.size() << _nThreads;
    int iDiv = (_voxelList.size()*_nThreads)/50 + 1;  // so the total number of outputs will be about 50
    for ( int jVoxel=0; jVoxel<_voxelList.size(); ++jVoxel)
    {
        _voxelList[jVoxel].binaryValue = 0;
        for (int jt=0; jt<_nTime; jt++)
        {
            if ( _includeVolume[jt] )
                _voxelList[jVoxel].binaryValue += static_cast<float>(_sourceFile.getStackValue(_voxelList[jVoxel],jt));
        }
        _voxelList[jVoxel].binaryValue /= static_cast<float>(_nAveraged);
        if ( jVoxel%iDiv == 0 ) emit progressVoxelAverager(0);
    } // loop over voxels
    FUNC_EXIT << _voxelList.size();
    emit finishedVoxelAverager(_voxelList);
}
*/
voxelSlicer::voxelSlicer(int firstSlice, int lastSlice, ImageIO imageOriginal, ImageIO imageReslice)
{
    _iFileType     = imageReslice.getCategory();
    _firstSlice    = firstSlice;
    _lastSlice     = lastSlice;
    _imageOriginal = imageOriginal;
    _imageReslice  = imageReslice;
}

void voxelSlicer::run()
{
    // _imageReslice contains the correct (resliced) header and dimensions already
    // _imageOriginal contains the original header and data
    // transformation: reslice space --> original space, then interpolate

    dPoint3D fwhm;
    fwhm.x = qCeil(_imageReslice.dx() / _imageOriginal.dx());
    fwhm.y = qCeil(_imageReslice.dy() / _imageOriginal.dy());
    fwhm.z = qCeil(_imageReslice.dz() / _imageOriginal.dz());
    iPoint3D iWidth;
    iWidth.x = qCeil(fwhm.x);
    iWidth.y = qCeil(fwhm.y);
    iWidth.z = qCeil(fwhm.z);
    iPoint3D kernelType, wrap;
    wrap.x = wrap.y = wrap.z = 0; // no wrapping
    enum kernels {Lanczos, Gaussian, Delta};
    if ( _iFileType == category_atlas )
        kernelType.x = kernelType.y = kernelType.z = Delta;  // delta
    else
    {  // Gaussian if fwhm > 1; otherwise Lanczos
        kernelType.x = (fwhm.x>1);  kernelType.y = (fwhm.y>1); kernelType.z = (fwhm.z>1);
    }
    //    kernelType.x = kernelType.y = kernelType.z = Gaussian;

    // Dummy statement to print info about interpolation (indicated by time point "-1")
    // Space_template.x = Space_template.y = Space_template.z = 0.;

    iPoint4D dim0 =  _imageReslice.getDimensions();
    bool inVolume;
    iPoint3D voxel;
    // Loop over voxels in the template space
    for (voxel.z=_firstSlice; voxel.z<=qMin(_lastSlice,dim0.z-1); voxel.z++)
    {
        for (voxel.y=0; voxel.y<_imageReslice.ny(); voxel.y++)
        {
            for (voxel.x=0; voxel.x<_imageReslice.nx(); voxel.x++)
            {
                // Convert from data indices --> space in the resliced space
                dPoint3D space = _imageReslice.convertVoxelToSpace(voxel);
                // Loop over all time points; hopefully this is usually just 1
                for (int jt=0; jt<_imageReslice.nt(); jt++)
                {
                    // Interpolate to the find the volume value at this spatial coordinate
                    double value = _imageOriginal.interpolateValueFromSpaceCoord(space,fwhm,iWidth,wrap,
                                                                                 kernelType,jt,inVolume);
                    _imageReslice.setStackValue(voxel,jt,value);
                } // x
            } // y
        } // t
        emit finishedSlice(voxel.z);
    } // z
    emit finished(_firstSlice, qMin(_lastSlice,dim0.z-1), _imageReslice);
}

overlayReader::overlayReader(OVLFile overlayFile, iPoint3D dimensions, bool currentSpaceShown)
{
    _file = overlayFile;
    _dimensions = dimensions;
    _currentSpaceShown = currentSpaceShown;
}

void overlayReader::run()
{
    FUNC_ENTER;
    // start with a clear voxel vector
    _file.voxelList.clear();

    // Read the overlay file
    QFile file(_file.path);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        QString errorText = "Error from overlayReader::run - attempting to open overlay file " + _file.path;
        emit errorReadingOverlay(errorText);
    }
    QTextStream in_stream(&file);

    fVoxel voxel;
    int iLine=1;
    QRegularExpression whiteSpaceComma("[,\\s]");// match a comma or a space
    while (!in_stream.atEnd())
    {
        QString line = in_stream.readLine();
        QString unCommented = line.left(line.indexOf("#"));
        if ( unCommented.isEmpty() )
            continue;
        QStringList valueList = unCommented.split(whiteSpaceComma,Qt::SkipEmptyParts); // clazy:exclude=use-static-qregularexpression
        if ( !valueList.at(0).compare("dimensions",Qt::CaseInsensitive) ||
             !valueList.at(0).compare("dimension",Qt::CaseInsensitive)  ||
             !valueList.at(0).compare("matrix",Qt::CaseInsensitive) )
        { // typically this should be on the 1st line, if used
            iPoint3D dim;
            dim.x = valueList.at(1).toInt();
            dim.y = valueList.at(2).toInt();
            dim.z = valueList.at(3).toInt();
            FUNC_INFO << _dimensions.x << _dimensions.y << _dimensions.z << dim.x << dim.y << dim.z;
            if ( _dimensions.x != 0 && (dim.x != _dimensions.x || dim.y != _dimensions.y || dim.z != _dimensions.z) )
            {
//                QString error = QString("Error: dimensions of overlay '%1' do not match dimensions of image volume.\n").arg(_file.name);
                QString error = "";
                emit errorReadingOverlay(error);
                emit successfullyReadOverlay(_file);
                return;
            }
        }
        else
        {
            voxel.x = valueList.at(0).toInt();
            voxel.y = valueList.at(1).toInt();
            voxel.z = valueList.at(2).toInt();
            if ( valueList.size() < 3 )
            {
                QString errorText = QString("Error: only read %1 values on line %2.").arg(valueList.size()).arg(iLine);
                emit errorReadingOverlay(errorText);
            }
            if ( valueList.size() == 3 )
                voxel.binaryValue = 1.;
            else
                voxel.binaryValue = valueList.at(3).toFloat();
            voxel.color = _file.color;  // save the color in the OVLFile structure
            _file.voxelList.append(voxel);
        }
        iLine++;
    }
    file.close();
    if ( (_dimensions.x * _dimensions.y * _dimensions.z == 0) && !_currentSpaceShown )
        emit errorReadingOverlay("");
    FUNC_EXIT << _dimensions.x << _dimensions.y << _dimensions.z;
    emit successfullyReadOverlay(_file);
}

voxelReader::voxelReader(ImageIO imageFile)
{
    _imageFile = imageFile;
}
/*
void voxelReader::run()
{
    imageHeader *hdr = _imageFile.getImageDataHeader();

    int swap_bytes  = utilIO::machine_byte_order() != hdr->byte_order;

    QFile file(hdr->fileName);
    if ( !file.open(QIODevice::ReadOnly) )
    {
        qInfo() << "Error: cannot open file " << hdr->fileName;
        return;
    }
    QDataStream in(&file);
    in.setByteOrder(QDataStream::LittleEndian);
    in.setFloatingPointPrecision(QDataStream::SinglePrecision);

    in.skipRawData(hdr->voxelOffset);
    if ( in.status() != QDataStream::Ok )
    {
        qInfo() << "Error reading the header for file" << hdr->fileName;
        if ( in.status() != QDataStream::ReadPastEnd )
            qInfo() << "read past end";
        else if ( in.status() != QDataStream::ReadCorruptData )
            qInfo() << "read corrupt data";
        return;
    }

    // Go to the first data voxel
    fMatrix fStack;
    fStack.resize(hdr->dim.t);

    // Read the array.
    for (int jt=0; jt<hdr->dim.t; jt++)
    {
        // Copy the integer array into the stack floating point array.
        fStack[jt].resize(hdr->points.z);

        if ( hdr->storageType == BFLOAT )
        {
            QVector<float> fData;
            fData.resize(hdr->points.z);
            for ( int j=0; j<fData.size(); j++)
                in >> fData[j];
            for (int jp=0; jp<hdr->points.z; jp++)
                fStack[jt][jp] = fData[jp];
        }
        else if ( hdr->storageType == BSHORT )
        {
            QVector<qint16> iData;
            iData.resize(hdr->points.z);
            for ( int j=0; j<iData.size(); j++)
                in >> iData[j];
            for (int jp=0; jp<hdr->points.z; jp++)
                fStack[jt][jp] = iData[jp];
        }
        else if ( hdr->storageType == BLONG )
        {
            QVector<qint32> iData;
            iData.resize(hdr->points.z);
            for ( int j=0; j<iData.size(); j++)
                in >> iData[j];
            for (int jp=0; jp<hdr->points.z; jp++)
                fStack[jt][jp] = iData[jp];
        }
        else if ( hdr->storageType == USHORT )
        {
            QVector<qint16> iData;
            iData.resize(hdr->points.z);
            for ( int j=0; j<iData.size(); j++)
                in >> iData[j];
            for (int jp=0; jp<hdr->points.z; jp++)
                fStack[jt][jp] = iData[jp];
        }
        else if ( hdr->storageType == UINT8 )
        {
            QVector<qint8> iData;
            iData.resize(hdr->points.z);
            for ( int j=0; j<iData.size(); j++)
                in >> iData[j];
            for (int jp=0; jp<hdr->points.z; jp++)
                fStack[jt][jp] = iData[jp];
        }
        if( in.status() != QDataStream::Ok )
        {
            qInfo() << "Error reading file" << hdr->fileName << "for time point" << jt;
            return;
        }
        emit readProgress(0);
    }

    // Scale the data if necessary.
    if ( hdr->scaleSlope != 0.f )
    {
        for (int jt=0; jt<hdr->dim.t; jt++)
        {
            for (int jp=0; jp<hdr->points.z; jp++)
            {
                fStack[jt][jp] *= hdr->scaleSlope;
                fStack[jt][jp] += hdr->scaleOffset;
            }
        }
        qInfo() << "scaled by" << hdr->scaleSlope;
    }
    file.close();

    fStack.swap(_imageFile._data.timePoints);
    emit finishedReadingVoxels(_imageFile);
}
void voxelReader::run()
{
    float *fData;
    short *iData;
    int *iiData;
    unsigned char *i8Data;
    unsigned short *iuData;

    FUNC_ENTER;

    imageHeader *hdr = _imageFile.getImageDataHeader();

    int swap_bytes  = utilIO::machine_byte_order() != hdr->byte_order;
    int output_info = hdr->points.t > 1000000 && hdr->dim.t > 10;

    // Open file file.
    QFile file(hdr->fileName);
    if ( !file.open(QIODevice::ReadOnly) )
    {
        qInfo() << "Error: cannot open file " << hdr->fileName;
        return;
    }
    // Go to the first data voxel
    if ( !file.seek(hdr->voxelOffset) )
//    ulong error = fseek(file, (long)(hdr->voxelOffset), SEEK_SET);
//    if ( error )
    {
        qInfo() << "Error stepping past header into data for file " << hdr->fileName;
        return;
    }

    // allocate one volume
    qint64 nAllocate = hdr->points.z;
    if ( hdr->dataType == 1 || hdr->dataType == 2 )
        nAllocate *= 2;  // complex
    else if ( hdr->dataType == 3 )
        nAllocate *= 3;  // 3-vector

    fData = allocate_vector<float>(nAllocate);
    if ( hdr->storageType == BSHORT ) // bshort
        iData = allocate_vector<short>(nAllocate);
    else if ( hdr->storageType == BLONG )  // blong
        iiData = allocate_vector<int>(nAllocate);
    else if ( hdr->storageType == USHORT ) // unsigned short
        iuData = allocate_vector<unsigned short>(nAllocate);
    else if ( hdr->storageType == UINT8 ) // unsigned char
        i8Data = allocate_vector<unsigned char>(nAllocate);

    qint64 nRead;
    // Read the array.
    for (int jt=0; jt<hdr->dim.t; jt++)
    {
        fVector dataForOneVolume;
        dataForOneVolume.resize(nAllocate);
        emit readProgress(0);
        //      if ( output_info && jt%(_data.hdr->dim.t/10) == 0 )
        //      ipercent = MIN(100,nearbyint((float)(jt*100)/(float)_data.hdr->dim.t));
        bool error;
        if ( hdr->storageType == BFLOAT )
        {
            qint64 nBytesToRead = nAllocate*sizeof(float);
            char *data = allocate_vector<char>(nBytesToRead);
            nRead = file.read(data,nBytesToRead);
            error = (nRead != nBytesToRead);
            nRead /= sizeof(float);  // convert to data words
            fData = reinterpret_cast<float *>(data);
//            error = fread((char *)fData,sizeof(float),nAllocate,file) != nAllocate;
            if ( swap_bytes )
            {
                for (qint64 j=0; j<nRead; j++)
                {
                    utilIO::swap( &(fData[j]) );
                    dataForOneVolume[j] = fData[j];
                }
            }
            else
                for (qint64 j=0; j<nRead; j++)
                    dataForOneVolume[j] = fData[j];
        }
        else if ( hdr->storageType == BSHORT )
        {
            qint64 nBytesToRead = nAllocate*sizeof(short);
            char *data = allocate_vector<char>(nBytesToRead);
            nRead = file.read(data,nBytesToRead);
            error = (nRead != nBytesToRead);
            nRead /= sizeof(short);  // convert to data words
            iData = reinterpret_cast<short *>(data);
//            error = fread((char *)iData,sizeof(short),nAllocate,file) != nAllocate;
            if ( swap_bytes )
            {
                for (qint64 j=0; j<nRead; j++)
                {
                    utilIO::swap( &(iData[j]) );
                    dataForOneVolume[j] = iData[j];
                }
            }
            else
                for (qint64 j=0; j<nRead; j++)
                    dataForOneVolume[j] = iData[j];
        }
        else if ( hdr->storageType == USHORT ) // unsigned short
        {
            qint64 nBytesToRead = nAllocate*sizeof(unsigned short);
            char *data = allocate_vector<char>(nBytesToRead);
            nRead = file.read(data,nBytesToRead);
            error = (nRead != nBytesToRead);
            nRead /= sizeof(unsigned short);  // convert to data words
            iuData = reinterpret_cast<unsigned short *>(data);
//            error = fread((char *)iuData,sizeof(unsigned short),nAllocate,file) != nAllocate;
            if ( swap_bytes )
            {
                for (qint64 j=0; j<nRead; j++)
                {
                    utilIO::swap( &(iuData[j]) );
                    dataForOneVolume[j] = iuData[j];
                }
            }
            else
                for (qint64 j=0; j<nRead; j++)
                    dataForOneVolume[j] = iuData[j];
        }
        else if ( hdr->storageType == UINT8 ) // unsigned chr
        {
            qint64 nBytesToRead = nAllocate*sizeof(unsigned char);
            char *data = allocate_vector<char>(nBytesToRead);
            nRead = file.read(data,nBytesToRead);
            error = (nRead != nBytesToRead);
            i8Data = reinterpret_cast<unsigned char *>(data);
//            error = fread((char *)i8Data,sizeof(unsigned char),nAllocate,file) != nAllocate;
            for (qint64 j=0; j<nAllocate; j++)
                dataForOneVolume[j] = i8Data[j];
        }
        else // if ( hdr->storageType == BLONG )  // blong
        {
            qint64 nBytesToRead = nAllocate*sizeof(int);
            char *data = allocate_vector<char>(nBytesToRead);
            nRead = file.read(data,nBytesToRead);
            error = (nRead != nBytesToRead);
            nRead /= sizeof(int);  // convert to data words
            iiData = reinterpret_cast<int *>(data);
//            error = fread((char *)iiData,sizeof(int),nAllocate,file) != nAllocate;
            if ( swap_bytes )
            {
                for (qint64 j=0; j<nRead; j++)
                {
                    utilIO::swap( &(iiData[j]) );
                    dataForOneVolume[j] = iiData[j];
                }
            }
            else
                for (qint64 j=0; j<nRead; j++)
                    dataForOneVolume[j] = iiData[j];
        }
        if ( error )
        {
            int nt = jt;
            qInfo() << "** Error reading time point " << jt+1 << " for file " << hdr->fileName;
            _imageFile._data.timePoints.resize(nt);
            _imageFile._data.hdr.dim.t = nt;
            _imageFile._data.hdr.points.t = _imageFile.nx() * _imageFile.ny() * _imageFile.nz() * nt;
            break;
        }
        else
        {
            // Scale the data if necessary.
            if ( hdr->scaleSlope != 0. )
            {
                for (qint64 j=0; j<nRead; j++)
                {
                    dataForOneVolume[j] *= hdr->scaleSlope;
                    dataForOneVolume[j] += hdr->scaleOffset;
                }
            }
            FUNC_INFO << "dataType" << hdr->dataType;
            // Copy the volume into the time series
            if ( hdr->dataType == 0 ) // scalar
                dataForOneVolume.swap(_imageFile._data.timePoints[jt].f1);
            else if ( hdr->dataType == 3 ) // 3-vector
            {
                for (qint64 j=0; j<nRead/3; j++)
                {
                    _imageFile._data.timePoints[jt].f3[j].x = dataForOneVolume[3*j];
                    _imageFile._data.timePoints[jt].f3[j].y = dataForOneVolume[3*j+1];
                    _imageFile._data.timePoints[jt].f3[j].z = dataForOneVolume[3*j+2];
                }
            }
            else // if ( hdr->dataType == 1 || hdr->dataType == 2 ) // complex
            {
                for (qint64 j=0; j<nRead/2; j++)
                {
                    _imageFile._data.timePoints[jt].fc[j].real = dataForOneVolume[2*j];
                    _imageFile._data.timePoints[jt].fc[j].imag = dataForOneVolume[2*j+1];
                }
            }
        }
    } // jt
    file.close();
//    fclose(file);

    if ( output_info )
    {
        //      printf("\r%s: 100 %c\n",FileName,37);
        //      fflush(stdout);
    }
    emit finishedReadingVoxels(_imageFile);
    FUNC_EXIT;
}
*/

void voxelReader::run()
{
    double *f64Data;
    float *fData;
    short *iData;
    int *iiData;
    unsigned char *i8Data;
    unsigned short *iuData;

    FUNC_ENTER;

    imageHeader *hdr = _imageFile.getImageDataHeader();
    FUNC_INFO << "dataType" << hdr->dataType;
    FUNC_INFO << "storageType" << hdr->storageType;

    int swap_bytes  = utilIO::machine_byte_order() != hdr->byte_order;
    int output_info = hdr->points.t > 1000000 && hdr->dim.t > 10;

    // Open file file.
    znzFile file;
    FUNC_INFO << "** read file" << hdr->fileName;
    int useCompression = hdr->fileType == NIFTI_NII_GZ;
    QByteArray array = hdr->fileName.toLocal8Bit();
    char* fileName = array.data();
    if ( !(file = znzopen(fileName,"r",useCompression)) )
    {
        qInfo() << "Error: cannot open file " << hdr->fileName;
        return;
    }
    // Go to the first data voxel
    FUNC_INFO << "seek" << hdr->voxelOffset;
    znzseek(file, static_cast<long>(hdr->voxelOffset), SEEK_SET) != hdr->voxelOffset;

    // allocate one volume
    unsigned long nAllocate = hdr->points.z;
    if ( hdr->dataType == 1 || hdr->dataType == 2 )
        nAllocate *= 2;  // complex
    else if ( hdr->dataType == 3 || hdr->dataType == 4 )
        nAllocate *= 3;  // 3-vector

    fData = allocate_vector<float>(nAllocate);
    if ( hdr->storageType == BSHORT ) // bshort
        iData = allocate_vector<short>(nAllocate);
    else if ( hdr->storageType == BLONG )  // blong
        iiData = allocate_vector<int>(nAllocate);
    else if ( hdr->storageType == USHORT ) // unsigned short
        iuData = allocate_vector<unsigned short>(nAllocate);
    else if ( hdr->storageType == UINT8 ) // unsigned char
        i8Data = allocate_vector<unsigned char>(nAllocate);
    else // if ( hdr->storageType == BDOUBLE ) // bytes
        f64Data = allocate_vector<double>(nAllocate);

    // Read the array.
    for (int jt=0; jt<hdr->dim.t; jt++)
    {
        emit readProgress(0);
        //      if ( output_info && jt%(_data.hdr->dim.t/10) == 0 )
        //      ipercent = MIN(100,nearbyint((float)(jt*100)/(float)_data.hdr->dim.t));
        bool error;
        if ( hdr->storageType == BFLOAT )
        {
            FUNC_INFO << "read float";
            error = znzread((char *)fData,sizeof(float),nAllocate,file) != nAllocate;
            if ( swap_bytes )
            {
                for (unsigned long j=0; j<nAllocate; j++)
                    utilIO::swap( &(fData[j]) );
            }
        }
        else if ( hdr->storageType == BSHORT )
        {
            error = znzread((char *)iData,sizeof(short),nAllocate,file) != nAllocate;
            if ( swap_bytes )
            {
                for (unsigned long j=0; j<nAllocate; j++)
                {
                    utilIO::swap( &(iData[j]) );
                    fData[j] = iData[j];
                }
            }
            else
                for (unsigned long j=0; j<nAllocate; j++)
                    fData[j] = iData[j];
        }
        else if ( hdr->storageType == USHORT ) // unsigned short
        {
            error = znzread((char *)iuData,sizeof(unsigned short),nAllocate,file) != nAllocate;
            if ( swap_bytes )
            {
                for (unsigned long j=0; j<nAllocate; j++)
                {
                    utilIO::swap( &(iuData[j]) );
                    fData[j] = iuData[j];
                }
            }
            else
                for (unsigned long j=0; j<nAllocate; j++)
                    fData[j] = iuData[j];
        }
        else if ( hdr->storageType == UINT8 ) // unsigned chr
        {
            error = znzread((char *)i8Data,sizeof(unsigned char),nAllocate,file) != nAllocate;
            for (unsigned long j=0; j<nAllocate; j++)
                fData[j] = i8Data[j];
        }
        else if ( hdr->storageType == BDOUBLE ) // unsigned char
        {
            error = znzread((char *)f64Data,sizeof(double),nAllocate,file) != nAllocate;
            for (unsigned long j=0; j<nAllocate; j++)
                fData[j] = f64Data[j];
        }
        else // if ( hdr->storageType == BLONG )  // blong
        {
            error = znzread((char *)iiData,sizeof(int),nAllocate,file) != nAllocate;
            if ( swap_bytes )
            {
                for (unsigned long j=0; j<nAllocate; j++)
                {
                    utilIO::swap( &(iiData[j]) );
                    fData[j] = iiData[j];
                }
            }
            else
                for (unsigned long j=0; j<nAllocate; j++)
                    fData[j] = iiData[j];
        }
        if ( error )
        {
            int nt = jt;
            qInfo() << "** Error reading time point " << jt+1 << " for file " << hdr->fileName;
            _imageFile._data.timePoints.resize(nt);
            _imageFile._data.hdr.dim.t = nt;
            _imageFile._data.hdr.points.t = _imageFile.nx() * _imageFile.ny() * _imageFile.nz() * nt;
            break;
        }
        else
        {
            // Scale the data if necessary.
            FUNC_INFO << "slope" << hdr->scaleSlope << "offset" << hdr->scaleOffset;
            if ( hdr->scaleSlope == 0. ) hdr->scaleSlope = 1.;
            for (unsigned long j=0; j<nAllocate; j++)
            {
                fData[j] *= hdr->scaleSlope;
                fData[j] += hdr->scaleOffset;
                // scrub data for unusual cases
                if ( qIsNaN(fData[j]) || qIsInf(fData[j]) ) fData[j] = 0.;
            }
            FUNC_INFO << "dataType" << hdr->dataType;
            // Copy the volume into the time series
            if ( hdr->dataType == 0 ) // scalar
            {
                FUNC_INFO << "assign scalar";
                // 4D data structures should already be allocate via readFileHeader
                for (unsigned long j=0; j<nAllocate; j++)
                {
                    _imageFile._data.timePoints[jt].f1[j] = fData[j];
                    if ( hdr->dim.t > 1 ) _imageFile._data.averageVolume[j] += fData[j];
                }
            }
            else if ( hdr->dataType == 3 || hdr->dataType == 4 ) // 3-vector
            {
                // 4D data structures should already be allocate via readFileHeader
                // _data.timePoints[jt].f3.resize(_data.hdr.points.z);
                for (unsigned long j=0; j<nAllocate/3; j++)
                {
                    _imageFile._data.timePoints[jt].f3[j].x = fData[3*j];
                    _imageFile._data.timePoints[jt].f3[j].y = fData[3*j+1];
                    _imageFile._data.timePoints[jt].f3[j].z = fData[3*j+2];
                    if ( hdr->dim.t > 1 )
                        _imageFile._data.averageVolume[j] +=
                                SQR(_imageFile._data.timePoints[jt].f3[j].x) +
                                SQR(_imageFile._data.timePoints[jt].f3[j].y) +
                                SQR(_imageFile._data.timePoints[jt].f3[j].z);
                }
            }
            else // if ( hdr->dataType == 1 || hdr->dataType == 2 ) // complex
            {
                // 4D data structures should already be allocate via readFileHeader
                // _data.timePoints[jt].fc.resize(_data.hdr.points.z);
                for (unsigned long j=0; j<nAllocate/2; j++)
                {
                    _imageFile._data.timePoints[jt].fc[j].real = fData[2*j];
                    _imageFile._data.timePoints[jt].fc[j].imag = fData[2*j+1];
                    if ( hdr->dim.t > 1 )
                        _imageFile._data.averageVolume[j] += SQR(_imageFile._data.timePoints[jt].fc[j].real)
                                +                            SQR(_imageFile._data.timePoints[jt].fc[j].imag);
                }
            }
        }
    } // jt
    if ( hdr->dim.t > 1 )
    {
        for (int j=0; j<hdr->points.z; j++)
            _imageFile._data.averageVolume[j] /= static_cast<double>(hdr->dim.t);
    }
    znzclose(file);
//    fclose(file);

    if ( output_info )
    {
        //      printf("\r%s: 100 %c\n",FileName,37);
        //      fflush(stdout);
    }
    bracketGoodSignal(hdr);
    emit finishedReadingVoxels(_imageFile);
    FUNC_EXIT;
}

void voxelReader::bracketGoodSignal(imageHeader *hdr)
{
    FUNC_ENTER;
    iPoint3D dim = _imageFile.getImageDimensions();
    iPoint3D line;
    // Average 1 line for all time points: x low
    bVector anyZero; anyZero.fill(false,12);
    dVector average; average.fill(0.,12);

    // sweep x
    for (int jLine=0; jLine<4; jLine++)
    {
        if ( jLine == 0 )
            line = {0,0,0};
        else if ( jLine == 1 )
            line = {0,0,dim.z-1};
        else if ( jLine == 2 )
            line = {0,dim.y-1,0};
        else
            line = {0,dim.y-1,dim.z-1};
        int iLine = jLine;
        for (line.x=0; line.x<dim.x; line.x++)
        {
            FUNC_INFO << "test" << line.x << line.y << line.z << "dim" << dim.x << dim.y <<dim.z;
            double value = SQR(_imageFile.getStackValue(line,0));
            anyZero[jLine] |= value == 0.;
            if ( anyZero[jLine] ) break;
            average[iLine] += value;
        } // line.x
        average[iLine] /= static_cast<double>(dim.x);
    } // jLine
    FUNC_INFO << "swept x" << average;

    // sweep y
    for (int jLine=0; jLine<4; jLine++)
    {
        if (jLine ==0)
            line = {0,0,0};
        else if (jLine ==1)
            line = {0,0,dim.z-1};
        else if (jLine ==2)
            line = {dim.x-1,0,0};
        else
            line = {dim.x-1,0,dim.z-1};
        int iLine=jLine+4;
//        for (int jt=0; jt<_imageFile.nt(); jt++)
        for (line.y=0; line.y<dim.y; line.y++)
        {
            double value = SQR(_imageFile.getStackValue(line,0));
            anyZero[iLine] |= value == 0.;
            if ( anyZero[iLine] ) break;
            average[iLine] += value;
        } // line.y
        average[iLine] /= static_cast<double>(dim.y);
    } // jLine
    FUNC_INFO << "swept y" << average;

    // sweep z
    for (int jLine=0; jLine<4; jLine++)
    {
        if (jLine ==0)
            line = {0,0,0};
        else if (jLine ==1)
            line = {0,dim.y-1,0};
        else if (jLine ==2)
            line = {dim.x-1,0,0};
        else
            line = {dim.x-1,dim.y-1,0};
        int iLine=jLine+8;
//        for (int jt=0; jt<_imageFile.nt(); jt++)
        for (line.z=0; line.z<dim.z; line.z++)
        {
            double value = SQR(_imageFile.getStackValue(line,0));
            anyZero[iLine] |= value == 0.;
            if ( anyZero[iLine] ) break;
            average[iLine] += value;
        } // line.
        average[iLine] /= static_cast<double>(dim.z);
    } // jLine
    FUNC_INFO << "swept z" << average;

    // test the lines
    double avMin=0.;
    for (int jLine=0; jLine<average.size(); jLine++)
    {
        FUNC_INFO << "line" << jLine << "has average" << average.at(jLine);
        if ( !anyZero.at(jLine) )
        {
            if ( avMin == 0. )
                avMin = average.at(jLine);
            else
                avMin = qMin(avMin,average.at(jLine));
        }
    }
    hdr->noiseLevel = qSqrt(avMin/2.);  // st dev of noise level

    FUNC_EXIT << "noiseLevel" << hdr->noiseLevel;
}

double ImageIO::calculateCostFunction(CostFunctionChoice costFunction, int iVolume,
                                      ImageIO *templateVolume, int iVolumeTemplate,
                                      bool includeTemplateZeroes)
{
    bitMap mask;
    iPoint3D dim = templateVolume->getImageDimensions();
    dPoint3D res = templateVolume->getImageResolution();
    mask.defineBitMap(dim,res);
    mask.setAllBits(); // include everything
    if ( !includeTemplateZeroes )
    { // exclude template zeroes
        for ( int j=0; j<templateVolume->voxelsPerVolume(); j++ )
        {
            double valueTemplate = templateVolume->getStackValue(j,iVolumeTemplate);
            if ( valueTemplate == 0. ) mask.setBitMapValue(j,0.,false);
        }
    }
    return calculateCostFunction(costFunction, iVolume, templateVolume, iVolumeTemplate, &mask);
}

double ImageIO::calculateCostFunction(CostFunctionChoice costFunction, int iVolume,
                                      ImageIO *templateVolume, int iVolumeTemplate,
                                      bitMap *mask)
{
    FUNC_ENTER;
    double cost = 0;
    if ( costFunction == CostFunction_MI || costFunction == CostFunction_NMI )
    { // Mutual information
        // Find the average value of non-zero voxels for the source and target
        double avTarget   = calculateAverage(iVolume, mask);
        double avTemplate = templateVolume->calculateAverage(iVolumeTemplate, mask);
        if ( avTarget == 0. ) avTarget = 1.;
        FUNC_INFO << "avTarget=" << avTarget << "avTemplate=" << avTemplate;
        // Find the min and max.
        iPoint3D dim = templateVolume->getImageDimensions();
        double min=1.e10;  double max=-1e10;
        for ( int j=0; j<templateVolume->voxelsPerVolume(); j++ )
        {
            if ( mask->isVoxelNonZero(j) )
            {
                // normalize all values by the respective average
                double valueTarget   = getStackValue(j,iVolume)                         / avTarget;
                double valueTemplate = templateVolume->getStackValue(j,iVolumeTemplate) / avTemplate;
                if ( valueTarget   < min ) min = valueTarget;
                if ( valueTemplate < min ) min = valueTemplate;
                if ( valueTarget   > max ) max = valueTarget;
                if ( valueTemplate > max ) max = valueTemplate;
            }
        }
        FUNC_INFO << "min1=" << min << "max1=" << max;
        // Initialize the histograms.
        dVector probI, probT;
        probI.fill(0.,NCOSTBINS);  probT.fill(0.,NCOSTBINS);
        dMatrix probIT;
        probIT.resize(NCOSTBINS);
        for (int jI=0; jI<NCOSTBINS; jI++)
            probIT[jI].fill(0.,NCOSTBINS);
        // Fill the histograms
        dim = templateVolume->getImageDimensions();
        double norm = 0.;
        for ( int j=0; j<templateVolume->voxelsPerVolume(); j++ )
        {
            if ( mask->isVoxelNonZero(j) )
            {
                double valueTarget   = getStackValue(j,iVolume)                         / avTarget;
                double valueTemplate = templateVolume->getStackValue(j,iVolumeTemplate) / avTemplate;

                iPoint2D iBinI, iBinT;
                dPoint2D weightBinI, weightBinT;
                utilMath::fuzzyBinning(valueTarget,    min,  max, NCOSTBINS, iBinI, weightBinI );
                utilMath::fuzzyBinning(valueTemplate,  min,  max, NCOSTBINS, iBinT, weightBinT );
                if ( iBinI.x >= 0 )
                    probI[iBinI.x] += weightBinI.x;
                if ( iBinI.y >= 0 )
                    probI[iBinI.y] += weightBinI.y;

                if ( iBinT.x >= 0 )
                {
                    probT[iBinT.x] += weightBinT.x;
                    if ( iBinI.x >= 0 )
                        probIT[iBinI.x][iBinT.x] += weightBinI.x*weightBinT.x;
                    if ( iBinI.y >= 0 )
                        probIT[iBinI.y][iBinT.x] += weightBinI.y*weightBinT.x;
                }
                if ( iBinT.y >= 0 )
                {
                    probT[iBinT.y] += weightBinT.y;
                    if ( iBinI.x >= 0 )
                        probIT[iBinI.x][iBinT.y] += weightBinI.x*weightBinT.y;
                    if ( iBinI.y >= 0 )
                        probIT[iBinI.y][iBinT.y] += weightBinI.y*weightBinT.y;
                }
                norm += 1.;
            } // nonZero
        } // loop over voxels
        FUNC_INFO << "norm = " << norm;
        // Normalize the probability functions.
        for (int jI=0; jI<NCOSTBINS; jI++)
        {
            probI[jI] /= norm;
            probT[jI] /= norm;
            for (int jT=0; jT<NCOSTBINS; jT++)
                probIT[jI][jT] /= norm;
        }
        if ( costFunction == CostFunction_MI )
        {
            // Calculate the cost function.
            for (int jI=0; jI<NCOSTBINS; jI++)
            { // sum over bins of target
                if ( probI[jI] != 0. )
                {
                    for (int jT=0; jT<NCOSTBINS; jT++)
                    { // sum over bins of target
                        if ( probT[jT] != 0. && probIT[jI][jT] != 0. )
                            cost += probIT[jI][jT] * qLn(probIT[jI][jT]/probI[jI]/probT[jT]);
                    }
                }
            }
        }
        else // normalized MI
        {
            double H_I=0.; double H_T=0.; double H_IT=0.;
            for (int jI=0; jI<NCOSTBINS; jI++)
            { // sum over bins of target
                if ( probI[jI] != 0. )
                    H_I -= probT[jI] * qLn(probI[jI]);
            }
            for (int jT=0; jT<NCOSTBINS; jT++)
            { // sum over bins of target
                if ( probT[jT] != 0. )
                    H_T -= probT[jT] * qLn((probT[jT]));
            }
            for (int jI=0; jI<NCOSTBINS; jI++)
            {
                for (int jT=0; jT<NCOSTBINS; jT++)
                { // sum over bins of target
                    if ( probIT[jI][jT] != 0. )
                        H_IT -= probIT[jI][jT] * qLn((probIT[jI][jT]));
                }
            }
            if ( H_IT != 0. ) cost = H_IT / (H_I + H_T);
        }
    }
    else if ( costFunction == CostFunction_WLS )
    {
        FUNC_INFO << "cost WLS";
        double sumWLS = 0.;
        iPoint3D dim = getImageDimensions();  iPoint3D voxel;
        for (voxel.z=0; voxel.z<dim.z; voxel.z++)
        {
            for (voxel.y=0; voxel.y<dim.y; voxel.y++)
            {
                for (voxel.x=0; voxel.x<dim.x; voxel.x++)
                {
                    if ( mask->isVoxelNonZero(voxel) )
                    {
                        double valueT = templateVolume->getStackValue(voxel,iVolumeTemplate);
                        double valueI = getStackValue(voxel,iVolume);
                        if ( valueI != 0. && valueT != 0. )
                        {
                            sumWLS += SQR(valueT - valueI);
                        }
                    }
                }
            }
        }
        cost = sumWLS;
    }
    else if ( costFunction == CostFunction_RV )
    { // same as WLS but weights = 1
        FUNC_INFO << "cost RV";
        double residue2 = 0.;
        iPoint3D dim = getImageDimensions();  iPoint3D voxel;
        for (voxel.z=0; voxel.z<dim.z; voxel.z++)
        {
            for (voxel.y=0; voxel.y<dim.y; voxel.y++)
            {
                for (voxel.x=0; voxel.x<dim.x; voxel.x++)
                {
                    if ( mask->isVoxelNonZero(voxel) )
                    {
                        double valueT = templateVolume->getStackValue(voxel,iVolumeTemplate);
                        double valueI = getStackValue(voxel,iVolume);
                        if ( valueI != 0. && valueT != 0. )
                            residue2 += SQR(valueI - valueT);
                    }
                }
            }
        }
        cost = residue2 / (dim.x*dim.y*dim.z);  // "cost" = residual variance, not accounting for model parameters
    }
    else if ( costFunction == CostFunction_NC )
    {
        double dotII = 0.; // target ("individual)
        double dotTT = 0.; // template
        double dotIT = 0.; // target . template
        iPoint3D dim = getImageDimensions();  iPoint3D voxel;
        FUNC_INFO << "cost NC, dim" << dim.x << dim.y << dim.z;
        for (voxel.z=0; voxel.z<dim.z; voxel.z++)
        {
            for (voxel.y=0; voxel.y<dim.y; voxel.y++)
            {
                for (voxel.x=0; voxel.x<dim.x; voxel.x++)
                {
                    if ( mask->isVoxelNonZero(voxel) )
                    {
                        double valueT = templateVolume->getStackValue(voxel,iVolumeTemplate);
                        double valueI = getStackValue(voxel,iVolume);
                        double weight = mask->getBitMapValue(voxel);
                        if ( valueI != 0. && valueT != 0. )
                        {
                            dotIT += weight * valueI * valueT;
                            dotII +=          valueI * valueI;
                            dotTT +=          valueT * valueT;
                        }
                    }
                }
            }
        }
        if ( dotII == 0. || dotTT == 0. )
            cost = 0.;
        else
            cost = dotIT / qSqrt(dotII) / qSqrt(dotTT);
        FUNC_INFO << cost;
    }
    else
        qFatal("Cost function unknown in ImageIO::costFunction");
    return(cost);
}
double ImageIO::calculateCostFunction(CostFunctionChoice costFunction, int iVolume,
                                      ImageIO *templateVolume, int iVolumeTemplate)
{
    FUNC_ENTER << iVolume;
    double cost = 0;
    if ( templateVolume->voxelsPerVolume() == 0 || voxelsPerVolume() == 0 ) return cost;

    if ( costFunction == CostFunction_MI || costFunction == CostFunction_NMI )
    { // Mutual information
        // Find the average value of non-zero voxels for the source and target
        double avTarget   = calculateAverage(iVolume);
        double avTemplate = templateVolume->calculateAverage(iVolumeTemplate);
        if ( avTarget == 0. ) avTarget = 1.;
        FUNC_INFO << "avTarget=" << avTarget << "avTemplate=" << avTemplate;
        // Find the min and max.
        iPoint3D dim = templateVolume->getImageDimensions();
        double min=1.e10;  double max=-1e10;
        FUNC_INFO << "nVoxels" << templateVolume->voxelsPerVolume() << voxelsPerVolume();
        for ( int j=0; j<templateVolume->voxelsPerVolume(); j++ )
        {
            // normalize all values by the respective average
            double valueTarget   = getStackValue(j,iVolume)                         / avTarget;
            double valueTemplate = templateVolume->getStackValue(j,iVolumeTemplate) / avTemplate;
            if ( valueTarget   < min ) min = valueTarget;
            if ( valueTemplate < min ) min = valueTemplate;
            if ( valueTarget   > max ) max = valueTarget;
            if ( valueTemplate > max ) max = valueTemplate;
        }
        FUNC_INFO << "min1=" << min << "max1=" << max;
        // Initialize the histograms.
        dVector probI, probT;
        probI.fill(0.,NCOSTBINS);  probT.fill(0.,NCOSTBINS);
        dMatrix probIT;
        probIT.resize(NCOSTBINS);
        for (int jI=0; jI<NCOSTBINS; jI++)
            probIT[jI].fill(0.,NCOSTBINS);
        // Fill the histograms
        dim = templateVolume->getImageDimensions();
        double norm = 0.;
        for ( int j=0; j<templateVolume->voxelsPerVolume(); j++ )
        {
            double valueTarget   = getStackValue(j,iVolume)                         / avTarget;
            double valueTemplate = templateVolume->getStackValue(j,iVolumeTemplate) / avTemplate;

            iPoint2D iBinI, iBinT;
            dPoint2D weightBinI, weightBinT;
            utilMath::fuzzyBinning(valueTarget,    min,  max, NCOSTBINS, iBinI, weightBinI );
            utilMath::fuzzyBinning(valueTemplate,  min,  max, NCOSTBINS, iBinT, weightBinT );
            if ( iBinI.x >= 0 )
                probI[iBinI.x] += weightBinI.x;
            if ( iBinI.y >= 0 )
                probI[iBinI.y] += weightBinI.y;

            if ( iBinT.x >= 0 )
            {
                probT[iBinT.x] += weightBinT.x;
                if ( iBinI.x >= 0 )
                    probIT[iBinI.x][iBinT.x] += weightBinI.x*weightBinT.x;
                if ( iBinI.y >= 0 )
                    probIT[iBinI.y][iBinT.x] += weightBinI.y*weightBinT.x;
            }
            if ( iBinT.y >= 0 )
            {
                probT[iBinT.y] += weightBinT.y;
                if ( iBinI.x >= 0 )
                    probIT[iBinI.x][iBinT.y] += weightBinI.x*weightBinT.y;
                if ( iBinI.y >= 0 )
                    probIT[iBinI.y][iBinT.y] += weightBinI.y*weightBinT.y;
            }
            norm += 1.;
        } // loop over voxels
        FUNC_INFO << "norm = " << norm;
        // Normalize the probability functions.
        for (int jI=0; jI<NCOSTBINS; jI++)
        {
            probI[jI] /= norm;
            probT[jI] /= norm;
            for (int jT=0; jT<NCOSTBINS; jT++)
                probIT[jI][jT] /= norm;
        }
        if ( costFunction == CostFunction_MI )
        {
            // Calculate the cost function.
            for (int jI=0; jI<NCOSTBINS; jI++)
            { // sum over bins of target
                if ( probI[jI] != 0. )
                {
                    for (int jT=0; jT<NCOSTBINS; jT++)
                    { // sum over bins of target
                        if ( probT[jT] != 0. && probIT[jI][jT] != 0. )
                            cost += probIT[jI][jT] * qLn(probIT[jI][jT]/probI[jI]/probT[jT]);
                    }
                }
            }
        }
        else // normalized MI
        {
            double H_I=0.; double H_T=0.; double H_IT=0.;
            for (int jI=0; jI<NCOSTBINS; jI++)
            { // sum over bins of target
                if ( probI[jI] != 0. )
                    H_I -= probT[jI] * qLn(probI[jI]);
            }
            for (int jT=0; jT<NCOSTBINS; jT++)
            { // sum over bins of target
                if ( probT[jT] != 0. )
                    H_T -= probT[jT] * qLn((probT[jT]));
            }
            for (int jI=0; jI<NCOSTBINS; jI++)
            {
                for (int jT=0; jT<NCOSTBINS; jT++)
                { // sum over bins of target
                    if ( probIT[jI][jT] != 0. )
                        H_IT -= probIT[jI][jT] * qLn((probIT[jI][jT]));
                }
            }
            if ( H_IT != 0. ) cost = H_IT / (H_I + H_T);
        }
    }
    else if ( costFunction == CostFunction_WLS )
    {
        FUNC_INFO << "cost WLS";
        double sumWLS = 0.;
        iPoint3D dim = getImageDimensions();  iPoint3D voxel;
        for (voxel.z=0; voxel.z<dim.z; voxel.z++)
        {
            for (voxel.y=0; voxel.y<dim.y; voxel.y++)
            {
                for (voxel.x=0; voxel.x<dim.x; voxel.x++)
                {
                    double valueT = templateVolume->getStackValue(voxel,iVolumeTemplate);
                    double valueI = getStackValue(voxel,iVolume);
                    if ( valueI != 0. && valueT != 0. )
                    {
                        sumWLS += SQR(valueT - valueI);
                    }
                }
            }
        }
        cost = sumWLS;
    }
    else if ( costFunction == CostFunction_RV )
    { // same as WLS but weights = 1
        FUNC_INFO << "cost RV";
        double residue2 = 0.;
        iPoint3D dim = getImageDimensions();  iPoint3D voxel;
        for (voxel.z=0; voxel.z<dim.z; voxel.z++)
        {
            for (voxel.y=0; voxel.y<dim.y; voxel.y++)
            {
                for (voxel.x=0; voxel.x<dim.x; voxel.x++)
                {
                    double valueT = templateVolume->getStackValue(voxel,iVolumeTemplate);
                    double valueI = getStackValue(voxel,iVolume);
                    if ( valueI != 0. && valueT != 0. )
                        residue2 += SQR(valueI - valueT);
                }
            }
        }
        cost = residue2 / (dim.x*dim.y*dim.z);  // "cost" = residual variance, not accounting for model parameters
    }
    else if ( costFunction == CostFunction_NC )
    {
        double dotII = 0.; // target ("individual)
        double dotTT = 0.; // template
        double dotIT = 0.; // target . template
        iPoint3D dim = templateVolume->getImageDimensions();  iPoint3D voxel;
        FUNC_INFO << "cost NC, dim" << dim.x << dim.y << dim.z;
        for (voxel.z=0; voxel.z<dim.z; voxel.z++)
        {
            for (voxel.y=0; voxel.y<dim.y; voxel.y++)
            {
                for (voxel.x=0; voxel.x<dim.x; voxel.x++)
                {
                    double valueT = templateVolume->getStackValue(voxel,iVolumeTemplate);
                    double valueI = getStackValue(voxel,iVolume);
                    if ( valueI != 0. && valueT != 0. )
                    {
                        dotIT += valueI * valueT;
                        dotII += valueI * valueI;
                        dotTT += valueT * valueT;
                    }
                }
            }
        }
        if ( dotII == 0. || dotTT == 0. )
            cost = 0.;
        else
            cost = dotIT / qSqrt(dotII) / qSqrt(dotTT);
        FUNC_INFO << cost;
    }
    else
        qFatal("Cost function unknown in ImageIO::costFunction");
    return(cost);
}
/*
double ImageIO::calculateCostFunction(CostFunctionChoice costFunction, int iVolume,
                                      ImageIO *templateVolume, int iVolumeTemplate,
                                      bool includeTemplateZeroes)
{
    FUNC_ENTER << includeTemplateZeroes;
    double cost = 0;
    if ( costFunction == CostFunction_MI || costFunction == CostFunction_NMI )
    { // Mutual information
        // Find the average value of non-zero voxels for the source and target
        double avTarget   = calculateAverage(iVolume, templateVolume, iVolumeTemplate,includeTemplateZeroes);
        double avTemplate = templateVolume->calculateAverage(iVolumeTemplate, templateVolume,
                                                             iVolumeTemplate, includeTemplateZeroes);
        if ( avTarget == 0. ) avTarget = 1.;
        FUNC_INFO << "avTarget=" << avTarget << "avTemplate=" << avTemplate;
        // Find the min and max.
        iPoint3D dim = templateVolume->getImageDimensions();
        double min=1.e10;  double max=-1e10;
        for ( int j=0; j<templateVolume->voxelsPerVolume(); j++ )
        {
            bool nonZero = (templateVolume->getStackValue(j,iVolumeTemplate) != 0.);
            if ( nonZero || includeTemplateZeroes )
            {
                // normalize all values by the respective average
                double valueTarget   = getStackValue(j,iVolume)                         / avTarget;
                double valueTemplate = templateVolume->getStackValue(j,iVolumeTemplate) / avTemplate;
                if ( valueTarget   < min ) min = valueTarget;
                if ( valueTemplate < min ) min = valueTemplate;
                if ( valueTarget   > max ) max = valueTarget;
                if ( valueTemplate > max ) max = valueTemplate;
            }
        }
        FUNC_INFO << "min =" << min << "max=" << max;
        // Initialize the histograms.
        dVector probI, probT;
        probI.fill(0.,NCOSTBINS);  probT.fill(0.,NCOSTBINS);
        dMatrix probIT;
        probIT.resize(NCOSTBINS);
        for (int jI=0; jI<NCOSTBINS; jI++)
            probIT[jI].fill(0.,NCOSTBINS);
        // Fill the histograms
        dim = templateVolume->getImageDimensions();
        double norm = 0.;
        for ( int j=0; j<templateVolume->voxelsPerVolume(); j++ )
        {
            bool nonZero = (templateVolume->getStackValue(j,iVolumeTemplate) != 0.);
            if ( nonZero || includeTemplateZeroes )
            {
                double valueTarget   = getStackValue(j,iVolume)                         / avTarget;
                double valueTemplate = templateVolume->getStackValue(j,iVolumeTemplate) / avTemplate;

                iPoint2D iBinI, iBinT;
                dPoint2D weightBinI, weightBinT;
                utilMath::fuzzyBinning(valueTarget,   min,  max, NCOSTBINS, iBinI, weightBinI );
                utilMath::fuzzyBinning(valueTemplate, min,  max, NCOSTBINS, iBinT, weightBinT );
                if ( iBinI.x >= 0 )
                    probI[iBinI.x] += weightBinI.x;
                if ( iBinI.y >= 0 )
                    probI[iBinI.y] += weightBinI.y;

                if ( iBinT.x >= 0 )
                {
                    probT[iBinT.x] += weightBinT.x;
                    if ( iBinI.x >= 0 )
                        probIT[iBinI.x][iBinT.x] += weightBinI.x*weightBinT.x;
                    if ( iBinI.y >= 0 )
                        probIT[iBinI.y][iBinT.x] += weightBinI.y*weightBinT.x;
                }
                if ( iBinT.y >= 0 )
                {
                    probT[iBinT.y] += weightBinT.y;
                    if ( iBinI.x >= 0 )
                        probIT[iBinI.x][iBinT.y] += weightBinI.x*weightBinT.y;
                    if ( iBinI.y >= 0 )
                        probIT[iBinI.y][iBinT.y] += weightBinI.y*weightBinT.y;
                }
                norm += 1.;
            } // nonZero || includeTemplateZeroes
        } // loop over voxels
        FUNC_INFO << "norm = " << norm;
        // Normalize the probability functions.
        for (int jI=0; jI<NCOSTBINS; jI++)
        {
            probI[jI] /= norm;
            probT[jI] /= norm;
            for (int jT=0; jT<NCOSTBINS; jT++)
                probIT[jI][jT] /= norm;
        }
        if ( costFunction == CostFunction_MI )
        {
            // Calculate the cost function.
            for (int jI=0; jI<NCOSTBINS; jI++)
            { // sum over bins of target
                if ( probI[jI] != 0. )
                {
                    for (int jT=0; jT<NCOSTBINS; jT++)
                    { // sum over bins of target
                        if ( probT[jT] != 0. && probIT[jI][jT] != 0. )
                            cost += probIT[jI][jT] * qLn(probIT[jI][jT]/probI[jI]/probT[jT]);
                    }
                }
            }
        }
        else // normalized MI
        {
            double H_I=0.; double H_T=0.; double H_IT=0.;
            for (int jI=0; jI<NCOSTBINS; jI++)
            { // sum over bins of target
                if ( probI[jI] != 0. )
                    H_I -= probT[jI] * qLn(probI[jI]);
            }
            for (int jT=0; jT<NCOSTBINS; jT++)
            { // sum over bins of target
                if ( probT[jT] != 0. )
                    H_T -= probT[jT] * qLn((probT[jT]));
            }
            for (int jI=0; jI<NCOSTBINS; jI++)
            {
                for (int jT=0; jT<NCOSTBINS; jT++)
                { // sum over bins of target
                    if ( probIT[jI][jT] != 0. )
                        H_IT -= probIT[jI][jT] * qLn((probIT[jI][jT]));
                }
            }
            if ( H_IT != 0. ) cost = H_IT / (H_I + H_T);
        }
    }
    else if ( costFunction == CostFunction_WLS )
    {
        FUNC_INFO << "cost WLS";
        double sumWLS = 0.;
        iPoint3D dim = getImageDimensions();  iPoint3D voxel;
        for (voxel.z=0; voxel.z<dim.z; voxel.z++)
        {
            for (voxel.y=0; voxel.y<dim.y; voxel.y++)
            {
                for (voxel.x=0; voxel.x<dim.x; voxel.x++)
                {
                    double valueT = templateVolume->getStackValue(voxel,iVolumeTemplate);
                    bool nonZero = (valueT != 0.);
                    if ( nonZero || includeTemplateZeroes )
                    {
                        double valueI = getStackValue(voxel,iVolume);
                        double weight = getMinDistanceToBoundaryWeight(voxel);
                        weight *= valueT;
                        sumWLS += weight * SQR(valueT - valueI);
                    }
                }
            }
        }
        cost = sumWLS;
    }
    else if ( costFunction == CostFunction_RV )
    { // same as WLS but weights = 1
        FUNC_INFO << "cost RV";
        double residue2 = 0.;
        iPoint3D dim = getImageDimensions();  iPoint3D voxel;
        for (voxel.z=0; voxel.z<dim.z; voxel.z++)
        {
            for (voxel.y=0; voxel.y<dim.y; voxel.y++)
            {
                for (voxel.x=0; voxel.x<dim.x; voxel.x++)
                {
                    double valueT = templateVolume->getStackValue(voxel,iVolumeTemplate);
                    bool nonZero = (valueT != 0.);
                    if ( nonZero || includeTemplateZeroes )
                    {
                        double valueT = templateVolume->getStackValue(voxel,iVolumeTemplate);
                        double valueI = getStackValue(voxel,iVolume);
                        residue2 += SQR(valueI - valueT);
                    }
                }
            }
        }
        cost = residue2 / (dim.x*dim.y*dim.z);  // "cost" = residual variance, not accounting for model parameters
    }
    else if ( costFunction == CostFunction_NC )
    {
        double dotII = 0.; // target ("individual)
        double dotTT = 0.; // template
        double dotIT = 0.; // target . template
        iPoint3D dim = getImageDimensions();
        FUNC_INFO << "cost NC, dims" << dim.x << dim.y << dim.z;
        iPoint3D voxel;
        for (voxel.z=0; voxel.z<dim.z; voxel.z++)
        {
            for (voxel.y=0; voxel.y<dim.y; voxel.y++)
            {
                for (voxel.x=0; voxel.x<dim.x; voxel.x++)
                {
                    double valueT = templateVolume->getStackValue(voxel,iVolumeTemplate);
                    double valueI = getStackValue(voxel,iVolume);
                    bool nonZero = (valueT != 0. && valueI != 0.);
                    if ( nonZero || includeTemplateZeroes )
                    {
//                        double weight = getMinDistanceToBoundaryWeight(voxel);
                        double weight = 1.;
                        dotIT += weight * valueI * valueT;
                        dotII +=          valueI * valueI;
                        dotTT +=          valueT * valueT;
                    }
                } // x
            } // y
        } // z
        if ( dotII == 0. || dotTT == 0. )
            cost = 0.;
        else
            cost = dotIT / qSqrt(dotII) / qSqrt(dotTT);
        FUNC_INFO << "costValue" << cost;
    }
    else
        qFatal("Cost function unknown in ImageIO::costFunction");
    return(cost);
}
*/
double ImageIO::getMinDistanceToBoundaryWeight(iPoint3D voxel)
{
    iPoint3D dim = getImageDimensions();
    dPoint3D res = getImageResolution();
    dPoint3D distance;
    distance.x = qMin(voxel.x+0.5,dim.x-voxel.x+0.5) * res.x;
    distance.y = qMin(voxel.y+0.5,dim.y-voxel.y+0.5) * res.y;
    distance.z = qMin(voxel.z+0.5,dim.z-voxel.z+0.5) * res.z;
    distance.x = qAbs(distance.x);
    distance.y = qAbs(distance.y);
    distance.z = qAbs(distance.z);

    // very minor cutoff
    dPoint3D FOV = getFOV();
    dPoint3D cutoff;
    cutoff.x = FOV.x / 10.;
    cutoff.y = FOV.y / 10.;
    cutoff.z = FOV.z / 10.;

    dPoint3D weight;
    weight.x = weight.y = weight.z = 1.;
    if ( distance.x < cutoff.x ) weight.x = distance.x / cutoff.x;
    if ( distance.y < cutoff.x ) weight.y = distance.y / cutoff.y;
    if ( distance.z < cutoff.z ) weight.z = distance.z / cutoff.z;

    double weightProduct = weight.x*weight.y*weight.z;
    return weightProduct;
}

double ImageIO::calculateAverage(int iTime)
{
    double average = 0.;
    int numberVoxels = 0;
    for ( int jVoxel=0; jVoxel<voxelsPerVolume(); ++jVoxel)
    {
        double value = getStackValue(jVoxel,iTime);
        average += getStackValue(jVoxel,iTime);
        if ( value != 0. ) numberVoxels++;
    }
    if ( numberVoxels != 0 )
        average /= static_cast<double>(numberVoxels);
    return average;
}

double ImageIO::calculateAverage(int iTime, bitMap *mask)
{
    double average = 0.;
    int numberVoxels = 0;
    for ( int jVoxel=0; jVoxel<voxelsPerVolume(); ++jVoxel)
    {
        if ( mask->isVoxelNonZero(jVoxel) )
        {
            average += getStackValue(jVoxel,iTime);
            numberVoxels++;
        }
    }
    if ( numberVoxels != 0 )
        average /= static_cast<double>(numberVoxels);
    return average;
}

double ImageIO::calculateAverage(int iTime, VoxelSet voxelList)
{
    double average = 0.;
    int numberVoxels = voxelList.size();
    if ( numberVoxels == 0 ) return 0.;
    for ( int jVoxel=0; jVoxel<numberVoxels; ++jVoxel)
    {
        fVoxel voxel = voxelList.at(jVoxel);
        average += getStackValue(voxel,iTime);
        numberVoxels++;
    }
    average /= static_cast<double>(numberVoxels);
    return average;
}

double ImageIO::calculateAverage(int iTime, OVLFile overlay)
{
    return calculateAverage(iTime, overlay.voxelList);
}

double ImageIO::calculateAverage(int iTime, ImageIO *templateFile, int iTimeTemplate, bool includeTemplateZeroes)
{
    double average = 0.;
    int numberVoxels = 0;
    for ( int jVoxel=0; jVoxel<voxelsPerVolume(); ++jVoxel)
    {
        double valueTemplate = templateFile->getStackValue(jVoxel,iTimeTemplate);
        bool nonZero = (valueTemplate != 0.);
        if ( nonZero || includeTemplateZeroes )
        {
            average += getStackValue(jVoxel,iTime);
            numberVoxels++;
        }
    }
    if ( numberVoxels != 0 )
        average /= static_cast<double>(numberVoxels);
    return average;
}

dPoint2D ImageIO::getThreshold(int iTime)
{
    dPoint2D extrema = getMinAndMax(iTime);
    int nBins = 200;  dVector histogram;  histogram.resize(nBins);
    fillHistogram(iTime,extrema.lower, extrema.upper, histogram);
    int iThresh = utilMath::OtsuThreshold(histogram);
    double binStep   = (extrema.upper - extrema.lower) / static_cast<double>(nBins);
    double OtsuThresh = extrema.lower + binStep * iThresh;

    double yMax=0.;  int iBinLow=0;  double sum=0.;
    for (int jBin=0; jBin<nBins; jBin++)
    {
        double yValue = histogram[jBin];
        sum += yValue;
        if ( yValue > yMax )
        {
            yMax = yValue;
            iBinLow = jBin;
        }
    }
    double fractionHigh = 0.99;
    double sumTargetHigh = fractionHigh * sum;
    sum=0.;  int iBinHigh=0;
    for (int jBin=0; jBin<nBins; jBin++)
    {
        double yValue = histogram[jBin];
        sum += yValue;
        if ( sum >= sumTargetHigh )
        {
            iBinHigh = jBin;
            break;
        }
    }
    extrema.upper = extrema.lower + iBinHigh * binStep;
    FUNC_INFO << "set upper to" << extrema.upper;
    extrema.lower = OtsuThresh;
    extrema.upper = qMax(extrema.lower,extrema.upper);
    if ( extrema.lower == extrema.upper ) extrema.upper *= 2.;

    return extrema;
}
dPoint2D ImageIO::getThreshold()
{
    dPoint2D extrema;
    extrema.lower =  1.e10;  // min
    extrema.upper = -1.e10;  // max
    for (int jTime=0; jTime<nt(); jTime++)
    {
        for (int jVoxel=0; jVoxel<voxelsPerVolume(); jVoxel++)
        {
            double value = getStackValue(jVoxel,jTime);
            if ( value < extrema.lower ) extrema.lower = value;
            if ( value > extrema.upper ) extrema.upper = value;
        }
    }
    return extrema;
}

d2Vector ImageIO::createProjectionInXYPlane(dPoint2D COM, int iZ, int angleDeg, double radiusStep, d2Vector &coordinates)
{ // lastOverlayCrossing is a (radius,signal) pairs
    d2Vector projection;
    double angle = angleDeg * PI/180.;  // deg

    dPoint3D res = getImageResolution();
    if ( radiusStep == 0. ) radiusStep = qMin(res.x,res.y);
    dPoint3D fwhm = res;
    iPoint3D iWidth;
    iWidth.x = qCeil(fwhm.x);
    iWidth.y = qCeil(fwhm.y);
    iWidth.z = qCeil(fwhm.z);
    iPoint3D wrap;          wrap.x = wrap.y = wrap.z = 0;
    enum kernels {Lanczos, Gaussian, Delta};
    iPoint3D kernelType;
    kernelType.x = kernelType.y = Lanczos;
    kernelType.z = Delta;

    // get z coordinate
    iPoint3D voxel = {0,0,iZ};
    dPoint3D point;  // point = (x,y,z)
    point.z = convertVoxelToSpace(voxel).z;
    dPoint2D pp;     // projection point = (radius,theta)
    double radius = 0.;
    bool inVolume=true;
    while ( inVolume )
    {
        point.x = COM.x + radius * qCos(angle);
        point.y = COM.y + radius * qSin(angle);
        double voxelValue = interpolateValueFromSpaceCoord( point, fwhm, iWidth, wrap,
                                                            kernelType, 0, inVolume);
        /*
        if ( angleDeg == 0 )
        {
            QINFO << "interpolate radius" << radius << "point" << point.x << point.y << point.z;
            QINFO << "involume" << inVolume << "value" << voxelValue;
        }
        */
        if ( inVolume )
        {
            // add the point to the projection
            pp.x = radius;      pp.y = voxelValue;
            projection.append(pp);
            // save the (x,y) pair
            pp.x = point.x;     pp.y = point.y;
            coordinates.append(pp);
        }
        radius += radiusStep;
    }

    return projection;
}

dcVector ImageIO::convertToComplexVolume(int iTime)
{
    iPoint3D dim = getImageDimensions();
    dcVector volume; volume.resize(dim.x*dim.y*dim.z);

    // Create a complex volume
    iPoint3D voxel;
    for (voxel.z=0; voxel.z<dim.z; voxel.z++)
    {
        for (voxel.y=0; voxel.y<dim.y; voxel.y++)
        {
            for (voxel.x=0; voxel.x<dim.x; voxel.x++)
            {
                volume[utilIO::index1d(dim, voxel)].real = getStackValue(voxel,iTime);
                volume[utilIO::index1d(dim, voxel)].imag = 0.;
            }
        }
    }
    return volume;
}

void ImageIO::logTransform(bool forward, bitMap mask)
{
    { // forward = ln; reverse = exp
        FUNC_ENTER;
        iPoint3D voxel;
        // Loop over all time points
        for (int jt=0; jt<nt(); jt++)
        {
            // Loop over voxels in the template space
            for (voxel.z=0; voxel.z<nz(); voxel.z++)
            {
                for (voxel.y=0; voxel.y<ny(); voxel.y++)
                {
                    for (voxel.x=0; voxel.x<nx(); voxel.x++)
                    {
                        if ( mask.isVoxelNonZero(voxel) )
                        {
                            double value = getStackValue(voxel,jt);
                            if ( forward )
                            { // ln
                                if ( value <= 0. ) // filter out non-pos values
                                    value = 0.;
                                else
                                    value = qLn(value);
                            }
                            else // exp: replace 0 values with 0, not 1
                                if ( value != 0. ) value = qExp(value);
                            // replace value
                            setStackValue(voxel,jt,value);
                        }
                    } // .x
                } // .y
            } // .z
        } // .t
    }
}

void ImageIO::logTransform(bool forward)
{ // forward = ln; reverse = exp
    FUNC_ENTER;
    iPoint4D voxel;
    // Loop over all time points
    for (voxel.t=0; voxel.t<nt(); voxel.t++)
    {
        // Loop over voxels in the template space
        for (voxel.z=0; voxel.z<nz(); voxel.z++)
        {
            for (voxel.y=0; voxel.y<ny(); voxel.y++)
            {
                for (voxel.x=0; voxel.x<nx(); voxel.x++)
                {
                    double value = getStackValue(voxel);
                    if ( forward )
                    { // ln
                        if ( value <= 0. ) // filter out non-pos values
                            value = 0.;
                        else
                            value = qLn(value);
                    }
                    else // exp: replace 0 values with 0, not 1
                        if ( value != 0. ) value = qExp(value);
                    // replace value
                    setStackValue(voxel,value);
                } // .x
            } // .y
        } // .z
    } // .t
}

double ImageIO::spatialPoly1D( int iPower, int iPos, iPoint2D range)
{ // legendre polynomial in space with exponential attenuation at edges
  double rVoxel, Basis;

  // Put the index on a scale of -1 to 1
  // rVoxel = -1 at range.lower
  // rVoxel =  1 at range.upper
  // slope = 2/(range.upper-range.lower)
  // 2/(range.upper-range.lower) * range.lower + b = -1
  // b = -1 - 2/(range.upper-range.lower)

  rVoxel = 2./static_cast<double>(range.upper-range.lower) * static_cast<double>(iPos-range.lower) - 1.;
  dPoint2D rRange; rRange.lower = -1.;  rRange.upper = 1.;

  double attenuation=1.;  double tau=0.3;
  if ( rVoxel > rRange.upper )
  {
      double t = rVoxel-rRange.upper;
      attenuation = qExp(-t/tau);
  }
  else if ( rVoxel < rRange.lower )
  {
      double t = rRange.lower-rVoxel;
      attenuation = qExp(-t/tau);
  }

  Basis = utilMath::polynomialLegendre(iPower, rVoxel) * attenuation;
  return(Basis);
}

void ImageIO::flattenByGLM(iPoint3D dimPoly, bitMap mask)
{
    iPoint3D dim = getImageDimensions();

    int nCoeff = dimPoly.x * dimPoly.y * dimPoly.z;      // 27 for 2nd order; 64 for 3rd order

    iPoint2D rangeX={{dim.x-1},{0}};
    iPoint2D rangeY={{dim.y-1},{0}};
    iPoint2D rangeZ={{dim.z-1},{0}};
    for (int jz=0; jz<nz(); jz++)
    {
        for (int jy=0; jy<ny(); jy++)
        {
            for (int jx=0; jx<nx(); jx++)
            {
                iPoint3D voxel={jx,jy,jz};
                if ( mask.isVoxelNonZero(voxel) )
                {
                    if ( jx < rangeX.lower ) rangeX.lower = jx;
                    if ( jy < rangeY.lower ) rangeY.lower = jy;
                    if ( jz < rangeZ.lower ) rangeZ.lower = jz;
                    if ( jx > rangeX.upper ) rangeX.upper = jx;
                    if ( jy > rangeY.upper ) rangeY.upper = jy;
                    if ( jz > rangeZ.upper ) rangeZ.upper = jz;
                }
            } // jx
        } // jy
    } // jz

    // Create polynomials in x, y, and z: dimension of matrices = [order][location]
    iPoint3D voxel;  dMatrix polyX;  polyX.resize(dimPoly.x);
    for (int jp=0; jp<dimPoly.x; jp++)
    {
        for (voxel.x=0; voxel.x<dim.x; voxel.x++)
            polyX[jp].append(spatialPoly1D(jp,voxel.x,rangeX));
    }
    dMatrix polyY;  polyY.resize(dimPoly.y);
    for (int jp=0; jp<dimPoly.y; jp++)
    {
        for (voxel.y=0; voxel.y<dim.y; voxel.y++)
            polyY[jp].append(spatialPoly1D(jp,voxel.y,rangeY));
    }
    dMatrix polyZ;  polyZ.resize(dimPoly.z);
    for (int jp=0; jp<dimPoly.z; jp++)
    {
        for (voxel.z=0; voxel.z<dim.z; voxel.z++)
            polyZ[jp].append(spatialPoly1D(jp,voxel.z,rangeZ));
    }

    iPoint3D iStepVoxel;
    iStepVoxel.x = qMax(1,dim.x/128);
    iStepVoxel.y = qMax(1,dim.y/128);
    iStepVoxel.z = qMax(1,dim.z/128);
//    iStepVoxel.x = iStepVoxel.y = iStepVoxel.z = 1;

    // Fit spatial voxels (v) signals Y with polynomial coefficients (c) using a design matrix X:
    // Y_v = X_vc * B_c
    // where
    // Y_v is a concatenation of all voxels inside the brain (0 outside),
    // X_vc has polynomial functions based upon (x,y,z). In each dimension, values are normalized (+-1)
    //      at edges of the FOV. Values are 0 where Y is zero.
    // B_c = {(X^t*X)^-1}_cc' * (X^t)_c'v Y_v
    //
    // To save memory, don't construct X explicitly as a matrix - just calculate components as needed.
    // Create matrix inverse_cc = Xt * X; Inverse = Xt_c'v * X_vc
    dMatrix inverse_cc;  inverse_cc.resize(nCoeff);
    for (int jCoeff=0; jCoeff<nCoeff; jCoeff++)
        inverse_cc[jCoeff].fill(0.,nCoeff);

    logTransform(true);

    for (voxel.z=0; voxel.z<dim.z; voxel.z+=iStepVoxel.z)
    {
        for (voxel.y=0; voxel.y<dim.y; voxel.y+=iStepVoxel.y)
        {
            for (voxel.x=0; voxel.x<dim.x; voxel.x+=iStepVoxel.x)
            {
                if ( mask.isVoxelNonZero(voxel) &&  getStackValue(voxel,0) != 0. )
                {
                    iPoint3D iPoly, iPoly1;  int iCoeff=0;
                    for (iPoly.z=0; iPoly.z<dimPoly.z; iPoly.z++)
                    {
                        for (iPoly.y=0; iPoly.y<dimPoly.y; iPoly.y++)
                        {
                            for (iPoly.x=0; iPoly.x<dimPoly.x; iPoly.x++, iCoeff++)
                            {
                                int iCoeff1 = 0;
                                for (iPoly1.z=0; iPoly1.z<dimPoly.z; iPoly1.z++)
                                {
                                    for (iPoly1.y=0; iPoly1.y<dimPoly.y; iPoly1.y++)
                                    {
                                        for (iPoly1.x=0; iPoly1.x<dimPoly.x; iPoly1.x++, iCoeff1++)
                                        {
                                            inverse_cc[iCoeff][iCoeff1] +=
                                                    polyX[iPoly.x][voxel.x] *polyY[iPoly.y][voxel.y] *polyZ[iPoly.z][voxel.z] *
                                                    polyX[iPoly1.x][voxel.x]*polyY[iPoly1.y][voxel.y]*polyZ[iPoly1.z][voxel.z];
                                        } // power x 1
                                    } // power y 1
                                } // power z 1
                            } // power x
                        } // power y
                    } // power z
                }
            } // voxel x
        } // voxel y
    } // voxel z
    // Invert Inverse = (X^tX)^-1.
    utilMath::invertSquareMatrix(inverse_cc);
    // Create matrix XtY_c = Xt_cv * Y_v

    dVector XtY_c;  XtY_c.fill(0.,nCoeff);
    for (voxel.z=0; voxel.z<dim.z; voxel.z+=iStepVoxel.z)
    {
        for (voxel.y=0; voxel.y<dim.y; voxel.y+=iStepVoxel.y)
        {
            for (voxel.x=0; voxel.x<dim.x; voxel.x+=iStepVoxel.x)
            {
                if ( mask.isVoxelNonZero(voxel) &&  getStackValue(voxel,0) != 0. )
                {
                    iPoint3D iPoly;  int iCoeff=0;
                    for (iPoly.z=0; iPoly.z<dimPoly.z; iPoly.z++)
                    {
                        for (iPoly.y=0; iPoly.y<dimPoly.y; iPoly.y++)
                        {
                            for (iPoly.x=0; iPoly.x<dimPoly.x; iPoly.x++, iCoeff++)
                                XtY_c[iCoeff] += polyX[iPoly.x][voxel.x]
                                        *        polyY[iPoly.y][voxel.y]
                                        *        polyZ[iPoly.z][voxel.z] * getStackValue(voxel,0);
                        } // power y
                    } // power z
                }
            } // index x
        } // index y
    } // index z


    // The hard work is done; now compute Beta. B_c = inverse_cc' * XtY_c'
    dVector beta_c;  beta_c.resize(nCoeff);
    for (int jc=0; jc<nCoeff; jc++)
    {
        beta_c[jc] = 0.;
        for (int jc1=0; jc1<nCoeff; jc1++)
            beta_c[jc] += inverse_cc[jc][jc1] * XtY_c[jc1];
    }
    dVector biasField; biasField.fill(0.,voxelsPerVolume());
    // Now compute the bias correction as Yfit_v = X_vc * B_c
    // Compute all terms except for the first basis function, which is power (0,0,0) = constant.
    double dcTerm;    double SOS = 0.;
    for (voxel.z=0; voxel.z<dim.z; voxel.z++)
    {
        for (voxel.y=0; voxel.y<dim.y; voxel.y++)
        {
            for (voxel.x=0; voxel.x<dim.x; voxel.x++)
            {
                int iVoxel = index3dFromVoxel(voxel);
                iPoint3D iPoly;  int iCoeff=0;
                for (iPoly.z=0; iPoly.z<dimPoly.z; iPoly.z++)
                {
                    for (iPoly.y=0; iPoly.y<dimPoly.y; iPoly.y++)
                    {
                        for (iPoly.x=0; iPoly.x<dimPoly.x; iPoly.x++, iCoeff++)
                        {
                            if ( (iPoly.x==0 && iPoly.y==0 && iPoly.z==0) )
                                dcTerm = beta_c[iCoeff];
                            else
                                biasField[iVoxel] += polyX[iPoly.x][voxel.x]
                                        *            polyY[iPoly.y][voxel.y]
                                        *            polyZ[iPoly.z][voxel.z] * beta_c[iCoeff];
                        } // power x
                    } // power y
                } // power z
            } // index x
        } // index y
    } // index z

    // Subtract off the bias correction.
    for (voxel.z=0; voxel.z<dim.z; voxel.z++)
    {
        for (voxel.y=0; voxel.y<dim.y; voxel.y++)
        {
            for (voxel.x=0; voxel.x<dim.x; voxel.x++)
            {
                int iVoxel = index3dFromVoxel(voxel);
                double value = getStackValue(iVoxel,0);
                if ( mask.isVoxelNonZero(voxel) )
                    SOS += SQR(value-biasField[iVoxel]);
                setStackValue(iVoxel,0,value - biasField[iVoxel]);
            } // .x
        } // .y
    } // .z
    /*
    for (int jVoxel=0; jVoxel<biasField.size(); jVoxel++)
    {
        // the dcTerm can be removed to give a value near 1
        double value = getStackValue(jVoxel,0);
//        setStackValue(jVoxel,0,value-dcTerm);
        setStackValue(jVoxel,0,value);
    }
    */
    logTransform(false);
}

void ImageIO::NLMeans(int iTime, bool dim3D, int patchVoxels, int searchVoxels)
{
    FUNC_ENTER;
    // iTime:   volume in file

    // rangePatch:  define patch extent using +- units: e.g., 1 --> 3 voxels
    // rangeSearch: define search range between patches: e.g., 3 --> 7 voxels
    // ensure that rangeSearch > rangePatch always to avoid boundary problems below.
    _timer.start();
    int nZPatch, nZSearch;
    if ( dim3D )
    {
        nZPatch  = patchVoxels;
        nZSearch = searchVoxels;
    }
    else
        nZPatch = nZSearch = 0;
    iPoint3D rangePatch  = {patchVoxels, patchVoxels, nZPatch};
    iPoint3D rangeSearch = {searchVoxels,searchVoxels,nZSearch};
    FUNC_INFO << "rangeSearch" << rangeSearch.x << rangeSearch.y << rangeSearch.z;

    double sigma2 = SQR(_data.hdr.noiseLevel);
    FUNC_INFO << "sigma2" << sigma2;
    double beta = 1.;  // should be close to 1
    double dofFactor;
    if ( dim3D )
        dofFactor = 7. /static_cast<double>(rangeSearch.x*rangeSearch.y*rangeSearch.z); // 7 parameters, 27 voxels
    else
        dofFactor = 7. /static_cast<double>(rangeSearch.x*rangeSearch.y); // 7 parameters, 27 voxels
    double heff2 = beta * sigma2 * dofFactor;

    // Fit spatial voxels (v) signals Y with polynomial coefficients (c) using a design matrix X:
    // Y_v = X_vc * B_c
    // where
    // Y_v is a concatenation of voxels in the region patch
    // X_vc has polynomial functions based upon (x,y,z). In each dimension, values are normalized (+-1) at edges of the FOV.
    // B_c = {(X^t*X)^-1}_cc' * (X^t)_c'v * Y_v = piX_cv * Y_v
    //
    ImageIO parImage;
    parImage.copyTemplateDimensions(this,4);  // s0, sx, sy, sz

    FUNC_INFO << "calculate parameter image";
    iPoint3D dim = getImageDimensions();
    int iVoxel=0;
    iPoint3D center, voxel;
    for (center.z=0, iVoxel=0; center.z<dim.z; center.z++)
    {
        for (center.y=0; center.y<dim.y; center.y++)
        {
            for (center.x=0; center.x<dim.x; center.x++, iVoxel++)
            {
                int nAverage = 0;
                double s0=0.;
                double sx=0.;   double x2=0.;
                double sy=0.;   double y2=0.;
                double sz=0.;   double z2=0.;
                iPoint3D voxelRel;
                for (voxelRel.z=-rangePatch.z; voxelRel.z<=rangePatch.z; voxelRel.z++)
                {
                    voxel.z = center.z + voxelRel.z;
                    for (voxelRel.y=-rangePatch.y; voxelRel.y<=rangePatch.y; voxelRel.y++)
                    {
                        voxel.y = center.y + voxelRel.y;
                        for (voxelRel.x=-rangePatch.x; voxelRel.x<=rangePatch.x; voxelRel.x++)
                        {
                            voxel.x = center.x + voxelRel.x;
                            if ( voxelInBounds(voxel) )
                            {
                                double signal = getStackValue(voxel,iTime);
                                s0 += signal;
                                sx += signal * voxelRel.x;      x2 += voxelRel.x * voxelRel.x;
                                sy += signal * voxelRel.y;      y2 += voxelRel.y * voxelRel.y;
                                if ( dim3D )
                                {
                                    sz += signal * voxelRel.z;
                                    z2 += voxelRel.z * voxelRel.z;
                                }
                                nAverage++;
                            }
                        } // voxelRel.x
                    } // voxelRel.y
                } // voxelRel.z
                double norm = 1.;
                if ( nAverage != 0 ) norm = static_cast<double>(nAverage);
                s0 /= norm;  x2 /= norm;  y2 /= norm; z2 /= norm;
                sx /= norm * x2;
                sy /= norm * y2;
                parImage.setStackValue(center,0,s0);
                parImage.setStackValue(center,1,sx);
                parImage.setStackValue(center,2,sy);
                if ( dim3D )
                {
                    sz /= norm * z2;
                    parImage.setStackValue(center,3,sz);
                }
            } // center.x
        } // center.y
    } // center.z
    utilIO::elapsedTime("Finding voxel parameters took",&_timer);

    FUNC_INFO << "fill NLM image";
    // Create a new volume that copies the old volume to start.
    // Then seach for similar voxels across a neighborhood of restricted rangeSearch (e.g., 3,3,3)
    ImageIO original;
    original.copyTemplateVolume(this,iTime);
    for (center.z=0; center.z<dim.z; center.z++)
    {
        for (center.y=0; center.y<dim.y; center.y++)
        {
            for (center.x=0; center.x<dim.x; center.x++)
            {
                iPoint3D voxelRel;
                double weightedSignal = 0.;
                double weightSum = 0.;
                for (voxelRel.z=-rangeSearch.z; voxelRel.z<=rangeSearch.z; voxelRel.z++)
                {
                    voxel.z = center.z + voxelRel.z;
                    for (voxelRel.y=-rangeSearch.y; voxelRel.y<=rangeSearch.y; voxelRel.y++)
                    {
                        voxel.y = center.y + voxelRel.y;
                        for (voxelRel.x=-rangeSearch.x; voxelRel.x<=rangeSearch.x; voxelRel.x++)
                        {
                            voxel.x = center.x + voxelRel.x;
                            if ( voxelInBounds(voxel) )
                            { // calculate distance between center patch and voxel patch
                                double signal = original.getStackValue(voxel,iTime);
                                double dist2 = SQR(parImage.getStackValue(center,0)-parImage.getStackValue(voxel,0))
                                        +      SQR(parImage.getStackValue(center,1)-parImage.getStackValue(voxel,1)) * .667
                                        +      SQR(parImage.getStackValue(center,2)-parImage.getStackValue(voxel,2)) * .667;
                                if ( dim3D ) dist2 += SQR(parImage.getStackValue(center,3)-parImage.getStackValue(voxel,3)) * .667;
                                double weight = qExp(-dist2/heff2);
                                // FUNC_INFO << "dist2 heff2 weight" << dist2 << heff2 << weight;
                                weightedSignal += weight * signal;
                                weightSum      += weight;
                            } // inbounds
                        } // voxelRel.x
                    } // voxelRel.y
                } // voxelRel.z
                // signal = original.getStackValue(center,iTime);
                // FUNC_INFO << "signal weightSum weightedSignal" << signal << weightSum << weightedSignal/weightSum;
                if ( weightSum != 0. ) weightedSignal /= weightSum;
                setStackValue(center,0,weightedSignal);
            } // center.x
        } // center.y
    } // center.z
    utilIO::elapsedTime("Finding weights took",&_timer);
    FUNC_EXIT;
}

#ifdef USE_FFTW
    void ImageIO::smoothVolumes(double fwhm)
    { // smooth by isotropic Gaussian with fwhm
        dPoint3D fwhm3d;
        fwhm3d.x = fwhm3d.y = fwhm3d.z = fwhm;
        smoothVolumes(fwhm3d);
    }
    void ImageIO::smoothVolumes(dPoint3D fwhm)
    { // smooth by isotropic Gaussian with fwhm
        for (int jt=0; jt<nt(); jt++)
            smoothVolume(jt, fwhm);
    }

    void ImageIO::smoothVolume(int iTime, double fwhm)
    { // smooth by isotropic Gaussian with fwhm
        FUNC_ENTER;
        if ( fwhm <= 0. ) return;
        dPoint3D fwhm3d;
        fwhm3d.x = fwhm3d.y = fwhm3d.z = fwhm;
        smoothVolume(iTime,fwhm3d);
    }

    void ImageIO::smoothVolume(int iTime, dPoint3D fwhm)
    { // set fwhm.z=0 for 2d smoothing
        /////////////////////
        /// This function is not thread-safe. Plans should be created in the calling thread.
        /// According to FFTW documentation, all functions but fftw_execute should be called from 1 thread at a time
        /////////////////////
        FUNC_ENTER << fwhm.x << fwhm.y << fwhm.z;
        if ( fwhm.x <= 0. || fwhm.y <= 0|| fwhm.z <= 0. ) return;
        iPoint3D dim = getImageDimensions();

        bool dim3D = fwhm.z != 0.;

        iPoint3D pad, padMax;
        padMax.x = dim.x/4;
        padMax.y = dim.y/4;
        padMax.z = dim.z/4;
        dPoint3D res = getImageResolution();
        pad.x = qMin(padMax.x,qCeil(static_cast<double>(2.*fwhm.x/res.x)));
        pad.y = qMin(padMax.y,qCeil(static_cast<double>(2.*fwhm.y/res.y)));
        pad.z = qMin(padMax.z,qCeil(static_cast<double>(2.*fwhm.z/res.z)));

        FUNC_INFO << "res" << res.x << res.y << res.z;
        FUNC_INFO << "fwhm" << fwhm.x << fwhm.y << fwhm.z;
        FUNC_INFO << "padMax" << padMax.x << padMax.y << padMax.z;
        FUNC_INFO << "pad" << pad.x << pad.y << pad.z;

        pad.x = pad.y = pad.z = 0;

        // Determine the 3D sigma_k
        // Multiply data by a k-space filter
        dPoint3D sigma2_k;
        sigma2_k.x = 2.355*2.355/4./PI/PI /fwhm.x/fwhm.x;
        sigma2_k.y = 2.355*2.355/4./PI/PI /fwhm.y/fwhm.y;
        if ( dim3D )
            sigma2_k.z = 2.355*2.355/4./PI/PI /fwhm.z/fwhm.z;
        else
            sigma2_k.z = 1.;

        dcVector volume = convertToComplexVolume(iTime);
        // inverse FT in x, y, and z
        utilMath::FFTW1D_Volume(dim, pad, volume, 0, false); // x
        utilMath::FFTW1D_Volume(dim, pad, volume, 1, false); // y
        if ( dim3D )
            utilMath::FFTW1D_Volume(dim, pad, volume, 2, false); // z

        dPoint3D resolution = getImageResolution();
        dPoint3D fov2;
        fov2.x = static_cast<double>(resolution.x*resolution.x) * dim.x*dim.x;
        fov2.y = static_cast<double>(resolution.y*resolution.y) * dim.y*dim.y;
        fov2.z = static_cast<double>(resolution.z*resolution.z) * dim.z*dim.z;
        dPoint3D sigma2;
        sigma2.x = sigma2_k.x * fov2.x; // sigma2 now dimensionless
        sigma2.y = sigma2_k.y * fov2.y;
        sigma2.z = sigma2_k.z * fov2.z;
        FUNC_INFO << "sigma2" << sigma2.x << sigma2.y << sigma2.z;
        iPoint3D voxel;
        for (voxel.z=0; voxel.z<dim.z; voxel.z++)
        {
            double z = voxel.z - dim.z/2 + 0.5;
            for (voxel.y=0; voxel.y<dim.y; voxel.y++)
            {
                double y = voxel.y - dim.y/2 + 0.5;
                for (voxel.x=0; voxel.x<dim.x; voxel.x++)
                {
                    double x = voxel.x - dim.x/2 + 0.5;
                    double scale = qExp(-x*x/2./sigma2.x)
                            *      qExp(-y*y/2./sigma2.y);
                    if ( dim3D )
                        scale *= qExp(-z*z/2./sigma2.z);
                    volume[utilIO::index1d(dim, voxel)].real *= scale;
                    volume[utilIO::index1d(dim, voxel)].imag *= scale;
                } // x
            } // y
        } // z

        // forward FT in x, y, and z
        utilMath::FFTW1D_Volume(dim, pad, volume, 0, true); // x
        utilMath::FFTW1D_Volume(dim, pad, volume, 1, true); // y
        if ( dim3D )
            utilMath::FFTW1D_Volume(dim, pad, volume, 2, true); // z
        // Put the result back into the original dataset.
        for (voxel.z=0; voxel.z<dim.z; voxel.z++)
        {
            for (voxel.y=0; voxel.y<dim.y; voxel.y++)
            {
                for (voxel.x=0; voxel.x<dim.x; voxel.x++)
                {
                    int index = utilIO::index1d(dim, voxel);
                    double value = qSqrt(volume[index].real*volume[index].real + volume[index].imag*volume[index].imag);
                    // double value = volume[index].real;
                    setStackValue(voxel,iTime,value);
                } // x
            } // y
        } // z
        FUNC_EXIT;
    }

    void ImageIO::smoothVolumeInThread(int iTime, dPoint3D fwhm, iPoint3D pad,
                                       dComplex *dcDataX, fftw_complex *fftwDataX, fftw_plan planBackwardX, fftw_plan planForwardX,
                                       dComplex *dcDataY, fftw_complex *fftwDataY, fftw_plan planBackwardY, fftw_plan planForwardY,
                                       dComplex *dcDataZ, fftw_complex *fftwDataZ, fftw_plan planBackwardZ, fftw_plan planForwardZ)
    { // set fwhm.z=0 for 2d smoothing
        /////////////////////
        /// Whenever smoothing happens in a thread separate from the main thread, plans should be created in the calling thread.
        /// According to FFTW documentation, all functions but fftw_execute should be called from 1 thread at a time
        /////////////////////
        FUNC_ENTER << fwhm.x << fwhm.y << fwhm.z;
        if ( fwhm.x <= 0. || fwhm.y <= 0|| fwhm.z <= 0. ) return;

        bool dim3D = fwhm.z != 0.;

        // Determine the 3D sigma_k
        // Multiply data by a k-space filter
        dPoint3D sigma2_k;
        sigma2_k.x = 2.355*2.355/4./PI/PI /fwhm.x/fwhm.x;
        sigma2_k.y = 2.355*2.355/4./PI/PI /fwhm.y/fwhm.y;
        if ( dim3D )
            sigma2_k.z = 2.355*2.355/4./PI/PI /fwhm.z/fwhm.z;
        else
            sigma2_k.z = 1.;

        dcVector volume = convertToComplexVolume(iTime);
        // inverse FT in x, y, and z
        iPoint3D dim = getImageDimensions();
        FUNC_INFO << "pre";
        utilMath::FFTW1DVolumeInThread(dim, volume, 0, false, pad,
                                       dcDataX, fftwDataX, planBackwardX, planForwardX,
                                       dcDataY, fftwDataY, planBackwardY, planForwardY,
                                       dcDataZ, fftwDataZ, planBackwardZ, planForwardZ); // x
        utilMath::FFTW1DVolumeInThread(dim, volume, 1, false, pad,
                                       dcDataX, fftwDataX, planBackwardX, planForwardX,
                                       dcDataY, fftwDataY, planBackwardY, planForwardY,
                                       dcDataZ, fftwDataZ, planBackwardZ, planForwardZ); // y
        if ( dim3D )
            utilMath::FFTW1DVolumeInThread(dim, volume, 2, false, pad,
                                           dcDataX, fftwDataX, planBackwardX, planForwardX,
                                           dcDataY, fftwDataY, planBackwardY, planForwardY,
                                           dcDataZ, fftwDataZ, planBackwardZ, planForwardZ); // z

        dPoint3D resolution = getImageResolution();
        dPoint3D fov2;
        fov2.x = static_cast<double>(resolution.x*resolution.x) * dim.x*dim.x;
        fov2.y = static_cast<double>(resolution.y*resolution.y) * dim.y*dim.y;
        fov2.z = static_cast<double>(resolution.z*resolution.z) * dim.z*dim.z;
        dPoint3D sigma2;
        sigma2.x = sigma2_k.x * fov2.x; // sigma2 now dimensionless
        sigma2.y = sigma2_k.y * fov2.y;
        sigma2.z = sigma2_k.z * fov2.z;
        FUNC_INFO << "sigma2" << sigma2.x << sigma2.y << sigma2.z;
        iPoint3D voxel;
        for (voxel.z=0; voxel.z<dim.z; voxel.z++)
        {
            double z = voxel.z - dim.z/2 + 0.5;
            for (voxel.y=0; voxel.y<dim.y; voxel.y++)
            {
                double y = voxel.y - dim.y/2 + 0.5;
                for (voxel.x=0; voxel.x<dim.x; voxel.x++)
                {
                    double x = voxel.x - dim.x/2 + 0.5;
                    double scale = qExp(-x*x/2./sigma2.x)
                            *      qExp(-y*y/2./sigma2.y);
                    if ( dim3D )
                        scale *= qExp(-z*z/2./sigma2.z);
                    volume[utilIO::index1d(dim, voxel)].real *= scale;
                    volume[utilIO::index1d(dim, voxel)].imag *= scale;
                } // x
            } // y
        } // z

        FUNC_INFO << "post";
        // forward FT in x, y, and z
        utilMath::FFTW1DVolumeInThread(dim, volume, 0, true, pad,
                                       dcDataX, fftwDataX, planBackwardX, planForwardX,
                                       dcDataY, fftwDataY, planBackwardY, planForwardY,
                                       dcDataZ, fftwDataZ, planBackwardZ, planForwardZ); // x
        utilMath::FFTW1DVolumeInThread(dim, volume, 1, true, pad,
                                       dcDataX, fftwDataX, planBackwardX, planForwardX,
                                       dcDataY, fftwDataY, planBackwardY, planForwardY,
                                       dcDataZ, fftwDataZ, planBackwardZ, planForwardZ); // y
        if ( dim3D )
            utilMath::FFTW1DVolumeInThread(dim, volume, 2, true, pad,
                                           dcDataX, fftwDataX, planBackwardX, planForwardX,
                                           dcDataY, fftwDataY, planBackwardY, planForwardY,
                                           dcDataZ, fftwDataZ, planBackwardZ, planForwardZ); // z
        // Put the result back into the original dataset.
        for (voxel.z=0; voxel.z<dim.z; voxel.z++)
        {
            for (voxel.y=0; voxel.y<dim.y; voxel.y++)
            {
                for (voxel.x=0; voxel.x<dim.x; voxel.x++)
                {
                    int index = utilIO::index1d(dim, voxel);
                    double value = qSqrt(volume[index].real*volume[index].real + volume[index].imag*volume[index].imag);
                    // double value = volume[index].real;
                    setStackValue(voxel,iTime,value);
                } // x
            } // y
        } // z
        FUNC_EXIT;
    }
#endif

void ImageIO::smoothVolumesNoFFT(double fwhm)
{ // smooth by isotropic Gaussian with fwhm
    dPoint3D fwhm3d;
    fwhm3d.x = fwhm3d.y = fwhm3d.z = fwhm;
    smoothVolumesNoFFT(fwhm3d);
}
void ImageIO::smoothVolumesNoFFT(dPoint3D fwhm)
{ // smooth by isotropic Gaussian with fwhm
    for (int jt=0; jt<nt(); jt++)
        smoothVolumeNoFFT(jt, fwhm);
}

void ImageIO::smoothVolumeNoFFT(int iTime, double fwhm)
{ // smooth by isotropic Gaussian with fwhm
    if ( fwhm <= 0. ) return;
    dPoint3D fwhm3d;
    fwhm3d.x = fwhm3d.y = fwhm3d.z = fwhm;
    smoothVolumeNoFFT(iTime,fwhm3d);
}
void ImageIO::smoothVolumeNoFFT(int iTime, dPoint3D fwhm)
{  // set fwhm.z=0 for 2d smoothing
    if ( fwhm.x <= 0. || fwhm.y <= 0|| fwhm.z <= 0. ) return;

    iPoint3D dim = getImageDimensions();
    dVector volume; volume.resize(dim.x*dim.y*dim.z);
    // copy the volume
    iPoint3D voxel;
    for (voxel.z=0; voxel.z<dim.z; voxel.z++)
    {
        for (voxel.y=0; voxel.y<dim.y; voxel.y++)
        {
            for (voxel.x=0; voxel.x<dim.x; voxel.x++)
                volume[utilIO::index1d(dim, voxel)] = getStackValue(voxel,iTime);
        }
    }

    dPoint3D distanceScale;
    distanceScale.x = fwhm.x;
    distanceScale.y = fwhm.y;
    distanceScale.z = fwhm.z;

    dPoint3D resolution = getImageResolution();
    iPoint3D voxelDist;

    // Define the limits in x, y, and z
    voxelDist.x = qMin(30,qCeil(static_cast<double>(distanceScale.x/resolution.x)));
    voxelDist.y = qMin(30,qCeil(static_cast<double>(distanceScale.y/resolution.y)));
    voxelDist.z = qMin(30,qCeil(static_cast<double>(distanceScale.z/resolution.z)));

    for (voxel.z=0; voxel.z<dim.z; voxel.z++)
    {
        for (voxel.y=0; voxel.y<dim.y; voxel.y++)
        {
            for (voxel.x=0; voxel.x<dim.x; voxel.x++)
            {
                double sum = 0.; double weightSum = 0.;
                iPoint3D voxelRel, neighbor;
                for (voxelRel.z=-voxelDist.z; voxelRel.z<=voxelDist.z; voxelRel.z++)
                {
                    neighbor.z = voxel.z + voxelRel.z;
                    dPoint3D weight3d;
                    weight3d.z = utilMath::GaussKernel(voxelRel.z*resolution.z,fwhm.z);
                    for (voxelRel.y=-voxelDist.y; voxelRel.y<=voxelDist.y; voxelRel.y++)
                    {
                        neighbor.y = voxel.y + voxelRel.y;
                        weight3d.y = utilMath::GaussKernel(voxelRel.y*resolution.y,fwhm.y);
                        for (voxelRel.x=-voxelDist.x; voxelRel.x<=voxelDist.x; voxelRel.x++)
                        {
                            neighbor.x = voxel.x + voxelRel.x;
                            weight3d.x = utilMath::GaussKernel(voxelRel.x*resolution.x,fwhm.x);
                            if ( neighbor.x >= 0 && neighbor.x < dim.x &&
                                 neighbor.y >= 0 && neighbor.y < dim.y &&
                                 neighbor.z >= 0 && neighbor.z < dim.z )
                            {
                                double weight = weight3d.x * weight3d.y * weight3d.z;
                                weightSum += weight;
                                sum += weight * volume[utilIO::index1d(dim, neighbor)];
                            }

                        } // voxelRel.x
                    }  // voxelRel.y
                }  // voxelRel.z
                if ( weightSum != 0. )
                    setStackValue(voxel,iTime,sum/weightSum);
            } // voxel.x
        } // voxel.y
    } // voxel.z
}

void ImageIO::dotProductWithKernelStride1(ImageIO *inputImage, int iTime, iPoint4D kernelSize, kernelSet kernel)
{
    FUNC_ENTER << iTime << inputImage->nt();
    // S(x,y,z) = S(x+i,y+j,z+k,d) * kernel(i,j,k,d)
    iPoint4D dim = getDimensions();
    FUNC_INFO << "conv dim" << dim.x << dim.y << dim.z << kernelSize.t;
    if ( inputImage->nt() != kernelSize.t )
    {
        qInfo() << "Error: kernel depth " << kernelSize.t << " != input image depth " << inputImage->nt();
        qFatal("Fatal error in ImageIO::dotProductWithKernelStride1");
    }
    else if ( iTime >= nt() )
    {
        qInfo() << "Error: depth index " << iTime << " >= image depth of stack" << nt();
        qFatal("Fatal error in ImageIO::dotProductWithKernelStride1");
    }

    // sum over image: stride 1
    iPoint3D voxel;
    for (voxel.z=0; voxel.z<dim.z; voxel.z++)
    {
        for (voxel.y=0; voxel.y<dim.y; voxel.y++)
        {
            for (voxel.x=0; voxel.x<dim.x; voxel.x++)
            {
                // for each voxel, sum over the kernel including the full depth
                double sum = 0.;
                iPoint4D neighbor;
                for (int jKernel=0; jKernel<kernel.size(); jKernel++)
                {
                    kernelVoxel k = kernel[jKernel];
                    neighbor.x = voxel.x + k.x;
                    neighbor.y = voxel.y + k.y;
                    neighbor.z = voxel.z + k.z;
                    neighbor.t = k.t;
                    bool inBounds = neighbor.x >= 0 && neighbor.x < dim.x
                            &&      neighbor.y >= 0 && neighbor.y < dim.y
                            &&      neighbor.z >= 0 && neighbor.z < dim.z;
                    if ( inBounds ) // save full 4D convolution in 1 volume specified by iTime (kernel(iTime))
                    {
//                        FUNC_INFO << "neighbor" << neighbor.x << neighbor.y << neighbor.z << neighbor.t;
                        sum += k.value * inputImage->getStackValue(neighbor);
                    }
                }
//                FUNC_INFO << "set stack value" << voxel.x << voxel.y << voxel.z << iTime;
                setStackValue(voxel,iTime,sum);
            } // voxel.x
        } // voxel.y
    } // voxel.z
    FUNC_EXIT;
}

void ImageIO::poolStride2(ImageIO *inputImage, int iTime, int iStride, kernelSet kernel, bool maxPool)
{
    FUNC_ENTER;
    iPoint3D dimIn  = inputImage->getImageDimensions();
    iPoint3D dimOut = getImageDimensions();
    int nTIn        = inputImage->nt();
    int nTOut       = nt();
    bool badDimMatch = (dimIn.x != iStride*dimOut.x) || (dimIn.y != iStride*dimOut.y)
            ||         (dimIn.z != iStride*dimOut.z) || (nTIn != nTOut );
    if ( badDimMatch )
    {
        QString error = QString("Error: input %1x%2x%3 (%4) not compatible without output %5x%6x%7 (%8) for stride %9")
                .arg(dimIn.x).arg(dimIn.y).arg(dimIn.z).arg(nTIn)
                .arg(dimOut.x).arg(dimOut.y).arg(dimOut.z).arg(nTOut).arg(iStride);
        qInfo() << error;
        qFatal("Output should be smaller by stride and depths should match in ImageIO::poolStride2");
    }
    else if ( iTime >= nt() )
    {
        qInfo() << "Error: depth index " << iTime << " >= image depth of stack" << nt();
        qFatal("Depths must always matchin ImageIO::poolStride2");
    }

    iPoint3D voxel;
    for (voxel.z=0; voxel.z<dimIn.z; voxel.z+=iStride)
    {
        for (voxel.y=0; voxel.y<dimIn.y; voxel.y+=iStride)
        {
            for (voxel.x=0; voxel.x<dimIn.x; voxel.x+=iStride)
            {
                // create a set of values across the 3D volume
                dVector setOfValues;
                iPoint4D neighbor;
                for (int jKernel=0; jKernel<kernel.size(); jKernel++)
                {
                    kernelVoxel k = kernel[jKernel];
                    neighbor.x = voxel.x + k.x;
                    neighbor.y = voxel.y + k.y;
                    neighbor.z = voxel.z + k.z;
                    neighbor.t = iTime;
                    if ( neighbor.x >= 0 && neighbor.x < dimIn.x &&
                         neighbor.y >= 0 && neighbor.y < dimIn.y &&
                         neighbor.z >= 0 && neighbor.z < dimIn.z )
                    {
                        double value = inputImage->getStackValue(neighbor);
                        setOfValues.append(value);
                    }
                } // kernel
                // For each t (depth), pool values and reduce dimension
                double value;
                if ( maxPool )
                {
                    value = -1.e10;
                    for (int jSet=0; jSet<setOfValues.size(); jSet++)
                        if ( setOfValues.at(jSet) < value ) value = setOfValues.at(jSet);
                }
                else
                {
                    value = 0.;
                    for (int jSet=0; jSet<setOfValues.size(); jSet++)
                        value += setOfValues.at(jSet);
                    value /= static_cast<double>(setOfValues.size());
                }
                iPoint3D voxelOut = {voxel.x/2,voxel.y/2,voxel.z/2};
                // save the 3D sum in a single time point = kernel(iTime)
                setStackValue(voxelOut,iTime,value);
            } // x
        } // y
    } // z
    FUNC_EXIT;
}

void ImageIO::cropUsingOverlay(bitMap mask)
{
    iPoint3D dimImage = getImageDimensions();
    iPoint3D dimMask  = mask.getDimensions();
    if ( (dimImage.x != dimMask.x) ||
         (dimImage.y != dimMask.y) ||
         (dimImage.z != dimMask.z) )
        qFatal("Error: image and mask dimensions don't match in ImageIO::cropUsingOverlay");
    for (int jt=0; jt<nt(); jt++)
    {
        for (int jVoxel=0; jVoxel<voxelsPerVolume(); jVoxel++)
        {
            if ( mask.isVoxelZero(jVoxel) )
                setStackValue(jVoxel,jt,0.);
        }
    }
}

void ImageIO::dilateOrErode(bool dilate, int iTime)
{
    // This could be modified to make a 2D alternative
    iPoint3D nRel;
    nRel.x = nRel.y = nRel.z = 1;

    iPoint3D dim = getImageDimensions();
    iPoint3D voxel;
    ImageIO original;
    original.copyTemplateVolume(this,iTime);
    for (voxel.z=0; voxel.z<dim.z; voxel.z++)
    {
        for (voxel.y=0; voxel.y<dim.y; voxel.y++)
        {
            for (voxel.x=0; voxel.x<dim.x; voxel.x++)
            {
                iPoint3D neighbor, rel;
                double extrema = 1.e10;
                if ( dilate ) extrema = -1.e10;
                // go through the local neighborhood
                for (rel.z=-nRel.z; rel.z<=nRel.z; rel.z++)
                {
                    neighbor.z = voxel.z + rel.z;
                    for (rel.y=-nRel.y; rel.y<=nRel.y; rel.y++)
                    {
                        neighbor.y = voxel.y + rel.y;
                        for (rel.x=-nRel.x; rel.x<=nRel.x; rel.x++)
                        {
                            neighbor.x = voxel.x + rel.x;
                            // Ignore indices that are out of bounds.
                            bool ignore = (neighbor.x < 0 || neighbor.x >= nx() )
                                    ||    (neighbor.y < 0 || neighbor.y >= ny() )
                                    ||    (neighbor.z < 0 || neighbor.z >= nz() );
                            if ( !ignore )
                            {
                                double value = original.getStackValue(neighbor,0);
                                if ( dilate )
                                    extrema = qMax(value,extrema);
                                else
                                    extrema = qMin(value,extrema);
                            }
                        } // jRel.x
                    } // jRel.y
                } // jRel.z
                // Assign the extrema to the local neighborhood
                for (rel.z=-nRel.z; rel.z<=nRel.z; rel.z++)
                {
                    neighbor.z = voxel.z + rel.z;
                    for (rel.y=-nRel.y; rel.y<=nRel.y; rel.y++)
                    {
                        neighbor.y = voxel.y + rel.y;
                        for (rel.x=-nRel.x; rel.x<=nRel.x; rel.x++)
                        {
                            neighbor.x = voxel.x + rel.x;
                            // Ignore indices that are out of bounds.
                            bool ignore = (neighbor.x < 0 || neighbor.x >= nx() )
                                    ||    (neighbor.y < 0 || neighbor.y >= ny() )
                                    ||    (neighbor.z < 0 || neighbor.z >= nz() );
                            if ( !ignore )
                                setStackValue(neighbor,iTime,extrema);
                        } // jRel.x
                    } // jRel.y
                } // jRel.z
            } // voxel.x
        } // voxel.y
    } // voxel.z
}

void ImageIO::setParallelCOMCoordinates(dPoint3D COM)
{
    iPoint3D dimension = getImageDimensions();
    dPoint3D resolution = getImageResolution();
    // dimensions: 3d to 4d
    _data.hdr.dim.x = dimension.x;
    _data.hdr.dim.y = dimension.y;
    _data.hdr.dim.z = dimension.z;
    // resolution: 3d to 4d
    _data.hdr.resolution.x = resolution.x;
    _data.hdr.resolution.y = resolution.y;
    _data.hdr.resolution.z = resolution.z;
    // origin
    dPoint3D origin;
    origin.x = -COM.x;
    origin.y = -COM.y;
    origin.z = -COM.z;
    // direction: align coordinates with indices
    iPoint3D direction;
    direction.x = direction.y = direction.z = 1;
    // transformation
    setParallelCoordinates(_data.hdr.resolution, origin, direction);
    // also make sure the time resolution is not 0
    if ( _data.hdr.resolution.t == 0. ) _data.hdr.resolution.t = 1.;
}

void ImageIO::setParallelCenteredCoordinates()
{
    iPoint3D dimension = getImageDimensions();
    dPoint3D resolution = getImageResolution();
    setParallelCenteredCoordinates(dimension, resolution);
}
void ImageIO::setParallelCenteredCoordinates(iPoint4D dimension, dPoint4D resolution)
{
    iPoint3D dim;  dim.x = dimension.x;  dim.y = dimension.y;  dim.z = dimension.z;
    dPoint3D res;  res.x = resolution.x;  res.y = resolution.y;  res.z = resolution.z;
    setParallelCenteredCoordinates(dim, res);
}
void ImageIO::setParallelCenteredCoordinates(iPoint4D dimension, dPoint3D resolution)
{
    iPoint3D dim;  dim.x = dimension.x;  dim.y = dimension.y;  dim.z = dimension.z;
    setParallelCenteredCoordinates(dim, resolution);
}
void ImageIO::setParallelCenteredCoordinates(iPoint3D dimension, dPoint3D resolution)
{ // dimension and resolution should be changed
    // dimensions: 3d to 4d
    _data.hdr.dim.x = dimension.x;
    _data.hdr.dim.y = dimension.y;
    _data.hdr.dim.z = dimension.z;

    // resolution: 3d to 4d
    _data.hdr.resolution.x = resolution.x;
    _data.hdr.resolution.y = resolution.y;
    _data.hdr.resolution.z = resolution.z;

    // origin
    dPoint3D origin;
    origin.x = -static_cast<double>(nx()-1) * _data.hdr.resolution.x / 2.;
    origin.y = -static_cast<double>(ny()-1) * _data.hdr.resolution.y / 2.;
    origin.z = -static_cast<double>(nz()-1) * _data.hdr.resolution.z / 2.;

    // direction: align coordinates with indices
    iPoint3D direction;
    direction.x = direction.y = direction.z = 1;

    // transformation
    setParallelCoordinates(_data.hdr.resolution, origin, direction);

    // also make sure the time resolution is not 0
    if ( _data.hdr.resolution.t == 0. ) _data.hdr.resolution.t = 1.;
}
iPoint3D ImageIO::determineDirectionMap(dMatrix ijk_to_xyz)
{
    iPoint3D directionMap={0,1,2}; // map of (ijk) onto (xyz)
    if ( qAbs(ijk_to_xyz[0][0]) > qAbs(ijk_to_xyz[1][0]) &&
         qAbs(ijk_to_xyz[0][0]) > qAbs(ijk_to_xyz[2][0]) )
        directionMap.x = 0;
    else if ( qAbs(ijk_to_xyz[1][0]) > qAbs(ijk_to_xyz[2][0]) )
        directionMap.x = 1;
    else
        directionMap.x = 2;
    // Y indices
    if ( qAbs(ijk_to_xyz[0][1]) > qAbs(ijk_to_xyz[1][1]) &&
         qAbs(ijk_to_xyz[0][1]) > qAbs(ijk_to_xyz[2][1]) )
        directionMap.y = 0;
    else if ( qAbs(ijk_to_xyz[1][1]) > qAbs(ijk_to_xyz[2][1]) )
        directionMap.y = 1;
    else
        directionMap.y = 2;
    // Z indices
    if ( qAbs(ijk_to_xyz[0][2]) > qAbs(ijk_to_xyz[1][2]) &&
         qAbs(ijk_to_xyz[0][2]) > qAbs(ijk_to_xyz[2][2]) )
        directionMap.z = 0;
    else if ( qAbs(ijk_to_xyz[1][2]) > qAbs(ijk_to_xyz[2][2]) )
        directionMap.z = 1;
    else
        directionMap.z = 2;
    return directionMap;
}

void ImageIO::initLoaded (int nThreads)
{ // not that "threads" could also be replaced by "volumes" in some cases
    _loaded.fill(false,nThreads);
//    FUNC_INFO << _loaded;
}
void ImageIO::setLoaded()
{
//    FUNC_ENTER << _loaded;
    int iThread = _loaded.lastIndexOf(false);
    if ( iThread < 0 )
    {
        QString errorText = QString("Error: attempted to access index %1 in setLoaded() with %2 threads associated with file %3").
                arg(iThread).arg(_loaded.size()).arg(getAbsoluteName());
        QMessageBox msgBox; msgBox.setText(errorText);  msgBox.setIcon(QMessageBox::Critical);  msgBox.exec();
        exit(1);
    }
    _loaded[iThread] = true;
//    FUNC_EXIT << _loaded;
}
bool ImageIO::noneLoaded()
{
    return !_loaded.contains(true);
}
bool ImageIO::isLoaded()
{
    bool loaded;
    if ( _loaded.size() == 0 )
        loaded = false;
    else
        loaded = !_loaded.contains(false);   // true if all values in the vector are true
    return loaded;
};

void ImageIO::setCurrentVolume(int iTime)
{
    _iTime = iTime;
    // If values fall outside the range (0 ..., nt-1), adjust using a cyclical pattern (assume +/- incrementing)
    if ( _iTime >= nt() )
        _iTime = 0;
    else if ( _iTime < 0 )
        _iTime = nt() - 1;
};

void ImageIO::resetDirName(QString newDirName)
{
    FUNC_ENTER << newDirName;
    QString localFileName = getFileName();
    FUNC_INFO << "local" << localFileName;
    QString fullName = newDirName + "/" + localFileName;
    FUNC_INFO << "full" << fullName;
    setFileName(fullName);
}

void ImageIO::resetFileName(ImageIO *templateImage, QString newName)
{
    FUNC_ENTER << newName;
    _imageFileName = newName;
    QString tempName  = templateImage->getAbsoluteName();
    QString dirName  = utilString::getDirectoryName(tempName);
    QString fullName = dirName + "/" + newName;
    setFileName(fullName);
    FUNC_EXIT;
}

void ImageIO::copyTemplate(ImageIO *templateImage, QString name)
{
    copyTemplate(templateImage);
    resetFileName(templateImage,name);
}
void ImageIO::copyTemplate(ImageIO *templateImage)
{// Copy the image template including data
    copyTemplateDimensions(templateImage, templateImage->nt());
    if ( is3Vector() )
    {
        if ( templateImage->getCategory() != category_vector )
            qFatal("Programming error: attempting to copy scale image into vector image");
        for ( int jt=0; jt<nt(); jt++ )
        {
            for ( int jp=0; jp<_data.hdr.points.z; ++jp )
                _data.timePoints[jt].f3[jp] = templateImage->_data.timePoints[jt].f3[jp];
        }
    }
    else
    {
        if ( templateImage->getCategory() == category_vector )
            qFatal("Programming error: attempting to copy vector image into scalar image");
        for ( int jt=0; jt<nt(); jt++ )
        {
            for ( int jp=0; jp<_data.hdr.points.z; ++jp )
                _data.timePoints[jt].f1[jp] = templateImage->_data.timePoints[jt].f1[jp];
        }
    }
    _data.empty = false;
}
void ImageIO::copyTemplateVolume(ImageIO *templateImage, int iVolume, QString name)
{
    copyTemplateVolume(templateImage,iVolume);
    resetFileName(templateImage,name);
}
void ImageIO::copyTemplateVolume(ImageIO *templateImage, int iVolume)
{// Copy the image template volume, including data
    FUNC_INFO << "ptr2 =" << this << templateImage;
    copyTemplateDimensions(templateImage, 1);
    if ( is3Vector() )
    {
        if ( templateImage->getCategory() != category_vector )
            qFatal("Programming error: attempting to copy scale image into vector image");
        for ( int jp=0; jp<_data.hdr.points.z; ++jp )
            _data.timePoints[0].f3[jp] = templateImage->_data.timePoints[iVolume].f3[jp];
    }
    else
    {
        if ( templateImage->getCategory() == category_vector )
            qFatal("Programming error: attempting to copy vector image into scalar image");
        for ( int jp=0; jp<_data.hdr.points.z; ++jp )
            _data.timePoints[0].f1[jp] = templateImage->_data.timePoints[iVolume].f1[jp];
    }
    _data.hdr.noiseLevel = templateImage->_data.hdr.noiseLevel;
    _data.empty = false;
}
void ImageIO::copyTemplateDimensions(ImageIO *templateImage, int nTime)
{ // Copy the image template, potentially change time points, and fill with zeroes
    FUNC_ENTER << nTime;
    resetFileName(templateImage,templateImage->getFileName());
    setDimensions(templateImage->getImageDimensions(), nTime);
    _data.hdr.resolution = templateImage->getResolution();
    setIjkToXyz(templateImage->getIjkToXyz());
    // copy the name and ID, as often a dataset is modified using this function to initialize it
    setFileID(templateImage->getFileID());
    FUNC_EXIT << templateImage->getFileID();
}

void ImageIO::setDimensions(const iPoint4D dim)
{
    setDimensions(dim.x, dim.y, dim.z, dim.t);
}
void ImageIO::setDimensions(const iPoint3D dim, const int nt)
{
    setDimensions(dim.x, dim.y, dim.z, nt);
}
void ImageIO::setDimensions(const iPoint3D dim)
{
    setDimensions(dim.x, dim.y, dim.z, nt());
}
void ImageIO::setDimensions(const int nX, const int nY, const int nZ, const int nT)
{
    FUNC_ENTER << nX << nY << nZ << nT;
    _data.hdr.dim.x = nX;  _data.hdr.dim.y = nY;  _data.hdr.dim.z = nZ;  _data.hdr.dim.t = nT;
    _data.hdr.points.x = nx();
    _data.hdr.points.y = nx() * ny();
    _data.hdr.points.z = nx() * ny() * nz();
    _data.hdr.points.t = nx() * ny() * nz() * nT;
    _data.color.reserve(_data.hdr.points.z);
    _data.color.resize(_data.hdr.points.z); // fill 3d color with undefined
    _data.color.fill(Color_Undefined,_data.hdr.points.z); // fill 3d color with undefined
    _data.timePoints.resize(nT);
    // this function might be used to create maps within the program, in which case the category defines the data type
    FUNC_INFO << "category" << _category << category_vector << category_RGB;
    if ( _category == category_vector )
        _data.hdr.dataType = 3;
    else if ( _category == category_RGB )
        _data.hdr.dataType = 4;
    else if ( _category == category_complex )
        _data.hdr.dataType = 1;
    else
        _data.hdr.dataType = 0;
    if ( _data.hdr.dataType == 0 ) // scalar
    {
        for ( int jt=0; jt<nT; jt++ )
        {
            _data.timePoints[jt].f1.reserve(_data.hdr.points.z);
            _data.timePoints[jt].f1.fill(0.,_data.hdr.points.z);       // fill all time points with 0
        }
    }
    else if ( _data.hdr.dataType == 3 || _data.hdr.dataType == 4 ) // vector
    {
        for ( int jt=0; jt<nT; jt++ )
            _data.timePoints[jt].f3.fill({0.,0.,0.},_data.hdr.points.z);       // fill all time points with 0
    }
    else // complex
    {
        for ( int jt=0; jt<nT; jt++ )
            _data.timePoints[jt].fc.fill({0.,0.},_data.hdr.points.z);       // fill all time points with 0
    }
    if ( nT > 1 )
    {
        _data.averageVolume.reserve(_data.hdr.points.z);
        _data.averageVolume.fill(0.,_data.hdr.points.z);
    }
    _data.empty = true;
    if ( _data.hdr.ijk_to_xyz.size() != 4 )
    {
        // Create a 4x4 identity matrix for the default transformation IF one has not been previously allocated
        _data.hdr.ijk_to_xyz.resize(4);
        _data.hdr.xyz_to_ijk.resize(4);
        for (int j=0; j<4; j++)
        {
            _data.hdr.ijk_to_xyz[j].fill(0.,4);
            _data.hdr.ijk_to_xyz[j][j] = 1.;
        }
        _data.hdr.xyz_to_ijk = _data.hdr.ijk_to_xyz;
    }
    setCurrentVolume(0);
    FUNC_EXIT;
}

dPoint3D ImageIO::getFOV()
{
    dPoint3D fov;
    fov.x = nx() * _data.hdr.resolution.x;
    fov.y = ny() * _data.hdr.resolution.y;
    fov.z = nz() * _data.hdr.resolution.z;
    return fov;
}
dPoint3D ImageIO::getOrigin()
{
    dPoint3D origin;
    origin.x = _data.hdr.ijk_to_xyz[0][3];
    origin.y = _data.hdr.ijk_to_xyz[1][3];
    origin.z = _data.hdr.ijk_to_xyz[2][3];
    return origin;
}
iPoint3D ImageIO::getDirection()
{
    iPoint3D direction;
    // Set the directions based upon the transformation matrix.
    if ( _data.hdr.ijk_to_xyz[0][0] > 0.  )
        direction.x = 1;
    else
        direction.x = -1;
    if ( _data.hdr.ijk_to_xyz[1][1] > 0.  )
        direction.y = 1;
    else
        direction.y = -1;
    if ( _data.hdr.ijk_to_xyz[2][2] > 0.  )
        direction.z = 1;
    else
        direction.z = -1;
    return direction;
}
void ImageIO::setOrigin(dPoint3D origin)
{
    _data.hdr.ijk_to_xyz[0][3] = origin.x;
    _data.hdr.ijk_to_xyz[1][3] = origin.y;
    _data.hdr.ijk_to_xyz[2][3] = origin.z;
    // update the transformations
    setIjkToXyz(_data.hdr.ijk_to_xyz);
}

void ImageIO::setResolution(dPoint4D resolution)
{
    // The resolution has changed, and so the coordinate transformation also must change.
    // The following calls set the resolution and defines centered coordinates.
    setParallelCenteredCoordinates(_data.hdr.dim,resolution);
}
void ImageIO::setResolution(dPoint3D resolution)
{
    // The resolution has changed, and so the coordinate transformation also must change.
    // The following calls set the resolution and defines centered coordinates.
    setParallelCenteredCoordinates(_data.hdr.dim,resolution);
}
void ImageIO::appendImageFileToTimeData(ImageIO imageFile)
{
    for (int jTime=0; jTime<imageFile.nt(); jTime++)
        _data.timePoints.append(imageFile._data.timePoints[jTime]);
    _data.hdr.dim.t = _data.timePoints.size();
    _data.hdr.points.t = nx() * ny() * nz() * nt();
}
void ImageIO::deleteFirstVolume()
{
    if ( _data.timePoints.size() == 0 )
        qFatal("ImageIO::removeFirstVolume programming error - Attempt to remove 1st volume of empty vector");
    _data.timePoints.removeFirst();
    _data.hdr.dim.t = _data.timePoints.size();
    _data.hdr.points.t = nx() * ny() * nz() * nt();
}
void ImageIO::setStackValue(const int i3d, const int iT, const Point3D value)
{ // this is for vectors
    CheckVoxelIndices(i3d, "setStackValue");
    QMutex mutex;
    mutex.lock();
    _data.timePoints[iT].f3[i3d] = value;
    mutex.unlock();
}
void ImageIO::setStackValue(const iPoint3D voxel, const int iT, const Point3D value)
{ // alternative for vectors
    CheckVoxelIndices(voxel,"setStackValue");
    int i3d = index3dFromVoxel(voxel);
    setStackValue(i3d,iT,value);
}
void ImageIO::setStackValue(const int i3d, const int iT, const double value)
{ // main entry for scalars
    CheckVoxelIndices(i3d, "setStackValue");
    QMutex mutex;
    mutex.lock();
    _data.timePoints[iT].f1[i3d] = static_cast<float>(value);
    mutex.unlock();
}
void ImageIO::setStackValue(const int iX, const int iY, const int iZ, const int iT, const double value)
{
    CheckVoxelIndices(iX, iY, iZ, "setStackValue");
    int i3d = index3dFromVoxel(iX, iY, iZ);
    setStackValue(i3d, iT, value);
}
void ImageIO::setStackValue(const fVoxel voxel, const int iT)
{
    CheckVoxelIndices(voxel,"setStackValue");
    int i3d = index3dFromVoxel(voxel);
    setStackValue(i3d,iT,voxel.binaryValue);
}

void ImageIO::setStackValue(const iPoint3D voxel, const int iT, const double value)
{
    CheckVoxelIndices(voxel,"setStackValue");
    int i3d = index3dFromVoxel(voxel);
    setStackValue(i3d,iT,value);
}
void ImageIO::setStackValue(const iPoint4D voxel, const double value)
{
    iPoint3D voxel3d;
    voxel3d.x = voxel.x;
    voxel3d.y = voxel.y;
    voxel3d.z = voxel.z;
    setStackValue(voxel3d, voxel.t, value);
}

void ImageIO::setParallelCoordinates(dPoint4D resolution,  dPoint3D origin, iPoint3D direction)
{ // call this function through setParallelCoordinates
    _data.hdr.resolution = resolution;
    dPoint3D res;  res.x  = resolution.x;  res.y  = resolution.y;  res.z = resolution.z;
    setParallelCoordinates(res, origin, direction);
}
void ImageIO::setParallelCoordinates(dPoint3D resolution,  dPoint3D origin, iPoint3D direction)
{ // call this function through setParallelCoordinates
    _data.hdr.resolution.x = resolution.x;
    _data.hdr.resolution.y = resolution.y;
    _data.hdr.resolution.z = resolution.z;

    dMatrix matrix;   matrix.resize(4);
    // Create a 4x4 identity matrix.
    for (int j=0; j<4; j++)
    {
        matrix[j].fill(0.,4);
        matrix[j][j] = 1.;
    }
    // Attach the pixel dimensions.
    matrix[0][0] *= resolution.x * direction.x;
    matrix[1][1] *= resolution.y * direction.y;
    matrix[2][2] *= resolution.z * direction.z;
    // Attach the pixel offsets.
    matrix[0][3] = origin.x;
    matrix[1][3] = origin.y;
    matrix[2][3] = origin.z;
    _data.hdr.ijk_to_xyz.swap(matrix);
    _data.hdr.xyz_to_ijk = _data.hdr.ijk_to_xyz;
    utilMath::invertSquareMatrix( _data.hdr.xyz_to_ijk );
}

void ImageIO::setIjkToXyz(dMatrix ijkToxyz)
{
    _data.hdr.ijk_to_xyz = ijkToxyz;
    _data.hdr.xyz_to_ijk = _data.hdr.ijk_to_xyz;
    utilMath::invertSquareMatrix( _data.hdr.xyz_to_ijk );
}

#undef  ASSIF                                 /* assign v to *p, if possible */
#define ASSIF(p,v) if( (p)!=NULL ) *(p) = (v)

void ImageIO::nifti_mat44_to_quatern( dMatrix R ,
                                      float *qb, float *qc, float *qd,
                                      float *qx, float *qy, float *qz,
                                      float *dx, float *dy, float *dz, float *qfac )
{
    double r11,r12,r13 , r21,r22,r23 , r31,r32,r33 ;
    double xd,yd,zd , a,b,c,d ;
    Mat33 P,Q ;

    /* offset outputs are read write out of input matrix  */

    ASSIF(qx,R[0][3]) ; ASSIF(qy,R[1][3]) ; ASSIF(qz,R[2][3]) ;

    /* load 3x3 matrix into local variables */

    r11 = R[0][0] ; r12 = R[0][1] ; r13 = R[0][2] ;
    r21 = R[1][0] ; r22 = R[1][1] ; r23 = R[1][2] ;
    r31 = R[2][0] ; r32 = R[2][1] ; r33 = R[2][2] ;

    /* compute lengths of each column; these determine grid spacings  */

    xd = qSqrt( r11*r11 + r21*r21 + r31*r31 ) ;
    yd = qSqrt( r12*r12 + r22*r22 + r32*r32 ) ;
    zd = qSqrt( r13*r13 + r23*r23 + r33*r33 ) ;

    /* if a column length is zero, patch the trouble */

    if( xd == 0.0l ){ r11 = 1.0l ; r21 = r31 = 0.0l ; xd = 1.0l ; }
    if( yd == 0.0l ){ r22 = 1.0l ; r12 = r32 = 0.0l ; yd = 1.0l ; }
    if( zd == 0.0l ){ r33 = 1.0l ; r13 = r23 = 0.0l ; zd = 1.0l ; }

    /* assign the output lengths */

    ASSIF(dx,xd) ; ASSIF(dy,yd) ; ASSIF(dz,zd) ;

    /* normalize the columns */

    r11 /= xd ; r21 /= xd ; r31 /= xd ;
    r12 /= yd ; r22 /= yd ; r32 /= yd ;
    r13 /= zd ; r23 /= zd ; r33 /= zd ;

    /* At this point, the matrix has normal columns, but we have to allow
      for the fact that the hideous user may not have given us a matrix
      with orthogonal columns.

      So, now find the orthogonal matrix closest to the current matrix.

      One reason for using the polar decomposition to get this
      orthogonal matrix, rather than just directly orthogonalizing
      the columns, is so that inputting the inverse matrix to R
      will result in the inverse orthogonal matrix at this point.
      If we just orthogonalized the columns, this wouldn't necessarily hold. */

    Q.m[0][0] = r11 ; Q.m[0][1] = r12 ; Q.m[0][2] = r13 ; /* load Q */
    Q.m[1][0] = r21 ; Q.m[1][1] = r22 ; Q.m[1][2] = r23 ;
    Q.m[2][0] = r31 ; Q.m[2][1] = r32 ; Q.m[2][2] = r33 ;

    P = nifti_mat33_polar(Q) ;  /* P is orthog matrix closest to Q */

    r11 = P.m[0][0] ; r12 = P.m[0][1] ; r13 = P.m[0][2] ; /* unload */
    r21 = P.m[1][0] ; r22 = P.m[1][1] ; r23 = P.m[1][2] ;
    r31 = P.m[2][0] ; r32 = P.m[2][1] ; r33 = P.m[2][2] ;

    /*                            [ r11 r12 r13 ]               */
    /* at this point, the matrix  [ r21 r22 r23 ] is orthogonal */
    /*                            [ r31 r32 r33 ]               */

    /* compute the determinant to determine if it is proper */

    zd = r11*r22*r33-r11*r32*r23-r21*r12*r33
            +r21*r32*r13+r31*r12*r23-r31*r22*r13 ;  /* should be -1 or 1 */

    if( zd > 0 ){             /* proper */
        ASSIF(qfac,1.0) ;
    } else {                  /* improper ==> flip 3rd column */
        ASSIF(qfac,-1.0) ;
        r13 = -r13 ; r23 = -r23 ; r33 = -r33 ;
    }

    /* now, compute quaternion parameters */

    a = r11 + r22 + r33 + 1.0l ;

    if( a > 0.5l ){                /* simplest case */
        a = 0.5l * qSqrt(a) ;
        b = 0.25l * (r32-r23) / a ;
        c = 0.25l * (r13-r31) / a ;
        d = 0.25l * (r21-r12) / a ;
    } else {                       /* trickier case */
        xd = 1.0 + r11 - (r22+r33) ;  /* 4*b*b */
        yd = 1.0 + r22 - (r11+r33) ;  /* 4*c*c */
        zd = 1.0 + r33 - (r11+r22) ;  /* 4*d*d */
        if( xd > 1.0 ){
            b = 0.5l * qSqrt(xd) ;
            c = 0.25l* (r12+r21) / b ;
            d = 0.25l* (r13+r31) / b ;
            a = 0.25l* (r32-r23) / b ;
        } else if( yd > 1.0 ){
            c = 0.5l * qSqrt(yd) ;
            b = 0.25l* (r12+r21) / c ;
            d = 0.25l* (r23+r32) / c ;
            a = 0.25l* (r13-r31) / c ;
        } else {
            d = 0.5l * qSqrt(zd) ;
            b = 0.25l* (r13+r31) / d ;
            c = 0.25l* (r23+r32) / d ;
            a = 0.25l* (r21-r12) / d ;
        }
        if( a < 0.0l ){ b=-b ; c=-c ; d=-d; a=-a; }
    }

    ASSIF(qb,b) ; ASSIF(qc,c) ; ASSIF(qd,d) ;
    return ;
}

/*---------------------------------------------------------------------------*/
/*! polar decomposition of a 3x3 matrix

   This finds the closest orthogonal matrix to input A
   (in both the Frobenius and L2 norms).

   Algorithm is that from NJ Higham, SIAM J Sci Stat Comput, 7:1160-1174.
*//*-------------------------------------------------------------------------*/
Mat33 ImageIO::nifti_mat33_polar( Mat33 A )
{
    Mat33 X , Y , Z ;
    float alp,bet,gam,gmi , dif=1.0 ;
    int k=0 ;

    X = A ;

    /* force matrix to be nonsingular */

    gam = nifti_mat33_determ(X) ;
    while( gam == 0.0 ){        /* perturb matrix */
        gam = 0.00001 * ( 0.001 + nifti_mat33_rownorm(X) ) ;
        X.m[0][0] += gam ; X.m[1][1] += gam ; X.m[2][2] += gam ;
        gam = nifti_mat33_determ(X) ;
    }

    while(1){
        Y = nifti_mat33_inverse(X) ;
        if( dif > 0.3 ){     /* far from convergence */
            alp = qSqrt( nifti_mat33_rownorm(X) * nifti_mat33_colnorm(X) ) ;
            bet = qSqrt( nifti_mat33_rownorm(Y) * nifti_mat33_colnorm(Y) ) ;
            gam = qSqrt( bet / alp ) ;
            gmi = 1.0 / gam ;
        } else {
            gam = gmi = 1.0 ;  /* close to convergence */
        }
        Z.m[0][0] = 0.5 * ( gam*X.m[0][0] + gmi*Y.m[0][0] ) ;
        Z.m[0][1] = 0.5 * ( gam*X.m[0][1] + gmi*Y.m[1][0] ) ;
        Z.m[0][2] = 0.5 * ( gam*X.m[0][2] + gmi*Y.m[2][0] ) ;
        Z.m[1][0] = 0.5 * ( gam*X.m[1][0] + gmi*Y.m[0][1] ) ;
        Z.m[1][1] = 0.5 * ( gam*X.m[1][1] + gmi*Y.m[1][1] ) ;
        Z.m[1][2] = 0.5 * ( gam*X.m[1][2] + gmi*Y.m[2][1] ) ;
        Z.m[2][0] = 0.5 * ( gam*X.m[2][0] + gmi*Y.m[0][2] ) ;
        Z.m[2][1] = 0.5 * ( gam*X.m[2][1] + gmi*Y.m[1][2] ) ;
        Z.m[2][2] = 0.5 * ( gam*X.m[2][2] + gmi*Y.m[2][2] ) ;

        dif = qFabs(Z.m[0][0]-X.m[0][0])+qFabs(Z.m[0][1]-X.m[0][1])
                +qFabs(Z.m[0][2]-X.m[0][2])+qFabs(Z.m[1][0]-X.m[1][0])
                +qFabs(Z.m[1][1]-X.m[1][1])+qFabs(Z.m[1][2]-X.m[1][2])
                +qFabs(Z.m[2][0]-X.m[2][0])+qFabs(Z.m[2][1]-X.m[2][1])
                +qFabs(Z.m[2][2]-X.m[2][2])                          ;

        k = k+1 ;
        if( k > 100 || dif < 3.e-6 ) break ;  /* convergence or exhaustion */
        X = Z ;
    }

    return Z ;
}

/*----------------------------------------------------------------------*/
/*! compute the max column norm of a 3x3 matrix
*//*--------------------------------------------------------------------*/
float ImageIO::nifti_mat33_colnorm( Mat33 A )  /* max column norm of 3x3 matrix */
{
    float r1,r2,r3 ;

    r1 = qFabs(A.m[0][0])+qFabs(A.m[1][0])+qFabs(A.m[2][0]) ;
    r2 = qFabs(A.m[0][1])+qFabs(A.m[1][1])+qFabs(A.m[2][1]) ;
    r3 = qFabs(A.m[0][2])+qFabs(A.m[1][2])+qFabs(A.m[2][2]) ;
    if( r1 < r2 ) r1 = r2 ;
    if( r1 < r3 ) r1 = r3 ;
    return r1 ;
}

/*----------------------------------------------------------------------*/
/*! compute the max row norm of a 3x3 matrix
*//*--------------------------------------------------------------------*/
float ImageIO::nifti_mat33_rownorm( Mat33 A )  /* max row norm of 3x3 matrix */
{
    float r1,r2,r3 ;

    r1 = qFabs(A.m[0][0])+qFabs(A.m[0][1])+qFabs(A.m[0][2]) ;
    r2 = qFabs(A.m[1][0])+qFabs(A.m[1][1])+qFabs(A.m[1][2]) ;
    r3 = qFabs(A.m[2][0])+qFabs(A.m[2][1])+qFabs(A.m[2][2]) ;
    if( r1 < r2 ) r1 = r2 ;
    if( r1 < r3 ) r1 = r3 ;
    return r1 ;
}


/*----------------------------------------------------------------------*/
/*! compute the inverse of a 3x3 matrix
*//*--------------------------------------------------------------------*/
Mat33 ImageIO::nifti_mat33_inverse( Mat33 R )   /* inverse of 3x3 matrix */
{
    double r11,r12,r13,r21,r22,r23,r31,r32,r33 , deti ;
    Mat33 Q ;
    /*  INPUT MATRIX:  */
    r11 = R.m[0][0]; r12 = R.m[0][1]; r13 = R.m[0][2];  /* [ r11 r12 r13 ] */
    r21 = R.m[1][0]; r22 = R.m[1][1]; r23 = R.m[1][2];  /* [ r21 r22 r23 ] */
    r31 = R.m[2][0]; r32 = R.m[2][1]; r33 = R.m[2][2];  /* [ r31 r32 r33 ] */

    deti = r11*r22*r33-r11*r32*r23-r21*r12*r33
            +r21*r32*r13+r31*r12*r23-r31*r22*r13 ;

    if( deti != 0.0l ) deti = 1.0l / deti ;

    Q.m[0][0] = deti*( r22*r33-r32*r23) ;
    Q.m[0][1] = deti*(-r12*r33+r32*r13) ;
    Q.m[0][2] = deti*( r12*r23-r22*r13) ;

    Q.m[1][0] = deti*(-r21*r33+r31*r23) ;
    Q.m[1][1] = deti*( r11*r33-r31*r13) ;
    Q.m[1][2] = deti*(-r11*r23+r21*r13) ;

    Q.m[2][0] = deti*( r21*r32-r31*r22) ;
    Q.m[2][1] = deti*(-r11*r32+r31*r12) ;
    Q.m[2][2] = deti*( r11*r22-r21*r12) ;

    return Q ;
}

/*----------------------------------------------------------------------*/
/*! compute the determinant of a 3x3 matrix
*//*--------------------------------------------------------------------*/
float ImageIO::nifti_mat33_determ( Mat33 R )   /* determinant of 3x3 matrix */
{
    double r11,r12,r13,r21,r22,r23,r31,r32,r33 ;
    /*  INPUT MATRIX:  */
    r11 = R.m[0][0]; r12 = R.m[0][1]; r13 = R.m[0][2];  /* [ r11 r12 r13 ] */
    r21 = R.m[1][0]; r22 = R.m[1][1]; r23 = R.m[1][2];  /* [ r21 r22 r23 ] */
    r31 = R.m[2][0]; r32 = R.m[2][1]; r33 = R.m[2][2];  /* [ r31 r32 r33 ] */

    return r11*r22*r33-r11*r32*r23-r21*r12*r33
            +r21*r32*r13+r31*r12*r23-r31*r22*r13 ;
}

int ImageIO::writeNiftiFile(QString dirName, bool useCompression)
{
    FUNC_ENTER;
    resetDirName(dirName);
    if ( writeNiftiHeader(getAbsoluteName(), useCompression) ) return(1);
    if ( writeNiftiData(getAbsoluteName(), useCompression) )   return(1);
    FUNC_EXIT;
    return(0);
}

int ImageIO::writeNiftiData(QString fileName, bool useCompression)
{
    FUNC_ENTER << fileName;

    fileStorageTypes type;
    if ( useCompression )
        type = NIFTI_NII_GZ;
    else
        type = NIFTI_NII;
    QString fullFileName = utilString::attachExtension(fileName,type);

    QByteArray array = fullFileName.toLocal8Bit();
    char* fName = array.data();
    // The file should exist and already have a header, so open in append mode.
    znzFile outputFile;
    FUNC_INFO << "append NIFTI data to file" << fName;
    if ( !(outputFile = znzopen(fName,"a",useCompression)) )
    {
        qInfo() << "Error: cannot open the output file" << fullFileName;
        return(1);
    }

    imageHeader hdr = _data.hdr;
    unsigned long nAllocate = voxelsPerVolume();
    if ( _category == category_vector )
        nAllocate *= 3;
    else if ( _category == category_complex )
        nAllocate *= 2;

    if ( hdr.storageType == BSHORT )
    {
        short *idata = allocate_vector<short>(nAllocate);
        for (int jt=0; jt<nt(); jt++)
        {
            if ( _category == category_vector )
            {
                for ( unsigned long jVoxel=0; jVoxel<nAllocate/3; jVoxel++ )
                {
                    idata[3*jVoxel]   = _data.timePoints[jt].f3[jVoxel].x;
                    idata[3*jVoxel+1] = _data.timePoints[jt].f3[jVoxel].y;
                    idata[3*jVoxel+2] = _data.timePoints[jt].f3[jVoxel].z;
                }
            }
            else if ( _category == category_RGB )
            {
                for ( unsigned long jVoxel=0; jVoxel<nAllocate; jVoxel++ )
                {
                    float x = _data.timePoints[jt].f3[jVoxel].x;
                    float y = _data.timePoints[jt].f3[jVoxel].y;
                    float z = _data.timePoints[jt].f3[jVoxel].z;
                    idata[jVoxel] = qSqrt(SQR(x)+SQR(y)+SQR(z));
                }
            }
            else if ( _category == category_complex )
            {
                for ( unsigned long jVoxel=0; jVoxel<nAllocate/2; jVoxel++ )
                {
                    idata[2*jVoxel]   = _data.timePoints[jt].fc[jVoxel].real;
                    idata[2*jVoxel+1] = _data.timePoints[jt].fc[jVoxel].imag;
                }
            }
            else
            {
                for ( unsigned long jVoxel=0; jVoxel<nAllocate; jVoxel++ )
                {
                    idata[jVoxel] = getStackValue(jVoxel,jt);
                    if ( idata[jVoxel] == -32768 ) idata[jVoxel] = 32767;
                }
            }
            char *data = reinterpret_cast<char *>(idata);
            unsigned long nWrite = znzwrite(data,sizeof(short),nAllocate,outputFile);
            if ( nWrite != nAllocate )
            {
                qInfo() << "Error writing time point" << jt << "for output file" << fullFileName << "\n";
                znzclose(outputFile);
                return(1);
            } // error
        } // jt
        delete_vector(idata);
    }
    else
    { // default = floats
        float *fdata = allocate_vector<float>(nAllocate);
        for (int jt=0; jt<nt(); jt++)
        {
            if ( _category == category_vector )
            {
                for ( unsigned long jVoxel=0; jVoxel<nAllocate/3; jVoxel++ )
                {
                    fdata[3*jVoxel]   = _data.timePoints[jt].f3[jVoxel].x;
                    fdata[3*jVoxel+1] = _data.timePoints[jt].f3[jVoxel].y;
                    fdata[3*jVoxel+2] = _data.timePoints[jt].f3[jVoxel].z;
                }
            }
            else if ( _category == category_RGB )
            {
                for ( unsigned long jVoxel=0; jVoxel<nAllocate; jVoxel++ )
                {
                    float x = _data.timePoints[jt].f3[jVoxel].x;
                    float y = _data.timePoints[jt].f3[jVoxel].y;
                    float z = _data.timePoints[jt].f3[jVoxel].z;
                    fdata[jVoxel] = qSqrt(SQR(x)+SQR(y)+SQR(z));
                }
            }
            else if ( _category == category_complex )
            {
                for ( unsigned long jVoxel=0; jVoxel<nAllocate/2; jVoxel++ )
                {
                    fdata[2*jVoxel]   = _data.timePoints[jt].fc[jVoxel].real;
                    fdata[2*jVoxel+1] = _data.timePoints[jt].fc[jVoxel].imag;
                }
            }
            else
            {
                for ( unsigned long jVoxel=0; jVoxel<nAllocate; jVoxel++ )
                    fdata[jVoxel] = getStackValue(jVoxel,jt);
            }
            char *data = reinterpret_cast<char *>(fdata);
            unsigned long nWrite = znzwrite(data,sizeof(float),nAllocate,outputFile);
            if ( nWrite != nAllocate )
            {
                qInfo() << "Error writing time point" << jt << "for output file" << fullFileName << "\n";
                znzclose(outputFile);
                return(1);
            } // error
        } // jt
        delete_vector(fdata);
    }
    znzclose(outputFile);
    FUNC_EXIT;
    return(0);
}

int ImageIO::writeNiftiHeader(QString fileName, bool useCompression)
{
    FUNC_ENTER << fileName;
    fileStorageTypes type;
    if ( useCompression )
        type = NIFTI_NII_GZ;
    else
        type = NIFTI_NII;
    QString fullFileName = utilString::attachExtension(fileName,type);
    FUNC_INFO << "fullFileName" << fullFileName;

    imageHeader hdr = _data.hdr;
    dMatrix mat = _data.hdr.ijk_to_xyz;
    float qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac;
    nifti_mat44_to_quatern( mat, &qb, &qc, &qd, &qx, &qy, &qz, &dx, &dy, &dz, &qfac );

    nifti_1_header nifti_header;
    nifti_header.quatern_b = qb;
    nifti_header.quatern_c = qc;
    nifti_header.quatern_d = qd;
    nifti_header.qoffset_x = qx;
    nifti_header.qoffset_y = qy;
    nifti_header.qoffset_z = qz;

    nifti_header.sizeof_hdr = 348;
    if ( _category == category_vector )
        nifti_header.dim[0] = 5;
    else if ( hdr.dim.t > 1 )
        nifti_header.dim[0] = 4;
    else
        nifti_header.dim[0] = 3;
    nifti_header.dim[1] = hdr.dim.x;
    nifti_header.dim[2] = hdr.dim.y;
    nifti_header.dim[3] = hdr.dim.z;
    nifti_header.dim[4] = hdr.dim.t;
    if ( _category == category_vector )
        nifti_header.dim[5] = 3;
    else
        nifti_header.dim[5] = 1;
    nifti_header.dim[6] = 1;
    nifti_header.dim[7] = 1;

    // Select an output type of short integers or floats
    FUNC_INFO << "determine shorts/floats" << getFileName();
    bool useShorts = true;
    int iTime = 0;
    while ( useShorts && iTime < nt() )
    {
        for ( int jVoxel=0; jVoxel<voxelsPerVolume(); jVoxel++ )
        {
            if ( _category == category_vector || _category == category_RGB )
            {
                useShorts &= (_data.timePoints[iTime].f3[jVoxel].x == static_cast<int>(_data.timePoints[iTime].f3[jVoxel].x));
                useShorts &= (_data.timePoints[iTime].f3[jVoxel].y == static_cast<int>(_data.timePoints[iTime].f3[jVoxel].y));
                useShorts &= (_data.timePoints[iTime].f3[jVoxel].z == static_cast<int>(_data.timePoints[iTime].f3[jVoxel].z));
            }
            else if ( _category == category_complex )
            {
                useShorts &= (_data.timePoints[iTime].fc[jVoxel].real == static_cast<int>(_data.timePoints[iTime].fc[jVoxel].real));
                useShorts &= (_data.timePoints[iTime].fc[jVoxel].imag == static_cast<int>(_data.timePoints[iTime].fc[jVoxel].imag));
            }
            else
            {
                double value = _data.timePoints[iTime].f1[jVoxel];
                useShorts &= (value == static_cast<int>(value)) && (value >= -32768) && (value <= 32767);
            }
        }
        iTime++;
    }
    FUNC_INFO << "useShorts" << useShorts;
    if ( useShorts )
        _data.hdr.storageType = BSHORT;
    else
        _data.hdr.storageType = BFLOAT;
    // end selection of shorts vs floats

    if ( _data.hdr.storageType == BSHORT )
    {
        nifti_header.datatype = NIFTI_TYPE_INT16;
        nifti_header.bitpix   = 16;
    }
    else
    { // default = floating point
        nifti_header.datatype = NIFTI_TYPE_FLOAT32;
        nifti_header.bitpix   = 32;
    }
    nifti_header.slice_start = 0;
    nifti_header.slice_end = hdr.dim.z-1;
    nifti_header.slice_code  = 0;
    nifti_header.slice_duration = 0.;

    nifti_header.pixdim[0]   = qfac;    // specifies parity of coordinate system
    nifti_header.pixdim[1]   = hdr.resolution.x;
    nifti_header.pixdim[2]   = hdr.resolution.y;
    nifti_header.pixdim[3]   = hdr.resolution.z;
    nifti_header.pixdim[4]   = hdr.resolution.t;
    nifti_header.pixdim[5]   = 0.;
    nifti_header.pixdim[6]   = 0.;
    nifti_header.pixdim[7]   = 0.;

    // temporary; ignore extended header info from original file
    if ( _data.hdr.DOF != 0. )
        nifti_header.vox_offset  = 368;  // 352 + 16-byte extension
    else
        nifti_header.vox_offset  = 352;  // 352 + 16-byte extension

    nifti_header.scl_slope   = 1.;
    nifti_header.scl_inter = 0;
    nifti_header.xyzt_units  = NIFTI_UNITS_MM | NIFTI_UNITS_SEC;
//    nifti_header.xyzt_units  = NIFTI_UNITS_MM;
    nifti_header.cal_min = nifti_header.cal_max = 0.;
    nifti_header.toffset = 0.;

    // Use the affine transformation (method 3).
    nifti_header.qform_code = 1;
    nifti_header.sform_code = 0;
    // Initialize all elements to 0.
    for (int j=0; j<4; j++)
        nifti_header.srow_x[j] = nifti_header.srow_y[j] = nifti_header.srow_z[j] = 0.;
    // Diagonal elements are resolution * direction.
    iPoint3D direction = getDirection();
    nifti_header.srow_x[0] = hdr.resolution.x * direction.x;
    nifti_header.srow_y[1] = hdr.resolution.y * direction.y;
    nifti_header.srow_z[2] = hdr.resolution.z * direction.z;
    // The last row is the origin.
    dPoint3D origin = getOrigin();
    nifti_header.srow_x[3] = origin.x;
    nifti_header.srow_y[3] = origin.y;
    nifti_header.srow_z[3] = origin.z;

    // Unused fields
    strcpy(nifti_header.data_type,"");
    strcpy(nifti_header.db_name,"");
    nifti_header.extents = 0;
    nifti_header.session_error = 0;
    nifti_header.glmin = nifti_header.glmax = 0;
    nifti_header.regular = 'r';
    nifti_header.dim_info = 0;

    if ( _category == category_vector )
      nifti_header.intent_code = NIFTI_INTENT_VECTOR;
    else
        nifti_header.intent_code = 0;
    strcpy(nifti_header.intent_name,"");
    nifti_header.intent_name[0] = 0;
    nifti_header.intent_p1 = nifti_header.intent_p2 = nifti_header.intent_p3 = 0.;

    strcpy(nifti_header.descrip,"");
    strcpy(nifti_header.aux_file,"");

    strcpy(nifti_header.magic,"n+1\0");

    znzFile outputFile;
    QByteArray array = fullFileName.toLocal8Bit();
    char* standardName = array.data();

    FUNC_INFO << "open file name" << standardName;
    int error = !(outputFile = znzopen(standardName,"w",useCompression));
    if ( error )
    {
        qInfo() << "Error: cannot open the header file" << fullFileName;
        return(1);
    }

    error = znzwrite((char *)(&nifti_header),sizeof(struct nifti_1_header),1,outputFile) != 1;
    if ( error )
    {
        qInfo() << "Error writing the NIFTI header into file" << fullFileName;
        return(1);
    }

    nifti1_extender pad;
    pad.extension[0] = pad.extension[1] = pad.extension[2] = pad.extension[3] = 0;
    if ( _data.hdr.DOF != 0. )
        pad.extension[0] = 1;  // indicates extension is present
    error = znzwrite(&pad, 4, 1, outputFile) != 1;
    if ( error )
    {
        qInfo() << "Error writing nifti header file extension";
        return(1);
    }
    if ( _data.hdr.DOF != 0. )
    {
        // Create a custom extension of 16 bytes as 4 integers.
        int extension[4];
        extension[0] = 16;             // number of bytes in extension
        extension[1] = 20;             // make up my own unregistered e-code
        extension[2] = static_cast<int>(_data.hdr.DOF); // degrees of freedom saved here
        extension[3] = 0;
        error = znzwrite(extension, sizeof(int), 4, outputFile) != 4;
        if ( error )
        {
            qInfo() << "Error writing nifti header extension for DOF";
            return(1);
        }
    }
    znzclose(outputFile);
    FUNC_EXIT;
    return(0);
}

QString ImageIO::getFileHeaderName()
{
    QString ext  = getExtension();
    FUNC_ENTER << ext;
    bool nifti = !ext.compare("nii") || !ext.compare("gz");
    if ( nifti )
    {
        FUNC_EXIT << getAbsoluteName();
        return getAbsoluteName();  // hdr and file are the same
    }
    else
    { // the hdr is separate from the data
        QString dir  = getDirectory();
        QString base = getBaseName();
        QString hdrName = dir + "/" + base + ".hdr";
        return hdrName;
    }
}

QString ImageIO::getBaseName()
{   // don't simply use completeBaseName, because the name could be "smooth3.0.nii.gz"
    QString baseName = _fileInfo.baseName();  // drop off last extension (and .)
    FUNC_INFO << "baseName" << baseName;
    if ( !baseName.contains("gz") )
    { // zipped, so remove one more (e.g., ".nii.gz")
        QFileInfo stripFile(baseName);
        baseName = stripFile.baseName();
    }
    FUNC_EXIT << baseName;
    return baseName;
}

void ImageIO::setFileName(QString rootName, int fileType)
{
    FUNC_ENTER << rootName << fileType;
    _data.hdr.fileType = fileType;
    QString fileName;
    if ( fileType == BSHORT )
        fileName = rootName + ".bshort";
    else if ( fileType == BFLOAT )
        fileName = rootName + ".bfloat";
    else if (fileType == NIFTI_NII )
        fileName = rootName + ".nii";
    else if (fileType == NIFTI_NII_GZ )
        fileName = rootName + ".nii.gz";
    else
        qFatal("Programming error: only allow types bshort, bfloat, and nii for output");
    _data.hdr.fileName = fileName;
    _fileInfo.setFile(fileName);
    FUNC_EXIT << fileName << _data.hdr.fileType << _data.hdr.fileName;
}

int ImageIO::readFileHeader( QString fileName, bool verbose )
{
    FUNC_ENTER << fileName;
    _data.hdr.fileName = fileName;
    _data.hdr.dim.x        = _data.hdr.dim.y        = _data.hdr.dim.z        = _data.hdr.dim.t        = 1;
    _data.hdr.resolution.x = _data.hdr.resolution.y = _data.hdr.resolution.z = _data.hdr.resolution.t = 1.;
    _data.hdr.points.x     = _data.hdr.points.y     = _data.hdr.points.z    = _data.hdr.points.t      = 0;
    // Set the default header
    _data.hdr.fileType     = NIFTI_NII;     // nifti
    _data.hdr.dataType     = 0;             // reals
    _data.hdr.storageType  = -1;            // floats
    _data.hdr.voxelOffset  = 0;
    _data.hdr.scaleOffset  = _data.hdr.DOF = 0.;
    _data.hdr.scaleSlope = 1.;
    _data.hdr.offDiagonal = false;
    _data.hdr.byte_order = utilIO::machine_byte_order();
    // set other parameters for the file
    _swap = false;

    setFileName(fileName);
    int fileType = _data.hdr.fileType = getFileType();
    FUNC_INFO << "fileName" << fileName << "type" << fileType << getFileHeaderName();
    if ( fileType == NIFTI_NII    || fileType == NIFTI_NII_GZ || fileType == NIFTI_IMG )
    {
        if ( fileType == NIFTI_IMG && verbose )
            qInfo() << "Files with .img extensions are assumed to be NIFTI, not ANALYZE";
        // The data storage type and other information will be read from the file header.
        if ( readNiftiHeader(verbose) )  return(1);
    }
    else
    {
        // Assign the storage type based upon the file name.
        _data.hdr.storageType = _data.hdr.fileType;
        // These formats use separate headers, so the voxel offset to data is 0.
        if ( read_xdisplay_header() ) return(1);
        // The storage type should be included in UNKNOWN files (e.g. Bruker)
        if ( _data.hdr.storageType < 0 )
        {
            qInfo() << "Error: the file type is not recognized,";
            qInfo() << "so it was assumed to be an xdisplay file without an extension.";
            qInfo() << "However, the storage type was not read from file" << getAbsoluteName() << "\n";
            qInfo() << "Valid types end in .nii, .img, .bshort, .bfloat, .blong, or _'";
            return(1);
        }
        _data.hdr.fileType = _data.hdr.storageType;  // in this was changed
    }

    if ( _data.hdr.dataType == 3 )
        _category = category_vector;
    else if ( _data.hdr.dataType == 4 )
        _category = category_RGB;
    else if ( _data.hdr.dataType == 1 || _data.hdr.dataType == 2 )
        _category = category_complex;
    // # of voxels in 2, 3, and 4d; also allocate vector
    setDimensions(_data.hdr.dim);

    if ( verbose )
    {
        if ( _data.hdr.dataType == 0 )
            qInfo() << "magnitude data, ";
        else if ( _data.hdr.dataType == 1 )
            qInfo() << "complex data, ";
        else if ( _data.hdr.dataType == 2 )
            qInfo() << "magnitude-phase data, ";
        else if ( _data.hdr.dataType == 3 )
            qInfo() << "vector data, ";
        else if ( _data.hdr.dataType == 4 )
            qInfo() << "RGB data, ";
        else
            qInfo() << "data type unknown:" << _data.hdr.dataType;

        if ( _data.hdr.storageType == BSHORT )
            qInfo() << "short ints, ";
        else if ( _data.hdr.storageType == BLONG )
            qInfo() << "long ints";
        else if ( _data.hdr.storageType == BFLOAT )
            qInfo() << "floating point";
        else if ( _data.hdr.storageType == USHORT )
            qInfo() << "unsigned shorts";
        else if ( _data.hdr.storageType == UINT8 )
            qInfo() << "unsigned char";
        else if ( _data.hdr.storageType == BDOUBLE )
            qInfo() << "double (64 bit)";
        else
            qInfo() << "data type unknown:" << _data.hdr.storageType;

        qInfo() << "voxels: " << nx() << ny() << nz() << _data.hdr.dim.t;
        qInfo() << "resolution: " << _data.hdr.resolution.x << _data.hdr.resolution.y << _data.hdr.resolution.z
                << _data.hdr.resolution.t;
    }

    FUNC_EXIT;
    return(0);
}

int ImageIO::read_xdisplay_header()
{
    FUNC_ENTER;
    QFile *file = new QFile();

    file->setFileName(getFileHeaderName());
    if (!file->open(QIODevice::ReadOnly | QIODevice::Text))
        return(1);
    QTextStream in_stream(file);

    int nLine = 1;
    while (!in_stream.atEnd())
    {
        FUNC_INFO << "read line" << nLine;
        QString line = in_stream.readLine();
        QString unCommented = line.left(line.indexOf("#"));
        if ( unCommented.isEmpty() ) continue;
        int iValue = decode_xdisplay_header_keyword(line);
        if ( iValue < 0 )
            break; // the keyword "info" means to ignore everything else in the file.
        else if ( iValue > 0 )
        { // This was an error.
            qInfo() << "Error reading line" << nLine << "of file '" << getFileHeaderName() << "'\n";
            return(1);
        }
        nLine++;
    }
    FUNC_INFO << "done reading";

    // Make sure that dimensions (x,y,z,t) have been read.
    if ( nx() == 0 )
    {
        qInfo() << "Error: Value for x dimension not found in header file " << getFileHeaderName();
        return(1);
    }
    else if ( ny() == 0 )
    {
        qInfo() << "Error: Value for y dimension not found in header file " << getFileHeaderName();
        return(1);
    }
    else if ( nz() == 0 )
    {
        qInfo() << "Error: Value for z dimension not found in header file " << getFileHeaderName();
        return(1);
    }
    else if ( _data.hdr.dim.t == 0 )
    {
        qInfo() << "Error: Value for t dimension not found in header file " << getFileHeaderName();
        return(1);
    }

    FUNC_INFO << "setParallelCoordinates";

    // Define the tranformation between voxel coordinates and space.
    setParallelCoordinates(_data.hdr.resolution, _xdOrigin, _xdDirection);

    FUNC_EXIT;
    return(0);
}

int ImageIO::decode_xdisplay_header_keyword( const QString line)
{
    QString unCommented = line.mid(0, line.indexOf("#"));
    // Split names by space or comma
    QRegularExpression whiteSpaceComma("[,\\s]");// match a comma or a space
    QStringList nameList = unCommented.split(whiteSpaceComma,Qt::SkipEmptyParts); // clazy:exclude=use-static-qregularexpression

    QString keyword = nameList.at(0);
    if ( nameList.size() < 2 ) return(1); // everything needs at least 1 argument
    int error = 0;

    if ( !keyword.compare("info") || !keyword.compare("info:")
         || !keyword.compare("information") || !keyword.compare("information:") )
        return(-1);
    else if ( !keyword.compare("x") )
        _data.hdr.dim.x = nameList.at(1).toFloat();
    else if ( !keyword.compare("y") )
        _data.hdr.dim.y = nameList.at(1).toFloat();
    else if ( !keyword.compare("z") )
        _data.hdr.dim.z = nameList.at(1).toFloat();
    else if ( !keyword.compare("t") )
        _data.hdr.dim.t = nameList.at(1).toFloat();
    else if ( !keyword.compare("matrix") )
    {
        if ( nameList.size() > 4 )
        {
            _data.hdr.dim.x = nameList.at(1).toFloat();
            _data.hdr.dim.y = nameList.at(2).toFloat();
            _data.hdr.dim.z = nameList.at(3).toFloat();
            _data.hdr.dim.t = nameList.at(4).toFloat();
        }
        else if ( nameList.size() > 3 )
        {
            _data.hdr.dim.x = nameList.at(1).toFloat();
            _data.hdr.dim.y = nameList.at(2).toFloat();
            _data.hdr.dim.z = nameList.at(3).toFloat();
            _data.hdr.dim.t = 1;
        }
        else
            error = 1;
    }
    else if ( !keyword.compare("byte-order") )
        _data.hdr.byte_order = nameList.at(1).toInt();
    else if ( !keyword.compare("resolution") )
    {
        if ( nameList.size() > 4 )
        {
            _data.hdr.resolution.x = nameList.at(1).toFloat();
            _data.hdr.resolution.y = nameList.at(2).toFloat();
            _data.hdr.resolution.z = nameList.at(3).toFloat();
            _data.hdr.resolution.t = nameList.at(4).toFloat();
        }
        else if ( nameList.size() > 3 )
        {
            _data.hdr.resolution.x = nameList.at(1).toFloat();
            _data.hdr.resolution.y = nameList.at(2).toFloat();
            _data.hdr.resolution.z = nameList.at(3).toFloat();
            _data.hdr.resolution.t = 1;
        }
        else
            error = 1;
    }
    else if ( !keyword.compare("origin") )
    {
        if ( nameList.size() < 4 )
            error = 1;
        else
        {
            _xdOrigin.x = nameList.at(1).toFloat();
            _xdOrigin.y = nameList.at(2).toFloat();
            _xdOrigin.z = nameList.at(3).toFloat();
        }
    }
    else if ( !keyword.compare("direction") )
    {
        if ( nameList.size() < 4 )
            error = 1;
        else
        {
            _xdDirection.x = nameList.at(1).toFloat();
            _xdDirection.y = nameList.at(2).toFloat();
            _xdDirection.z = nameList.at(3).toFloat();
            if ( qAbs(_xdDirection.x) != 1. ||
                 qAbs(_xdDirection.y) != 1. ||
                 qAbs(_xdDirection.z) != 1. )
            {
                qInfo() << "Error: allowed values for direction (%g,%g,%g) are 1 and -1.";
                error = 1;
            }
        }
    }
    else if ( !keyword.compare("data-slope") || !keyword.compare("data-scale") )
    {
        if ( nameList.size() < 2 )
            error = 1;
        else
            _data.hdr.scaleSlope = nameList.at(1).toFloat();
    }
    else if ( !keyword.compare("data-type") )
    {
        keyword = nameList.at(1);
        _data.hdr.dataType = 0;   // default is magnitude only
        if ( !keyword.compare("magnitude") )
            _data.hdr.dataType = 0; // magnitude
        else if ( !keyword.compare("complex") )
            _data.hdr.dataType = 1; // real-imaginary
        else if ( !keyword.compare("real-imaginary") )
            _data.hdr.dataType = 1; // real-imaginary
        else if ( !keyword.compare("magnitude-phase") )
            _data.hdr.dataType = 2; // magnitude-phase
        else if ( !keyword.compare("vector") )
            _data.hdr.dataType = 3; // magnitude-phase
        else
        {
            qInfo() << "Error: arguments for data-type are";
            qInfo() << "1) magnitude (this is the default),";
            qInfo() << "2) complex (or real-imaginary)";
            qInfo() << "3) magnitude-phase";
            error = 1;
        }
    }
    else if ( !keyword.compare("storage-type") )
    {
        keyword = nameList.at(1);
        if ( !keyword.compare("short") )
            _data.hdr.storageType = BSHORT; // short
        else if ( !keyword.compare("float") )
            _data.hdr.storageType = BFLOAT; // float
        else if ( !keyword.compare("long") )
            _data.hdr.storageType = BLONG; // long
        else if ( !keyword.compare("unsigned-short") )
            _data.hdr.storageType = USHORT; // unsigned-short
        else if ( !keyword.compare("unsigned-char") )
            _data.hdr.storageType = UINT8; // unsigned-short
        else if ( !keyword.compare("double") )
            _data.hdr.storageType = BDOUBLE; // unsigned-short
        else
        {
            qInfo() << "Error: arguments for storage-type are ...";
            qInfo() << "    short, float, long, unsigned-short, unsigned-char";
            error = 1;
        }
    }
    else if ( !keyword.compare("degrees-of-freedom") )
        _data.hdr.DOF = nameList.at(1).toDouble();
    else // not found
        error = 1;

    return(error);
}

int ImageIO::readNiftiHeader( bool verbose )
{
    QString hdrName = getFileHeaderName();
    FUNC_ENTER << hdrName;

    QFile testFile(hdrName);
    if (!testFile.open(QIODevice::ReadOnly | QIODevice::Text))
    {  //  the default didn't work, so check the other NIFTI type (.nii , .nii.gz)
        QString extensionTest;  int fileType;
        if ( _data.hdr.fileType == NIFTI_NII )
        {
            extensionTest = ".nii.gz";
            fileType = NIFTI_NII_GZ;
        }
        else
        {
            extensionTest = ".nii";
            fileType = NIFTI_NII;
        }
        QString baseName = utilString::getFileNameWithoutExtension(getAbsoluteName());
        hdrName  = baseName + extensionTest;
        FUNC_INFO << "try" << baseName << extensionTest << hdrName;
        testFile.setFileName(hdrName);
        if (!testFile.open(QIODevice::ReadOnly | QIODevice::Text))
        {  // neither worked
            FUNC_INFO << "could not open";
            qInfo() << "Error: cannot open Nifti file " << hdrName;
            return(1);
        }
        else
        {
            FUNC_INFO << "could open";
            setFileName(baseName,fileType);
        }
        testFile.close();
    } // testFile
    FUNC_INFO << "now read" << getFileHeaderName();

    znzFile file;
    QByteArray array = getFileHeaderName().toLocal8Bit();
    char* fileName = array.data();
    int useCompression = _data.hdr.fileType == NIFTI_NII_GZ;
    if ( !(file = znzopen(fileName,"r",useCompression)) )
    {
        qInfo() << "Error: cannot open Nifti file " << getFileHeaderName();
        return(1);
    }

    nifti_1_header nifti_header;
    nifti1_extender pad;
    int nread = znzread(&nifti_header, sizeof(nifti_header), 1, file);
    if (nread != 1)
    {
        qWarning() << "Error reading header file" << getFileHeaderName();
        exit(1);
    }
    // Read the pad at the end of the header.
    nread = znzread(&pad, 4, 1, file);
    if ( (nread!=0) && (pad.extension[0]!=0) )
    {
        int esize;
        if ( znzread(&esize, 4, 1, file) != 1 )
        {
            qInfo() << "Error: A NIFTI1 extension was indicated, but it is missing.";
            znzclose(file);
            return(1);
        }
        int ecode;
        if ( znzread(&ecode, 4, 1, file) != 1 )
        {
            qInfo() << "Error: A NIFTI1 extension was indicated, but the ecode is missing.";
            znzclose(file);
            return(1);
        }
        if ( ecode == 20 ) // my unregistered code
        {
            int lDOF;
            if ( znzread(&lDOF, 4, 1, file) != 1 )
            {
                qInfo() << "Error read degrees of freedom.";
                znzclose(file);
                return(1);
            }
            _data.hdr.DOF = static_cast<double>(lDOF);
        }
    }
    znzclose(file);

    // Examine dim[0] to see if it is in the range 1->7
    bool swap_bytes = (nifti_header.dim[0]<1 || nifti_header.dim[0]>7);
    if (swap_bytes )
        SwapNiftiHeader(&nifti_header);
    // Assign values to the x,y,z,t dimensions.
    _data.hdr.dim.x    = nifti_header.dim[1];
    _data.hdr.dim.y    = nifti_header.dim[2];
    _data.hdr.dim.z    = nifti_header.dim[3];
    _data.hdr.dim.t    = nifti_header.dim[4];
    if ( nifti_header.dim[0] < 4 )
        _data.hdr.dim.t = 1;
    else if ( nifti_header.dim[0] > 5 )
    {
        qInfo() << "Warning: fastmap does not currently handle data of dimensionality > 5 (" << nifti_header.dim[0] << ")";
        qInfo() << "Proceeding using the 1st 5 dimensions...";
    }
    // Allow 3-vectors at each voxel
    bool vectorType = nifti_header.dim[0] == 5 && nifti_header.dim[5] == 3;
    bool rgbType = nifti_header.datatype == NIFTI_TYPE_RGB24;
    FUNC_INFO << "test types" << vectorType << rgbType;
    if ( vectorType )
        _data.hdr.dataType = 3;   // vector
    else if ( rgbType )
        _data.hdr.dataType = 4;   // vector
    else
        _data.hdr.dataType = 0;   // default is magnitude only
    FUNC_INFO << "dataType1" << _data.hdr.dataType;

    _data.hdr.resolution.x = nifti_header.pixdim[1];
    _data.hdr.resolution.y = nifti_header.pixdim[2];
    _data.hdr.resolution.z = nifti_header.pixdim[3];
    _data.hdr.resolution.t = nifti_header.pixdim[4];
    if ( _data.hdr.resolution.t == 0. ) _data.hdr.resolution.t = 1.;  // don't allow a time step of 0
    // The data can contain a slope and offset.
    _data.hdr.scaleSlope  = nifti_header.scl_slope;
    _data.hdr.scaleOffset = nifti_header.scl_inter;
    //  header->qfac         = nifti_header.pixdim[0];
    //  header->quaternion.x = nifti_header.quatern_b;
    //  header->quaternion.y = nifti_header.quatern_c;
    //  header->quaternion.z = nifti_header.quatern_d;
    // Define the data storage type.
    // Only a few data types are supported.

    if ( _data.hdr.dataType == 4 )
        _data.hdr.storageType = UINT8;
    else if ( nifti_header.datatype == NIFTI_TYPE_INT16 )
        _data.hdr.storageType = BSHORT;
    else if ( nifti_header.datatype == NIFTI_TYPE_UINT16 )
        _data.hdr.storageType = USHORT;
    else if ( nifti_header.datatype == NIFTI_TYPE_INT32 )
        _data.hdr.storageType = BLONG;
    else if ( nifti_header.datatype == NIFTI_TYPE_FLOAT32 )
        _data.hdr.storageType = BFLOAT;
    else if ( nifti_header.datatype == NIFTI_TYPE_UINT8 )
        _data.hdr.storageType = UINT8;
    else if ( nifti_header.datatype == NIFTI_TYPE_FLOAT64 )
        _data.hdr.storageType= BDOUBLE;
    else if ( nifti_header.datatype == NIFTI_TYPE_COMPLEX64 )
    {
        _data.hdr.storageType = BFLOAT;
        _data.hdr.dataType    = 1;    // complex
    }
    else
    {
        qInfo() << "Error: the Nifti data storage type" <<  nifti_header.datatype << "is not supported";
        return(1);
    }
    FUNC_INFO << "dataType2" << _data.hdr.dataType;

    // define the voxel offset
    _data.hdr.voxelOffset = nifti_header.vox_offset;

    // Make sure that dimensions (x,y,z,t) have been read and make sense.
    if ( nx() <= 0 )
    {
        qInfo() << "Error: Value for x dimension not found in header file " << getFileHeaderName();
        return(1);
    }
    else if ( ny() <= 0 )
    {
        qInfo() << "Error: Value for y dimension not found in header file " << getFileHeaderName();
        return(1);
    }
    else if ( nz() <= 0 )
    {
        qInfo() << "Error: Value for z dimension not found in header file " << getFileHeaderName();
        return(1);
    }
    else if ( _data.hdr.dim.t <= 0 )
    {
        qInfo() << "Error: Value for t dimension not found in header file " << getFileHeaderName();
        return(1);
    }

    // Define the tranformation between voxel coordinates and space.
    if ( nifti_header.qform_code == 0 ) // "old" Analyze 7.5 method
    { // Setting the origin to zero seems wrong, but this follows the nifti1.h documentation.
        if ( verbose )
            qInfo() << "The NIFTI header for file" << fileName << "uses Method 1: the old analyze format";
        // origin
        dPoint3D origin;  origin.x = origin.y = origin.z = 0.;
        // direction: align coordinates with indices
        iPoint3D direction; direction.x = direction.y = direction.z = 1;
        setParallelCoordinates(_data.hdr.resolution,origin,direction);
    }
    else // qform_code > 0
    {
        if ( nifti_header.sform_code == 0 )
        {
            if ( verbose )
                qInfo() << "The NIFTI header for file" << fileName << "uses Method 2: quaternion rotations";
            _data.hdr.ijk_to_xyz = nifti_quatern_to_mat44(nifti_header.quatern_b,   // qb
                                                          nifti_header.quatern_c,   // qc
                                                          nifti_header.quatern_d,   // qd
                                                          nifti_header.qoffset_x,   // offset_x
                                                          nifti_header.qoffset_y,   // offset_y
                                                          nifti_header.qoffset_z,   // offset_z
                                                          nifti_header.pixdim[1],   // resolution_x
                    nifti_header.pixdim[2],   // resolution_y
                    nifti_header.pixdim[3],   // resolution_z
                    nifti_header.pixdim[0]);  // qfac
            _data.hdr.xyz_to_ijk = _data.hdr.ijk_to_xyz;
            utilMath::invertSquareMatrix( _data.hdr.xyz_to_ijk );
        }
        else
        {
            if ( verbose )
                qInfo() << "The NIFTI header for file" << fileName << "uses Method 3: affine transformation of voxel indices";
            _data.hdr.ijk_to_xyz = general_affine_to_mat44( nifti_header.srow_x,nifti_header.srow_y,nifti_header.srow_z);
            _data.hdr.xyz_to_ijk = _data.hdr.ijk_to_xyz;
            utilMath::invertSquareMatrix( _data.hdr.xyz_to_ijk );
        }
    }
    // Issue a warning for non-zero off-diagnonal elements.
    _data.hdr.offDiagonal = false;
    for (int j1=0; j1<3; j1++)
    {
        for (int j2=0; j2<3; j2++)
        {
            if ( j1 != j2 && _data.hdr.ijk_to_xyz[j1][j2] != 0. )
                _data.hdr.offDiagonal = true;
        }
    }
    if ( _data.hdr.offDiagonal && verbose )
    {
        qInfo() << "Warning: voxel-to-coordinate matrix has non-zero off-diagonal elements.";
        for (int j1=0; j1<4; j1++)
        {
            for (int j2=0; j2<4; j2++)
                qInfo() << "#1 ijk_to_xyz[" << j1 << "][" << j2 << "] = " <<  _data.hdr.ijk_to_xyz[j1][j2];
        }
    }
    return(0);
}

/*! Swap siz bytes at a time from the given array of n sets of bytes
 *
 *  Declared void * so that the fields from the headers can be passed through
 *  without casting.
 */
void ImageIO::SwapBytes(size_t n, int size, void *bytes)
{
    size_t i;
    char *cp0 = (char *)bytes, *cp1, *cp2;
    char swap;

    for(i=0; i < n; i++) {
        cp1 = cp0;
        cp2 = cp0 + (size-1);
        while (cp2 > cp1)
        {
            swap = *cp1; *cp1 = *cp2; *cp2 = swap;
            cp1++; cp2--;
        }
        cp0 += size;
    }
}

/*
 *  Byte swap the individual fields of a NIFTI-1 header struct.
 */
void ImageIO::SwapNiftiHeader(struct nifti_1_header *h)
{
    SwapBytes(1, 4, &h->sizeof_hdr);
    SwapBytes(1, 4, &h->extents);
    SwapBytes(1, 2, &h->session_error);

    SwapBytes(8, 2, h->dim);
    SwapBytes(1, 4, &h->intent_p1);
    SwapBytes(1, 4, &h->intent_p2);
    SwapBytes(1, 4, &h->intent_p3);

    SwapBytes(1, 2, &h->intent_code);
    SwapBytes(1, 2, &h->datatype);
    SwapBytes(1, 2, &h->bitpix);
    SwapBytes(1, 2, &h->slice_start);

    SwapBytes(8, 4, h->pixdim);

    SwapBytes(1, 4, &h->vox_offset);
    SwapBytes(1, 4, &h->scl_slope);
    SwapBytes(1, 4, &h->scl_inter);
    SwapBytes(1, 2, &h->slice_end);

    SwapBytes(1, 4, &h->cal_max);
    SwapBytes(1, 4, &h->cal_min);
    SwapBytes(1, 4, &h->slice_duration);
    SwapBytes(1, 4, &h->toffset);
    SwapBytes(1, 4, &h->glmax);
    SwapBytes(1, 4, &h->glmin);

    SwapBytes(1, 2, &h->qform_code);
    SwapBytes(1, 2, &h->sform_code);

    SwapBytes(1, 4, &h->quatern_b);
    SwapBytes(1, 4, &h->quatern_c);
    SwapBytes(1, 4, &h->quatern_d);
    SwapBytes(1, 4, &h->qoffset_x);
    SwapBytes(1, 4, &h->qoffset_y);
    SwapBytes(1, 4, &h->qoffset_z);

    SwapBytes(4, 4, h->srow_x);
    SwapBytes(4, 4, h->srow_y);
    SwapBytes(4, 4, h->srow_z);

    return ;
}

double ImageIO::getStackValue(int i3d, int iTime)
{
    double value;
    QMutex mutex;
    mutex.lock();
    if ( _category == category_vector || _category == category_RGB )
    {
        Point3D vecData = _data.timePoints[iTime].f3[i3d];
        if ( _vectorScalar == vectorScalar_Magnitude )
            value = qSqrt(SQR(static_cast<double>(vecData.x)) + SQR(static_cast<double>(vecData.y)) + SQR(static_cast<double>(vecData.z)));
        else if ( _vectorScalar == vectorScalar_x )
            value = static_cast<double>(vecData.x);
        else if ( _vectorScalar == vectorScalar_y )
            value = static_cast<double>(vecData.y);
        else // if ( _vectorScalar == vectorScalar_z )
            value = static_cast<double>(vecData.z);
    }
    else
        value = static_cast<double>(_data.timePoints[iTime].f1[i3d]);
    mutex.unlock();
    return value;
}

double ImageIO::getStackValue(int i2d, int iZ, int it)
{
    int i3d = iZ * _data.hdr.points.y + i2d;
    return getStackValue(i3d, it);
}
double ImageIO::getStackValue(int iX, int iY, int iZ, int it)
{
    int i3d = index3dFromVoxel(iX, iY, iZ);
    return getStackValue(i3d, it);
}
double ImageIO::getStackValue(iPoint3D voxel, int it)
{
    int i3d = index3dFromVoxel(voxel);
    return getStackValue(i3d, it);
}
double ImageIO::getStackValue(fVoxel voxel, int it)
{
    int i3d = index3dFromVoxel(voxel);
    return getStackValue(i3d, it);
}
double ImageIO::getStackValue(fVoxelLight voxelLight, int it)
{
    int i3d = index3dFromVoxel(voxelLight);
    return getStackValue(i3d, it);
}
double ImageIO::getStackValue(iPoint4D voxel4D)
{
    int i3d = voxel4D.z * _data.hdr.points.y + voxel4D.y * _data.hdr.points.x + voxel4D.x;
    return getStackValue(i3d, voxel4D.t);
}
Point3D ImageIO::getStackVectorValue(fVoxel voxel, int iTime)
{
    iPoint3D voxel3d; voxel3d.x = voxel.x; voxel3d.y = voxel.y; voxel3d.z = voxel.z;
    return getStackVectorValue(voxel3d, iTime);
}

Point3D ImageIO::getStackVectorValue(iPoint3D voxel, int iTime)
{
    if ( _category != category_RGB )
        qFatal("Programming error in getStackVectorValue: wrong data type for this file.");
    QMutex mutex;
    mutex.lock();
    if ( iTime >= _data.timePoints.size() ) qFatal("_iTime out of range in getVolumeValue");
    int i3d = voxel.z * _data.hdr.points.y + voxel.y * _data.hdr.points.x + voxel.x;
    Point3D value = _data.timePoints[iTime].f3[i3d];
    mutex.unlock();
    return value;
}
Point3D ImageIO::getVolumeVectorValue(iPoint3D voxel)
{
    int i3d = index3dFromVoxel(voxel);
    return getVolumeVectorValue(i3d);
}
Point3D ImageIO::getVolumeVectorValue(int i3d)
{
    if ( _category != category_vector )
        qFatal("Programming error in getVolumeVectorValue: wrong data type for this file.");
    QMutex mutex;
    mutex.lock();
    if ( _iTime >= _data.timePoints.size() ) qFatal("_iTime out of range in getVolumeValue");
    Point3D value = _data.timePoints[_iTime].f3[i3d];
    mutex.unlock();
    return value;
}
double ImageIO::getVolumeValue(int i3d)
{
//    FUNC_ENTER << "category" << _category;
    float value;
    if ( _category == category_vector || _category == category_RGB )
    {
        QMutex mutex;
        mutex.lock();
        if ( _iTime >= _data.timePoints.size() ) qFatal("_iTime out of range in getVolumeValue");
        Point3D vecData = _data.timePoints[_iTime].f3[i3d];
        mutex.unlock();
        if ( _vectorScalar == vectorScalar_Magnitude )
            value = static_cast<float>(qSqrt(SQR(static_cast<double>(vecData.x)) + SQR(static_cast<double>(vecData.y)) + SQR(static_cast<double>(vecData.z))));
        else if ( _vectorScalar == vectorScalar_x )
            value = vecData.x;
        else if ( _vectorScalar == vectorScalar_y )
            value = vecData.y;
        else // vectorScalar_z
            value = vecData.z;
    }
    else
    {
        QMutex mutex;
        mutex.lock();
        if ( _iTime >= _data.timePoints.size() ) qFatal("_iTime out of range in getVolumeValue");
        if ( i3d >= _data.timePoints[_iTime].f1.size() )
        {
            FUNC_INFO << "i3d" << i3d << _data.timePoints[_iTime].f1.size();
            FUNC_INFO << nx() << ny() << nz() << nt();
            exit(1);
        }
        value = _data.timePoints[_iTime].f1[i3d];
        mutex.unlock();
    }
//    FUNC_EXIT << value;
    return static_cast<double>(value);
}
double ImageIO::getVolumeValue(int iX, int iY, int iZ)
{
    int i3d = index3dFromVoxel(iX, iY, iZ);
    return getVolumeValue(i3d);
}
double ImageIO::getVolumeValue(iPoint3D voxel)
{
    int i3d = index3dFromVoxel(voxel);
    return getVolumeValue(i3d);
}
double ImageIO::getVolumeValue(fVoxel voxel)
{
    int i3d = index3dFromVoxel(voxel);
    return getVolumeValue(i3d);
}

bool ImageIO::voxelInBounds(iPoint3D voxel)
{
    bool inBounds = voxel.x >= 0 && voxel.x < nx();
    inBounds     &= voxel.y >= 0 && voxel.y < ny();
    inBounds     &= voxel.z >= 0 && voxel.z < nz();
    return inBounds;
}
bool ImageIO::voxelInBounds(iPoint4D voxel)
{
    bool inBounds = voxel.x >= 0 && voxel.x < nx();
    inBounds     &= voxel.y >= 0 && voxel.y < ny();
    inBounds     &= voxel.z >= 0 && voxel.z < nz();
    inBounds     &= voxel.t >= 0 && voxel.t < nt();
    return inBounds;
}
bool ImageIO::voxelInBounds(fVoxel voxel)
{
    bool inBounds = voxel.x >= 0 && voxel.x < nx();
    inBounds     &= voxel.y >= 0 && voxel.y < ny();
    inBounds     &= voxel.z >= 0 && voxel.z < nz();
    return inBounds;
}
bool ImageIO::voxelInBounds(fVoxelLight voxel)
{
    bool inBounds = voxel.x >= 0 && voxel.x < nx();
    inBounds     &= voxel.y >= 0 && voxel.y < ny();
    inBounds     &= voxel.z >= 0 && voxel.z < nz();
    return inBounds;
}

void ImageIO::CheckVoxelIndices( fVoxel voxel, const QString callingFunction )
{
    if ( !voxelInBounds(voxel) )
    {
        qWarning() << "file name" << getAbsoluteName();
        qWarning() << "voxel indices " << voxel.x << voxel.y << voxel.z << " out of bounds " <<
                      nx() << ny() << nz();
        qWarning() << "calling function: " << callingFunction;
        exit(1);
    }
}
void ImageIO::CheckVoxelIndices( iPoint3D voxel, const QString callingFunction )
{
    if ( !voxelInBounds(voxel) )
    {
        qWarning() << "calling function: " << callingFunction;
        qWarning() << "file name" << getAbsoluteName();
        qWarning() << "voxel indices " << voxel.x << voxel.y << voxel.z << " out of bounds "
                   << nx() << ny() << nz();
        exit(1);
    }
}
void ImageIO::CheckVoxelIndices( int iX, int iY, int iZ, const QString callingFunction )
{
    if ( iX < 0 || iX >= nx() ||
         iY < 0 || iY >= ny() ||
         iZ < 0 || iZ >= nz() )
    {
        qWarning() << "voxel out of bounds; file name =" << getAbsoluteName();
        qWarning() << "voxel indices " << iX << iY << iZ << " out of bounds " << nx() << ny() << nz();
        qWarning() << "calling function: " << callingFunction;
        exit(1);
    }
}
void ImageIO::CheckVoxelIndices( int i3d, const QString callingFunction )
{
    if ( i3d < 0 || i3d >= _data.hdr.points.z )
    {
        qInfo() << "file name" << getAbsoluteName();
        qInfo() << "voxel indices " << i3d << " out of bounds " << _data.hdr.points.z;
        qInfo() << "calling function: " << callingFunction;
        exit(1);
    }
}

dPoint3D ImageIO::coordinateFromVoxel(iPoint3D voxel)
{
    double DataIndex[4];
    dPoint3D Coord={0.,0.,0.};
    CheckVoxelIndices(voxel,"coordinateFromVoxel (iPoint3D)");

    DataIndex[0] = voxel.x;
    DataIndex[1] = voxel.y;
    DataIndex[2] = voxel.z;
    DataIndex[3] = 1.;
    /*
    for (int j1=0; j1<4; j1++)
    {
        for (int j2=0; j2<4; j2++)
            FUNC_INFO << "ijk_to_xyz[" << j1 << "][" << j2 << "] = " <<  _data.hdr.ijk_to_xyz[j1][j2];
    }
    */
    for (int j=0; j<4; j++)
        Coord.x += _data.hdr.ijk_to_xyz[0][j] * DataIndex[j];
    for (int j=0; j<4; j++)
        Coord.y += _data.hdr.ijk_to_xyz[1][j] * DataIndex[j];
    for (int j=0; j<4; j++)
        Coord.z += _data.hdr.ijk_to_xyz[2][j] * DataIndex[j];
    return Coord;
}
dPoint3D ImageIO::coordinateFromVoxel(int iX, int iY, int iZ)
{
    double DataIndex[4];
    dPoint3D Coord={0.,0.,0.};
    CheckVoxelIndices(iX,iY,iZ,"coordinateFromVoxel (3 integers)");

    DataIndex[0] = iX;
    DataIndex[1] = iY;
    DataIndex[2] = iZ;
    DataIndex[3] = 1.;

    for (int j=0; j<4; j++)
        Coord.x += _data.hdr.ijk_to_xyz[0][j] * DataIndex[j];
    for (int j=0; j<4; j++)
        Coord.y += _data.hdr.ijk_to_xyz[1][j] * DataIndex[j];
    for (int j=0; j<4; j++)
        Coord.z += _data.hdr.ijk_to_xyz[2][j] * DataIndex[j];
    return Coord;
}

dMatrix ImageIO::nifti_quatern_to_mat44( float qb, float qc, float qd,
                                         float qx, float qy, float qz,
                                         float dx, float dy, float dz, float qfac )
{
    dMatrix R;
    R.resize(4);
    for (int j=0; j<4; j++)
        R[j].resize(4);

    qreal a,b=qb,c=qc,d=qd , xd,yd,zd ;

    /* last row is always [ 0 0 0 1 ] */

    R[3][0]=R[3][1]=R[3][2] = 0.0 ; R[3][3]= 1.0 ;

    /* compute a parameter from b,c,d */

    a = 1.0l - (b*b + c*c + d*d) ;
    if( a < 1.e-7l )
    {                   /* special case */
        a = 1.0l / qSqrt(b*b+c*c+d*d) ;
        b *= a ; c *= a ; d *= a ;        /* normalize (b,c,d) vector */
        a = 0.0l ;                        /* a = 0 ==> 180 degree rotation */
    } else{
        a = qSqrt(a) ;                     /* angle = 2*arccos(a) */
    }

    /* load rotation matrix, including scaling factors for voxel sizes */

    xd = (dx > 0.0) ? dx : 1.0l ;       /* make sure are positive */
    yd = (dy > 0.0) ? dy : 1.0l ;
    zd = (dz > 0.0) ? dz : 1.0l ;

    if( qfac < 0.0 ) zd = -zd ;         /* left handedness? */

    R[0][0] =        (a*a+b*b-c*c-d*d) * xd ;
    R[0][1] = 2.0l * (b*c-a*d        ) * yd ;
    R[0][2] = 2.0l * (b*d+a*c        ) * zd ;
    R[1][0] = 2.0l * (b*c+a*d        ) * xd ;
    R[1][1] =        (a*a+c*c-b*b-d*d) * yd ;
    R[1][2] = 2.0l * (c*d-a*b        ) * zd ;
    R[2][0] = 2.0l * (b*d-a*c        ) * xd ;
    R[2][1] = 2.0l * (c*d+a*b        ) * yd ;
    R[2][2] =        (a*a+d*d-c*c-b*b) * zd ;

    /* load offsets */

    R[0][3] = qx ; R[1][3] = qy ; R[2][3] = qz ;

    return R ;
}

dMatrix ImageIO::general_affine_to_mat44( float *srow_x, float *srow_y, float *srow_z )
{
    dMatrix R;
    R.resize(4);
    for (int j=0; j<4; j++)
        R[j].resize(4);

    R[0][0] = srow_x[0];
    R[0][1] = srow_x[1];
    R[0][2] = srow_x[2];
    R[0][3] = srow_x[3];

    R[1][0] = srow_y[0];
    R[1][1] = srow_y[1];
    R[1][2] = srow_y[2];
    R[1][3] = srow_y[3];

    R[2][0] = srow_z[0];
    R[2][1] = srow_z[1];
    R[2][2] = srow_z[2];
    R[2][3] = srow_z[3];

    // last row is always [ 0 0 0 1 ]
    R[3][0]=R[3][1]=R[3][2]=0.0;  R[3][3]=1.0;

    bool AllZero[3];
    for (int j2=0; j2<3; j2++)
    {
        AllZero[j2] = true;
        for (int j1=0; j1<3; j1++)
            AllZero[j2] &= R[j2][j1] == 0.;
    }
    if ( AllZero[0] || AllZero[1] || AllZero[2] )
    {
        qInfo() << "Invalid Affine transformation matrix: resetting to identity matrix.\n";
        for (int j2=0; j2<3; j2++)
        {
            for (int j1=0; j1<3; j1++)
            {
                if ( j1 == j2 )
                    R[j1][j2] = 1.;
                else
                    R[j1][j2] = 0.;
            }
        }
    }
    /*
    for (int j2=0; j2<4; j2++)
    {
        for (int j1=0; j1<4; j1++)
            qDebug() << "R[" << j1 << "][" << j2 << "] =" << R[j1][j2];
    }
    */
    return R ;
}
int ImageIO::readImageData()
{
    QProgressBar *progressBar=nullptr;
    return readImageData(progressBar);
}

int ImageIO::readImageData(QProgressBar *progressBar)
// This function is replaced by voxelReader::readImageData(), which is suited to multiple threads
{
    double *f64Data;
    float *fData;
    short *iData;
    int *iiData;
    unsigned short *iuData;
    unsigned char *i8Data;

    int swap_bytes  = utilIO::machine_byte_order() != _data.hdr.byte_order;
    int output_info = _data.hdr.points.t > 1000000 && _data.hdr.dim.t > 10;
    if ( progressBar && output_info )
    {
        qInfo() << "set progressbar min and max" << (_data.hdr.dim.t-1)*(nz()-1);
        progressBar->setMinimum(0);
        progressBar->setMaximum((_data.hdr.dim.t-1)*(nz()-1));
    }

    FUNC_INFO << "set absolute?" << getAbsoluteName();
    // set the absolute name
    if ( _fileInfo.exists() && _fileInfo.isFile() )
        FUNC_INFO << "absolute file name" << getAbsoluteName();

    FUNC_INFO << "**** read file" << getAbsoluteName();
    // Open file file.
    znzFile file;
    QByteArray array = getAbsoluteName().toLocal8Bit();
    char* fileName = array.data();
    int useCompression = _data.hdr.fileType == NIFTI_NII_GZ;
    FUNC_INFO << "useCompression" << useCompression << "on file" << fileName;
    if ( !(file = znzopen(fileName,"r",useCompression)) )
    {
        qInfo() << "Error: cannot open file" << getFileHeaderName();
        return(1);
    }
    // Go to the first data voxel
    znzseek(file, (long)(_data.hdr.voxelOffset), SEEK_SET);

    unsigned long nAllocate = _data.hdr.points.z;
    if ( _data.hdr.dataType == 1 || _data.hdr.dataType == 2 )
        nAllocate *= 2;  // complex
    else if ( _data.hdr.dataType == 3 || _data.hdr.dataType == 4 )
        nAllocate *= 3;  // 3-vector

    fData = allocate_vector<float>(nAllocate);
    if ( _data.hdr.storageType == BSHORT ) // bshort
        iData = allocate_vector<short>(nAllocate);
    else if ( _data.hdr.storageType == BLONG )  // blong
        iiData = allocate_vector<int>(nAllocate);
    else if ( _data.hdr.storageType == USHORT ) // unsigned short
        iuData = allocate_vector<unsigned short>(nAllocate);
    else if ( _data.hdr.storageType == UINT8 ) // bytes
        i8Data = allocate_vector<unsigned char>(nAllocate);
    else if ( _data.hdr.storageType == BDOUBLE ) // bytes
        f64Data = allocate_vector<double>(nAllocate);

    // 4D data structures should already be allocate via readFileHeader
//    _data.timePoints.resize(_data.hdr.dim.t);
//    _data.color.fill(Color_Undefined,_data.hdr.points.z);  // fill 3d color with undefined

    int error;
    // Read the array.
    for (int jt=0; jt<_data.hdr.dim.t; jt++)
    {
        if ( progressBar && output_info )
        {
            qInfo() << "set progress to " << jt;
            progressBar->setValue(jt);
        }
        //      if ( output_info && jt%(_data.hdr.dim.t/10) == 0 )
        //      ipercent = MIN(100,nearbyint((float)(jt*100)/(float)_data.hdr.dim.t));

        // Read the array, swap bytes if necessary, and place results into fData
        if ( _data.hdr.storageType == BFLOAT )      // bfloat
        {
            FUNC_INFO << "read float";
            error = znzread((char *)fData,sizeof(float),nAllocate,file) != nAllocate;
            if ( swap_bytes )
            {
                for (unsigned long j=0; j<nAllocate; j++)
                    utilIO::swap( &(fData[j]) );
            }
        }
        else if ( _data.hdr.storageType == BSHORT ) // bshort
        {
            error = znzread((char *)iData,sizeof(short),nAllocate,file) != nAllocate;
            if ( swap_bytes )
            {
                for (unsigned long j=0; j<nAllocate; j++)
                {
                    utilIO::swap( &(iData[j]) );
                    fData[j] = iData[j];
                }
            }
            else
                for (unsigned long j=0; j<nAllocate; j++)
                    fData[j] = iData[j];
        }
        else if ( _data.hdr.storageType == BLONG )  // blong
        {
            error = znzread((char *)iiData,sizeof(int),nAllocate,file) != nAllocate;
            if ( swap_bytes )
            {
                for (unsigned long j=0; j<nAllocate; j++)
                {
                    utilIO::swap( &(iiData[j]) );
                    fData[j] = iiData[j];
                }
            }
            else
                for (unsigned long j=0; j<nAllocate; j++)
                    fData[j] = iiData[j];
        }
        else if ( _data.hdr.storageType == USHORT ) // unsigned short
        {
            error = znzread((char *)iuData,sizeof(unsigned short),nAllocate,file) != nAllocate;
            if ( swap_bytes )
            {
                for (unsigned long j=0; j<nAllocate; j++)
                {
                    utilIO::swap( &(iuData[j]) );
                    fData[j] = iuData[j];
                }
            }
            else
                for (unsigned long j=0; j<nAllocate; j++)
                    fData[j] = iuData[j];
        }
        else if ( _data.hdr.storageType == UINT8 ) // unsigned char
        {
            error = znzread((char *)i8Data,sizeof(unsigned char),nAllocate,file) != nAllocate;
            for (unsigned long j=0; j<nAllocate; j++)
                fData[j] = i8Data[j];
        }
        else if ( _data.hdr.storageType == BDOUBLE ) // unsigned char
        {
            error = znzread((char *)f64Data,sizeof(double),nAllocate,file) != nAllocate;
            for (unsigned long j=0; j<nAllocate; j++)
                fData[j] = f64Data[j];
        }
        // errors?
        if ( error )
        {
            qInfo() << "* Error reading time point " << jt+1 << " for file " << fileName;
            return(1);
        }

        // Scale the data if necessary.
        FUNC_INFO << "slope" << _data.hdr.scaleSlope;
        if ( _data.hdr.scaleSlope == 0. ) _data.hdr.scaleSlope = 1.;
        for (unsigned long j=0; j<nAllocate; j++)
        {
            fData[j] *= _data.hdr.scaleSlope;
            fData[j] += _data.hdr.scaleOffset;
            // scrub data for unusual cases
            if ( qIsNaN(fData[j]) || qIsInf(fData[j]) ) fData[j] = 0.;
        }
        // Copy the volume into the time series
        if ( _data.hdr.dataType == 0 ) // scalar
        {
            FUNC_INFO << "assign scalar";
            // 4D data structures should already be allocate via readFileHeader
            // _data.timePoints[jt].f1.resize(_data.hdr.points.z);
            for (unsigned long j=0; j<nAllocate; j++)
            {
                _data.timePoints[jt].f1[j] = fData[j];
                if ( _data.hdr.dim.t > 1 ) _data.averageVolume[j] += fData[j];
            }
        }
        else if ( _data.hdr.dataType == 3 || _data.hdr.dataType == 4 ) // 3-vector
        {
            // 4D data structures should already be allocate via readFileHeader
            // _data.timePoints[jt].f3.resize(_data.hdr.points.z);
            for (unsigned long j=0; j<nAllocate/3; j++)
            {
                _data.timePoints[jt].f3[j].x = fData[3*j];
                _data.timePoints[jt].f3[j].y = fData[3*j+1];
                _data.timePoints[jt].f3[j].z = fData[3*j+2];
                if ( _data.hdr.dim.t > 1 )
                    _data.averageVolume[j] += SQR(_data.timePoints[jt].f3[j].x)
                            +                 SQR(_data.timePoints[jt].f3[j].y)
                            +                 SQR(_data.timePoints[jt].f3[j].z);
            }
        }
        else // if ( _data.hdr.dataType == 1 || _data.hdr.dataType == 2 ) // complex
        {
            // 4D data structures should already be allocate via readFileHeader
            // _data.timePoints[jt].fc.resize(_data.hdr.points.z);
            for (unsigned long j=0; j<nAllocate/2; j++)
            {
                _data.timePoints[jt].fc[j].real = fData[2*j];
                _data.timePoints[jt].fc[j].imag = fData[2*j+1];
                if ( _data.hdr.dim.t > 1 )
                    _data.averageVolume[j] += SQR(_data.timePoints[jt].fc[j].real) + SQR(_data.timePoints[jt].fc[j].imag);
            }
        }
    } // jt
    if ( _data.hdr.dim.t > 1 )
    {
        for (int j=0; j<_data.hdr.points.z; j++)
            _data.averageVolume[j] /= static_cast<double>(_data.hdr.dim.t);
    }

    znzclose(file);

    if ( output_info )
    {
        //      printf("\r%s: 100 %c\n",FileName,37);
        //      fflush(stdout);
    }
    setCurrentVolume(0);
    return(0);
}

void ImageIO::resliceImageData( ImageIO *templateImage )
{
    QProgressBar *progressBar=0;
    resliceImageData( templateImage, progressBar );
}
void ImageIO::resliceImageData( ImageIO *templateImage, QProgressBar *progressBar )
{
    // First make a copy of the data
    ImageIO copyStack = *this;

    imageHeader *templateHeader = templateImage->getImageDataHeader();

    // Copy everything first (coordinate transformations, etc)
    int nt_saved = nt();
    _data = templateImage->_data;
    // Reset image dimensions and allocate (resize) data
    setDimensions(nx(),ny(),nz(),nt_saved);

    dPoint3D fwhm;
    fwhm.x = qCeil(templateHeader->resolution.x / _data.hdr.resolution.x);
    fwhm.y = qCeil(templateHeader->resolution.y / _data.hdr.resolution.y);
    fwhm.z = qCeil(templateHeader->resolution.z / _data.hdr.resolution.z);
    iPoint3D iWidth;
    iWidth.x = qCeil(fwhm.x);
    iWidth.y = qCeil(fwhm.y);
    iWidth.z = qCeil(fwhm.z);
    iPoint3D wrap, kernelType;
    wrap.x = wrap.y = wrap.z = 0; // no wrapping
    enum kernels {Lanczos, Gaussian, Delta};
    // Gaussian if fwhm > 1; otherwise Lanczos
    kernelType.x = (fwhm.x>1);  kernelType.y = (fwhm.y>1); kernelType.z = (fwhm.z>1);

    // Dummy statement to print info about interpolation (indicated by time point "-1")
    // Space_template.x = Space_template.y = Space_template.z = 0.;
    // VolumeValue = modify->interpolateValueFromSpaceCoord( Space_template, fwhm, wrap, kernelType, -1);

    if ( progressBar )
    {
        progressBar->setMinimum(0);
        progressBar->setMaximum((nt()-1)*(nz()-1));
    }
    bool inVolume;
    iPoint3D voxel;
    // Loop over all time points
    for (int jt=0; jt<nt(); jt++)
    {
        // Loop over voxels in the template space
        for (voxel.z=0; voxel.z<nz(); voxel.z++)
        {
            int iProgress = voxel.z * jt;
            if ( progressBar ) progressBar->setValue(iProgress);
            for (voxel.y=0; voxel.y<ny(); voxel.y++)
            {
                for (voxel.x=0; voxel.x<nx(); voxel.x++)
                {
                    // Convert from data indices --> space in the registered system.
                    dPoint3D space = templateImage->convertVoxelToSpace(voxel);
                    // Interpolate to the find the volume value at this spatial coordinate
                    double value = copyStack.interpolateValueFromSpaceCoord( space, fwhm, iWidth, wrap,
                                                                             kernelType, jt, inVolume);
                    setStackValue(voxel,jt,value);
                } // indexTemplate.x
            } // indexTemplate.y
            progressBar->setValue(voxel.z);
        } // indexTemplate.z
    } // time loop
}

dPoint3D ImageIO::convertVoxelToSpace( fVoxel voxel)
{
    dPoint3D rVoxel;
    rVoxel.x = voxel.x;
    rVoxel.y = voxel.y;
    rVoxel.z = voxel.z;
    return convertVoxelToSpace(rVoxel);
}

dPoint3D ImageIO::convertVoxelToSpace( iPoint3D iVoxel)
{
    dPoint3D rVoxel;
    rVoxel.x = iVoxel.x;
    rVoxel.y = iVoxel.y;
    rVoxel.z = iVoxel.z;
    return convertVoxelToSpace(rVoxel);
}

dPoint3D ImageIO::convertVoxelToSpace( dPoint3D rVoxel)
{
    double dataVector[4];
    dataVector[0] = rVoxel.x;
    dataVector[1] = rVoxel.y;
    dataVector[2] = rVoxel.z;
    dataVector[3] = 1;

    dPoint3D space={0.,0.,0.};
    for ( int j=0; j<4; j++)
        space.x += _data.hdr.ijk_to_xyz[0][j] * dataVector[j];
    for ( int j=0; j<4; j++)
        space.y += _data.hdr.ijk_to_xyz[1][j] * dataVector[j];
    for ( int j=0; j<4; j++)
        space.z += _data.hdr.ijk_to_xyz[2][j] * dataVector[j];
    return space;
}

iPoint3D ImageIO::convertSpaceToVoxel(dPoint3D space)
{
    dPoint3D  rVoxel;
    iPoint3D iVoxel;
    rVoxel = convertSpaceToVoxelFraction(space);
    iVoxel.x = qRound(rVoxel.x);
    iVoxel.y = qRound(rVoxel.y);
    iVoxel.z = qRound(rVoxel.z);
    return iVoxel;
}
dPoint3D ImageIO::convertSpaceToVoxelFraction(dPoint3D space)
{
    dVector spaceToVoxel;
    spaceToVoxel.append(space.x);
    spaceToVoxel.append(space.y);
    spaceToVoxel.append(space.z);
    spaceToVoxel.append(1.);
    utilMath::matrixTransform(_data.hdr.xyz_to_ijk,spaceToVoxel);
    dPoint3D rVoxel;
    rVoxel.x = spaceToVoxel.at(0);
    rVoxel.y = spaceToVoxel.at(1);
    rVoxel.z = spaceToVoxel.at(2);
    return(rVoxel);
}
/*
dPoint3D ImageIO::convertSpaceToVoxelFraction(dPoint3D space)
{
    double spaceVector[4];
    spaceVector[0] = space.x;
    spaceVector[1] = space.y;
    spaceVector[2] = space.z;
    spaceVector[3] = 1;

    dPoint3D rVoxel={0.,0.,0.};    
    for ( int j=0; j<4; j++)
        rVoxel.x += _data.hdr.xyz_to_ijk[0][j] * spaceVector[j];
    for ( int j=0; j<4; j++)
        rVoxel.y += _data.hdr.xyz_to_ijk[1][j] * spaceVector[j];
    for ( int j=0; j<4; j++)
        rVoxel.z += _data.hdr.xyz_to_ijk[2][j] * spaceVector[j];

    return(rVoxel);
}
*/

double ImageIO::interpolateValueFromSpaceCoord(dPoint3D space, dPoint3D fwhm, iPoint3D iWidth, iPoint3D wrap,
                                               iPoint3D kernelType, int iTime)
// Give the coordinate "Space", interpolate to find the value in the stack at time = iTime
{
    // Convert the space coordinate to a voxel coordinate
    bool inVolume;
    dPoint3D rIndex = convertSpaceToVoxelFraction(space);
    return interpolateValueFromrIndex(rIndex, fwhm, iWidth, wrap, kernelType, iTime, inVolume);
}

double ImageIO::interpolateValueFromSpaceCoord(dPoint3D space, dPoint3D fwhm, iPoint3D iWidth, iPoint3D wrap,
                                               iPoint3D kernelType, int iTime, bool &inVolume)
// Give the coordinate "Space", interpolate to find the value in the stack at time = iTime
{
    // Convert the space coordinate to a voxel coordinate
    dPoint3D rIndex = convertSpaceToVoxelFraction(space);
    return interpolateValueFromrIndex(rIndex, fwhm, iWidth, wrap, kernelType, iTime, inVolume);
}

double ImageIO::interpolateValueFromSpaceCoordBilinear(dPoint3D space, int iTime)
// Give the coordinate "Space", interpolate to find the value in the stack at time = iTime
{
    // Convert the space coordinate to a voxel coordinate
    dPoint3D rIndex = convertSpaceToVoxelFraction(space);
    return interpolateValueFromrIndexBilinear(rIndex, iTime);
}

double ImageIO::interpolateValueFromrIndex(dPoint4D rIndex, dPoint3D fwhm, iPoint3D iWidth, iPoint3D wrap,
                                           iPoint3D kernelType, int iTime, bool &inVolume)
{
    dPoint3D rIndex3d;  rIndex3d.x = rIndex.x;  rIndex3d.y = rIndex.y;  rIndex3d.z = rIndex.z;
    return interpolateValueFromrIndex(rIndex3d, fwhm, iWidth, wrap, kernelType, iTime, inVolume);
}
double ImageIO::interpolateValueFromrIndex(dPoint4D rIndex, dPoint3D fwhm, iPoint3D iWidth, iPoint3D wrap,
                                           iPoint3D kernelType, int iTime)
{
    bool inVolume;
    dPoint3D rIndex3d;  rIndex3d.x = rIndex.x;  rIndex3d.y = rIndex.y;  rIndex3d.z = rIndex.z;
    return interpolateValueFromrIndex(rIndex3d, fwhm, iWidth, wrap, kernelType, iTime, inVolume);
}
double ImageIO::interpolateValueFromrIndex(dPoint3D rIndex, dPoint3D fwhm, iPoint3D iWidth,
                                           iPoint3D wrap, iPoint3D kernelType, int iTime)
{
    bool inVolume;
    return interpolateValueFromrIndex(rIndex, fwhm, iWidth, wrap, kernelType, iTime, inVolume);
}
double ImageIO::interpolateValueFromrIndex(dPoint3D rIndex, dPoint3D fwhm, iPoint3D iWidth, iPoint3D wrap,
                                           iPoint3D kernelType, int iTime, bool &inVolume)
// Give the index coordinate "rIndex", interpolate to find the value in the stack at time = iTime
//
// rIndex    : 3-dimensional index
// wrap      : controls wrapping in x, y, z
// kernelType: 0 = Lanczos, 1 = Gaussian, 2 = Delta
// fwhm      : fwhm of kernel in voxel units; usually 1 for interpolation or relWidth for down-sampling
// iWidth    : include 2*iWidth points, so iWidth=1 uses just 2 points (neighbors)
//             usually iWidth = (int)fwhm
//
// iTime     : time point in stack to interpolate
{
    enum kernels {Lanczos, Gaussian, Delta};

    if ( kernelType.x == Delta ) iWidth.x = 0;
    if ( kernelType.y == Delta ) iWidth.y = 0;
    if ( kernelType.z == Delta ) iWidth.z = 0;

    if ( iTime == -1 )
    { // just print and return
        QString infoString = QString("Kernel width = (%1, %2, %3)\n").arg(iWidth.x).arg(iWidth.y).arg(iWidth.z);
        qInfo() << infoString;
        infoString = "Kernel type  = (";
        switch ( kernelType.x )
        {
        case Lanczos:  infoString += "Lanczos, ";  break;
        case Gaussian: infoString += "Gaussian, "; break;
        case Delta:    infoString += "Delta, ";    break;
        }
        switch ( kernelType.y )
        {
        case Lanczos:  infoString += "Lanczos, ";  break;
        case Gaussian: infoString += "Gaussian, "; break;
        case Delta:    infoString += "Delta, ";    break;
        }
        switch ( kernelType.z )
        {
        case Lanczos:  infoString += "Lanczos)";  break;
        case Gaussian: infoString += "Gaussian)"; break;
        case Delta:    infoString += "Delta)";    break;
        }
        return(0.);
    }

    dVector coeffX, coeffY, coeffZ;
    coeffX.resize(2*iWidth.x);   coeffY.resize(2*iWidth.y);   coeffZ.resize(2*iWidth.z);

    iPoint4D dim = getDimensions();

    // index0 is the nearest integer location in the array
    iPoint3D index0;
    index0.x = qRound(rIndex.x);   index0.y = qRound(rIndex.y);   index0.z = qRound(rIndex.z);

    if ( wrap.x )
    { // wrap x
        if ( index0.x < 0 )
            index0.x += dim.x;
        else if ( index0.x >= dim.x )
            index0.x -= dim.x;
    }
    if ( wrap.y )
    { // wrap y
        if ( index0.y < 0 )
            index0.y += dim.y;
        else if ( index0.y >= dim.y )
            index0.y -= dim.y;
    }
    if ( wrap.z )
    { // wrap z
        if ( index0.z < 0 )
            index0.z += dim.z;
        else if ( index0.z >= dim.z )
            index0.z -= dim.z;
    }
    inVolume = ( index0.x>=0 && index0.x<dim.x &&
                 index0.y>=0 && index0.y<dim.y &&
                 index0.z>=0 && index0.z<dim.z );

    double VolumeValue = 0.;
    if ( inVolume )
    {
        // define low/high limits
        iPoint3D iLow, iHigh, jIndex, jRel;
        if ( kernelType.x == Delta )
            iLow.x = iHigh.x = qRound(rIndex.x);
        else
        {
            iLow.x  = qFloor((double)rIndex.x);
            iHigh.x = qCeil((double)rIndex.x);
        }
        if ( kernelType.y == Delta )
            iLow.y = iHigh.y = qRound(rIndex.y);
        else
        {
            iLow.y  = qFloor((double)rIndex.y);
            iHigh.y = qCeil((double)rIndex.y);
        }
        if ( kernelType.z == Delta )
            iLow.z = iHigh.z = qRound(rIndex.z);
        else
        {
            iLow.z  = qFloor((double)rIndex.z);
            iHigh.z = qCeil((double)rIndex.z);
        }
        // For intermediate points, at least 2 neighbors, plus addition specified by iWidth
        // For on-grid points,
        iPoint3D step;
        step.x = qMax(1,iHigh.x - iLow.x + (iWidth.x-1));
        step.y = qMax(1,iHigh.y - iLow.y + (iWidth.y-1));
        step.z = qMax(1,iHigh.z - iLow.z + (iWidth.z-1));

        // Create interpolation coefficients in each dimension
        for (jIndex.x=iLow.x-iWidth.x, jRel.x=0;
             jIndex.x<iHigh.x+iWidth.x; jIndex.x++, jRel.x++)
        {
            if ( kernelType.x == Gaussian )
                coeffX[jRel.x] = utilMath::GaussKernel((double)(rIndex.x-jIndex.x),fwhm.x);
            else if ( kernelType.x == Delta )
                coeffX[jRel.x] = 1.;
            else
                coeffX[jRel.x] = utilMath::LanczosKernel((double)(rIndex.x-jIndex.x),iWidth.x);
        }
        for (jIndex.y=iLow.y-iWidth.y, jRel.y=0;
             jIndex.y<iHigh.y+iWidth.y; jIndex.y++, jRel.y++)
        {
            if ( kernelType.y == Gaussian )
                coeffY[jRel.y] = utilMath::GaussKernel((double)(rIndex.y-jIndex.y),fwhm.y);
            else if ( kernelType.y == Delta )
                coeffY[jRel.y] = 1.;
            else
                coeffY[jRel.y] = utilMath::LanczosKernel((double)(rIndex.y-jIndex.y),iWidth.y);
        }
        for (jIndex.z=iLow.z-iWidth.z, jRel.z=0;
             jIndex.z<iHigh.z+iWidth.z; jIndex.z++, jRel.z++)
        {
            if ( kernelType.z == Gaussian )
                coeffZ[jRel.z] = utilMath::GaussKernel((double)(rIndex.z-jIndex.z),fwhm.z);
            else if ( kernelType.z == Delta )
                coeffZ[jRel.z] = 1.;
            else
                coeffZ[jRel.z] = utilMath::LanczosKernel((double)(rIndex.z-jIndex.z),iWidth.z);
        }

        double sum, sum_weight;
        sum = sum_weight = 0.;
        // Now loop over a finite range and interpolate
        for (jIndex.z=iLow.z-iWidth.z, jRel.z=0;
             jIndex.z<iHigh.z+iWidth.z; jIndex.z++, jRel.z++)
        {
            index0.z = jIndex.z;  // We can re-use index0 because we have the low/high range
            if ( wrap.z )
            { // wrap z
                if ( index0.z < 0 )
                    index0.z += dim.z;
                else if ( index0.z >= dim.z )
                    index0.z -= dim.z;
            }
            for (jIndex.y=iLow.y-iWidth.y, jRel.y=0;
                 jIndex.y<iHigh.y+iWidth.y; jIndex.y++, jRel.y++)
            {
                // wrap y
                index0.y = jIndex.y;
                if ( wrap.y )
                { // wrap y
                    if ( index0.y < 0 )
                        index0.y += dim.y;
                    else if ( index0.y >= dim.y )
                        index0.y -= dim.y;
                }
                for (jIndex.x=iLow.x-iWidth.x, jRel.x=0;
                     jIndex.x<iHigh.x+iWidth.x; jIndex.x++, jRel.x++)
                {
                    // wrap x
                    index0.x = jIndex.x;
                    if ( wrap.x )
                    { // wrap x
                        if ( index0.x < 0 )
                            index0.x += dim.x;
                        else if ( index0.x >= dim.x )
                            index0.x -= dim.x;
                    }
                    inVolume = ( index0.x>=0 && index0.x<dim.x &&
                                 index0.y>=0 && index0.y<dim.y &&
                                 index0.z>=0 && index0.z<dim.z );
                    if ( inVolume )
                    {
                        // sum over neighbors in x
                        double weight = coeffX[jRel.x] * coeffY[jRel.y] * coeffZ[jRel.z];
                        sum_weight += weight;
                        sum        += weight * getStackValue(index0,iTime);
                    } // in volume
                } // jIndex.x
            } // jIndex.y
        } // jIndex.z
        if ( sum_weight == 0. ) sum_weight = 1.;  // prevent division by zero
        VolumeValue = sum / sum_weight;
    } // if inVolume
    return(VolumeValue);
}

double ImageIO::interpolateValueFromrIndexBilinear(dPoint3D rIndex, int iTime)
// Give the index coordinate "rIndex", interpolate to find the value in the stack at time = iTime
//
// rIndex    : 3-dimensional index
// iTime     : time point in stack to interpolate
{

    iPoint3D iLow, iHigh;
    iLow.x  = qFloor((double)rIndex.x);
    iHigh.x = qCeil((double)rIndex.x);
    iLow.y  = qFloor((double)rIndex.y);
    iHigh.y = qCeil((double)rIndex.y);
    iLow.z  = qFloor((double)rIndex.z);
    iHigh.z = qCeil((double)rIndex.z);

    dVector coeffX, coeffY, coeffZ;
    coeffX.resize(2);   coeffY.resize(2);   coeffZ.resize(2);

    iPoint4D dim = getDimensions();
    double sum, sum_weight;
    sum = sum_weight = 0.;
    iPoint3D jIndex, jRel;
    for (jIndex.x=iLow.x-1, jRel.x=0;
         jIndex.x<iHigh.x+1; jIndex.x++, jRel.x++)
    {
        coeffX[jRel.x] = utilMath::GaussKernel((double)(rIndex.x-jIndex.x),1.);
        for (jIndex.y=iLow.y-1, jRel.y=0;
             jIndex.y<iHigh.y+1; jIndex.y++, jRel.y++)
        {
            coeffY[jRel.y] = utilMath::GaussKernel((double)(rIndex.y-jIndex.y),1.);
            for (jIndex.z=iLow.z-1, jRel.z=0;
                 jIndex.z<iHigh.z+1; jIndex.z++, jRel.z++)
            {
                coeffZ[jRel.z] = utilMath::GaussKernel((double)(rIndex.z-jIndex.z),1.);
                bool inVolume = ( jIndex.x>=0 && jIndex.x<dim.x &&
                                  jIndex.y>=0 && jIndex.y<dim.y &&
                                  jIndex.z>=0 && jIndex.z<dim.z );
                if ( inVolume )
                {
                    // sum over neighbors in x
                    double weight = coeffX[jRel.x] * coeffY[jRel.y] * coeffZ[jRel.z];
                    sum_weight += weight;
                    sum        += weight * getStackValue(jIndex,iTime);
                } // in volume
            } // z
        } // x
    } // y
    if ( sum_weight == 0. ) sum_weight = 1.;  // prevent division by zero
    double volumeValue = sum / sum_weight;
    return(volumeValue);
}

void ImageIO::resample( iPoint3D dimensionNew, dPoint3D resolutionNew )
{ // Resample a data_stack into DimensionNew. Assume centered coordinates.

    // Make a copy of the original data as the source.
    ImageIO original;
    original.copyTemplate(this);
    // Use centered coordinates for the source
    original.setParallelCenteredCoordinates();

    // Change the dimension and reset the input stack, which becomes the target.
    setDimensions(dimensionNew);   // this does not set any resolution, etc.
    setResolution(resolutionNew);  // this defines centered coordinates for the target

    // Changing this transformation does not affect origin & direction, so the transformation can be remade at the end
    //    setParallelCenteredCoordinates();

    dPoint3D resolution = getImageResolution();
    dPoint3D relRes;
    relRes.x = qCeil(resolutionNew.x/resolution.x);
    relRes.y = qCeil(resolutionNew.y/resolution.y);
    relRes.z = qCeil(resolutionNew.z/resolution.z);
    iPoint3D iWidth;
    iWidth.x = qCeil(relRes.x);
    iWidth.y = qCeil(relRes.y);
    iWidth.z = qCeil(relRes.z);
    iPoint3D wrap;  wrap.x = wrap.y = wrap.z = 0; // no wrapping
    iPoint3D GaussianKernel;
    GaussianKernel.x = (relRes.x>1);
    GaussianKernel.y = (relRes.y>1);
    GaussianKernel.z = (relRes.z>1);

    bool inVolume;
    for (int jt=0; jt<nt(); jt++)
    {
        iPoint3D indexReg;
        iPoint3D dim = getImageDimensions();  // new dimension
        // Loop over voxels in the downsampled resolution space.
        for (indexReg.z=0; indexReg.z<dim.z; indexReg.z++)
        {
            for (indexReg.y=0; indexReg.y<dim.y; indexReg.y++)
            {
                for (indexReg.x=0; indexReg.x<dim.x; indexReg.x++)
                {
                    // Convert from data indices --> space in the registered system.
                    dPoint3D space = convertVoxelToSpace( indexReg );
                    // Interpolate to find the volume value at this spatial coordinate
                    double volumeValue = original.interpolateValueFromSpaceCoord(space,relRes,iWidth,wrap,
                                                                                 GaussianKernel,jt,inVolume);
                    setStackValue(indexReg,jt,volumeValue);
                }
            }
        }
    } // time loop
    //    setParallelCoordinates();
}

void ImageIO::upSampleOrDownSample( iPoint3D dimensionNew )
{ // Down-sample a data_stack into DimensionLow. Assume centered coordinates in the low-res stack.
    iPoint3D dimension  = getImageDimensions();
    dPoint3D resolution = getImageResolution();
    FUNC_INFO << dimension.x << dimension.y << dimension.z;
    FUNC_INFO << dimensionNew.x << dimensionNew.y << dimensionNew.z;
    FUNC_INFO << resolution.x << resolution.y << resolution.z;
    dPoint3D  resolutionNew;
    resolutionNew.x = resolution.x * static_cast<double>(dimension.x) / static_cast<double>(dimensionNew.x);
    resolutionNew.y = resolution.y * static_cast<double>(dimension.y) / static_cast<double>(dimensionNew.y);
    resolutionNew.z = resolution.z * static_cast<double>(dimension.z) / static_cast<double>(dimensionNew.z);
    FUNC_INFO << resolutionNew.x << resolutionNew.y << resolutionNew.z;
    resample(dimensionNew, resolutionNew);
}

void ImageIO::edgeFilter(int iTime)
{ // Resample a data_stack into DimensionNew. Assume centered coordinates.
    FUNC_ENTER << iTime;
    // Make a copy of the original data as the source.
    FUNC_INFO << 0.5;
    ImageIO original;
    FUNC_INFO << 0.75;
    original.copyTemplateVolume(this,iTime);
    FUNC_INFO << 1;

    double sobel[3][3] = { {3. , 10. ,  3.},
                           {0. ,  0. ,  0.},
                           {-3., -10., -3.} };

    iPoint3D voxel;
    iPoint3D dim = getImageDimensions();  // new dimension
    FUNC_INFO << "dim" << dim.x << dim.y << dim.z;
    // Loop over voxels in the downsampled resolution space.
    for (voxel.z=0; voxel.z<dim.z; voxel.z++)
    {
        for (voxel.y=0; voxel.y<dim.y; voxel.y++)
        {
            for (voxel.x=0; voxel.x<dim.x; voxel.x++)
            {
                double sumX=0.;  double sumY=0.;  double ampSum=0.;
                for (int jXRel=-1; jXRel<=1; jXRel+=2)
                {
                    int iX = voxel.x + jXRel;
                    if ( iX < 0 ) iX = dim.x-1;
                    if ( iX >= dim.x ) iX = 0;
                    for (int jYRel=-1; jYRel<=1; jYRel+=2)
                    {
                        int iY = voxel.y + jYRel;
                        if ( iY < 0 ) iY = dim.y-1;
                        if ( iY >= dim.y ) iY = 0;
                        double amp = original.getStackValue(iX,iY,voxel.z,0);
                        sumX += amp * sobel[jXRel+1][jYRel+1];
                        sumY += amp * sobel[jYRel+1][jXRel+1];
                        ampSum += amp;
                    } // jYRel
                } // jXRel
                ampSum += original.getStackValue(voxel.x,voxel.y,voxel.z,0);
                if ( ampSum == 0. )
                    setStackValue(voxel,0,0.);
                else // normalizing by local amp seems to be important to keep uniform intensities
                    setStackValue(voxel,0,qSqrt(sumX*sumX + sumY*sumY)/ampSum);
            } // x
        } // y
    } // z
    FUNC_EXIT;
}

void ImageIO::edgeFilter3D(int iTime)
{ // Resample a data_stack into DimensionNew. Assume centered coordinates.
    FUNC_ENTER;
    // Make a copy of the original data as the source.
    ImageIO original;
    original.copyTemplateVolume(this,iTime);

    double sobelMinus[3][3] = {  {1., 2., 1.},
                                 {2., 4., 2.},
                                 {1., 2., 1.} };
    double sobelPlus[3][3]  = {  {-1., -2., -1.},
                                 {-2., -4., -2.},
                                 {-1., -2., -1.} };
    double sobel2D[3][3]    = { {1. , 2. ,  1.},
                                {0. ,  0.,  0.},
                                {-1., -2., -1.} };

    iPoint3D voxel;
    iPoint3D dim = getImageDimensions();  // new dimension
    // Loop over voxels in the downsampled resolution space.
    for (voxel.z=0; voxel.z<dim.z; voxel.z++)
    {
        for (voxel.y=0; voxel.y<dim.y; voxel.y++)
        {
            for (voxel.x=0; voxel.x<dim.x; voxel.x++)
            {
                // z
                double sumZ=0.;     double ampSum = 0.;
                for (int jXRel=-1; jXRel<=1; jXRel++)
                {
                    int iX = voxel.x + jXRel;  if ( iX < 0 ) iX = dim.x-1;  if ( iX >= dim.x ) iX = 0;
                    for (int jYRel=-1; jYRel<=1; jYRel++)
                    {
                        int iY = voxel.y + jYRel;   if ( iY < 0 ) iY = dim.y-1;  if ( iY >= dim.y ) iY = 0;
                        if ( voxel.z == 0 || voxel.z == dim.z-1 )
                        {
                            double amp = original.getStackValue(iX,iY,voxel.z,0);
                            sumZ += amp * sobel2D[jXRel+1][jYRel+1];
                            ampSum += amp;
                        }
                        else
                        {
                            int iZMinus  = voxel.z - 1;
                            double amp = original.getStackValue(iX,iY,iZMinus,0);
                            sumZ += amp * sobelMinus[jXRel+1][jYRel+1];
                            ampSum += amp;
                            int iZPlus = voxel.z + 1;
                            amp = original.getStackValue(iX, iY, iZPlus, 0);
                            sumZ += amp * sobelPlus[jXRel+1][jYRel+1];
                            ampSum += amp;
                        } // end slices
                    } // jYRel
                } // jXRel

                // x: average in y,z and differential in x
                double sumX=0.;
                for (int jYRel=-1; jYRel<=1; jYRel++)
                {
                    int iY = voxel.y + jYRel;  if ( iY < 0 ) iY = dim.y-1;  if ( iY >= dim.y ) iY = 0;
                    for (int jZRel=-1; jZRel<=1; jZRel++)
                    {
                        int iZ = voxel.z + jZRel;   if ( iZ < 0 ) iZ = dim.z-1;  if ( iZ >= dim.z ) iZ = 0;
                        int iXMinus  = voxel.x - 1; if ( iXMinus < 0 ) iXMinus = dim.x-1;
                        double amp = original.getStackValue(iXMinus,iY,iZ,0);
                        sumX += amp * sobelMinus[jYRel+1][jZRel+1];
                        int iXPlus = voxel.x + 1;  if ( iXPlus >= dim.x ) iXPlus = 0;
                        amp = original.getStackValue(iXPlus,iY,iZ,0);
                        sumX += amp * sobelPlus[jYRel+1][jZRel+1];
                    } // jZRel
                } // jYRel
                // y
                double sumY=0.;
                for (int jXRel=-1; jXRel<=1; jXRel++)
                {
                    int iX = voxel.x + jXRel;  if ( iX < 0 ) iX = dim.x-1;  if ( iX >= dim.x ) iX = 0;
                    for (int jZRel=-1; jZRel<=1; jZRel++)
                    {
                        int iZ = voxel.z + jZRel;   if ( iZ < 0 ) iZ = dim.z-1;  if ( iZ >= dim.z ) iZ = 0;
                        int iYMinus  = voxel.y - 1; if ( iYMinus < 0 ) iYMinus = dim.y-1;
                        double amp = original.getStackValue(iX,iYMinus,iZ,0);
                        sumY += amp * sobelMinus[jXRel+1][jZRel+1];
                        int iYPlus = voxel.y + 1;  if ( iYPlus >= dim.y ) iYPlus = 0;
                        amp = original.getStackValue(iX, iYPlus,iZ,0);
                        sumY += amp * sobelPlus[jXRel+1][jZRel+1];
                    } // jZRel
                } // jXRel

                if ( ampSum == 0. )
                    setStackValue(voxel,0,0.);
                else // normalizing by local amp seems to be important to keep uniform intensities across RF bias
                    setStackValue(voxel,0,qSqrt(sumX*sumX + sumY*sumY + sumZ*sumZ)/ampSum);
            } // x
        } // y
    } /// z

    FUNC_EXIT;
}

void ImageIO::edgeFilter2D(int iTime)
{ // Resample a data_stack into DimensionNew. Assume centered coordinates.
    FUNC_ENTER;
    // Make a copy of the original data as the source.
    ImageIO original;
    original.copyTemplateVolume(this,iTime);

    double sobel2D[3][3]    = { {1. , 2. ,  1.},
                                {0. ,  0.,  0.},
                                {-1., -2., -1.} };

    iPoint3D voxel;
    iPoint3D dim = getImageDimensions();  // new dimension
    // Loop over voxels in the downsampled resolution space.
    for (voxel.z=0; voxel.z<dim.z; voxel.z++)
    {
        for (voxel.y=0; voxel.y<dim.y; voxel.y++)
        {
            for (voxel.x=0; voxel.x<dim.x; voxel.x++)
            {
                // z
                double sum=0.;     double ampSum = 0.;
                for (int jXRel=-1; jXRel<=1; jXRel++)
                {
                    int iX = voxel.x + jXRel;  if ( iX < 0 ) iX = dim.x-1;  if ( iX >= dim.x ) iX = 0;
                    for (int jYRel=-1; jYRel<=1; jYRel++)
                    {
                        int iY = voxel.y + jYRel;   if ( iY < 0 ) iY = dim.y-1;  if ( iY >= dim.y ) iY = 0;
                        double amp = original.getStackValue(iX,iY,voxel.z,0);
                        sum += amp * sobel2D[jXRel+1][jYRel+1];
                        ampSum += amp;
                    } // jYRel
                } // jXRel
                if ( ampSum == 0. )
                    setStackValue(voxel,0,0.);
                else // normalizing by local amp seems to be important to keep uniform intensities across RF bias
                    setStackValue(voxel,0,sum/ampSum);
            } // x
        } // y
    } /// z

    FUNC_EXIT;
}

void ImageIO::fillHistogram(int iTime, double lowerLimit, double upperLimit, dVector &histogram)
{ // histogram bins must already be allocated
    FUNC_ENTER << lowerLimit << upperLimit;
    bitMap overlay;
    fillHistogram(iTime, lowerLimit, upperLimit, overlay, histogram);
}

void ImageIO::fillHistogram(int iTime, double lowerLimit, double upperLimit, bitMap overlay, dVector &histogram)
{ // histogram bins must already be allocated
    FUNC_ENTER << lowerLimit << upperLimit;
    iPoint3D dim = getImageDimensions();
    int nBins = histogram.size();
    double binStep   = (upperLimit - lowerLimit) / static_cast<double>(nBins);
    histogram.fill(0.,nBins);

    if ( overlay.someBitsAreSet() )
    {
        iPoint3D voxel;
        for ( voxel.z=0; voxel.z<dim.z; voxel.z++ )
        {
            for ( voxel.y=0; voxel.y<dim.y; voxel.y++ )
            {
                for ( voxel.x=0; voxel.x<dim.x; voxel.x++ )
                {
                    double value = getStackValue(voxel,iTime);
                    if ( overlay.getBitMapValue(voxel) != 0. )
                    {
                        int iBin = (value - lowerLimit) / binStep;
                        iBin = qMax(0,(qMin(iBin,nBins-1)));  // set low values to 0 and high values to nBins-1
                        histogram[iBin]++;
                    }
                } // x
            } // y
        } // z
    } // bitmap is defined
    else
    { // use all voxels
        for (int jVoxel=0; jVoxel<voxelsPerVolume(); jVoxel++)
        {
            double value = getStackValue(jVoxel,iTime);
            int iBin = (value - lowerLimit) / binStep;
            iBin = qMax(0,(qMin(iBin,nBins-1)));  // set low values to 0 and high values to nBins-1
            histogram[iBin]++;
        }
    }
}

dPoint2D ImageIO::getMinAndMax(int iTime)
{
    dPoint2D minMax = {{1.e30},{-1.e30}};
    for (int jVoxel=0; jVoxel<voxelsPerVolume(); jVoxel++)
    {
        double value = getStackValue(jVoxel,iTime);
        if ( !qIsNaN(value) && !qIsInf(value) )
        {
            minMax.x = qMin(value,minMax.x);
            minMax.y = qMax(value,minMax.y);
        }
    }
    FUNC_EXIT << "*** minMax:" << minMax.x << minMax.y;
    return minMax;
}

double ImageIO::getAverageSignalAfterThreshold(int iTime, double peakFraction)
{
    dPoint2D minMax = getMinAndMax(iTime);
    dVector histogram;  histogram.resize(100);
    bitMap overlay;  // needs to be passed here, but does not need to be defined
    fillHistogram(0, minMax.x, minMax.y, overlay, histogram);
    FUNC_INFO << "histogram:" << histogram;

    //////////////////////////////////////////////////////////////////
    int nBins = histogram.size();
    double yMax=0.;  int iBinMax=0;
    for (int jBin=0; jBin<nBins-1; jBin++)  // do not include the last "overflow" bin
    {
        double yValue = histogram[jBin];
        if ( yValue > yMax )
        {
            yMax = yValue;
            iBinMax = jBin;
        }
    }
    int iBinSelect = iBinMax;
    FUNC_INFO << "*** iBinMax =" << iBinMax;
    for (int jBin=iBinMax; jBin<nBins-1; jBin++)
    {
        if ( histogram[jBin] < yMax * peakFraction )
        {
            iBinSelect = jBin;
            break;
        }
    }
    double binStep   = (minMax.y - minMax.x) / static_cast<double>(nBins);
    // sum over voxels above the threshold
    double average = 0.;
    int nVoxelsAboveThreshold=0;
    for (int jBin=iBinSelect; jBin<nBins-1; jBin++)
    {
        double signal = binStep * jBin + minMax.x;
        average += signal * histogram[jBin];
        nVoxelsAboveThreshold += histogram[jBin];
    }
    if ( nVoxelsAboveThreshold != 0 )
        average /= static_cast<double>(nVoxelsAboveThreshold);
    FUNC_INFO << "*** average =" << average;
    return average;
}

void ImageIO::multiplyConstant(double multiplier)
{
    for (int jTime=0; jTime<nt(); jTime++)
    {
        for (int jVoxel=0; jVoxel<voxelsPerVolume(); jVoxel++)
        {
            double value = getStackValue(jVoxel,jTime);
            setStackValue(jVoxel,jTime,value*multiplier);
        }
    }
}

void ImageIO::multiplyConstant(int iTime, double multiplier)
{
    for (int jVoxel=0; jVoxel<voxelsPerVolume(); jVoxel++)
    {
        double value = getStackValue(jVoxel,iTime);
        setStackValue(jVoxel,iTime,value*multiplier);
    }
}

void ImageIO::addImageData(ImageIO *addFile, int iTime, int iTimeAdd)
{
    iPoint3D dim = getImageDimensions();
    iPoint3D dimAdd = addFile->getImageDimensions();
    if ( dimAdd.x != dim.x || dimAdd.y != dim.y || dimAdd.z != dim.z ) return;
    if ( iTime >= nt() || iTimeAdd >= addFile->nt() ) return;

    for (int jVoxel=0; jVoxel<voxelsPerVolume(); jVoxel++)
    {
        double value    = getStackValue(jVoxel,iTime);
        double valueAdd = addFile->getStackValue(jVoxel,iTimeAdd);
        double newValue = qAbs(value + valueAdd);
        setStackValue(jVoxel,iTime,newValue);
    }
}

void bitMap::setRectangularOverlay( iPoint3D center, iPoint3D width, bool clear)
{
    if ( clear )
        clearBitMap();

    _voxelSet.clear();
    iPoint3D offset;
    offset.x = center.x - width.x/2;
    offset.y = center.y - width.y/2;
    offset.z = center.z - width.z/2;
    iPoint3D voxel;
    iPoint3D dim = getImageDimensions();
    for (voxel.z = qMax(offset.z,0); voxel.z < qMin(offset.z+width.z,dim.z); voxel.z++)
    {
        for (voxel.y = qMax(offset.y,0); voxel.y < qMin(offset.y+width.y,dim.y); voxel.y++)
        {
            for (voxel.x = qMax(offset.x,0); voxel.x < qMin(offset.x+width.x,dim.x); voxel.x++)
                setBitMapValue(voxel,1.,Color_Undefined,false);
        }
    }
    updateBitMapVoxelSet();  // don't change the color, which should be carried by the input voxelSet
}

void bitMap::setOrAppendBitMap( bitMap overlay, bool clear)
{
    if ( clear )
        clearBitMap();

    _voxelSet.resize(overlay._voxelSet.size());
    for (int jVoxel=0; jVoxel<_voxelSet.size(); jVoxel++)
        setBitMapValue(overlay._voxelSet[jVoxel],false);
    updateBitMapVoxelSet();  // don't change the color, which should be carried by the input voxelSet
    setFileName("interactive.ovl");
}

void bitMap::setOverlayData( OVLFile ovl, bool clear)
{
    FUNC_ENTER << ovl.name << clear;
    if ( clear )
        clearBitMap();

    _voxelSet.resize(ovl.voxelList.size());
    for (int jVoxel=0; jVoxel<ovl.voxelList.size(); jVoxel++)
    {
        double bitValue = getBitMapValue(ovl.voxelList[jVoxel]);
        if ( ovl.voxelList[jVoxel].binaryValue > bitValue )
            setBitMapValue(ovl.voxelList[jVoxel],false);
    }
    updateBitMapVoxelSet();  // don't change the color, which should be carried by the input voxelSet
    setFileName(ovl.path);
}

int bitMap::readBitmapFile(QString fileName, overlayColor color )
{
    FUNC_ENTER;
    if ( readFileHeader(fileName, false) )  return(1);
    if ( nt() != 1 ) return (1);
    if ( readImageData() ) return(1);
    updateBitMapVoxelSet(color);
    FUNC_EXIT;
    return(0);
}

int bitMap::readOverlayFile(QString fileName, overlayColor color )
{
    FUNC_ENTER << "color" << color << fileName << _overlayInputLogic;
    if ( nt() != 1 )
        qFatal("Error: in readOverlayFile, # time points != 1");
    CheckVoxelIndices(0,0,0,"readOverlayFile");

    QFile file(fileName);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return(1);
    QTextStream in_stream(&file);

    bool preExistingOverlay = someBitsAreSet();
    bitMap original;
    if ( preExistingOverlay ) original = *this;

    clearBitMap();
    _voxelSet.resize(0);
    fVoxel voxel;
    int iLine=1;
    QRegularExpression whiteSpaceComma("[,\\s]");// match a comma or a space
    while (!in_stream.atEnd())
    {
        QString line = in_stream.readLine();
        QString unCommented = line.left(line.indexOf("#"));
        if ( unCommented.isEmpty() )
            continue;
        QStringList valueList = unCommented.split(whiteSpaceComma,Qt::SkipEmptyParts);
        if ( !valueList.at(0).compare("dimensions",Qt::CaseInsensitive) )
        { // typically this should be on the 1st line, if used
            iPoint3D dim;
            dim.x = valueList.at(1).toInt();
            dim.y = valueList.at(2).toInt();
            dim.z = valueList.at(3).toInt();
            iPoint3D dim0 = getImageDimensions();
            if ( dim.x != dim0.x || dim.y != dim0.y || dim.z != dim0.z )
            {
//                qInfo() << "Error: dimensions of overlay do not match dimensions of image volume in bitMap::readOverlayFile.";
                return(1);
            }
        }
        else
        {
            voxel.x = valueList.at(0).toInt();
            voxel.y = valueList.at(1).toInt();
            voxel.z = valueList.at(2).toInt();
            if ( valueList.size() < 3 )
                qInfo() << "Error: only read " << valueList.size() << " values on line " << iLine;
            else if ( valueList.size() == 3 )
                voxel.binaryValue = 1.;
            else
                voxel.binaryValue = valueList.at(3).toFloat();
            voxel.color = color;
            setBitMapValue(voxel,false);
            _data.empty = false;
        }
        iLine++;
    }

    // Apply overlay logic with a pre-existing overlay
    if ( preExistingOverlay )
    { // loop over voxels in new overlay, and test the same voxel in the old overlay
        for (int jVoxel=0; jVoxel<voxelsPerVolume(); jVoxel++)
        {
            double value   = getBitMapValue(jVoxel);
            double valueOld = original.getBitMapValue(jVoxel);
            if ( _overlayInputLogic == 0 )        // OR = +
                value = qMin(1.,qMax(0.,valueOld+value));
            else if ( _overlayInputLogic == 1 )   // AND = *
                value = qMin(1.,qMax(0.,valueOld*value));
            else                            // NOT = -
                value = qMin(1.,qMax(0.,valueOld-value));
            setStackValue(jVoxel,0,value);
        }
    }
    updateBitMapVoxelSet(color);
    file.close();
    setFileName(fileName);
    FUNC_EXIT << _voxelSet.size();
    return(0);
}

int bitMap::writeOverlayFile(QString fileName)
{
    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        qInfo() << "Error opening file " << fileName;
        return(1);
    }
    QTextStream out_stream(&file);

    out_stream << "dimensions " << nx() << " " << ny() << " " << nz() << Qt::endl;
    for ( int jVoxel=0; jVoxel<_voxelSet.size(); jVoxel++)
    {
        if ( _voxelSet[jVoxel].binaryValue == 1. )
            out_stream << _voxelSet[jVoxel].x << " " << _voxelSet[jVoxel].y << " " << _voxelSet[jVoxel].z << Qt::endl;
        else if ( _voxelSet[jVoxel].binaryValue > 0.01 )
            out_stream << _voxelSet[jVoxel].x << " " << _voxelSet[jVoxel].y << " " << _voxelSet[jVoxel].z << " "
                       << _voxelSet[jVoxel].binaryValue << Qt::endl;
    }
    file.close();
    return(0);
}
int bitMap::writeBitmapFile(QString fileName, bool useCompression)
{
    if ( writeNiftiHeader(fileName, useCompression) ) return(1);
    if ( writeNiftiData(fileName, useCompression) )   return(1);
    return(0);
}

void bitMap::moveOverlay(iPoint3D move)
{
    if ( _data.empty )
    {
        qInfo() << "Empty overlay in mirror function: programming error??";
        return;
    }

    iPoint3D dim = getImageDimensions();
    VoxelSet voxelMove;
    voxelMove.resize(_voxelSet.size());
    for ( int jVoxel=0; jVoxel<_voxelSet.size(); jVoxel++)
    {
        voxelMove[jVoxel] = _voxelSet[jVoxel];
        voxelMove[jVoxel].x += move.x;
        voxelMove[jVoxel].y += move.y;
        voxelMove[jVoxel].z += move.z;
        if ( voxelMove[jVoxel].x < 0 ) voxelMove[jVoxel].x = dim.x + voxelMove[jVoxel].x;
        if ( voxelMove[jVoxel].y < 0 ) voxelMove[jVoxel].y = dim.y + voxelMove[jVoxel].y;
        if ( voxelMove[jVoxel].z < 0 ) voxelMove[jVoxel].z = dim.z + voxelMove[jVoxel].z;
        if ( voxelMove[jVoxel].x >= dim.x ) voxelMove[jVoxel].x = voxelMove[jVoxel].x - dim.x;
        if ( voxelMove[jVoxel].y >= dim.y ) voxelMove[jVoxel].y = voxelMove[jVoxel].y - dim.y;
        if ( voxelMove[jVoxel].z >= dim.z ) voxelMove[jVoxel].z = voxelMove[jVoxel].z - dim.z;
    }

    clearBitMap(); // this set _voxelSet size to 0
    for ( int jVoxel=0; jVoxel<voxelMove.size(); jVoxel++)
        setBitMapValue(voxelMove[jVoxel],false);
    updateBitMapVoxelSet();  // don't change the color, which should be conserved
}

void bitMap::mirrorOverlay()
{
    if ( _data.empty )
    {
        qInfo() << "Empty overlay in mirror function: programming error??";
        return;
    }

    VoxelSet voxelMirror;
    voxelMirror.resize(_voxelSet.size());
    for ( int jVoxel=0; jVoxel<_voxelSet.size(); jVoxel++)
    {
        voxelMirror[jVoxel] = _voxelSet[jVoxel];
        voxelMirror[jVoxel].x = (nx()-1) - _voxelSet[jVoxel].x;
    }

    clearBitMap(); // this set _voxelSet size to 0
    for ( int jVoxel=0; jVoxel<voxelMirror.size(); jVoxel++)
        setBitMapValue(voxelMirror[jVoxel],false);
    updateBitMapVoxelSet();  // don't change the color, which should be conserved
}

void bitMap::setOverlayDimensions(const iPoint3D voxel)
{
    setOverlayDimensions(voxel.x, voxel.y, voxel.z);
}
void bitMap::setOverlayDimensions(const int nx, const int ny, const int nz)
{
    setDimensions(nx, ny, nz, 1);
    _voxelSet.resize(0);
}

void bitMap::defineBitMap(iPoint3D dim)
{ // a "bitmap" with value from 0. to 1.
    setOverlayDimensions(dim);
    _data.hdr.resolution.x = _data.hdr.resolution.y = _data.hdr.resolution.z = 1;
}
void bitMap::defineBitMap(iPoint3D dim, dPoint3D res)
{ // a "bitmap" with value from 0. to 1.
    setOverlayDimensions(dim);
    _data.hdr.resolution.x = res.x;
    _data.hdr.resolution.y = res.y;
    _data.hdr.resolution.z = res.z;
}
void bitMap::defineBitMap(iPoint3D dim, dMatrix ijkToXyz, dMatrix xyzToIjk)
{ // a "bitmap" with value from 0. to 1.
    setOverlayDimensions(dim);
    _data.hdr.resolution.x = qAbs(ijkToXyz[0][0]);
    _data.hdr.resolution.y = qAbs(ijkToXyz[1][1]);
    _data.hdr.resolution.z = qAbs(ijkToXyz[2][2]);
    _data.hdr.ijk_to_xyz = ijkToXyz;
    _data.hdr.xyz_to_ijk = xyzToIjk;
}
void bitMap::destroyBitMap()
{
    iPoint3D dim={0,0,0};
    dPoint3D res={0.,0.,0.};
    defineBitMap(dim,res);
    _voxelSet.clear();
}
void bitMap::clearBitMap()
{
    if ( nt() == 0 ) return;
    if ( nt() != 1 )
        qFatal("Error: in defineBitMap, # time points != 1:");
    else
    {
        CheckVoxelIndices(0,0,0,"clearBitMap");
        for (int j=0; j<_data.hdr.points.z; j++)
        {
//            setStackValue(j,0,0.);
            _data.timePoints[0].f1[j] = 0.;
            _data.color[j]   = Color_Undefined;
        }
    }
    _data.empty = true;
    _voxelSet.resize(0);
}

void bitMap::setAllBits()
{
    if ( nt() == 0 ) return;
    if ( nt() != 1 )
        qFatal("Error: in defineBitMap, # time points != 1:");
    else
    {
        CheckVoxelIndices(0,0,0,"setAllBits");
        for (int j=0; j<_data.hdr.points.z; j++)
        {
            _data.timePoints[0].f1[j] = 1.;
            _data.color[j]   = Color_Undefined;
        }
    }
    updateBitMapVoxelSet();
}

iPoint3D bitMap::getNearestNonzeroVoxel(const iPoint3D seedVoxel)
{ // starting at the seed voxel, cluster the overlay field
    // check that the seed voxel is set
    FUNC_ENTER << seedVoxel.x << seedVoxel.y << seedVoxel.z << "value" << isVoxelNonZero(seedVoxel);
    if ( isVoxelNonZero(seedVoxel) ) return seedVoxel;
    FUNC_INFO << "# voxels set" << getNumberVoxelsInOverlay();

    int expand=1;
    iPoint3D jRel; iPoint3D voxel;
    for (jRel.z=-expand; jRel.z<=expand; jRel.z++)
    {
        voxel.z = seedVoxel.z + jRel.z;
        for (jRel.y=-expand; jRel.y<=expand; jRel.y++)
        {
            voxel.y = seedVoxel.y + jRel.y;
            for (jRel.x=-expand; jRel.x<=expand; jRel.x++)
            {
                voxel.x = seedVoxel.x + jRel.x;
                // Ignore indices that are out of bounds.
                bool ignore = (voxel.x < 0 || voxel.x >= nx() )
                        ||    (voxel.y < 0 || voxel.y >= ny() )
                        ||    (voxel.z < 0 || voxel.z >= nz() );
                // Ignore center voxel.
                ignore |= (jRel.x==0 && jRel.y==0 && jRel.z==0);
                // Ignore voxels that have already been identified, or that are not set in the overlay
                if ( !ignore )
                {
                    if ( isVoxelNonZero(voxel) )
                         return voxel;
                }
            } // x
        } // y
    } // z
    return seedVoxel;  // don't change anything if not found
}

void bitMap::clusterOverlayField(const iPoint3D seedVoxel)
{ // starting at the seed voxel, cluster the overlay field
    FUNC_ENTER << seedVoxel.x << seedVoxel.y << seedVoxel.z;
    FUNC_INFO << "original" << getBitMapValue(seedVoxel);
    // Make a copy of the overlay
    bitMap originalOverlay;
    originalOverlay.copyTemplateDimensions(this,1);
    originalOverlay.clearBitMap();
    for (int jVoxel=0; jVoxel<_voxelSet.size(); jVoxel++)
        originalOverlay.setBitMapValue(_voxelSet[jVoxel],false);
    clearBitMap(); // Clear the original overlay

    // check that the seed voxel is set
    FUNC_INFO << "new" << originalOverlay.getBitMapValue(seedVoxel);
    if ( originalOverlay.getBitMapValue(seedVoxel) == 0. ) // the seed itself is not contained in the overlay
        return;
    else
        setBitMapValue(seedVoxel,1.,getBitMapColor(seedVoxel),true);

    int nvoxelsInField = 1;     int nvoxelsInFieldLast=0;   int expand=0;
    while ( nvoxelsInField > nvoxelsInFieldLast )
    {
        nvoxelsInFieldLast = nvoxelsInField;
        // Increment the expansion range.
        expand++;
        iPoint3D jRel; fVoxel voxel;
        for (jRel.z=-expand; jRel.z<=expand; jRel.z++)
        {
            voxel.z = seedVoxel.z + jRel.z;
            for (jRel.y=-expand; jRel.y<=expand; jRel.y++)
            {
                voxel.y = seedVoxel.y + jRel.y;
                for (jRel.x=-expand; jRel.x<=expand; jRel.x++)
                {
                    voxel.x = seedVoxel.x + jRel.x;
                    // Ignore indices that are out of bounds.
                    bool ignore = (voxel.x < 0 || voxel.x >= nx() )
                            ||    (voxel.y < 0 || voxel.y >= ny() )
                            ||    (voxel.z < 0 || voxel.z >= nz() );
                    // Ignore center voxel.
                    ignore |= (jRel.x==0 && jRel.y==0 && jRel.z==0);
                    // Ignore voxels that have already been identified, or that are not set in the overlay
                    if ( !ignore ) ignore |= getBitMapValue(voxel) != 0. || originalOverlay.getVolumeValue(voxel) == 0.;
                    if ( !ignore )
                    {
                        voxel.binaryValue = originalOverlay.getBitMapValue(voxel);
                        voxel.color = originalOverlay.getBitMapColor(voxel);
                        setBitMapVoxelIfTouchingConnectedNeighbor(voxel);
                    }
                } // x
            } // y
        } // z
        nvoxelsInField = getBitsSetInBitMap();
    }

    // At this point, the new overlay is binary. Go back and replace values (potentially non-binary using the copy)
    iPoint3D dim = getImageDimensions();  iPoint3D voxel;
    for (voxel.z=0; voxel.z<dim.z; voxel.z++)
    {
        for (voxel.y=0; voxel.y<dim.y; voxel.y++)
        {
            for (voxel.x=0; voxel.x<dim.x; voxel.x++)
                if ( getBitMapValue(voxel) != 0. ) setBitMapValue(voxel,
                                                                  originalOverlay.getBitMapValue(voxel),
                                                                  originalOverlay.getBitMapColor(voxel),
                                                                  true);
        }
    }
    FUNC_EXIT;
}

void bitMap::applyThresholds(bool reset, ImageIO *sourceFile, int iTime, double lowThresh, double highThresh)
{ // reset=true: clear the bitmap and start over; reset=false: only look at voxel previously set in overlay
    if ( reset ) clearBitMap();
    iPoint3D dim = getImageDimensions();
    iPoint3D voxel;
    for (voxel.z=0; voxel.z<dim.z; voxel.z++)
    {
        for (voxel.y=0; voxel.y<dim.y; voxel.y++)
        {
            for (voxel.x=0; voxel.x<dim.x; voxel.x++)
            {
                if ( reset || isVoxelNonZero(voxel) )
                {
                    double imageValue = sourceFile->getStackValue(voxel,iTime);
                    bool inRange = (imageValue >= lowThresh) && (imageValue <= highThresh);
                    if ( inRange )
                        setBitMapValue(voxel,1.,Color_Undefined,false);
                    else
                        setBitMapValue(voxel,0.,Color_Undefined,false); // reset to 0 in case it was previously set
                }
            }
        }
    }
    updateBitMapVoxelSet(Color_Undefined);
}

void bitMap::growOverlayField(const iPoint3D seedVoxel, ImageIO *colorFile, double threshold)
{ // starting at the seed voxel, grow the overlay field using values above threshold
    FUNC_INFO;

    // make a copy of the current bitmap and then clear this
    bitMap priorBitMap;
    priorBitMap.copyTemplateVolume(this,0);
    priorBitMap._voxelSet = _voxelSet;
    clearBitMap();
    // check that the seed voxel is above theshold; also check sign
    double seedValue = colorFile->getVolumeValue(seedVoxel);
    double sign;
    if ( seedValue > threshold )
        sign = 1.;
    else if ( -seedValue > threshold )
        sign = -1.;
    else  // the seed itself is not above the threshold
        return;

    setBitMapValue(seedVoxel,1.,getBitMapColor(seedVoxel),true);

    int nvoxelsInField = 1;     int nvoxelsInFieldLast=0;   int expand=0;
    while ( nvoxelsInField > nvoxelsInFieldLast )
    {
        nvoxelsInFieldLast = nvoxelsInField;
        // Increment the expansion range.
        expand++;
        iPoint3D jRel; fVoxel voxel;
        for (jRel.z=-expand; jRel.z<=expand; jRel.z++)
        {
            voxel.z = seedVoxel.z + jRel.z;
            for (jRel.y=-expand; jRel.y<=expand; jRel.y++)
            {
                voxel.y = seedVoxel.y + jRel.y;
                for (jRel.x=-expand; jRel.x<=expand; jRel.x++)
                {
                    voxel.x = seedVoxel.x + jRel.x;
                    // Ignore indices that are out of bounds.
                    bool ignore = (voxel.x < 0 || voxel.x >= nx() )
                            ||    (voxel.y < 0 || voxel.y >= ny() )
                            ||    (voxel.z < 0 || voxel.z >= nz() );
                    // Ignore center voxel.
                    ignore |= (jRel.x==0 && jRel.y==0 && jRel.z==0);
                    // Ignore voxels that have already been identified, or that are not above threshold
                    if ( !ignore ) ignore |= getBitMapValue(voxel) != 0. || sign * colorFile->getVolumeValue(voxel) < threshold;
                    if ( !ignore )
                    {
                        voxel.binaryValue = 1.;
                        voxel.color = Color_Undefined;
                        setBitMapVoxelIfTouchingConnectedNeighbor(voxel);
                    }
                } // x
            } // y
        } // z
        nvoxelsInField = getBitsSetInBitMap();
    }

    // add the prior bitmap
    for (int jVoxel=0; jVoxel<priorBitMap.getNumberVoxelsInOverlay(); jVoxel++)
    {
        fVoxel voxel = priorBitMap.getBitMapVoxel(jVoxel);
        setBitMapValue(voxel,false);
    }
    updateBitMapVoxelSet();
    FUNC_EXIT;
}

void bitMap::setBitMapVoxelIfTouchingConnectedNeighbor(const fVoxel testVoxel)
{
    iPoint3D jRel;
    for (jRel.z=-1; jRel.z<=1; jRel.z++)
    {
        iPoint3D neighbor;
        neighbor.z = testVoxel.z + jRel.z;
        for (jRel.y=-1; jRel.y<=1; jRel.y++)
        {
            neighbor.y = testVoxel.y + jRel.y;
            for (jRel.x=-1; jRel.x<=1; jRel.x++)
            {
                neighbor.x = testVoxel.x + jRel.x;
                // Ignore center voxel.
                //bool ignore = (jRel.x==0 && jRel.y==0 && jRel.z==0);
                // Ignore indices that are out of bounds.
                bool ignore = (neighbor.x < 0 || neighbor.x >= nx() )
                        ||    (neighbor.y < 0 || neighbor.y >= ny() )
                        ||    (neighbor.z < 0 || neighbor.z >= nz() );
                if ( !ignore )
                {
                    // Has the neighbor voxel been identified as selected previously?
                    bool touch = getVolumeValue(neighbor) != 0.;
                    bool corner = (jRel.x * jRel.y * jRel.z != 0);
                    if ( touch && !corner )
                        setBitMapValue(testVoxel,true);
                } // !ignore
            } // x
        } // y
    } // z
}

double bitMap::getMedianSignal(ImageIO *file, bool aboveThreshold, double thresholdRatio)
{ // a threshold of zero circumvents value tests
    dVector signalVector;
    if ( thresholdRatio == 0. )
    {
        for ( int jVoxel=0; jVoxel<_voxelSet.size(); jVoxel++)
            signalVector.append(file->getVolumeValue(_voxelSet[jVoxel]));
    }
    else
    {
        for ( int jVoxel=0; jVoxel<_voxelSet.size(); jVoxel++)
        {
            double signal = _voxelSet[jVoxel].binaryValue;
            bool ok = signal > thresholdRatio * thresholdRatio;
            if ( !aboveThreshold ) ok = signal < thresholdRatio * thresholdRatio;
            if ( ok ) signalVector.append(file->getVolumeValue(_voxelSet[jVoxel]));
        }
    }
    double medianSignal = utilMath::medianValue(signalVector);
    return medianSignal;
}

dPoint3D ImageIO::computeCOM(int iVolume)
{
    dPoint3D origin = computeCOMIndex(iVolume, false);
    dPoint3D res = getImageResolution();
    origin.x *= res.x;
    origin.y *= res.y;
    origin.z *= res.z;
    return origin;
}

dPoint3D ImageIO::computeCOMIndex(int iVolume)
{
    return computeCOMIndex(iVolume, false);
}
dPoint3D ImageIO::computeCOMIndex(int iVolume, bool inverse)
{ // center of mass ( weight ~ signal^2)
    dPoint3D COM={0.,0.,0.};
    double amp2Sum=0.;
    iPoint3D voxel;
    for (voxel.z=0; voxel.z<nz(); voxel.z++)
    {
        for (voxel.y=0; voxel.y<ny(); voxel.y++)
        {
            for (voxel.x=0; voxel.x<nx(); voxel.x++)
            {
                double amp2 = SQR(getStackValue(voxel,iVolume));
                if ( inverse )
                    if ( amp2 != 0. ) amp2 = 1./amp2;
                if  ( amp2 < 0. || qIsNaN(amp2) || qIsInf(amp2 ) ) amp2 = 0.;
                amp2Sum += amp2;
                COM.x += amp2 * voxel.x;
                COM.y += amp2 * voxel.y;
                COM.z += amp2 * voxel.z;
            } // x
        } // y
    } // z
    if ( amp2Sum != 0. )
    {
        COM.x /= amp2Sum;
        COM.y /= amp2Sum;
        COM.z /= amp2Sum;
    }
    return COM;
}

dPoint3D ImageIO::computeCOMStdev(int iVolume, bool inverse, dPoint3D COM)
{
    dPoint3D COMStdev={0.,0.,0.};
    double ampSum=0.;
    iPoint3D voxel;
    for (voxel.z=0; voxel.z<nz(); voxel.z++)
    {
        for (voxel.y=0; voxel.y<ny(); voxel.y++)
        {
            for (voxel.x=0; voxel.x<nx(); voxel.x++)
            {
                // Convert from data indices --> space in the registered system.
                double amp = getStackValue(voxel,iVolume);  // typical case
                if ( inverse )
                    if ( amp != 0. ) amp = 1./amp;
                ampSum += amp;
                COMStdev.x += amp * SQR(voxel.x-COM.x);
                COMStdev.y += amp * SQR(voxel.y-COM.y);
                COMStdev.z += amp * SQR(voxel.z-COM.z);
            } // x
        } // y
    } // z
    if ( ampSum != 0. )
    {
        COMStdev.x = qSqrt(COMStdev.x/ampSum);
        COMStdev.y = qSqrt(COMStdev.y/ampSum);
        COMStdev.z = qSqrt(COMStdev.z/ampSum);
    }
    return COMStdev;
}

iPoint3D bitMap::computeCOMIndex()
{
    iPoint3D iCentroid={0,0,0};
    if ( !someBitsAreSet() ) return iCentroid;

    double weight = 0.;
    dPoint3D dCentroid={0.,0.,0.};
    for ( int jVoxel=0; jVoxel<_voxelSet.size(); jVoxel++)
    {
        fVoxel voxel = _voxelSet[jVoxel];
        dCentroid.x += static_cast<double>(voxel.x * voxel.binaryValue);
        dCentroid.y += static_cast<double>(voxel.y * voxel.binaryValue);
        dCentroid.z += static_cast<double>(voxel.z * voxel.binaryValue);
        weight += static_cast<double>(voxel.binaryValue);
    } // jVoxel
    dCentroid.x /= weight;
    dCentroid.y /= weight;
    dCentroid.z /= weight;
    // return ints
    iCentroid.x = static_cast<int>(dCentroid.x);
    iCentroid.y = static_cast<int>(dCentroid.y);
    iCentroid.z = static_cast<int>(dCentroid.z);
    FUNC_EXIT << "COM" << iCentroid.x << iCentroid.y << iCentroid.z;
    return iCentroid;
}

d3Vector bitMap::computeCOMs()
{ // return 3D COM for each slice in voxel units; if 3rd dimension < 0, no COM on slice
    d3Vector COM; COM.fill({0.,0.,0.},nz());
    iPoint3D voxel;
    for (voxel.z=0; voxel.z<nz(); voxel.z++)
    {
        COM[voxel.z].x = 0.;
        COM[voxel.z].y = 0.;
        COM[voxel.z].z = -1.;
        double ampSum=0.;
        for (voxel.y=0; voxel.y<ny(); voxel.y++)
        {
            for (voxel.x=0; voxel.x<nx(); voxel.x++)
            {
                // Convert from data indices --> space in the registered system.
                double amp = getBitMapValue(voxel);
                ampSum += amp;
                COM[voxel.z].x += amp * voxel.x;
                COM[voxel.z].y += amp * voxel.y;
            } // x
        } // y
        if ( ampSum != 0. )
        {
            COM[voxel.z].x /= ampSum;
            COM[voxel.z].y /= ampSum;
            COM[voxel.z].z = voxel.z;
        }
    } // z
    return COM;
}

void bitMap::dilateOverlayUsingThreshold(ImageIO *file, double thresholdRatio,
                                         bool aboveThreshold, bool twoDimensional, bool fillSurroundedVoxels)
{
    FUNC_ENTER;
    if ( !someBitsAreSet() ) return;

    double medianSignal;
    if ( fillSurroundedVoxels )
        medianSignal = getMedianSignal(file, aboveThreshold, thresholdRatio);
    else
        medianSignal = getMedianSignal(file, aboveThreshold, 0.);  // small speed increase if it doesn't have to check threshold
//    medianSignal = 0.15;
    FUNC_INFO << "median =" << medianSignal;

    iPoint3D nRel;
    nRel.x = nRel.y = nRel.z = 1;
    if ( twoDimensional ) nRel.z = 0;

    // make a copy of the current bitmap and then clear this
    bitMap priorBitMap;
    priorBitMap.copyTemplateVolume(this,0);
    priorBitMap._voxelSet = _voxelSet;
    clearBitMap();

    iPoint2D rangeX = {{nx()},{-1}};
    iPoint2D rangeY = {{ny()},{-1}};
    iPoint2D rangeZ = {{nz()},{-1}};

    int nVoxels = priorBitMap.getBitsSetInBitMap();
    for ( int jVoxel=0; jVoxel<nVoxels; jVoxel++)
    {
        fVoxel voxel = priorBitMap.getBitMapVoxel(jVoxel);
        fVoxel neighbor, rel;
        // go through the local neighborhood
        for (rel.z=-nRel.z; rel.z<=nRel.z; rel.z++)
        {
            neighbor.z = voxel.z + rel.z;
            for (rel.y=-nRel.y; rel.y<=nRel.y; rel.y++)
            {
                neighbor.y = voxel.y + rel.y;
                for (rel.x=-nRel.x; rel.x<=nRel.x; rel.x++)
                {
                    neighbor.x = voxel.x + rel.x;
                    // Ignore indices that are out of bounds.
                    bool ignore = (neighbor.x < 0 || neighbor.x >= nx() )
                            ||    (neighbor.y < 0 || neighbor.y >= ny() )
                            ||    (neighbor.z < 0 || neighbor.z >= nz() );
                    if ( !ignore )
                    {
                        double neighborSignal = file->getVolumeValue(neighbor);
                        bool ok = neighborSignal > medianSignal / thresholdRatio;
                        if ( !aboveThreshold ) ok = neighborSignal < medianSignal * thresholdRatio;
                        if ( !ok && fillSurroundedVoxels ) ok |= voxelSurrounded(neighbor);
                        // FUNC_INFO << "neighborSignal, threshold, ok" << neighborSignal << medianSignal * thresholdRatio << ok;
                        if ( ok )
                        {
                            neighbor.binaryValue = 1.;  neighbor.color = voxel.color;
                            setBitMapValue(neighbor,false);
                            rangeX.lower = qMin(rangeX.lower,neighbor.x);
                            rangeY.lower = qMin(rangeY.lower,neighbor.y);
                            rangeZ.lower = qMin(rangeZ.lower,neighbor.z);
                            rangeX.upper = qMax(rangeX.upper,neighbor.x);
                            rangeY.upper = qMax(rangeY.upper,neighbor.y);
                            rangeZ.upper = qMax(rangeZ.upper,neighbor.z);
                        }
                    }
                } // jRel.x
            } // jRel.y
        } // jRel.z
    } // jVoxel
    updateBitMapVoxelSet(rangeX, rangeY, rangeZ);
//    updateBitMapVoxelSet();
    FUNC_EXIT << getBitsSetInBitMap();
}

void bitMap::fillSurroundedVoxels()
{
    FUNC_ENTER << getBitsSetInBitMap();
    iPoint2D rangeX = {{nx()},{-1}};
    iPoint2D rangeY = {{ny()},{-1}};
    iPoint2D rangeZ = {{nz()},{-1}};

    int nVoxels = getBitsSetInBitMap();
    for ( int jVoxel=0; jVoxel<nVoxels; jVoxel++)
    {
        iPoint3D voxel = getBitMapVoxel3D(jVoxel);
        rangeX.lower = qMin(rangeX.lower,voxel.x);
        rangeY.lower = qMin(rangeY.lower,voxel.y);
        rangeZ.lower = qMin(rangeZ.lower,voxel.z);
        rangeX.upper = qMax(rangeX.upper,voxel.x);
        rangeY.upper = qMax(rangeY.upper,voxel.y);
        rangeZ.upper = qMax(rangeZ.upper,voxel.z);
    } // jVoxel

    int nVoxelsNew = 1;
    iPoint3D voxel;
    while ( nVoxelsNew > 0 )
    {
        nVoxelsNew = 0;
        for (voxel.z=rangeZ.lower; voxel.z<=rangeZ.upper; voxel.z++)
        {
            for (voxel.y=rangeY.lower; voxel.y<=rangeY.upper; voxel.y++)
            {
                for (voxel.x=rangeX.lower; voxel.x<=rangeX.upper; voxel.x++)
                {
                    if ( isVoxelZero(voxel) && voxelSurrounded(voxel) )
                    {
                        setBitMapValue(voxel,1.,false);
                        nVoxelsNew++;
                    }
                } // x
            } // y
        } // z
        updateBitMapVoxelSet(rangeX, rangeY, rangeZ);
    } // still increasing
    FUNC_EXIT << getBitsSetInBitMap();
}

bool bitMap::voxelSurrounded(iPoint3D voxel)
{
    iPoint3D dim = getDimensions();
    bool edge = (voxel.x == 0) || (voxel.x == dim.x-1);
    edge     |= (voxel.y == 0) || (voxel.y == dim.y-1);
    edge     |= (voxel.z == 0) || (voxel.z == dim.z-1);
    if ( edge ) return false;
    iPoint3D testVoxel;

    // +x
    bool edgeXPlus  = true;  bool edgeXMinus = true;
    testVoxel = voxel;
    for (testVoxel.x=voxel.x+1; testVoxel.x<dim.x; testVoxel.x++)
    {
        // it is an edge voxel if all voxels are zero to the edge
        edgeXPlus &= isVoxelZero(testVoxel);
        if ( !edgeXPlus ) break;
    }
    if ( !edgeXPlus )
    {
        // -x
        for (testVoxel.x=voxel.x-1; testVoxel.x>=0; testVoxel.x--)
        {
            // it is an edge voxel if all voxels are zero to the edge
            edgeXMinus &= isVoxelZero(testVoxel);
            if ( !edgeXMinus ) break;
        }
    }
    bool edgeX = edgeXPlus | edgeXMinus;

    // +y
    testVoxel = voxel;
    bool edgeYPlus = true;  bool edgeYMinus = true;
    for (testVoxel.y=voxel.y+1; testVoxel.y<dim.y; testVoxel.y++)
    {
        // it is an edge voxel if all voxels are zero to the edge
        edgeYPlus &= isVoxelZero(testVoxel);
        if ( !edgeYPlus ) break;
    }
    if ( !edgeYPlus )
    {
        // -y
        for (testVoxel.y=voxel.y-1; testVoxel.y>=0; testVoxel.y--)
        {
            // it is an edge voxel if all voxels are zero to the edge
            edgeYMinus &= isVoxelZero(testVoxel);
            if ( !edgeYMinus ) break;
        }
    }
    bool edgeY = edgeYPlus | edgeYMinus;

    // +z
    testVoxel = voxel;
    bool edgeZPlus = true;  bool edgeZMinus = true;
    for (testVoxel.z=voxel.z+1; testVoxel.z<dim.z; testVoxel.z++)
    {
        // it is an edge voxel if all voxels are zero to the edge
        edgeZPlus &= isVoxelZero(testVoxel);
        if ( !edgeZPlus ) break;
    }
    if ( !edgeZPlus )
    {
        // -x
        for (testVoxel.z=voxel.z-1; testVoxel.z>=0; testVoxel.z--)
        {
            // it is an edge voxel if all voxels are zero to the edge
            edgeZMinus &= isVoxelZero(testVoxel);
            if ( !edgeZMinus ) break;
        }
    }
    bool edgeZ = edgeZPlus | edgeZMinus;

    bool surroundedXY  = !(edgeX || edgeY);
    bool surroundedXZ  = !(edgeX || edgeZ);
    bool surroundedYZ  = !(edgeY || edgeZ);
    return surroundedXY || surroundedXZ || surroundedYZ;
//    bool surroundedXYZ = !(edgeX || edgeY || edgeZ);
//    return surroundedXYZ;
}

void bitMap::medianFilterOverlay(bool twoDimensional, int iHalfWidth)
{
    if ( !someBitsAreSet() ) return;

    bitMap priorBitMap;
    priorBitMap.copyTemplateVolume(this,0);
    priorBitMap._voxelSet = _voxelSet;
    clearBitMap();

    // This could be modified to make a 2D alternative
    iPoint3D nRel;
    nRel.x = nRel.y = nRel.z = iHalfWidth;
    if ( twoDimensional ) nRel.z = 0;

    iPoint2D rangeX = {{nx()},{-1}};
    iPoint2D rangeY = {{ny()},{-1}};
    iPoint2D rangeZ = {{nz()},{-1}};

    int nVoxels = priorBitMap.getBitsSetInBitMap();
    for ( int jVoxel=0; jVoxel<nVoxels; jVoxel++)
    {
        fVoxel voxel = priorBitMap.getBitMapVoxel(jVoxel);
        fVoxel neighbor, rel;  neighbor.binaryValue = voxel.binaryValue;  // copy value for "contains" below
        // go through the local neighborhood
        dVector neighborVector;
        for (rel.z=-nRel.z; rel.z<=nRel.z; rel.z++)
        {
            neighbor.z = voxel.z + rel.z;
            for (rel.y=-nRel.y; rel.y<=nRel.y; rel.y++)
            {
                neighbor.y = voxel.y + rel.y;
                for (rel.x=-nRel.x; rel.x<=nRel.x; rel.x++)
                {
                    neighbor.x = voxel.x + rel.x;
                    // Ignore indices that are out of bounds.
                    if ( voxelInBounds(neighbor) )
                    {
                        double neighborValue=priorBitMap.getBitMapValue(neighbor);
                        rangeX.lower = qMin(rangeX.lower,neighbor.x);
                        rangeY.lower = qMin(rangeY.lower,neighbor.y);
                        rangeZ.lower = qMin(rangeZ.lower,neighbor.z);
                        rangeX.upper = qMax(rangeX.upper,neighbor.x);
                        rangeY.upper = qMax(rangeY.upper,neighbor.y);
                        rangeZ.upper = qMax(rangeZ.upper,neighbor.z);
                        neighborVector.append(neighborValue);
                    }
                } // jRel.x
            } // jRel.y
        } // jRel.z
        double median = utilMath::medianValue(neighborVector);
        // Assign the max to the local neighborhood
        // go through the local neighborhood
        for (rel.z=-nRel.z; rel.z<=nRel.z; rel.z++)
        {
            neighbor.z = voxel.z + rel.z;
            for (rel.y=-nRel.y; rel.y<=nRel.y; rel.y++)
            {
                neighbor.y = voxel.y + rel.y;
                for (rel.x=-nRel.x; rel.x<=nRel.x; rel.x++)
                {
                    neighbor.x = voxel.x + rel.x;
                    if ( voxelInBounds(neighbor) )
                    {
                        neighbor.binaryValue = median;
                        neighbor.color = voxel.color;
                        setBitMapValue(neighbor,false);
                    }
                } // jRel.x
            } // jRel.y
        } // jRel.z
    } // jVoxel
    updateBitMapVoxelSet(rangeX, rangeY, rangeZ);
    FUNC_EXIT << someBitsAreSet();
//    updateBitMapVoxelSet();
}
void bitMap::dilateOverlay(bool twoDimensional)
{
    if ( !someBitsAreSet() ) return;

    bitMap priorBitMap;
    priorBitMap.copyTemplateVolume(this,0);
    priorBitMap._voxelSet = _voxelSet;
    clearBitMap();

    // This could be modified to make a 2D alternative
    iPoint3D nRel;
    nRel.x = nRel.y = nRel.z = 1;
    if ( twoDimensional ) nRel.z = 0;

    iPoint2D rangeX = {{nx()},{-1}};
    iPoint2D rangeY = {{ny()},{-1}};
    iPoint2D rangeZ = {{nz()},{-1}};

    int nVoxels = priorBitMap.getBitsSetInBitMap();
    for ( int jVoxel=0; jVoxel<nVoxels; jVoxel++)
    {
        fVoxel voxel = priorBitMap.getBitMapVoxel(jVoxel);
        double fMax = -1.e10;
        fVoxel neighbor, rel;  neighbor.binaryValue = voxel.binaryValue;  // copy value for "contains" below
        // go through the local neighborhood
        for (rel.z=-nRel.z; rel.z<=nRel.z; rel.z++)
        {
            neighbor.z = voxel.z + rel.z;
            for (rel.y=-nRel.y; rel.y<=nRel.y; rel.y++)
            {
                neighbor.y = voxel.y + rel.y;
                for (rel.x=-nRel.x; rel.x<=nRel.x; rel.x++)
                {
                    neighbor.x = voxel.x + rel.x;
                    // Ignore indices that are out of bounds.
                    if ( voxelInBounds(neighbor) )
                    {
                        double neighborValue=priorBitMap.getBitMapValue(neighbor);
                        if ( neighborValue > fMax ) fMax = neighborValue;
                        rangeX.lower = qMin(rangeX.lower,neighbor.x);
                        rangeY.lower = qMin(rangeY.lower,neighbor.y);
                        rangeZ.lower = qMin(rangeZ.lower,neighbor.z);
                        rangeX.upper = qMax(rangeX.upper,neighbor.x);
                        rangeY.upper = qMax(rangeY.upper,neighbor.y);
                        rangeZ.upper = qMax(rangeZ.upper,neighbor.z);
                    }
                } // jRel.x
            } // jRel.y
        } // jRel.z
        // Assign the max to the local neighborhood
        // go through the local neighborhood
        for (rel.z=-nRel.z; rel.z<=nRel.z; rel.z++)
        {
            neighbor.z = voxel.z + rel.z;
            for (rel.y=-nRel.y; rel.y<=nRel.y; rel.y++)
            {
                neighbor.y = voxel.y + rel.y;
                for (rel.x=-nRel.x; rel.x<=nRel.x; rel.x++)
                {
                    neighbor.x = voxel.x + rel.x;
                    if ( voxelInBounds(neighbor) )
                    {
                        neighbor.binaryValue = fMax;
                        neighbor.color = voxel.color;
                        setBitMapValue(neighbor,false);
                    }
                } // jRel.x
            } // jRel.y
        } // jRel.z
    } // jVoxel
    updateBitMapVoxelSet(rangeX, rangeY, rangeZ);
    FUNC_EXIT << someBitsAreSet();
//    updateBitMapVoxelSet();
}
void bitMap::erodeOverlay(bool medianFilter, bool twoDimensional)
{
    if ( !someBitsAreSet() ) return;

    bitMap priorBitMap;
    priorBitMap.copyTemplateVolume(this,0);
    priorBitMap._voxelSet = _voxelSet;
    clearBitMap();

    // This could be modified to make a 2D alternative
    iPoint3D nRel;
    nRel.x = nRel.y = nRel.z = 1;
    if ( twoDimensional ) nRel.z = 0;

    iPoint2D rangeX = {{nx()},{-1}};
    iPoint2D rangeY = {{ny()},{-1}};
    iPoint2D rangeZ = {{nz()},{-1}};

    int nVoxels = priorBitMap.getBitsSetInBitMap();
    for ( int jVoxel=0; jVoxel<nVoxels; jVoxel++)
    {
        fVoxel voxel = priorBitMap.getBitMapVoxel(jVoxel);
        double fMin = 1.e10;
        int countNeighbors=0;
        double sumNeighbors=0.;
        fVoxel neighbor, rel;  neighbor.binaryValue = voxel.binaryValue;  // copy value for "contains" below
        // go through the local neighborhood
        for (rel.z=-nRel.z; rel.z<=nRel.z; rel.z++)
        {
            neighbor.z = voxel.z + rel.z;
            for (rel.y=-nRel.y; rel.y<=nRel.y; rel.y++)
            {
                neighbor.y = voxel.y + rel.y;
                for (rel.x=-nRel.x; rel.x<=nRel.x; rel.x++)
                {
                    neighbor.x = voxel.x + rel.x;
                    if ( voxelInBounds(neighbor) )
                    {
                        double neighborValue=priorBitMap.getBitMapValue(neighbor);
                        if ( neighborValue < fMin ) fMin = neighborValue;
                        countNeighbors++;
                        sumNeighbors += neighborValue;
                        rangeX.lower = qMin(rangeX.lower,neighbor.x);
                        rangeY.lower = qMin(rangeY.lower,neighbor.y);
                        rangeZ.lower = qMin(rangeZ.lower,neighbor.z);
                        rangeX.upper = qMax(rangeX.upper,neighbor.x);
                        rangeY.upper = qMax(rangeY.upper,neighbor.y);
                        rangeZ.upper = qMax(rangeZ.upper,neighbor.z);
                    }
                } // jRel.x
            } // jRel.y
        } // jRel.z
        if ( medianFilter )
        {
            voxel.binaryValue = 1.;
            if ( countNeighbors != 0 )
            {
                double testValue = sumNeighbors / static_cast<double>(countNeighbors);
                if ( testValue <= 0.5 ) voxel.binaryValue = 0.;  // value could be sliding to use variable degrees of erosion
            }
        }
        else
            voxel.binaryValue = fMin;
        setBitMapValue(voxel,false);
    } // jVoxel
    updateBitMapVoxelSet(rangeX, rangeY, rangeZ);
//    updateBitMapVoxelSet();
}

void bitMap::outerShell(bitMap source)
{ // this only works if the source overlay is solid; create a shell equal to the current boundary
    FUNC_ENTER << source.getBitsSetInBitMap();
    *this = source;
    FUNC_INFO << "this" << getBitsSetInBitMap();
    erodeOverlay();          // smaller than source
    FUNC_INFO << "bits set" << getBitsSetInBitMap() << source.getBitsSetInBitMap();
    // do a subtraction: source - this  ; set value=qMax(0,value)
    for ( int jVoxel=0; jVoxel<source.getBitsSetInBitMap(); jVoxel++)  // sum over larger bitmap
    {
        fVoxel voxel       = source.getBitMapVoxel(jVoxel);
        double value       = source.getBitMapValue(voxel);
        double valueErode  = getBitMapValue(voxel);
        FUNC_INFO << jVoxel << value << valueErode;
        value -= valueErode;
        voxel.binaryValue = qMax(0.,value);  // but it should be always > 0
        setBitMapValue(voxel,false);
    }
    updateBitMapVoxelSet();
    FUNC_EXIT << getNumberVoxelsInOverlay();
}

bitMap bitMap::findPerimeter()
{
    // a voxel is a perimeter voxel if it is nonzero but connected to a zero
    bitMap perimeter = *this;
    perimeter.clearBitMap();
    iPoint3D dim = getImageDimensions();
    iPoint3D voxel;
    for (voxel.z=0; voxel.z<dim.z; voxel.z++)
    {
        for (voxel.y=0; voxel.y<dim.y; voxel.y++)
        {
            for (voxel.x=0; voxel.x<dim.x; voxel.x++)
            {
                /*
                if ( isVoxelUnity(voxel) )
                {
                    overlayColor color = getBitMapColor(voxel);
                    // Four cases: +-x and +-y
                    iPoint3D neighbor = voxel;
                    bool edge = false;
                    neighbor.x += 1;
                    if (neighbor.x < dim.x) edge |= !isVoxelUnity(neighbor);
                    neighbor.x -= 2;
                    if ( neighbor.x >= 0 )  edge |= !isVoxelUnity(neighbor);
                    neighbor = voxel;
                    neighbor.y += 1;
                    if (neighbor.y < dim.y) edge |= !isVoxelUnity(neighbor);
                    neighbor.y -= 2;
                    if ( neighbor.y >= 0 )  edge |= !isVoxelUnity(neighbor);
                    if ( edge )
                        perimeter.setBitMapValue(voxel,1.,color,false);
                    else
                        perimeter.setBitMapValue(voxel,0.,color,false);
                }
            }
            */
                if ( isVoxelNonZero(voxel) )
                {
                    overlayColor color = getBitMapColor(voxel);
                    // Four cases: +-x and +-y
                    iPoint3D neighbor = voxel;
                    bool edge = false;
                    neighbor.x += 1;
                    if (neighbor.x < dim.x) edge |= isVoxelZero(neighbor);
                    neighbor.x -= 2;
                    if ( neighbor.x >= 0 )  edge |= isVoxelZero(neighbor);
                    neighbor = voxel;
                    neighbor.y += 1;
                    if (neighbor.y < dim.y) edge |= isVoxelZero(neighbor);
                    neighbor.y -= 2;
                    if ( neighbor.y >= 0 )  edge |= isVoxelZero(neighbor);
                    if ( edge )
                        perimeter.setBitMapValue(voxel,1.,color,false);
                    else
                        perimeter.setBitMapValue(voxel,0.,color,false);
                }
            }
        }
    }
    perimeter.updateBitMapVoxelSet();
    return perimeter;
}

bitMap bitMap::constraintField(double cutOffVoxels)
{
    iPoint3D dim = getDimensions();
    dPoint3D res = getResolution();
    bitMap constraintField;
    constraintField.defineBitMap(dim,res);
    dPoint3D distanceCutoff = {cutOffVoxels*res.x,cutOffVoxels*res.y,cutOffVoxels*res.z};
    for (int j=0; j<voxelsPerVolume(); j++)
    {
        iPoint3D voxel = utilIO::index3d(dim,j);
        if ( isVoxelZero(voxel) )
            constraintField.setBitMapValue(voxel,0.,false);
        else // if ( voxel.x == 14 && voxel.y == 14 && voxel.z == 14 )
        {
            dPoint3D minDist={1.e10,1.e10,1.e10};
            iPoint3D neighbor = voxel;
            // -x
            for (neighbor.x=voxel.x; neighbor.x>=0; neighbor.x--)
                if ( isVoxelZero(neighbor) ) break;
            neighbor.x = qMax(0,neighbor.x);
            double dist = qAbs(voxel.x-neighbor.x)*res.x;
            if ( dist < minDist.x ) minDist.x = dist;
            // +x
            for (neighbor.x=voxel.x; neighbor.x<dim.x; neighbor.x++)
                if ( isVoxelZero(neighbor) ) break;
            neighbor.x = qMin(dim.x-1,neighbor.x);
            dist = qAbs(voxel.x-neighbor.x)*res.x;
            if ( dist < minDist.x ) minDist.x = dist;

            // -y
            neighbor = voxel;
            for (neighbor.y=voxel.y; neighbor.y>=0; neighbor.y--)
                if ( isVoxelZero(neighbor) ) break;
            neighbor.y = qMax(0,neighbor.y);
            dist = qAbs(voxel.y-neighbor.y)*res.y;
            if ( dist < minDist.y ) minDist.y = dist;
            // +y
            for (neighbor.y=voxel.y; neighbor.y<dim.y; neighbor.y++)
                if ( isVoxelZero(neighbor) ) break;
            neighbor.y = qMin(dim.y-1,neighbor.y);
            dist = qAbs(voxel.y-neighbor.y)*res.y;
            if ( dist < minDist.y ) minDist.y = dist;

            // -z
            neighbor = voxel;
            for (neighbor.z=voxel.z; neighbor.z>=0; neighbor.z--)
                if ( isVoxelZero(neighbor) ) break;
            neighbor.z = qMax(0,neighbor.z);
            dist = qAbs(voxel.z-neighbor.z)*res.z;
            if ( dist < minDist.z ) minDist.z = dist;
            // +z
            for (neighbor.z=voxel.z; neighbor.z<dim.z; neighbor.z++)
                if ( isVoxelZero(neighbor) ) break;
            neighbor.z = qMin(dim.z-1,neighbor.z);
            dist = qAbs(voxel.z-neighbor.z)*res.z;
            if ( dist < minDist.z ) minDist.z = dist;

            dPoint3D constraint = {1.,1.,1.};
            if ( minDist.x < distanceCutoff.x ) constraint.x = minDist.x / distanceCutoff.x;
            if ( minDist.y < distanceCutoff.y ) constraint.y = minDist.y / distanceCutoff.y;
            if ( minDist.z < distanceCutoff.z ) constraint.z = minDist.z / distanceCutoff.z;
            double constraintValue = qSqrt(constraint.x) * qSqrt(constraint.y) * qSqrt(constraint.z);
            constraintField.setBitMapValue(voxel,constraintValue,false);
        }
    }
    constraintField.updateBitMapVoxelSet();
    return constraintField;
}

void bitMap::selectOverlayFromThresholdGrayScale(ImageIO *signalFile, bool above, double value, bool binary)
{
    bool existingOverlay = _voxelSet.size() != 0;
    if  ( existingOverlay )
    {   // only to loop over existing voxels
        for (int jVoxel=0; jVoxel<_voxelSet.size(); jVoxel++)
        {   // Trim values of the existing overlay that are not above or below the threshold
            fVoxel voxel = _voxelSet[jVoxel];
            double signal = signalFile->getVolumeValue(voxel);
            if ( (above && signal < value) || (!above && signal > value) )
            {
                voxel.binaryValue = 0.;
                setBitMapValue(voxel,false);
            }
        }
    }
    else
    { // there is no existing overlay, just apply the threshold to everything
        clearBitMap();
        iPoint4D dim = signalFile->getDimensions();
        iPoint3D voxel;
        for ( voxel.z=0; voxel.z<dim.z; voxel.z++ )
        {
            for ( voxel.y=0; voxel.y<dim.y; voxel.y++ )
            {
                for ( voxel.x=0; voxel.x<dim.x; voxel.x++ )
                {
                    double signal = signalFile->getVolumeValue(voxel);
                    if ( (above && signal >= value) || (!above && signal <= value) )
                    {
                        double value = 1.;
                        if ( !binary ) value = signal;
                        setBitMapValue(voxel,value,Color_Undefined,true);
                    }
                }
            }
        }
    }
    updateBitMapVoxelSet();  // don't change the color, which should be carried by the input voxelSet
}

void bitMap::selectOverlayFromThresholdColorScale(ImageIO *signalFile, bool above, double value)
{
    bool existingOverlay = _voxelSet.size() != 0;
    if  ( existingOverlay )
    {   // only to loop over existing voxels
        for (int jVoxel=0; jVoxel<_voxelSet.size(); jVoxel++)
        {   // Trim values of the existing overlay that are not above or below the threshold
            fVoxel voxel = _voxelSet[jVoxel];
            double signal = signalFile->getVolumeValue(voxel);
            if ( (above && signal < value) || (!above && signal > value) )
            {
                voxel.binaryValue = 0.;
                setBitMapValue(voxel,false);
            }
        }
    }
    else
    { // there is no existing overlay, just apply the threshold to everything
        clearBitMap();
        iPoint4D dim = signalFile->getDimensions();
        iPoint3D voxel;
        for ( voxel.z=0; voxel.z<dim.z; voxel.z++ )
        {
            for ( voxel.y=0; voxel.y<dim.y; voxel.y++ )
            {
                for ( voxel.x=0; voxel.x<dim.x; voxel.x++ )
                {
                    double signal = signalFile->getVolumeValue(voxel);
                    if ( (above && signal >= value) || (!above && signal <= value) )
                        setBitMapValue(voxel,1.,Color_Undefined,true);
                }
            }
        }
    }
    updateBitMapVoxelSet();  // don't change the color, which should be carried by the input voxelSet
}

double bitMap::getNumberVoxelsInOverlay()
{
    double nVoxelsInOverlay = 0.;
    for (int jVoxel=0; jVoxel<_voxelSet.size(); jVoxel++)
    {
        double weight = _voxelSet[jVoxel].binaryValue;
        nVoxelsInOverlay += weight;
    }
    return nVoxelsInOverlay;
}

double bitMap::getVolumeInOverlay()
{
    double volumeInOverlay = 0.;
    for (int jVoxel=0; jVoxel<_voxelSet.size(); jVoxel++)
    {
        double weight = _voxelSet[jVoxel].binaryValue;
        volumeInOverlay  += weight * _data.hdr.resolution.x * _data.hdr.resolution.y * _data.hdr.resolution.z;
    }
    return volumeInOverlay;
}

void bitMap::updateBitMapVoxelSet(overlayColor inputColor)
{
    iPoint2D rangeX = {{0},{nx()-1}};
    iPoint2D rangeY = {{0},{ny()-1}};
    iPoint2D rangeZ = {{0},{nz()-1}};
    updateBitMapVoxelSet(inputColor,rangeX,rangeY,rangeZ);
}

void bitMap::updateBitMapVoxelSet(overlayColor inputColor,
                                  iPoint2D rangeX, iPoint2D rangeY, iPoint2D rangeZ)
{   // convert a bitmap (image) to a voxel list (vector)
    // this function can be slow when values are set more than once because
    // it needs to search the list _voxelSet
    FUNC_ENTER;
    _voxelSet.clear();
    for (int jColor=0; jColor<NCOLORS; jColor++)
        _colorInVoxelSet[jColor] = false;

    fVoxel voxel;
    for (voxel.z=rangeZ.lower; voxel.z<=rangeZ.upper; voxel.z++)
    {
        for (voxel.y=rangeY.lower; voxel.y<=rangeY.upper; voxel.y++)
        {
            for (voxel.x=rangeX.lower; voxel.x<=rangeX.upper; voxel.x++)
            {
                int i3d = index3dFromVoxel(voxel);
                voxel.index3d = i3d;
                voxel.binaryValue = qMax(0.,qMin(1.,getBitMapValue(i3d)));
                if ( voxel.binaryValue > 0.01 )  // transformations (smoothing, registration) can leave low values
                {   // Add the voxel; If the voxel color is unknown and a known color is passed, update the color
                    voxel.color = static_cast<overlayColor>(_data.color[i3d]);
                    if ( voxel.color == Color_Undefined )
                    {
                        if (inputColor != Color_Undefined )
                        {  // recolor undefined colors whenever the input color is defined
                            voxel.color = inputColor;
                            _data.color[i3d] = inputColor;
                            _colorInVoxelSet[(int)inputColor] = true;
                        }
                    }
                    else // if ( voxel.color >= Color_Cyan &&  voxel.color < Color_Undefined ) // voxel.color != Color_Undefined
                    {
                        int iColor = qMax(0,qMin(NCOLORS-1,(int)voxel.color));
                        _colorInVoxelSet[iColor] = true;
                    }
                    _voxelSet.append(voxel);
                }
            } // x
        } // y
    } // z
    _data.empty = _voxelSet.size() == 0;
    FUNC_EXIT;
}

void bitMap::setBitMapValue(const int i3d, const float value, bool updateVoxelSet)
{
    fVoxel voxelWithValue;
    iPoint3D voxel = utilIO::index3d(getDimensions(),i3d);
    voxelWithValue.x = voxel.x; voxelWithValue.y = voxel.y; voxelWithValue.z = voxel.z;
    voxelWithValue.binaryValue = value;
    voxelWithValue.color = Color_Undefined;
    setBitMapValue(voxelWithValue,updateVoxelSet);
}

void bitMap::setBitMapValue(const iPoint3D voxel, const float value, bool updateVoxelSet)
{
    fVoxel voxelWithValue;
    voxelWithValue.x = voxel.x; voxelWithValue.y = voxel.y; voxelWithValue.z = voxel.z;
    voxelWithValue.binaryValue = value;
    voxelWithValue.color = Color_Undefined;
    setBitMapValue(voxelWithValue,updateVoxelSet);
}

void bitMap::setBitMapValueBroadBrush(fVoxel voxel, bool updateVoxelSet)
{ // this function can be slow when values are set more than once because it needs to search the list _voxelSet
    CheckVoxelIndices(voxel,"setBitMapValue");
    if ( nt() != 1 )
        qFatal("Error: in setBitMapValue, # time points != 1:");
    else
    {
        fVoxel neighborVoxel = voxel;
        iPoint3D jRel;
        for ( jRel.z=-1; jRel.z<=1; jRel.z++)
        {
            neighborVoxel.z = voxel.z + jRel.z;
            for ( jRel.y=-1; jRel.y<=1; jRel.y++)
            {
                neighborVoxel.y = voxel.y + jRel.y;
                for ( jRel.x=-1; jRel.x<=1; jRel.x++)
                {
                    neighborVoxel.x = voxel.x + jRel.x;
                    int i3d = index3dFromVoxel(neighborVoxel);
                    if ( i3d >=0 && i3d < voxelsPerVolume() )
                    {
                        if ( updateVoxelSet )
                        {
                            if ( getBitMapValue(i3d) == 0. )
                                _voxelSet.append(neighborVoxel);    // only include the value in the vector if it's not already there
                            else
                            { // if the voxel is already set, we need to find it
                                for ( int jVoxel=0; jVoxel<_voxelSet.size(); jVoxel++ )
                                {
                                    if ( _voxelSet[jVoxel].x == neighborVoxel.x && _voxelSet[jVoxel].y == neighborVoxel.y && _voxelSet[jVoxel].z == neighborVoxel.z )
                                        _voxelSet[jVoxel].binaryValue = neighborVoxel.binaryValue;
                                }
                            }
                        }
                        setStackValue(i3d,0,static_cast<double>(neighborVoxel.binaryValue));
                        _data.color[i3d] = voxel.color;
                    }
                }
            }
        }
    }
    _data.empty = false;
}
void bitMap::setBitMapValueBroadBrush(const iPoint3D voxel3d, const float value, overlayColor color, bool updateVoxelSet)
{
    fVoxel voxel;
    voxel.x = voxel3d.x;  voxel.y = voxel3d.y;  voxel.z = voxel3d.z;
    voxel.binaryValue = value;
    voxel.color = color;
    setBitMapValueBroadBrush(voxel,updateVoxelSet);
}
void bitMap::setBitMapValue(fVoxelLight voxelLight, bool updateVoxelSet)
{
    fVoxel voxel;
    voxel.x = voxelLight.x;
    voxel.y = voxelLight.y;
    voxel.z = voxelLight.z;
    voxel.binaryValue = voxelLight.binaryValue;
    voxel.color = Color_Undefined;
    voxel.index3d = index3dFromVoxel(voxel);
    setBitMapValue(voxel,updateVoxelSet);
}
void bitMap::setBitMapValue(fVoxel voxel, bool updateVoxelSet)
{ // this function can be slow when values are set more than once because it needs to search the list _voxelSet
    if ( nt() != 1 )
        qFatal("Error: in setBitMapValue, # time points != 1:");
    else
    {
        CheckVoxelIndices(voxel,"setBitMapValue");
        int i3d = index3dFromVoxel(voxel);
        voxel.index3d = i3d;
        if ( updateVoxelSet )
        {
            if ( getBitMapValue(i3d) == 0. )
                _voxelSet.append(voxel);    // only include the value in the vector if it's not already there
            else
            { // if the voxel is already set, we need to find it
                for ( int jVoxel=0; jVoxel<_voxelSet.size(); jVoxel++ )
                {
                    if ( _voxelSet[jVoxel].x == voxel.x && _voxelSet[jVoxel].y == voxel.y && _voxelSet[jVoxel].z == voxel.z )
                    {
                        _voxelSet[jVoxel].binaryValue = voxel.binaryValue;
                        _voxelSet[jVoxel].index3d = i3d;
                    }
                }
            }
        }
        setStackValue(i3d,0,voxel.binaryValue);
        _data.color[i3d] = voxel.color;
    }
    _data.empty = false;
}
void bitMap::setBitMapValue(const iPoint3D voxel3d, const float value, overlayColor color, bool updateVoxelSet)
{
    fVoxel voxel;
    voxel.x = voxel3d.x;  voxel.y = voxel3d.y;  voxel.z = voxel3d.z;  voxel.binaryValue = value;  voxel.color = color;
    setBitMapValue(voxel,updateVoxelSet);
}
void bitMap::setBitMapValue(const int iX, const int iY, const int iZ, const float value, overlayColor color, bool updateVoxelSet)
{
    fVoxel voxel;
    voxel.x = iX;  voxel.y = iY;  voxel.z = iZ;  voxel.binaryValue = value;  voxel.color = color;
    setBitMapValue(voxel,updateVoxelSet);
}

double bitMap::getBitMapValue(int i3d)
{
    if ( nt() != 1 )
    {
        qInfo() << "# time points in bitmap = " << nt();
        qFatal("Error: in getBitMapValue, # time points != 1");
    }
    else
    {
        CheckVoxelIndices(i3d, "getBitMapValue");
        return getVolumeValue(i3d);
    }
}

double bitMap::getBitMapValue(int iX, int iY, int iZ)
{
    if ( nt() != 1 )
    {
        qInfo() << "# time points in bitmap = " << nt();
        qFatal("Error: in getBitMapValue, # time points != 1");
    }
    else
    {
        CheckVoxelIndices(iX,iY,iZ,"getBitMapValue");
        int i3d = index3dFromVoxel(iX, iY, iZ);
        return getVolumeValue(i3d);
    }
}
double bitMap::getBitMapValue(iPoint3D voxel)
{
    if ( nt() != 1 )
    {
        qInfo() << "# time points in bitmap = " << nt();
        qFatal("Error: in getBitMapValue, # time points != 1");
    }
    else
    {
        CheckVoxelIndices(voxel,"getBitMapValue");
        int i3d = index3dFromVoxel(voxel);
        return getVolumeValue(i3d);
    }
}
double bitMap::getBitMapValue(fVoxel voxel)
{
    if ( nt() != 1 )
        qFatal("Error: in getBitMapValue, # time points != 1");
    else
    {
        CheckVoxelIndices(voxel,"getBitMapValue");
        int i3d = index3dFromVoxel(voxel);
        return getVolumeValue(i3d);
    }
}

overlayColor bitMap::getBitMapColor(int iX, int iY, int iZ)
{
    if ( nt() != 1 )
        qFatal("Error: in getBitMapColor, # time points != 1");
    else
    {
        CheckVoxelIndices(iX,iY,iZ,"getBitMapColor");
        int i3d = index3dFromVoxel(iX, iY, iZ);
        overlayColor color = (overlayColor)_data.color[i3d];
        return color;
    }
}
overlayColor bitMap::getBitMapColor(iPoint3D voxel)
{
    if ( nt() != 1 )
        qFatal("Error: in getBitMapColor, # time points != 1");
    else
    {
        CheckVoxelIndices(voxel,"getBitMapColor");
        int i3d = index3dFromVoxel(voxel);
        overlayColor color = (overlayColor)_data.color[i3d];
        return color;
    }
}
overlayColor bitMap::getBitMapColor(fVoxel voxel)
{
    if ( nt() != 1 )
        qFatal("Error: in getBitMapColor, # time points != 1");
    else
    {
        CheckVoxelIndices(voxel,"getBitMapColor");
        int i3d = index3dFromVoxel(voxel);
        overlayColor color = (overlayColor)_data.color[i3d];
        return color;
    }
}
