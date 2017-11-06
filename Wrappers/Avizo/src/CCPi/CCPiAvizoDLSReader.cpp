/*
 *  Template of a read routine
 */

#include <hxcore/HxData.h>
#include <hxcore/HxMessage.h>

#include "CCPiNexusReader.h"
#include "CCPiNexusWidgetDialog.h"

#include <hxfield/HxUniformScalarField3.h>
#include <hxfield/HxUniformCoord3.h>
#include <hxfield/HxRegScalarField3.h>
#include <hxfield/internal/HxRectilinearCoord3.h>
#include "api.h"

typedef struct{
	void *data;
	int  ndims;
	int  *dims;
	CCPiNexusReader::DATATYPE dataType;
	double* axisData;
} TomoData;

void CCPiRegisterDLSDataset(std::string name, TomoData* data, TomoData* angles, TomoData* image_key, TomoData* pixel_size_x,TomoData* pixel_size_y);

TomoData* readData(CCPiNexusReader& reader, std::string name)
{
	TomoData* tData = new TomoData();
	reader.ReadCompleteData(name, &tData->data, &tData->ndims, &tData->dims, &tData->dataType,&tData->axisData);
	return tData;
}

CCPI_API
int CCPiAvizoDLSReader(const char* filename)
{

	std::string data_entry = "/entry1/tomo_entry/data/data";
	std::string rotation_angle = "/entry1/tomo_entry/data/rotation_angle";
	std::string image_key = "/entry1/tomo_entry/instrument/detector/image_key";
	std::string pixel_size_x = "/entry1/tomo_entry/instrument/detector/x_pixel_size";
	std::string pixel_size_y = "/entry1/tomo_entry/instrument/detector/y_pixel_size";

	theMsg->stream() <<"Started Reading the Datasets"<<std::endl;

	CCPiNexusReader reader(filename);
	TomoData *data = readData(reader, data_entry);
	TomoData *angles = readData(reader, rotation_angle);
	TomoData *keys = readData(reader, image_key);
	TomoData *pixel_x = readData(reader, pixel_size_x);
	TomoData *pixel_y = readData(reader, pixel_size_y);

	std::cout<<"Number of Dimensions"<<data->ndims<<std::endl;
	CCPiRegisterDLSDataset(data_entry, data, angles, keys, pixel_x, pixel_y);
    return 1;
}


template <class T>
HxRegScalarField3* CCPiRegisterDLSAvizoDataset(std::string name, T *data, int ndims, int *dims, McPrimType::Type dataType,double *axisData)
{
	int newDims[3];
	if(ndims==3){ //Row major , column major
		newDims[0]=dims[2];
		newDims[1]=dims[1];
		newDims[2]=dims[0];
	} else  if(ndims==2) {
		newDims[0]=dims[1];
		newDims[1]=dims[0];
		newDims[2]=1;
	} else {
		return NULL;
	}
	theMsg->stream() <<"Dataset Dimensions "<<dims[0]<<" "<<dims[1]<<" "<<dims[2]<<std::endl;
	HxRegScalarField3 * field;

	if(axisData!=NULL) {
		field = new HxRegScalarField3((newDims), dataType, C_RECTILINEAR);
	}else{
		field = new HxUniformScalarField3((newDims), dataType);
	}
	T *rawout=(T *)field->lattice().dataPtr();
	hsize_t totalsize = 1;
	for(int idx =0; idx < ndims;idx++) totalsize *=dims[idx];
	for(long idx=0;idx<totalsize;idx++)
		rawout[idx]=((T*)data)[idx];

	if(axisData!=NULL) {
		HxRectilinearCoord3 * coords =(HxRectilinearCoord3 *) field->lattice().coords();

		float *coordValues = coords->coordX();
		for(int i=0;i<newDims[0];i++)
		{
			*(coordValues+i) = *(axisData+newDims[2]+newDims[1]+i);
		}

		coordValues = coords->coordY();
		for(int i=0;i<newDims[1];i++)
		{
			*(coordValues+i) = *(axisData+newDims[2]+i);
		}

		coordValues = coords->coordZ();
		for(int i=0;i<newDims[2];i++)
		{
			*(coordValues+i) = *(axisData+i);
		}
	}else{
		HxUniformCoord3 * coords =(HxUniformCoord3 *) field->lattice().coords();
		McBox3f bx=coords->getBoundingBox();
		bx[0]=0;
		if(newDims[0]>1)
			bx[1]=(float)newDims[0]-1;
		else
			bx[1]=1;
		bx[2]=0;
		if(newDims[1]>1)
			bx[3]=(float)newDims[1]-1;
		else
			bx[3]=1;
		bx[4]=0;
		if(ndims==3 || newDims[2]==1)
			bx[5]=(float)newDims[2]-1;
		else
			bx[5] = 1;
	}
	return field;
}


void CCPiRegisterDLSDataset(std::string name, TomoData* data, TomoData* angles, TomoData* image_key, TomoData* pixel_size_x,TomoData* pixel_size_y)
{
	HxRegScalarField3 * field;
	switch(data->dataType)
	{
	case CCPiNexusReader::CHAR:
		field = CCPiRegisterDLSAvizoDataset(name, (char*) data->data,data->ndims,data->dims,McPrimType::MC_INT8,data->axisData);
		break;
	case CCPiNexusReader::UCHAR:
		field = CCPiRegisterDLSAvizoDataset(name, (unsigned char*) data->data,data->ndims,data->dims,McPrimType::MC_UINT8,data->axisData);
		break;
	case CCPiNexusReader::SHORT:
		field = CCPiRegisterDLSAvizoDataset(name, (short*) data->data,data->ndims,data->dims,McPrimType::MC_INT16,data->axisData);
		break;
	case CCPiNexusReader::USHORT:
		field = CCPiRegisterDLSAvizoDataset(name, (unsigned short*) data->data,data->ndims,data->dims,McPrimType::MC_UINT16,data->axisData);
		break;
	case CCPiNexusReader::INT:
		field = CCPiRegisterDLSAvizoDataset(name, (int*) data->data,data->ndims,data->dims,McPrimType::MC_INT32,data->axisData);
		break;
	case CCPiNexusReader::USINT:
		field = CCPiRegisterDLSAvizoDataset(name, (unsigned int*) data->data,data->ndims,data->dims,McPrimType::MC_UINT32,data->axisData);
		break;
	case CCPiNexusReader::LONG:
		field = CCPiRegisterDLSAvizoDataset(name, (long*) data->data,data->ndims,data->dims,McPrimType::MC_INT32,data->axisData);
		break;
	case CCPiNexusReader::ULONG:
		field = CCPiRegisterDLSAvizoDataset(name, (unsigned long*) data->data,data->ndims,data->dims,McPrimType::MC_UINT32,data->axisData);
		break;
	case CCPiNexusReader::LLONG:
		field = CCPiRegisterDLSAvizoDataset(name, (long long*) data->data,data->ndims,data->dims,McPrimType::MC_INT64,data->axisData);
		break;
	case CCPiNexusReader::ULLONG:
		field = CCPiRegisterDLSAvizoDataset(name, (unsigned long long*) data->data,data->ndims,data->dims,McPrimType::MC_UINT64,data->axisData);
		break;
	case CCPiNexusReader::FLOAT:
		field = CCPiRegisterDLSAvizoDataset(name, (float*) data->data,data->ndims,data->dims,McPrimType::MC_FLOAT,data->axisData);
		break;
	case CCPiNexusReader::DOUBLE:
		field = CCPiRegisterDLSAvizoDataset(name, (double*) data->data,data->ndims,data->dims,McPrimType::MC_DOUBLE,data->axisData);
		break;
	}
	/* Set angles, image key and pixel size as parameters*/
	/*TODO:: Change the angles and image key data to be generic instead of double*/
	field->parameters.set("Angles", angles->dims[0], (double*)angles->data); 
	field->parameters.set("ImageKey", image_key->dims[0], (double*)image_key->data);
	double pixel_size[2];
	if(pixel_size_x->dataType != CCPiNexusReader::DOUBLE) {
		pixel_size[0]=1;
		pixel_size[1]=1;
	}else{
		pixel_size[0]=((double*)pixel_size_x->data)[0];
		pixel_size[1]=((double*)pixel_size_y->data)[0];
	}
	field->parameters.set("PixelSize", 2, pixel_size);
	HxData::registerData(field, name.c_str());
}



