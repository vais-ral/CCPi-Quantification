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
#include <hxfield/HxRectilinearCoord3.h>
#include "api.h"

void CCPiRegisterUniformDataset(std::string name, void *data, int ndims, int *dims, NexusReader::DATATYPE dataType);
void CCPiRegisterRegularDataset(std::string name, void *data, int ndims, int *dims, NexusReader::DATATYPE dataType, double* axisData);
CCPI_API
int CCPiAvizoNexusReader(const char* filename)
{

	CCPiNexusWidgetDialog* nexusDlg=new CCPiNexusWidgetDialog(filename);
	nexusDlg->exec();

	std::cout<<"Started Reading the Datasets"<<std::endl;
	std::vector<std::string> output = nexusDlg->GetSelectedDataSetList();
	for(std::vector<std::string>::iterator itr = output.begin(); itr!=output.end();itr++)
	{
		CCPiNexusReader reader(filename);
		int ndims;
		int *dims;
		NexusReader::DATATYPE dataType;
		void* data;
		double* axisData=NULL;
		reader.ReadCompleteData(*(itr), &data, &ndims, &dims, &dataType,&axisData);
		std::cout<<"Number of Dimensions"<<ndims<<std::endl;
		if(axisData==NULL)
		{
			CCPiRegisterUniformDataset(*(itr), data, ndims, dims, dataType);
		}else{
			CCPiRegisterRegularDataset(*(itr), data, ndims, dims, dataType,axisData);
		}
	}

    return 1;
}




template <class T>
void CCPiRegisterAvizoUniformDataset(std::string name, T *data, int ndims, int *dims, McPrimType::Type dataType)
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
		return;
	}
	HxUniformScalarField3 * field;

	field = new HxUniformScalarField3((newDims), dataType);  
	T *rawout=(T *)field->lattice.dataPtr();
	hsize_t totalsize = 1;
	for(int idx =0; idx < ndims;idx++) totalsize *=dims[idx];
	for(long idx=0;idx<totalsize;idx++)
		rawout[idx]=((T*)data)[idx];

	HxUniformCoord3 * coords =(HxUniformCoord3 *) field->lattice.coords();
	float * bx;
	bx=coords->bbox();
	bx[0]=0;
	if(newDims[0]!=1)
		bx[1]=(float)newDims[0]-1;
	else
		bx[1]=1;
	bx[2]=0;
	if(newDims[1]!=1)
		bx[3]=(float)newDims[1]-1;
	else
		bx[3]=1;
	bx[4]=0;
	if(ndims==3 || newDims[2]==1)
		bx[5]=(float)newDims[2]-1;
	else
		bx[5] = 1;
	HxData::registerData(field, name.c_str());
}

void CCPiRegisterUniformDataset(std::string name, void *data, int ndims, int *dims, NexusReader::DATATYPE dataType)
{
	switch(dataType)
	{
	case NexusReader::CHAR:
		CCPiRegisterAvizoUniformDataset(name, (char*) data,ndims,dims,McPrimType::mc_int8);
		break;
	case NexusReader::UCHAR:
		CCPiRegisterAvizoUniformDataset(name, (unsigned char*) data,ndims,dims,McPrimType::mc_uint8);
		break;
	case NexusReader::SHORT:
		CCPiRegisterAvizoUniformDataset(name, (short*) data,ndims,dims,McPrimType::mc_int16);
		break;
	case NexusReader::USHORT:
		CCPiRegisterAvizoUniformDataset(name, (unsigned short*) data,ndims,dims,McPrimType::mc_uint16);
		break;
	case NexusReader::INT:
		CCPiRegisterAvizoUniformDataset(name, (int*) data,ndims,dims,McPrimType::mc_int32);
		break;
	case NexusReader::UINT:
		CCPiRegisterAvizoUniformDataset(name, (unsigned int*) data,ndims,dims,McPrimType::mc_uint32);
		break;
	case NexusReader::LONG:
		CCPiRegisterAvizoUniformDataset(name, (long*) data,ndims,dims,McPrimType::mc_int32);
		break;
	case NexusReader::ULONG:
		CCPiRegisterAvizoUniformDataset(name, (unsigned long*) data,ndims,dims,McPrimType::mc_uint32);
		break;
	case NexusReader::LLONG:
		CCPiRegisterAvizoUniformDataset(name, (long long*) data,ndims,dims,McPrimType::mc_int64);
		break;
	case NexusReader::ULLONG:
		CCPiRegisterAvizoUniformDataset(name, (unsigned long long*) data,ndims,dims,McPrimType::mc_uint64);
		break;
	case NexusReader::FLOAT:
		CCPiRegisterAvizoUniformDataset(name, (float*) data,ndims,dims,McPrimType::mc_float);
		break;
	case NexusReader::DOUBLE:
		CCPiRegisterAvizoUniformDataset(name, (double*) data,ndims,dims,McPrimType::mc_double);
		break;
	}
}



template <class T>
void CCPiRegisterAvizoRegularDataset(std::string name, T *data, int ndims, int *dims, McPrimType::Type dataType,double *axisData)
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
		return;
	}
	HxRegScalarField3 * field;

	field = new HxRegScalarField3((newDims), dataType,c_rectilinear);  
	T *rawout=(T *)field->lattice.dataPtr();
	hsize_t totalsize = 1;
	for(int idx =0; idx < ndims;idx++) totalsize *=dims[idx];
	for(long idx=0;idx<totalsize;idx++)
		rawout[idx]=((T*)data)[idx];

	HxRectilinearCoord3 * coords =(HxRectilinearCoord3 *) field->lattice.coords();

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
	HxData::registerData(field, name.c_str());
}

void CCPiRegisterRegularDataset(std::string name, void *data, int ndims, int *dims, NexusReader::DATATYPE dataType, double* axisData)
{
	switch(dataType)
	{
	case NexusReader::CHAR:
		CCPiRegisterAvizoRegularDataset(name, (char*) data,ndims,dims,McPrimType::mc_int8, axisData);
		break;
	case NexusReader::UCHAR:
		CCPiRegisterAvizoRegularDataset(name, (unsigned char*) data,ndims,dims,McPrimType::mc_uint8, axisData);
		break;
	case NexusReader::SHORT:
		CCPiRegisterAvizoRegularDataset(name, (short*) data,ndims,dims,McPrimType::mc_int16, axisData);
		break;
	case NexusReader::USHORT:
		CCPiRegisterAvizoRegularDataset(name, (unsigned short*) data,ndims,dims,McPrimType::mc_uint16, axisData);
		break;
	case NexusReader::INT:
		CCPiRegisterAvizoRegularDataset(name, (int*) data,ndims,dims,McPrimType::mc_int32, axisData);
		break;
	case NexusReader::UINT:
		CCPiRegisterAvizoRegularDataset(name, (unsigned int*) data,ndims,dims,McPrimType::mc_uint32, axisData);
		break;
	case NexusReader::LONG:
		CCPiRegisterAvizoRegularDataset(name, (long*) data,ndims,dims,McPrimType::mc_int32, axisData);
		break;
	case NexusReader::ULONG:
		CCPiRegisterAvizoRegularDataset(name, (unsigned long*) data,ndims,dims,McPrimType::mc_uint32, axisData);
		break;
	case NexusReader::LLONG:
		CCPiRegisterAvizoRegularDataset(name, (long long*) data,ndims,dims,McPrimType::mc_int64, axisData);
		break;
	case NexusReader::ULLONG:
		CCPiRegisterAvizoRegularDataset(name, (unsigned long long*) data,ndims,dims,McPrimType::mc_uint64, axisData);
		break;
	case NexusReader::FLOAT:
		CCPiRegisterAvizoRegularDataset(name, (float*) data,ndims,dims,McPrimType::mc_float, axisData);
		break;
	case NexusReader::DOUBLE:
		CCPiRegisterAvizoRegularDataset(name, (double*) data,ndims,dims,McPrimType::mc_double, axisData);
		break;
	}
}