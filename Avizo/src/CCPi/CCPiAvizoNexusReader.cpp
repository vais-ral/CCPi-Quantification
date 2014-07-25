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

void CCPiRegisterUniformDataset(std::string name, void *data, int ndims, int *dims, CCPiNexusReader::DATATYPE dataType);
void CCPiRegisterRegularDataset(std::string name, void *data, int ndims, int *dims, CCPiNexusReader::DATATYPE dataType, double* axisData);
CCPI_API
int CCPiAvizoNexusReader(const char* filename)
{

	CCPiNexusWidgetDialog* nexusDlg=new CCPiNexusWidgetDialog(filename);
	nexusDlg->exec();


	theMsg->stream() <<"Started Reading the Datasets"<<std::endl;

	//std::vector<std::string> output = nexusDlg->GetSelectedDataSetList();
	for(int index=0; index < nexusDlg->GetSelectedDataSetCount();index++)
	{
		std::string datasetname(nexusDlg->GetSelectedDataSet(index));
		CCPiNexusReader reader(filename);
		int ndims;
		int *dims;
		CCPiNexusReader::DATATYPE dataType;
		void* data;
		double* axisData=NULL;
		reader.ReadCompleteData(datasetname, &data, &ndims, &dims, &dataType,&axisData);
		std::cout<<"Number of Dimensions"<<ndims<<std::endl;
		if(axisData==NULL)
		{
			CCPiRegisterUniformDataset(datasetname, data, ndims, dims, dataType);
		}else{
			CCPiRegisterRegularDataset(datasetname, data, ndims, dims, dataType,axisData);
		}
	}
	delete nexusDlg;
	//std::vector<std::string> output = nexusDlg->GetSelectedDataSetList();
	//for(std::vector<std::string>::iterator itr = output.begin(); itr!=output.end();itr++)
	//{
	//	CCPiNexusReader reader(filename);
	//	int ndims;
	//	int *dims;
	//	CCPiNexusReader::DATATYPE dataType;
	//	void* data;
	//	double* axisData=NULL;
	//	reader.ReadCompleteData(*(itr), &data, &ndims, &dims, &dataType,&axisData);
	//	std::cout<<"Number of Dimensions"<<ndims<<std::endl;
	//	if(axisData==NULL)
	//	{
	//		CCPiRegisterUniformDataset(*(itr), data, ndims, dims, dataType);
	//	}else{
	//		CCPiRegisterRegularDataset(*(itr), data, ndims, dims, dataType,axisData);
	//	}
	//}

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

void CCPiRegisterUniformDataset(std::string name, void *data, int ndims, int *dims, CCPiNexusReader::DATATYPE dataType)
{
	switch(dataType)
	{
	case CCPiNexusReader::CHAR:
		CCPiRegisterAvizoUniformDataset(name, (char*) data,ndims,dims,McPrimType::mc_int8);
		break;
	case CCPiNexusReader::UCHAR:
		CCPiRegisterAvizoUniformDataset(name, (unsigned char*) data,ndims,dims,McPrimType::mc_uint8);
		break;
	case CCPiNexusReader::SHORT:
		CCPiRegisterAvizoUniformDataset(name, (short*) data,ndims,dims,McPrimType::mc_int16);
		break;
	case CCPiNexusReader::USHORT:
		CCPiRegisterAvizoUniformDataset(name, (unsigned short*) data,ndims,dims,McPrimType::mc_uint16);
		break;
	case CCPiNexusReader::INT:
		CCPiRegisterAvizoUniformDataset(name, (int*) data,ndims,dims,McPrimType::mc_int32);
		break;
	case CCPiNexusReader::UINT:
		CCPiRegisterAvizoUniformDataset(name, (unsigned int*) data,ndims,dims,McPrimType::mc_uint32);
		break;
	case CCPiNexusReader::LONG:
		CCPiRegisterAvizoUniformDataset(name, (long*) data,ndims,dims,McPrimType::mc_int32);
		break;
	case CCPiNexusReader::ULONG:
		CCPiRegisterAvizoUniformDataset(name, (unsigned long*) data,ndims,dims,McPrimType::mc_uint32);
		break;
	case CCPiNexusReader::LLONG:
		CCPiRegisterAvizoUniformDataset(name, (long long*) data,ndims,dims,McPrimType::mc_int64);
		break;
	case CCPiNexusReader::ULLONG:
		CCPiRegisterAvizoUniformDataset(name, (unsigned long long*) data,ndims,dims,McPrimType::mc_uint64);
		break;
	case CCPiNexusReader::FLOAT:
		CCPiRegisterAvizoUniformDataset(name, (float*) data,ndims,dims,McPrimType::mc_float);
		break;
	case CCPiNexusReader::DOUBLE:
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

void CCPiRegisterRegularDataset(std::string name, void *data, int ndims, int *dims, CCPiNexusReader::DATATYPE dataType, double* axisData)
{
	switch(dataType)
	{
	case CCPiNexusReader::CHAR:
		CCPiRegisterAvizoRegularDataset(name, (char*) data,ndims,dims,McPrimType::mc_int8, axisData);
		break;
	case CCPiNexusReader::UCHAR:
		CCPiRegisterAvizoRegularDataset(name, (unsigned char*) data,ndims,dims,McPrimType::mc_uint8, axisData);
		break;
	case CCPiNexusReader::SHORT:
		CCPiRegisterAvizoRegularDataset(name, (short*) data,ndims,dims,McPrimType::mc_int16, axisData);
		break;
	case CCPiNexusReader::USHORT:
		CCPiRegisterAvizoRegularDataset(name, (unsigned short*) data,ndims,dims,McPrimType::mc_uint16, axisData);
		break;
	case CCPiNexusReader::INT:
		CCPiRegisterAvizoRegularDataset(name, (int*) data,ndims,dims,McPrimType::mc_int32, axisData);
		break;
	case CCPiNexusReader::UINT:
		CCPiRegisterAvizoRegularDataset(name, (unsigned int*) data,ndims,dims,McPrimType::mc_uint32, axisData);
		break;
	case CCPiNexusReader::LONG:
		CCPiRegisterAvizoRegularDataset(name, (long*) data,ndims,dims,McPrimType::mc_int32, axisData);
		break;
	case CCPiNexusReader::ULONG:
		CCPiRegisterAvizoRegularDataset(name, (unsigned long*) data,ndims,dims,McPrimType::mc_uint32, axisData);
		break;
	case CCPiNexusReader::LLONG:
		CCPiRegisterAvizoRegularDataset(name, (long long*) data,ndims,dims,McPrimType::mc_int64, axisData);
		break;
	case CCPiNexusReader::ULLONG:
		CCPiRegisterAvizoRegularDataset(name, (unsigned long long*) data,ndims,dims,McPrimType::mc_uint64, axisData);
		break;
	case CCPiNexusReader::FLOAT:
		CCPiRegisterAvizoRegularDataset(name, (float*) data,ndims,dims,McPrimType::mc_float, axisData);
		break;
	case CCPiNexusReader::DOUBLE:
		CCPiRegisterAvizoRegularDataset(name, (double*) data,ndims,dims,McPrimType::mc_double, axisData);
		break;
	}
}