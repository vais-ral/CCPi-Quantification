#include "CCPiNexusReader.h"
#include <iostream>
#include <map>
#include "qmutex.h"

CCPiNexusReader::CCPiNexusReader(std::string filename)
{
	Filename = filename;

	FileId = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
}

CCPiNexusReader::~CCPiNexusReader()
{
	H5Fclose(FileId);
}

void CCPiNexusReader::ReadCompleteData(std::string datasetPath, void** data,int *ndims, int** dims, DATATYPE* dataType, double** axisData)
{
	hid_t dataset_id = H5Dopen(FileId,datasetPath.c_str(),H5P_DEFAULT);
	hid_t dataspace = H5Dget_space(dataset_id);
	int rank = H5Sget_simple_extent_ndims(dataspace);
	hsize_t *dims_out = new hsize_t[rank];
	hsize_t *maxdims =  new hsize_t[rank];
	H5Sget_simple_extent_dims(dataspace,dims_out,maxdims);

	hid_t datasetDataType = H5Dget_type(dataset_id);
	hid_t nativeDataType = H5Tget_native_type(datasetDataType,H5T_DIR_DEFAULT);

	*ndims = rank;
	*dims = new int[rank];
	for(int idx=0;idx<rank;idx++) (*dims)[idx] = dims_out[idx];
	*dataType = GetDataType(nativeDataType);

	void* allocatedData = AllocateMemory(datasetDataType, *ndims, dims_out);
	int status = H5Dread(dataset_id, nativeDataType,H5S_ALL,H5S_ALL,H5P_DEFAULT,allocatedData);
	*data = allocatedData;
	if(status>0) std::cout<<"Successfully read the data"<<std::endl;

	//Check to read if there are axis information
	if(isSignalData(dataset_id))
	{
		std::cout<<"It's a Signal Data"<<std::endl;
		bool status = ReadAxisData(datasetPath, *ndims, dims_out, (void**)axisData);
		if(!status)
		{
			delete[] ((double*)*axisData);
			*axisData=NULL;
		}
	}
	H5Tclose(nativeDataType);
	H5Tclose(datasetDataType);
	H5Dclose(dataset_id);
	delete[] dims_out;
	delete[] maxdims;
}

void CCPiNexusReader::ReadCompleteDataNoAllocation(std::string datasetPath, void* data)
{
	hid_t dataset_id = H5Dopen(FileId,datasetPath.c_str(),H5P_DEFAULT);
	hid_t dataspace = H5Dget_space(dataset_id);
	int rank = H5Sget_simple_extent_ndims(dataspace);
	hsize_t *dims_out = new hsize_t[rank];
	hsize_t *maxdims =  new hsize_t[rank];
	H5Sget_simple_extent_dims(dataspace,dims_out,maxdims);

	hid_t datasetDataType = H5Dget_type(dataset_id);
	hid_t nativeDataType = H5Tget_native_type(datasetDataType,H5T_DIR_DEFAULT);

	int status = H5Dread(dataset_id, nativeDataType,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);

	if(status>0) std::cout<<"Successfully read the data"<<std::endl;

	H5Tclose(nativeDataType);
	H5Tclose(datasetDataType);
	H5Dclose(dataset_id);
	delete[] dims_out;
	delete[] maxdims;
}

void* CCPiNexusReader::AllocateMemory(hid_t datatype, int ndims, hsize_t *dims)
{
	hsize_t totalsize = 1;
	for(int idx =0; idx < ndims;idx++) totalsize *=dims[idx];
	std::cout<<"Allocating memory of "<<totalsize<<std::endl;
	return AllocateMemory(datatype,totalsize);
}

void* CCPiNexusReader::AllocateMemory(hid_t datatype, hsize_t totalsize)
{
	std::cout<<"Allocating memory of "<<totalsize<<std::endl;
	if(H5Tget_class(datatype) == H5T_INTEGER && H5Tget_size(datatype) == H5Tget_size(H5T_NATIVE_CHAR))
	{
		std::cout<<"It's a Char type"<<std::endl;
		char *data = new char[totalsize];
		return data;
	} else if(H5Tget_class(datatype) == H5T_INTEGER && H5Tget_size(datatype) ==  H5Tget_size(H5T_NATIVE_SHORT)) {
		short *data = new short[totalsize];
		std::cout<<"It's a short type"<<std::endl;
		return data;
	} else if(H5Tget_class(datatype) == H5T_INTEGER && H5Tget_size(datatype) ==  H5Tget_size(H5T_NATIVE_INT)) {
		int *data = new int[totalsize];
		std::cout<<"It's a int type"<<std::endl;
		return data;
	} else if(H5Tget_class(datatype) == H5T_INTEGER && H5Tget_size(datatype) ==  H5Tget_size(H5T_NATIVE_LONG)) {
		long *data = new long[totalsize];
		std::cout<<"It's a long type"<<std::endl;
		return data;
	} else if(H5Tget_class(datatype) == H5T_INTEGER && H5Tget_size(datatype) == H5Tget_size(H5T_NATIVE_ULLONG)) {
		long long *data = new long long[totalsize];
		std::cout<<"It's a long long type"<<std::endl;
		return data;
	} else if(H5Tget_class(datatype) == H5T_FLOAT && H5Tget_size(datatype) == H5Tget_size(H5T_NATIVE_FLOAT)) {
		float *data = new float[totalsize];
		std::cout<<"It's a float type"<<std::endl;
		return data;
	} else if(H5Tget_class(datatype) == H5T_FLOAT && H5Tget_size(datatype) == H5Tget_size(H5T_NATIVE_DOUBLE)) {
		std::cout<<"It's a Double type"<<std::endl;
		double *data = new double[totalsize];
		return data;
	} else if(H5Tget_class(datatype) == H5T_FLOAT) {
		long double *data = new long double[totalsize];
		return data;
	} else {
		return NULL;
	}
}

CCPiNexusReader::DATATYPE CCPiNexusReader::GetDataType(hid_t datatype)
{

	if(H5Tget_class(datatype) == H5T_INTEGER && H5Tget_size(datatype) ==  H5Tget_size(H5T_NATIVE_CHAR) && H5Tget_sign(datatype) == H5T_SGN_2 )
	{
		return CHAR;
	} else if(H5Tget_class(datatype) == H5T_INTEGER && H5Tget_size(datatype) ==  H5Tget_size(H5T_NATIVE_UCHAR) && H5Tget_sign(datatype) == H5T_SGN_NONE ){
		return UCHAR;
	} else if(H5Tget_class(datatype) == H5T_INTEGER && H5Tget_size(datatype) ==  H5Tget_size(H5T_NATIVE_SHORT) && H5Tget_sign(datatype) == H5T_SGN_2 ){
		return SHORT;
	} else if(H5Tget_class(datatype) == H5T_INTEGER && H5Tget_size(datatype) ==  H5Tget_size(H5T_NATIVE_SHORT) && H5Tget_sign(datatype) == H5T_SGN_NONE ) {
		return USHORT;
	} else if(H5Tget_class(datatype) == H5T_INTEGER && H5Tget_size(datatype) ==  H5Tget_size(H5T_NATIVE_INT) && H5Tget_sign(datatype) == H5T_SGN_2 ) {
		return INT;
	} else if(H5Tget_class(datatype) == H5T_INTEGER && H5Tget_size(datatype) ==  H5Tget_size(H5T_NATIVE_INT) && H5Tget_sign(datatype) == H5T_SGN_NONE ) {
		return UINT;
	} else if(H5Tget_class(datatype) == H5T_INTEGER && H5Tget_size(datatype) ==  H5Tget_size(H5T_NATIVE_LONG) && H5Tget_sign(datatype) == H5T_SGN_2 ) {
		return LONG;
	} else if(H5Tget_class(datatype) == H5T_INTEGER && H5Tget_size(datatype) ==  H5Tget_size(H5T_NATIVE_LONG) && H5Tget_sign(datatype) == H5T_SGN_NONE ) {
		return ULONG;
	} else if(H5Tget_class(datatype) == H5T_INTEGER && H5Tget_size(datatype) ==  H5Tget_size(H5T_NATIVE_LLONG) && H5Tget_sign(datatype) == H5T_SGN_2 ) {
		return LLONG;
	} else if(H5Tget_class(datatype) == H5T_INTEGER && H5Tget_size(datatype) ==  H5Tget_size(H5T_NATIVE_LLONG) && H5Tget_sign(datatype) == H5T_SGN_NONE ) {
		return ULLONG;
	} else if(H5Tget_class(datatype) == H5T_FLOAT && H5Tget_size(datatype) ==  H5Tget_size(H5T_NATIVE_FLOAT)) {
		return FLOAT;
	} else if(H5Tget_class(datatype) == H5T_FLOAT && H5Tget_size(datatype) ==  H5Tget_size(H5T_NATIVE_DOUBLE)) {
		return DOUBLE;
	} else if(H5Tget_class(datatype) == H5T_FLOAT && H5Tget_size(datatype) ==  H5Tget_size(H5T_NATIVE_LDOUBLE)) {
		return LDOUBLE;
	} else {
		return UNKNOWN;
	}

}

int  CCPiNexusReader::GetDataNumberOfDimensions(std::string datasetPath)
{
	hid_t dataset_id = H5Dopen(FileId,datasetPath.c_str(),H5P_DEFAULT);
	hid_t dataspace = H5Dget_space(dataset_id);
	int rank = H5Sget_simple_extent_ndims(dataspace);
	H5Dclose(dataset_id);
	return rank;
}

void CCPiNexusReader::GetDataDimensions(std::string datasetPath, int *dims)
{
	hid_t dataset_id = H5Dopen(FileId,datasetPath.c_str(),H5P_DEFAULT);
	hid_t dataspace = H5Dget_space(dataset_id);
	int rank = H5Sget_simple_extent_ndims(dataspace);
	hsize_t *dims_out = new hsize_t[rank];
	hsize_t *maxdims =  new hsize_t[rank];
	H5Sget_simple_extent_dims(dataspace,dims_out,maxdims);
	for(int i=0;i<rank;i++)
		dims[i]=dims_out[i];
	H5Dclose(dataset_id);
	delete[] dims_out;
	delete[] maxdims;
}

CCPiNexusReader::DATATYPE CCPiNexusReader::GetDataType(std::string datasetPath)
{
	hid_t dataset_id = H5Dopen(FileId,datasetPath.c_str(),H5P_DEFAULT);
	hid_t datasetDataType = H5Dget_type(dataset_id);
	hid_t nativeDataType = H5Tget_native_type(datasetDataType,H5T_DIR_DEFAULT);
	DATATYPE dataType = GetDataType(nativeDataType);
	H5Dclose(dataset_id);
	return dataType;
}

bool CCPiNexusReader::isSignalData(hid_t dataset)
{
	hid_t signal_id = H5Aopen_name(dataset,"signal");
	if(signal_id == -1) return false;
	return true;
}

bool ReadAttributeInDataset(hid_t dataset_id, std::string attributeName, int *value);
std::map<int,std::string> HDFAxisDataInGroup;
QMutex HDFAxisDataMutex;
herr_t AxisDatasetsInGroup(hid_t loc_id, const char *name, const H5L_info_t *info, void *opdata)
{
	herr_t          status;
	H5O_info_t      infobuf;
	status = H5Oget_info_by_name (loc_id, name, &infobuf, H5P_DEFAULT);
	if (name[0] == '.')         /* Root group, do not print '.' */
		printf ("  (Group)\n");
	else
		switch (infobuf.type) {
		case H5O_TYPE_GROUP:
			//				HDFTreeDataMap.insert(std::pair<std::string,std::string>(name,name));
			printf ("%s  (Group)\n", name);
			break;
		case H5O_TYPE_DATASET:
			{
				printf ("%s  (Group)\n", name);
				hid_t dataset_id = H5Dopen(loc_id,name,H5P_DEFAULT);
				hid_t dataspace = H5Dget_space(dataset_id);
				int rank = H5Sget_simple_extent_ndims(dataspace);
				if(rank==1){ //could be an axis
					int axis_id, primary;
					bool success;
					success = ReadAttributeInDataset(dataset_id,"axis",&axis_id);
					if(success)
					{
						success = ReadAttributeInDataset(dataset_id,"primary",&primary);

						if(primary==1 && success)
							HDFAxisDataInGroup[axis_id]=name;
						else
							HDFAxisDataInGroup.insert(std::pair<int,std::string>(axis_id,name));
					}
				}

				H5Sclose(dataspace);
				H5Dclose(dataset_id);
			}
			break;
		case H5O_TYPE_NAMED_DATATYPE:
			//				HDFTreeDataMap.insert(std::pair<std::string,std::string>(name,name));
			printf ("%s  (Datatype)\n", name);
			break;
		default:
			//				HDFTreeDataMap.insert(std::pair<std::string,std::string>(name,name));
			printf ("%s  (Unknown)\n", name);
	}
	return 0;
}

bool ReadAttributeInDataset(hid_t dataset_id, std::string attributeName, int *value)
{
	bool success = false;
	hid_t attribute_id = H5Aopen_name(dataset_id,attributeName.c_str());
	if(attribute_id != -1) { //Found a attribute
		hid_t attribute_space = H5Aget_space(attribute_id);
		hid_t attribute_type = H5Aget_type(attribute_id);
		hid_t atype  = H5Tcopy(H5T_C_S1);
		if(H5Tis_variable_str(attribute_type) && H5T_STRING == H5Tget_class(attribute_type)) ///HEY we need to deal with variable string seperately
		{
			H5Tset_size(atype, H5T_VARIABLE);
			char *string_out=NULL;
			H5Aread(attribute_id, atype, &string_out);		
			std::cout<<"Read string "<<string_out<<std::endl;
			//Now add that to the Axis list
			*(value) = atoi(string_out);
			success=true;
			//free(string_out); //TODO:: Memory leak here need to check how to free hdf allocated memory
		} else if( H5T_STRING == H5Tget_class(attribute_type) && !H5Tis_variable_str(attribute_type))
		{
			H5Tset_size(atype, H5Tget_size(attribute_type)+1);
			char *string_out= (char*) malloc(sizeof(char)*(H5Tget_size(attribute_type)+1));
			H5Aread(attribute_id, atype, string_out);		
			std::cout<<"Read string "<<string_out<<std::endl;
			//Now add that to the Axis list
			*(value) = atoi(string_out);
			free(string_out);
			success=true;
		} else if(H5T_INTEGER == H5Tget_class(attribute_type))
		{
			H5Aread(attribute_id,attribute_type,value);
			success=true;			
		}
		H5Tclose(atype);
		H5Aclose(attribute_id);
		return success;
	}
}

bool CCPiNexusReader::ReadAxisData(std::string datasetPath, int ndims, hsize_t *dims,void** axisData)
{
	//strip the first level to go one level up in the tree
	size_t pos = datasetPath.find_last_of("/");
	if(pos==std::string::npos)
		return false; //Did find the correct dataset
	std::string datasetParent = datasetPath.substr(0,pos);

	hid_t group_id = H5Gopen(FileId,datasetParent.c_str(),H5P_DEFAULT);
	std::map<int,std::string> axisList;

	HDFAxisDataMutex.lock(); // lock the mutext so no one can change when its being used
	HDFAxisDataInGroup.clear();
	H5Literate(group_id, H5_INDEX_NAME, H5_ITER_INC, NULL, AxisDatasetsInGroup,  NULL);
	for(std::map<int,std::string>::iterator itr=HDFAxisDataInGroup.begin();itr!=HDFAxisDataInGroup.end();itr++)
		axisList[itr->first] = itr->second;
	HDFAxisDataMutex.unlock();

	bool success = true;
	if(axisList.size()<ndims)
		success=false;
	//Now we have axis information lets read the data.
	//AllocateMemory for the axis data
	hsize_t totalSize = 0;
	for(int id =0; id<ndims;id++) totalSize += dims[id];
	*axisData = AllocateMemory(H5T_NATIVE_DOUBLE, totalSize);

	for(std::map<int,std::string>::iterator itr=axisList.begin(); itr != axisList.end();itr++)
	{
		std::cout<<"Axis Id:"<<itr->first-1<<" Axis Name:"<<itr->second<<std::endl;
		bool status = ReadOneAxisDataAndSetInOutput(itr->first-1,datasetParent+"/"+itr->second, ndims, dims, (double*)*axisData); //Axis index start from 1
		if(!status)
		{
			success=false;
			break;
		}
	}
	H5Gclose(group_id);
	return success;
}

bool CCPiNexusReader::ReadOneAxisDataAndSetInOutput(int axisId, std::string datasetPath, int axisNDims, hsize_t *axisDims, double* axisData)
{
	void* data;
	int ndims;
	int *dims;
	DATATYPE dataType;
	ReadCompleteData(datasetPath, &data, &ndims, &dims, &dataType,NULL);
	if(axisDims[axisId]!=dims[0]||ndims!=1){
		std::cout<<"Might be a problem with the dimension"<<std::endl;
		delete[] dims;
		return false;
	}
	std::cout<<"Reading AxisId "<<axisId<<std::endl;
	//Copy the data 	
	hsize_t offset = 0;
	for(int i=0;i<axisId;i++)
		offset += axisDims[i];
	CopyAndDeleteData(dataType, dims[0],data,axisData+offset);
	delete[] dims;
	return true;
}

bool CCPiNexusReader::CopyAndDeleteData(DATATYPE dataType,int num, void* data, double *axisData)
{
	switch(dataType)
	{
	case CHAR:
		CopyAndDeleteDataTemplate(num,(char*) data,axisData);
		break;
	case UCHAR:
		CopyAndDeleteDataTemplate(num,(unsigned char*) data,axisData);
		break;
	case SHORT:
		CopyAndDeleteDataTemplate(num,(short*) data,axisData);
		break;
	case USHORT:
		CopyAndDeleteDataTemplate(num,(unsigned short*) data,axisData);
		break;
	case INT:
		CopyAndDeleteDataTemplate(num,(int*) data,axisData);
		break;
	case UINT:
		CopyAndDeleteDataTemplate(num,(unsigned int*) data,axisData);
		break;
	case LONG:
		CopyAndDeleteDataTemplate(num,(long*) data,axisData);
		break;
	case ULONG:
		CopyAndDeleteDataTemplate(num,(unsigned long*) data,axisData);
		break;
	case LLONG:
		CopyAndDeleteDataTemplate(num,(long long*) data,axisData);
		break;
	case ULLONG:
		CopyAndDeleteDataTemplate(num,(unsigned long long*) data,axisData);
		break;
	case FLOAT:
		CopyAndDeleteDataTemplate(num,(float*) data,axisData);
		break;
	case DOUBLE:
		CopyAndDeleteDataTemplate(num,(double*) data,axisData);
		break;
	case LDOUBLE:
		CopyAndDeleteDataTemplate(num,(long double*) data,axisData);
		break;
	default:
		return false;
	}
	return true;
}



template<class T>
void CCPiNexusReader::CopyAndDeleteDataTemplate(int num, T* data, double *axisData)
{
	for(hsize_t id=0;id<num;id++)
	{
		*(axisData+id) = *(((T*)data)+id);
	}
	delete[] ((T*)data);
}

