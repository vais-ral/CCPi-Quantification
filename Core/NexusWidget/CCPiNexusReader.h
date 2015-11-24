/**
 * This reads the data from nexus files given the filename and the data path
 * Author : Mr. Srikanth Nagella
 */

#ifndef CCPINEXUSREADER_H
#define CCPINEXUSREADER_H

#include "hdf5.h"
#include <string>
#include <vector>

#include "CCPiDefines.h"
class CCPI_EXPORT CCPiNexusReader
{
public:
	enum DATATYPE { UCHAR, CHAR, SHORT, USHORT, INT, UINT, LONG, ULONG, LLONG, ULLONG, FLOAT, DOUBLE, LDOUBLE, UNKNOWN } ;
	CCPiNexusReader(std::string filename);
	~CCPiNexusReader();
	void ReadCompleteData(std::string datasetPath, void** data,int *ndims, int** dims, DATATYPE* dataType, double** axisData,bool isAxis=false);
	void ReadCompleteDataNoAllocation(std::string datasetPath, void* data);
	int  GetDataNumberOfDimensions(std::string datasetPath);
	void GetDataDimensions(std::string datasetPath, int *dims);
	std::vector<std::string> GetVariableNames();
	DATATYPE GetDataType(std::string datasetPath);
private:
	std::string Filename;
	hid_t		FileId;

	void* AllocateMemory(hid_t datatype, int ndims, hsize_t *dims);
	void* AllocateMemory(hid_t datatype, hsize_t totalsize);
	DATATYPE GetDataType(hid_t datatype);
	bool isSignalData(hid_t dataset);
	std::string CCPiNexusReader::getParentDatasetName(std::string datasetPath);
	bool ReadAxisData(std::string datasetPath, int ndims, hsize_t *dims, void** axisData);
	bool ReadOneAxisDataAndSetInOutput(int axisId,std::string datasetPath, int axisNDims, hsize_t *axisDims,  double* axisData);
	bool ReadAxisDataNxsV2(std::string datasetPath, int ndims, hsize_t *dims,void** axisData);
	bool CopyAndDeleteData(DATATYPE dataType,int num, void* data, double *axisData);
	template<class T>
	void CopyAndDeleteDataTemplate(int num, T* data, double *axisData);
};

#endif
