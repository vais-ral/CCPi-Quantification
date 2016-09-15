#include <QApplication>
#include <hxcore/HxMessage.h>
#include <CCPiNexusReader.h>
#include "CCPiAvizoNexusCmdReader.h"


HX_INIT_CLASS(CCPiAvizoNexusCmdReader, HxData)

extern void CCPiRegisterUniformDataset(std::string name, void *data, int ndims, int *dims, CCPiNexusReader::DATATYPE dataType);
extern void CCPiRegisterRegularDataset(std::string name, void *data, int ndims, int *dims, CCPiNexusReader::DATATYPE dataType, double* axisData);

CCPiAvizoNexusCmdReader::CCPiAvizoNexusCmdReader():
portAction(this, "action", QApplication::translate("CCPiAvizoNexusCmdReader", "Action")),
portFilename(this, "nexusFilePath", QApplication::translate("CCPiAvizoNexusCmdReader", "NexusFilePath")),
portDatasetPath(this, "datasetPath", QApplication::translate("CCPiAvizoNexusCmdReader", "Dataset Path")),
portDimensionOne(this, "dimensionOne", QApplication::translate("CCPiAvizoNexusCmdReader", "Dimension 1 (Start,Stride,Num)"), 3),
portDimensionTwo(this, "dimensionTwo", QApplication::translate("CCPiAvizoNexusCmdReader", "Dimension 2 (Start,Stride,Num)"), 3),
portDimensionThree(this, "dimensionThree", QApplication::translate("CCPiAvizoNexusCmdReader", "Dimension 3 (Start,Stride,Num)"),3)
{


}

CCPiAvizoNexusCmdReader::~CCPiAvizoNexusCmdReader()
{
}

void CCPiAvizoNexusCmdReader::compute()
{
	theMsg->stream() << " Started reading using Nexus Reader" << std::endl;
	if (!portAction.wasHit())
	{
		if (!portFilename.getFilename().isEmpty())
		{ //TODO: Check if the extension is nxs
			CCPiNexusReader reader(portFilename.getFilename().toUtf8().data());
			std::vector<std::string> variableNames = reader.GetVariableNames();
			portDatasetPath.setNum((int)variableNames.size());
			for (int i = 0; i < variableNames.size(); i++){
				portDatasetPath.setLabel(i, QString(variableNames[i].c_str()));
			}
			if (!portDatasetPath.getLabel().isEmpty())
			{
				size_t nDims = reader.GetDataNumberOfDimensions(portDatasetPath.getLabel(portDatasetPath.getValue()).toUtf8().data());
				//TODO:Check if the dimension is more than three
				int *dims = new int[nDims];
				reader.GetDataDimensions(portDatasetPath.getLabel(portDatasetPath.getValue()).toUtf8().data(), dims);
				resetDimensionWidgets();
				if (nDims >= 1)
				{
					portDimensionOne.setValue(2, dims[0]);
					portDimensionOne.setValue(1, 1);
					portDimensionOne.setValue(0, 0);
				}
				if (nDims >= 2)
				{
					portDimensionTwo.setValue(2, dims[1]);
					portDimensionTwo.setValue(1, 1);
					portDimensionTwo.setValue(0, 0);
				}
				if (nDims >= 3)
				{
					portDimensionThree.setValue(2, dims[2]);
					portDimensionThree.setValue(1, 1);
					portDimensionThree.setValue(0, 0);
				}
				delete[] dims;
			}
		}
		return;
	}
	else {
		if (portFilename.getFilename().isEmpty() || portDatasetPath.getLabel().isEmpty() || portDimensionOne.getValue(2) == 0 ||
			portDimensionTwo.getValue(2) == 0)
		{
			theMsg->stream() << "Please select the filename, dataset paths and dimensions properly" << std::endl;
			return;
		}
		std::string datasetname = portDatasetPath.getLabel(portDatasetPath.getValue()).toUtf8().data();
		CCPiNexusReader reader(portFilename.getFilename().toUtf8().data());
		int ndims = reader.GetDataNumberOfDimensions(datasetname);
		int *dims = new int[ndims];
		reader.GetDataDimensions(datasetname, dims);
		CCPiNexusReader::DATATYPE dataType;
		void* data;
		double* axisData = NULL;
		reader.ReadPartialData(datasetname, getStartValues(), getCountValues(), getStrideValues(),
			&data, &dataType, &axisData, false);
		for (int i = 0; i<ndims; i++)
			dims[i] = getCountValues()[i];

		if (axisData == NULL)
		{
			CCPiRegisterUniformDataset(datasetname, data, ndims, dims, dataType);
		}
		else{
			CCPiRegisterRegularDataset(datasetname, data, ndims, dims, dataType, axisData);
		}
		delete[] dims;
	}

}

void CCPiAvizoNexusCmdReader::resetDimensionWidgets()
{

	for (int i = 0; i < 3;i++)
	{
		portDimensionOne.setValue(i, 0);
		portDimensionTwo.setValue(i, 0);
		portDimensionThree.setValue(i, 0);
	}
}

std::vector<long> CCPiAvizoNexusCmdReader::getStartValues()
{
	std::vector<long> result;
	result.push_back(portDimensionOne.getValue(0));
	result.push_back(portDimensionTwo.getValue(0));
	result.push_back(portDimensionThree.getValue(0));
	return result;
}

std::vector<long> CCPiAvizoNexusCmdReader::getStrideValues()
{
	std::vector<long> result;
	result.push_back(portDimensionOne.getValue(1));
	result.push_back(portDimensionTwo.getValue(1));
	result.push_back(portDimensionThree.getValue(1));
	return result;
}

std::vector<long> CCPiAvizoNexusCmdReader::getCountValues()
{
	std::vector<long> result;
	result.push_back(portDimensionOne.getValue(2));
	result.push_back(portDimensionTwo.getValue(2));
	result.push_back(portDimensionThree.getValue(2));
	return result;
}