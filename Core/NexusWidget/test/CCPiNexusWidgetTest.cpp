/****************************************************************************
** This is test program for Nexus Widget
* Author: Mr. Srikanth Nagella
* Date: 24.07.2014
****************************************************************************/

#include <QtGui>

#include "CCPiNexusTreeModel.h"
#include "CCPiNexusWidgetDialog.h"
#include "CCPiNexusSubsetWidget.h"
#include "CCPiNexusReader.h"
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

	std::string filename(argv[1]);

	CCPiNexusWidgetDialog* nexusDlg=new CCPiNexusWidgetDialog(filename);
	nexusDlg->exec();

	std::cout<<"Started Reading the Datasets"<<std::endl;
	std::vector<std::string> output = nexusDlg->GetSelectedDataSetList();
	for(std::vector<std::string>::iterator itr = output.begin(); itr!=output.end();itr++)
	{
		CCPiNexusReader reader(filename);
		int ndims;
		int *dims;
		CCPiNexusReader::DATATYPE dataType;
		void* data;
		double* axisData;
		ndims = reader.GetDataNumberOfDimensions(*(itr));
		dims = new int[ndims];
		reader.GetDataDimensions(*(itr), dims);

		//create a vector of long dimensions
		std::vector<long> dimVector;
		for(int i=0;i<ndims;i++)
			dimVector.push_back(dims[i]);
		//create a selection widget
		CCPiNexusSubsetWidget *subsetWidget = new CCPiNexusSubsetWidget(0, std::vector<std::string>(), dimVector);
		subsetWidget->exec();
		std::vector<long> newdims = subsetWidget->getSelectedCountValues();
		std::cout<<"selected count"<<newdims[0]<<"X"<<newdims[1]<<"X"<<newdims[2]<<std::endl;
		//working
		//data = new double[newdims[0]*newdims[1]*newdims[2]];
		//reader.ReadPartialDataNoAllocation(*(itr), data, subsetWidget->getSelectedStartValues(), subsetWidget->getSelectedCountValues(), subsetWidget->getSelectedStrideValues());
	    reader.ReadPartialData(*(itr), subsetWidget->getSelectedStartValues(), subsetWidget->getSelectedCountValues(), subsetWidget->getSelectedStrideValues(),
			                    &data, &dataType, &axisData, false); 
		for(int i=0;i<subsetWidget->getSelectedCountValues()[0];i++)
		{
			std::cout<<" "<<axisData[i];
		}
		std::cout<<std::endl;
//		reader.ReadCompleteData(*(itr), &data, &ndims, &dims, &dataType, &axisData);
		//std::cout<<"Number of Dimensions"<<ndims<<std::endl;
		//for(int i=0;i<dims[0]+dims[1]+dims[2];i++)
		//	std::cout<<axisData[i]<<std::endl;
	}
	std::cout<<"Completed Reading the Datasets"<<std::endl;
    return app.exec();
}
