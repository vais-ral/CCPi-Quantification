/****************************************************************************
** This is test program for Nexus Widget
* Author: Mr. Srikanth Nagella
* Date: 24.07.2014
****************************************************************************/

#include <QtGui>

#include "CCPiNexusTreeModel.h"
#include "CCPiNexusWidgetDialog.h"
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
		reader.ReadCompleteData(*(itr), &data, &ndims, &dims, &dataType, &axisData);
		//std::cout<<"Number of Dimensions"<<ndims<<std::endl;
		//for(int i=0;i<dims[0]+dims[1]+dims[2];i++)
		//	std::cout<<axisData[i]<<std::endl;
	}
	std::cout<<"Completed Reading the Datasets"<<std::endl;
    return app.exec();
}
