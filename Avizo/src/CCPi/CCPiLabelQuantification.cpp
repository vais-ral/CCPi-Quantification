/** 
 * @file Implementation file for module to perform quantification calculations on a
 * labelled volume.
 * @author David Worth, Scientifc Computing Department, STFC
 * @date April 2013
 */

#include <QApplication>

#include "omp.h"

#include "Quan3D.hpp"       // Class that control quantification calculation
#include "QuanWorker.hpp"   // Worker class to does calculation for one label

#include <hxcore/HxMessage.h>
#include <hxcore/internal/HxWorkArea.h>
#include <hxfield/HxUniformScalarField3.h>

#include "CCPiImageData.h"
#include "CCPiAvizoUserInterface.h"
#include "CCPiLabelQuantification.h"
#include "CCPiLabelQuantificationITKImpl.h"
#include "CCPiLabelQuantificationResult.h"



HX_INIT_CLASS(CCPiLabelQuantification,HxCompModule)

CCPiLabelQuantification::CCPiLabelQuantification() :
    HxCompModule(HxUniformScalarField3::getClassTypeId()),
    portAction(this,"action",
        QApplication::translate("CCPiLabelQuantification", "Action")),
    portMinSize(this,"minSize",
        QApplication::translate("CCPiLabelQuantification", "Minimum feature size")),
    portVoxelSize(this,"voxelSize",
        QApplication::translate("CCPiLabelQuantification", "Voxel size"),1),
    portOutputFile(this,"outputFile",
        QApplication::translate("CCPiLabelQuantification", "Results file")),
    m_field(NULL)
{
    portAction.setLabel(0,"DoIt");
    portMinSize.setMinMax(0,4000);
    portMinSize.setValue(100.0);
    portVoxelSize.setValue(1.0);
    portOutputFile.setValue("3D_data.csv");
}

CCPiLabelQuantification::~CCPiLabelQuantification()
{
}

/**
 * Method to do the calculation and show the results.
 * Uses a template function to treat different data types.
 */
void CCPiLabelQuantification::compute()
{
    // Check whether the action port button was clicked
    if (!portAction.wasHit()) return;

    // Access the input data object. The member portData, which is of type
    // HxConnection, is inherited from HxModule.
    m_field = (HxUniformScalarField3*) portData.getSource();
    
    // Turn application into busy state but don't activiate the Stop button
    theWorkArea->startWorkingNoStop(
        QApplication::translate("CCPiLabelQuantification", "Computing quantites for the image"));

    // Call the tempalted function to do calculation for integer types
    switch (m_field->primType()) {
    
      case McPrimType::MC_UINT8:
        runQuantification((unsigned char*)m_field->lattice().dataPtr(), 3);
        break;
      case McPrimType::MC_INT16:
        runQuantification((short*)m_field->lattice().dataPtr(), 4);
        break;
      case McPrimType::MC_UINT16:
        runQuantification((unsigned short*)m_field->lattice().dataPtr(), 5);
        break;
      case McPrimType::MC_INT32:
        runQuantification((int*)m_field->lattice().dataPtr(), 6);
        break;
      default:
        theMsg->stream() << "The input data has a datatype that this module" <<
                            " cannot handle" << std::endl;
   }
    
    // Stop progress bar
    theWorkArea->stopWorking();
}

/**
 * Template function to run calculation
 * @param data          Raw data from the input image
 * @param vtkDataType   Value to indicate IT type to VTK classes.
 *                      Can be found in vtkType.h
 * 
 * @return @TODO
 */
template <class IT>
void CCPiLabelQuantification::runQuantification(IT *data, int vtkDataType)
{
    // Create the controller class for calculations
    CCPiQuantification3D<IT> quan3D;

    theWorkArea->setProgressValue( 0.01 );
    theWorkArea->setProgressInfo("Initialising...");

    // Get the origin of the data (lower left position of bounding box)
    float origin[3];
    origin[0] = m_field->getBoundingBox()[0];
	origin[1] = m_field->getBoundingBox()[2];
	origin[2] = m_field->getBoundingBox()[4];
    float voxelSize[3];
	voxelSize[0] = m_field->getVoxelSize().getValue()[0];
	voxelSize[1] = m_field->getVoxelSize().getValue()[1];
	voxelSize[2] = m_field->getVoxelSize().getValue()[2];
    // Get the min and max values of input data
    float min = 0.0, max = 0.0;
    m_field->getRange(min, max);
    
	long dims[3];
	dims[0]= m_field->lattice().getDims()[0];
	dims[1]= m_field->lattice().getDims()[1];
	dims[2]= m_field->lattice().getDims()[2];
	CCPiImageData<IT>* image = new CCPiImageData<IT>(data, dims, false);
	CCPiAvizoUserInterface *ui = new CCPiAvizoUserInterface();
	CCPiLabelQuantificationITKImpl<IT> quantification(image, (CCPiUserApplicationInterface*)ui, origin, dims, voxelSize,  min, max, portMinSize.getValue(), vtkDataType);
	quantification.Compute();

	//Register spreadsheet volume fraction
	setResult(createSpreadsheetOutput(m_field->getName(),quantification.GetOutput()));

    theWorkArea->setProgressValue(1.0);
    theWorkArea->setProgressInfo("Processing Complete");
    delete image;
	delete ui;
}

/**
 * Creates the volume path as a spreadsheet and puts in the work area
 * @param name prefix for the module
 * @param volpathMap is map of sphere diameter and its volume path
 * @return the output in a spreadsheet.
 */
HxSpreadSheet* CCPiLabelQuantification::createSpreadsheetOutput(std::string prefix,CCPiLabelQuantificationResult* quantResult)
{

/*	HxSpreadSheet *output = new HxSpreadSheet();
//	output->addTable(("Accessible Volume("+prefix+")").c_str());
	output->setLabel(("LabelQuantification("+prefix+")").c_str());
	int tableId = 0;
	std::vector<std::string> columnNames = quantResult->GetQuantityNames();
	std::list<int> labelIndexes = quantResult->GetLabelIndexes();
	for(std::vector<std::string>::iterator itr = columnNames.begin(); itr!=columnNames.end();itr++)
	{
		int columnId = output->addColumn((*itr).c_str(),HxSpreadSheet::Column::FLOAT,tableId);
	}

	output->setNumRows(labelIndexes.size(),tableId);
	output->setCurrentTableIndex(tableId);

	int columnId=0;
	for(std::vector<std::string>::iterator column_itr = columnNames.begin();column_itr!=columnNames.end();column_itr++, columnId++)
	{
		int rowId=0;
		for(std::list<int>::iterator row_itr=labelIndexes.begin();row_itr!=labelIndexes.end();row_itr++, rowId++)
		{
			double value = quantResult->GetValue(*column_itr, *row_itr);
			int columnId = output->findColumn((*column_itr).c_str(), HxSpreadSheet::Column::FLOAT,tableId);
			HxSpreadSheet::Column *column = output->column(columnId,tableId);
			column->setValue(rowId,value);
		}
	}
    return output;
*/
    return NULL;
}
