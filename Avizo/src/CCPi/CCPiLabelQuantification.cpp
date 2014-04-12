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
#include <hxcore/HxWorkArea.h>
#include <hxfield/HxUniformScalarField3.h>

#include "CCPiLabelQuantification.h"



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
    m_field = (HxUniformScalarField3*) portData.source();
    
    // Turn application into busy state but don't activiate the Stop button
    theWorkArea->startWorkingNoStop(
        QApplication::translate("CCPiLabelQuantification", "Computing quantites for the image"));

    // Call the tempalted function to do calculation for integer types
    switch (m_field->primType()) {
    
      case McPrimType::mc_uint8:
        runQuantification((unsigned char*)m_field->lattice.dataPtr(), 3);
        break;
      case McPrimType::mc_int16:
        runQuantification((short*)m_field->lattice.dataPtr(), 4);
        break;
      case McPrimType::mc_uint16:
        runQuantification((unsigned short*)m_field->lattice.dataPtr(), 5);
        break;
      case McPrimType::mc_int32:
        runQuantification((int*)m_field->lattice.dataPtr(), 6);
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
    origin[0] = m_field->bbox()[0];
    origin[1] = m_field->bbox()[2];
    origin[2] = m_field->bbox()[4];
    
    // Get the min and max values of input data
    float min = 0.0, max = 0.0;
    m_field->getRange(min, max);
    
    // Initialise the controller
    quan3D.Initialise(data, m_field->lattice.dims(), 1,
      vtkDataType, origin, m_field->getVoxelSize().getValue(), 
      min, max);

    quan3D.SetMinFeatureSize(portMinSize.getValue());
    quan3D.CreateVoxelIndexList();

    quan3D.PrepareForQuantification();

    quan3D.PrintSummaryData();  

    quan3D.WriteCSVData(portOutputFile.getValue());

    int totalVoxels = quan3D.GetNumVoxelValues(), n = 0;

    theWorkArea->setProgressValue(0.05);
    theWorkArea->setProgressInfo("Processing...");

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < totalVoxels; i++) {

        CCPiQuantificationWorker<IT> *worker = NULL;

        // Do the real work
        #pragma omp critical(nextworker)
        {
            worker = quan3D.GetNextWorker();
        }
        if (worker != NULL) {

            if (0 == worker->Run()) {
                #pragma omp critical(writefile)
                {
                    worker->WriteCSVData(portOutputFile.getValue());
                }
            }
            delete worker;
        }
        #pragma omp atomic
        n++;
        if (omp_get_thread_num() == 0) {
            theWorkArea->setProgressValue((float)n/(float)totalVoxels);
            theWorkArea->setProgressInfo("Quantification underway");
        }
    }
    theMsg->stream() << "Quantification complete" << std::endl;

    theWorkArea->setProgressValue(1.0);
    theWorkArea->setProgressInfo("Processing Complete");
    
}
