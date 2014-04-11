/*
 *  Template of a compute module
 */

/** 
 * @file Implementation file for a simple module used for learning about extending Avizo.
 * @author David Worth, Scientifc Computing Department, STFC
 * @date April 2013
 */

#include <QApplication>

#include <hxcore/HxMessage.h>               // For output in Avizo console
#include <hxcore/HxWorkArea.h>              // Busy-cursor and progress bar
#include <hxfield/HxUniformScalarField3.h>  // Class representing 3D images

#include "CCPiComputeThreshold.h"

// Required macro added by XPand Development Wizard
HX_INIT_CLASS(CCPiComputeThreshold,HxCompModule)

/**
 * Default constructor calls superclass with the 3D data class and constructs
 * the action and range ports. The range port has 2 float fields.
 */
CCPiComputeThreshold::CCPiComputeThreshold() :
    HxCompModule(HxUniformScalarField3::getClassTypeId()),
    portAction(this,"action",QApplication::translate("CCPiComputeThreshold", "Action")),
    portRange(this,"range",QApplication::translate("CCPiComputeThreshold", "Range"),2) 
{
    portAction.setLabel(0,QApplication::translate("CCPiComputeThreshold", "DoIt"));
}

/**
 * Detructor. 
 */
CCPiComputeThreshold::~CCPiComputeThreshold()
{
}

/**
 * Method to do the thresholding calculation and show the results
 */
void CCPiComputeThreshold::compute()
{
    // Check whether the action port button was clicked
    if (!portAction.wasHit()) return;
    
    // Access the input data object. The member portData, which is of type
    // HxConnection, is inherited from HxModule.
    HxUniformScalarField3 *field = (HxUniformScalarField3*) portData.source();
    
    // Check whether the input port is connected
    if (!field) return;
    
    // Get the input parameters from the user interface
    float minValue = portRange.getValue(0);
    float maxValue = portRange.getValue(1);
    
    // Access size of data volume
    const int *dims = field->lattice.dims();
    
    // Create an output with same size and primitive data type as input
    HxUniformScalarField3 *output = createOutput(field);
    
    // Output shall have same bounding box as input
    output->coords()->setBoundingBox(field->bbox());
    
    // Turn application into busy state but don't activiate the Stop button
    theWorkArea->startWorkingNoStop(
        QApplication::translate("CCPiComputeThreshold", "Computing threshold"));
    
    // Loop through whole field and count the pixels below min and above max. Set
    // their values in output to zero. Copy other values to output.
    int belowCnt = 0, aboveCnt = 0;
    for (int k = 0; k < dims[2]; k++) {
        // Set progress bar, the argument ranges between 0 and 1
        theWorkArea->setProgressValue( (float)(k+1)/dims[2] );
        for (int j = 0; j < dims[1]; j++) {
            for (int i = 0; i < dims[0]; i++) {
            
                // Function evalReg returns the value at the specific grid node.
                // It implicitly casts the result to float if necessary.
                // Check 
                float value = field->evalReg(i,j,k);
                float newValue = 0;
                if (value < minValue) belowCnt++;
                else if (value > maxValue) aboveCnt++;
                else newValue = value;
                
                output->set(i,j,k,newValue);
                
            }
        }
    }
    
    // Register result - adds data object to project view in fot already present.
    // Also connects object's master port to compute module
    setResult(output);
    
    // Stop progress bar
    theWorkArea->stopWorking();
    
    // Finally print the result
    theMsg->printf("%d voxels < %g, %d voxels > %g\n",
        belowCnt, minValue, aboveCnt, maxValue);
}

/**
 * Create an output with same size and primitive data type as input.
 * Check if there is a result that can be re-used. If so check type, size and
 * datatype match current input.
 * @param field The input field
 * @return The output to use for showing data in application
 */
HxUniformScalarField3* CCPiComputeThreshold::createOutput(HxUniformScalarField3 *field)
{
    // Check if there is a result which we can reuse
    HxUniformScalarField3 *output = (HxUniformScalarField3*) getResult();
    
    // Check for proper type
    if ( output && !output->isOfType(HxUniformScalarField3::getClassTypeId()) )
        output = NULL;
    
    // Check if size and primitive type still match current input
    const int *dims = field->lattice.dims();
    if (output) {
        const int *outdims = output->lattice.dims();
        if ( dims[0] != outdims[0] || dims[1] != outdims[1] ||
            dims[2] != outdims[2] || field->primType() != output->primType() )
            
            output = NULL;
    }
    
    // If necessary create a new result data set
    if (!output) {
        output = new HxUniformScalarField3(dims, field->primType());
        output->composeLabel(field->getName(), "masked");
    }
    
    return output;
}

