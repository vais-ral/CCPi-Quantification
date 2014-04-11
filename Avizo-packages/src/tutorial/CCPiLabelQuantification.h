/** 
 * @file Header file for module to perform quantification calculations on a
 * labelled volume.
 * @author David Worth, Scientifc Computing Department, STFC
 * @date April 2013
 */

#ifndef CCPILABELQUANTIFICATION_H
#define CCPILABELQUANTIFICATION_H

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxPortFilename.h>
#include <hxcore/HxPortFloatTextN.h>
#include <hxcore/HxPortIntSlider.h>

#include "tutorialAPI.h"

class HxUniformScalarField3; // Forward declaration

/** 
 * A module to calculate several characteristics from a labelled image. 
 * It uses the class CCPiQuantification3D to do all the work.
 *
 * The following characteristics are calculated:
 * \li Volume by voxel counts
 * \li Equivalent sphere diameter by voxel counts
 * \li Bounding box diagonal
 * \li Principal Component Analysis
 * \li Ellipsoid fitting by PCA
 * \li Equivalent circle diameter by PCA
 * \li Isosurface by marching cube
 * \li Surface area
 * \li Surface volume
 * \li Equivalent sphere diameter from surface volume
 * \li Sphercity
 * \li Normalised surface area to volume ratio (Radius*Sa/Vol)
 * 
 * It has user interface to set the minimum size of features to include in
 * the calculation.
 */
class TUTORIAL_API CCPiLabelQuantification : public HxCompModule
{
    HX_HEADER(CCPiLabelQuantification);

  public:

    CCPiLabelQuantification();
    ~CCPiLabelQuantification();

    HxPortDoIt portAction;
    
    /** Port providing int text input field with slider for minimum feature size */
    HxPortIntSlider portMinSize;
    /** Port providing float valued input field for voxel size */
    HxPortFloatTextN portVoxelSize;
    /** Port to allow specification of results output file */
    HxPortFilename portOutputFile;

    /** Perform the calculation in this module. Called by Avizo. */
    virtual void compute();
    
  private:
  
    /// The input data object
    HxUniformScalarField3 *m_field;

    /**
     * Template function to run calculation
     * 
     * @return @TODO
     */
    template<class IT> void runQuantification(IT *data, int vtkDataType);
};

#endif // CCPILABELQUANTIFICATION_H
