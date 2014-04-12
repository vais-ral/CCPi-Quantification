/** 
 * @file Implementation file for module using ITK k-means filter.
 * @author David Worth, Scientifc Computing Department, STFC
 * @date April 2013
 */

#include <QApplication>

// Classes from ITK
#include "itkImage.h"
#include "itkImportImageFilter.h"
#include "itkScalarImageKmeansImageFilter.h"


#include <hxcore/HxMessage.h>               // For output in Avizo console
#include <hxcore/HxWorkArea.h>              // Busy-cursor and progress bar
#include <hxfield/HxUniformScalarField3.h>  // Class representing 3D images

#include "CCPiKMeansFilterITK.h"

HX_INIT_CLASS(CCPiKMeansFilterITK,HxCompModule)

CCPiKMeansFilterITK::CCPiKMeansFilterITK() :
    HxCompModule(HxUniformScalarField3::getClassTypeId()),
    portAction(this,"action",QApplication::translate("CCPiKMeansFilterITK", "Action")),
    portClasses(this,"classes",QApplication::translate("CCPiKMeansFilterITK", "Number of classes"))
{
    portAction.setLabel(0,"DoIt");
    portClasses.setMinMax(2,255);
}

CCPiKMeansFilterITK::~CCPiKMeansFilterITK()
{
}

/**
 * Method to do the calculation and show the results.
 * Uses a template function to treat different data types.
 */
void CCPiKMeansFilterITK::compute()
{
    // Check whether the action port button was clicked
    if (!portAction.wasHit()) return;

    // Access the input data object. The member portData, which is of type
    // HxConnection, is inherited from HxModule.
    HxUniformScalarField3 *field = (HxUniformScalarField3*) portData.source();
    
    // Get the origin of the data (lower left position of bounding box)
    float origin[3];
    origin[0] = field->bbox()[0];
    origin[1] = field->bbox()[2];
    origin[2] = field->bbox()[4];
    
    // Create an output with same size as input. Data type must be unsigned char
    // for ITK filter we are using.
    HxUniformScalarField3 *output = createOutput(field);
    
    // Output shall have same bounding box as input
    output->coords()->setBoundingBox(field->bbox());
    
    // Turn application into busy state but don't activiate the Stop button
    theWorkArea->startWorkingNoStop(
        QApplication::translate("CCPiKMeansFilterITK", "Computing k-means thresholds"));

    std::vector<float> estimatedMeans;
    switch (field->primType()) {
    
      case McPrimType::mc_uint8:
        estimatedMeans = runKMeansFilter((unsigned char*)field->lattice.dataPtr(),
                                         field->lattice.dims(),
                                         field->getVoxelSize().getValue(), origin, 
                                         (unsigned char*)output->lattice.dataPtr());
        break;
      case McPrimType::mc_int16:
         estimatedMeans = runKMeansFilter((short*)field->lattice.dataPtr(),
                                         field->lattice.dims(),
                                         field->getVoxelSize().getValue(), origin, 
                                         (unsigned char*)output->lattice.dataPtr());
        break;
      case McPrimType::mc_uint16:
        estimatedMeans = runKMeansFilter((unsigned short*)field->lattice.dataPtr(),
                                         field->lattice.dims(),
                                         field->getVoxelSize().getValue(), origin, 
                                         (unsigned char*)output->lattice.dataPtr());
        break;
      case McPrimType::mc_int32:
        estimatedMeans = runKMeansFilter((int*)field->lattice.dataPtr(),
                                         field->lattice.dims(),
                                         field->getVoxelSize().getValue(), origin, 
                                         (unsigned char*)output->lattice.dataPtr());
        break;
      case McPrimType::mc_float:
        estimatedMeans = runKMeansFilter((float*)field->lattice.dataPtr(),
                                         field->lattice.dims(),
                                         field->getVoxelSize().getValue(), origin, 
                                         (unsigned char*)output->lattice.dataPtr());
        break;
      case McPrimType::mc_double:
        estimatedMeans = runKMeansFilter((double*)field->lattice.dataPtr(),
                                         field->lattice.dims(),
                                         field->getVoxelSize().getValue(), origin, 
                                         (unsigned char*)output->lattice.dataPtr());
         break;
      default:
        theMsg->stream() << "The input data has a datatype that this module" <<
                            " cannot handle" << std::endl;
   }
    
    // Register result - adds data object to project view in fot already present.
    // Also connects object's master port to compute module
    setResult(output);

    // Stop progress bar
    theWorkArea->stopWorking();

    // Finally print the result
    for(unsigned int i = 0 ; i < estimatedMeans.size() ; ++i)
    {
        theMsg->stream() <<  "cluster[" << i << "] ";
        theMsg->stream() << "    estimated mean : " << estimatedMeans[i] << std::endl;
    }
}

/**
 * Template function to run the ITK k-means image filter
 * @param data The raw data from the input image
 * @param dims The dimensions of the image, (i,j,k)
 * @param voxelSize Size of a voxel in the image
 * @param origin    Origin position of the image
 * @param output    Raw data for the output image. Set to NULL if no image
 *                  output required.
 * @return The estimated means of the threshold classes.
 */
template <class IT>
std::vector<float> CCPiKMeansFilterITK::runKMeansFilter(
    IT *data, const int *dims, const float *voxelSize, const float *origin,
    unsigned char *output)
{
    // Typedefs for image types
    typedef itk::Image< IT, 3 >     ImageType;
    typedef itk::Image< IT, 3 >     OutputImageType;
    // Typedef for importing image data from VolView data
    typedef itk::ImportImageFilter< IT, 3 > ImportFilterType;

    // Typedef for filter
    typedef itk::ScalarImageKmeansImageFilter< ImageType > KMeansFilterType;

    // Size of image
    int imgSize = dims[0]*dims[1]*dims[2];

    // Create the filters
    typename ImportFilterType::Pointer importFilter = ImportFilterType::New();
    typename KMeansFilterType::Pointer kmeansFilter = KMeansFilterType::New();

    kmeansFilter->SetUseNonContiguousLabels(true);

    // Output of import filter is input to k-means filter
    kmeansFilter->SetInput( importFilter->GetOutput() );

    // Data required for importing the image
    typename ImportFilterType::IndexType start;
    start.Fill(0);
    typename ImportFilterType::SizeType size;
    size[0] = dims[0];
    size[1] = dims[1];
    size[2] = dims[2];
    typename ImportFilterType::RegionType region;
    region.SetIndex(start);
    region.SetSize(size);
    importFilter->SetSpacing(voxelSize);
    importFilter->SetOrigin(origin);
    importFilter->SetRegion(region);

    // Set data for full image import filter
    importFilter->SetImportPointer(data, imgSize, false);
    
    // Set raw data for output of the filter
    kmeansFilter->GetOutput()->GetPixelContainer()->SetImportPointer(
        output, imgSize, false);


    // Set up classes for k-means analysis
    const unsigned int numberOfClasses = portClasses.getValue();
    for (unsigned int ic = 0; ic < numberOfClasses; ic++) {
        kmeansFilter->AddClassWithInitialMean(255*ic/numberOfClasses);
    }
    
    // Run the filter
    kmeansFilter->Update();
    
    // Retrieve the class means and store them in the return object
    typename KMeansFilterType::ParametersType estimatedMeans = kmeansFilter->GetFinalMeans();
    std::vector<float> means(estimatedMeans.Size());
    for (unsigned int i = 0; i < estimatedMeans.Size(); i++) means[i] = estimatedMeans[i];
    
    return means;

}

/**
 * Create an output with same size as input. The ITK class used for k-means
 * analysis means type must be unsigned char.
 * Check if there is a result that can be re-used. If so check type and size
 * match current input.
 * @param field The input field
 * @return The output to use for showing data in application
 */
HxUniformScalarField3* CCPiKMeansFilterITK::createOutput(HxUniformScalarField3 *field)
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
            dims[2] != outdims[2] || output->primType() != McPrimType::mc_uint8 )
            
            output = NULL;
    }
    
    // If necessary create a new result data set
    if (!output) {
        output = new HxUniformScalarField3(dims, McPrimType::mc_uint8);
        output->composeLabel(field->getName(), "thresholded");
    }
    
    return output;
}
