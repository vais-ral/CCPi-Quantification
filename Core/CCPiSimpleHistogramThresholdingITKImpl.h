/**
 * This plugin performs a simple thresholding of an image. It is designed for
 * images which have histogram with two distinct peaks and sets the threshold
 * value to be the average of the two peak intensity values.
 * It will not work well for other images!
 *
 * Author: Mr. Srikanth Nagella
 * Date  : 16.05.2014
 */

#ifndef CCPISIMPLEHISTOGRAMTHRESHOLDINGITKIMPL_H
#define CCPISIMPLEHISTOGRAMTHRESHOLDINGITKIMPL_H

#include "itkImage.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkBinaryThresholdImageFilter.h"

#include <vector>

template <class T>
class CCPiSimpleHistogramThresholdingITKImpl
{
  // Typedefs for image types
  typedef itk::Image< T, 3 >     ImageType;
  typedef itk::Image< unsigned char, 3 >     OutputImageType;

  // Typedef for histogram generator
  typedef itk::Statistics::ScalarImageToHistogramGenerator< ImageType > HistogramGeneratorType;
  // Typedef for histogram
  typedef typename HistogramGeneratorType::HistogramType  HistogramType;
  // Typedef for thresholding filter
  typedef itk::BinaryThresholdImageFilter< ImageType, OutputImageType > ThresholdFilterType;
	

public:
	CCPiSimpleHistogramThresholdingITKImpl(T* inputImage, const int *volumeDims, const float *voxelSize, const float *origin,float minIntensity, float maxIntensity);
	void Compute();
	
private:
	typename ImageType::Pointer		 Image;
	float							 MinPixelIntensity;
	float							 MaxPixelIntensity;
	OutputImageType::Pointer		 OutputImage;
	float							 Peaks[2];

	typename ImageType::Pointer			CreateImageFromRawPointer(T* inputImage, const int *volumeDims, const float *voxelSize, const float *origin);
	typename HistogramType::Pointer		ComputeHistogram(typename ImageType::Pointer inputImage, float minIntensity, float maxIntensity);
	std::vector<float>					FindTwoPeaksInHistogram(typename HistogramType::Pointer histogram);
	typename OutputImageType::Pointer	ThresholdAnImage(typename ImageType::Pointer inputImage, double thresholdValue);
};

#endif
