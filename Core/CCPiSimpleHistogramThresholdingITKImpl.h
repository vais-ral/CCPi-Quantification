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
public:
	// Typedefs for image types
	typedef itk::Image< T, 3 >     ImageType;
	typedef itk::Image< unsigned char, 3 >     OutputImageType;

	// Typedef for histogram generator
	typedef itk::Statistics::ScalarImageToHistogramGenerator< ImageType > HistogramGeneratorType;
	// Typedef for histogram
	typedef typename HistogramGeneratorType::HistogramType  HistogramType;
	// Typedef for thresholding filter
	typedef itk::BinaryThresholdImageFilter< ImageType, OutputImageType > ThresholdFilterType;

	/**
	 * Constructor which takes in input image data (raw pointer) and its structure information
	 * and min,max instensity to be used in thresholding
	 * @inputImage : input image raw pointer
	 * @volumeDims : 3D dimensions of the input image
	 * @voxelSize  : voxel sizes
	 * @origin	   : origin coordinates of the iamge
	 * @minIntensity : minimum intensity to be used by the thresholding
	 * @maxIntensiyt : maximum intensity to be used by the thresholding
	 */
	CCPiSimpleHistogramThresholdingITKImpl(T* inputImage, const int *volumeDims, const float *voxelSize, const float *origin,float minIntensity, float maxIntensity);
	/**
	 * Here is the entry point for the computation 
	 */
	void Compute();
	/** 
	 * @return output image after thresholding
	 */
	typename OutputImageType::Pointer	GetOutputImage();
	/**
	 * @return returns the two peak values in the histogram
	 */
	std::vector<float>					GetPeaks();

private:
	typename ImageType::Pointer		 Image;
	float							 MinPixelIntensity;
	float							 MaxPixelIntensity;
	OutputImageType::Pointer		 OutputImage;
	std::vector<float>				 Peaks;

	/**
	 * creates a ITK image object from the raw image data
	 * @inputImage : input image raw data pointer
	 * @volumeDims : image dimensions
	 * @voxelSize  : voxel sizes in image
	 * @origin	   : origin values of the image
	 * @return ITK image object 
	 */
	typename ImageType::Pointer					CreateImageFromRawPointer(T* inputImage, const int *volumeDims, const float *voxelSize, const float *origin);
	/**
	 * Computes Histogram from image
	 * @inputImage : input image to which histogram has to be calculated
	 * @minIntensity : lower threshold for the histogram
	 * @maxIntensity : upper threshold for the histogram
	 * @return histogram.
	 */
	typename HistogramGeneratorType::Pointer	ComputeHistogram(typename ImageType::Pointer inputImage, float minIntensity, float maxIntensity);
	/**
	 * Computes and returns two peaks in the histogram this is a very basic searching which
	 * finds the two consecutive peaks.
	 * @histogram : histogram
	 * @minPixelIntensity: minimum pixel intensity
	 * @maxPixelIntensity: maximum pixel intensity
	 * @return two peak intensity values
	 */
	std::vector<float>							FindTwoPeaksInHistogram(const typename HistogramType *histogram,float minPixelIntensity, float maxPixelIntensity);
	/**
	 * Thresholds the input image with a threshold value 
	 * @inputImage : input image
	 * @thresholdValue : threshold value
	 * @return : return result image after thesholding
	 */
	typename OutputImageType::Pointer			ThresholdAnImage(typename ImageType::Pointer inputImage, double thresholdValue);
};


template <class T>
CCPiSimpleHistogramThresholdingITKImpl<T>::CCPiSimpleHistogramThresholdingITKImpl(T* inputImage, const int *volumeDims, const float *voxelSize, const float *origin,float minIntensity, float maxIntensity)
{
	Image = CreateImageFromRawPointer(inputImage, volumeDims, voxelSize, origin);
	MinPixelIntensity = minIntensity;
	MaxPixelIntensity = maxIntensity;
}

template <class T>
typename CCPiSimpleHistogramThresholdingITKImpl<T>::ImageType::Pointer 
	CCPiSimpleHistogramThresholdingITKImpl<T>::CreateImageFromRawPointer(T* inputImage, const int *volumeDims, const float *voxelSize, const float *origin)
{
	// Typedef for importing image data from VolView data
	typedef itk::ImportImageFilter< T, 3 > ImportFilterType;
	// Create import filter
	typename ImportFilterType::Pointer importFilter = ImportFilterType::New();

	typename ImportFilterType::IndexType start;
	start.Fill(0);
	typename ImportFilterType::SizeType size;
	size[0] = volumeDims[0];
	size[1] = volumeDims[1];
	size[2] = volumeDims[2];
	typename ImportFilterType::RegionType region;
	region.SetIndex(start);
	region.SetSize(size);
	importFilter->SetSpacing(voxelSize);
	importFilter->SetOrigin(origin);
	importFilter->SetRegion(region);

	// Size of image
	long imgSize = ((long)volumeDims[0])*volumeDims[1]*volumeDims[2];

	// Set data for full image import filter
	importFilter->SetImportPointer(inputImage, imgSize, false);
	importFilter->Update();

	return importFilter->GetOutput();
}

template <class T>
void CCPiSimpleHistogramThresholdingITKImpl<T>::Compute()
{

	CCPiSimpleHistogramThresholdingITKImpl<T>::HistogramGeneratorType::Pointer histogramGenerator = ComputeHistogram(Image, MinPixelIntensity, MaxPixelIntensity);
	const typename CCPiSimpleHistogramThresholdingITKImpl<T>::HistogramType *histogram = histogramGenerator->GetOutput();
	Peaks = FindTwoPeaksInHistogram(histogram,MinPixelIntensity, MaxPixelIntensity);
	double thresholdValue = (Peaks[0]+Peaks[1])/2.0;
	OutputImage = ThresholdAnImage(Image, thresholdValue);
}

template <class T>
typename CCPiSimpleHistogramThresholdingITKImpl<T>::HistogramGeneratorType::Pointer
	CCPiSimpleHistogramThresholdingITKImpl<T>::ComputeHistogram(typename CCPiSimpleHistogramThresholdingITKImpl<T>::ImageType::Pointer inputImage, float minIntensity, float maxIntensity)
{
	typename CCPiSimpleHistogramThresholdingITKImpl<T>::HistogramGeneratorType::Pointer histogramGenerator = HistogramGeneratorType::New();
	// Output of import filter is input to distance map filter
	histogramGenerator->SetInput( inputImage );

	// Number of bins = range of input image so 1 bin covers 1 intensity
	const unsigned int numBins = maxIntensity-minIntensity;
	histogramGenerator->SetNumberOfBins(numBins);
	histogramGenerator->SetMarginalScale( 10.0 );

	histogramGenerator->SetHistogramMin(minIntensity-0.5);
	histogramGenerator->SetHistogramMax(maxIntensity+0.5);

	// Generate the histogram
	histogramGenerator->Compute();  
	return histogramGenerator;
}


template <class T>
std::vector<float> CCPiSimpleHistogramThresholdingITKImpl<T>::FindTwoPeaksInHistogram(const typename CCPiSimpleHistogramThresholdingITKImpl<T>::HistogramType *histogram,float minPixelIntensity, float maxPixelIntensity)
{
	const unsigned int histogramSize = histogram->Size();

	// Want to find first and last peak so have array size 2 to store bin numbers
	bool goingUp = false;
	std::vector<float> peaks(2); // A copy will be returned
	int peakId = 0;
	std::vector<int> highCount(2,0); // Two highest peaks of histogram in ascending order
	std::vector<float> highs(2); // Associated grey levels of peaks in highCount

	for( unsigned int bin=0; bin < histogramSize-1; bin++ ) {
		// Check for a peak
		if ( histogram->GetFrequency( bin+1, 0 ) < histogram->GetFrequency( bin, 0 ) ) {
			if (goingUp == true ) {
				peaks[peakId] = bin;
				peakId = 1;
				goingUp = false;
			}
		}
		else {
			goingUp = true;
		}
		// Check if this value is bigger than the 2 current highest
		if (histogram->GetFrequency(bin, 0) > highCount[0]) {
			if (histogram->GetFrequency(bin, 0) > highCount[1]) {
				highs[0] = highs[1];
				highs[1] = bin;
				highCount[1] = histogram->GetFrequency(bin, 0);
			}
			else {
				highs[0] = bin;
				highCount[0] = histogram->GetFrequency(bin, 0);
			}
		}
	}

	// Transform peaks from bin number to actual range
	peaks[0] += minPixelIntensity;
	peaks[1] += minPixelIntensity;
	// Transform highs from bin number to actual range
	highs[0] += minPixelIntensity;
	highs[1] += minPixelIntensity;

	return peaks;
}

template<class T>
typename CCPiSimpleHistogramThresholdingITKImpl<T>::OutputImageType::Pointer
	CCPiSimpleHistogramThresholdingITKImpl<T>::ThresholdAnImage(typename CCPiSimpleHistogramThresholdingITKImpl<T>::ImageType::Pointer image, double thresholdValue)
{
	typename CCPiSimpleHistogramThresholdingITKImpl<T>::ThresholdFilterType::Pointer thresholdFilter = CCPiSimpleHistogramThresholdingITKImpl<T>::ThresholdFilterType::New();
	thresholdFilter->SetInput( image );
	thresholdFilter->SetOutsideValue(1); // Value >= threshold
	thresholdFilter->SetInsideValue(0);  // Value < threshold
	thresholdFilter->SetLowerThreshold(0.0);
	thresholdFilter->SetUpperThreshold(thresholdValue);  
	thresholdFilter->Update();

	return thresholdFilter->GetOutput();
}


template<class T>
typename CCPiSimpleHistogramThresholdingITKImpl<T>::OutputImageType::Pointer	CCPiSimpleHistogramThresholdingITKImpl<T>::GetOutputImage()
{
	return OutputImage;
}


template<class T>
std::vector<float> 	CCPiSimpleHistogramThresholdingITKImpl<T>::GetPeaks()
{
	return Peaks;
}
#endif
