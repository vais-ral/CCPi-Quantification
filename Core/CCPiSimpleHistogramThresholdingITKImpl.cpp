#include "CCPiSimpleHistogramThresholdingITKImpl.h"	

#include "itkImportImageFilter.h"


template <class T>
CCPiSimpleHistogramThresholdingITKImpl<T>::CCPiSimpleHistogramThresholdingITKImpl(T* inputImage, const int *volumeDims, const float *voxelSize, const float *origin,float minIntensity, float maxIntensity)
{
	Image = CreateImageFromRawPointer(inputImage, volumeDims, voxelSize, origin);
	MinIntensity = minIntensity;
	MaxIntensity = maxIntensity;
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
    importFilter->SetImportPointer(data, imgSize, false);
    importFilter->Update();

	return importFilter->GetOutput();
}

template <class T>
void CCPiSimpleHistogramThresholdingITKImpl<T>::Compute()
{

	CCPiSimpleHistogramThresholdingITKImpl<T>::HistogramType::Pointer histogram = ComputeHistogram(Image, MinIntensity, MaxIntensity);
	std::vector<float> peaks = FindTwoPeaksInHistogram(histogram);
	double thresholdValue = (peaks[0]+peaks[1])/2.0;
	OutputImage = ThresholdAnImage(image, thresholdValue);
}

template <class T>
typename CCPiSimpleHistogramThresholdingITKImpl<T>::HistogramType::Pointer 
	CCPiSimpleHistogramThresholdingITKImpl<T>::ComputeHistogram(typename CCPiSimpleHistogramThresholdingITKImpl<T>::ImageType::Pointer inputImage, float minIntensity, float maxIntensity)
{
    typename HistogramGeneratorType::Pointer histogramGenerator = HistogramGeneratorType::New();

    // Output of import filter is input to distance map filter
    histogramGenerator->SetInput( inputImage );

	// Number of bins = range of input image so 1 bin covers 1 intensity
    const unsigned int numBins = maxIntensity-minIntensity;
    histogramGenerator->SetNumberOfBins(numBins);
    histogramGenerator->SetMarginalScale( 10.0 );

    histogramGenerator->SetHistogramMin(min-0.5);
    histogramGenerator->SetHistogramMax(max+0.5);

	// Generate the histogram
    histogramGenerator->Compute();  
	return histogramGenerator->GetOutput();
}


template <class T>
std::vector<float> CCPiSimpleHistogramThresholdingITKImpl<T>::FindTwoPeaksInHistogram(typename CCPiSimpleHistogramThresholdingITKImpl<T>::HistogramType::Pointer histogram)
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
	peaks[0] += min;
	peaks[1] += min;
	// Transform highs from bin number to actual range
	highs[0] += min;
	highs[1] += min;

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