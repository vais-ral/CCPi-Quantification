/**
 * This is an Label Quantification ITK/VTK implementation. This is an wrapper which calls
 * CCPiQuantification3D class to perform the quantification. This call will manage the 
 * Threads to process all the data.
 *
 * Author: Mr. Srikanth Nagella
 * Date:   10.07.2014
 */

#ifndef CCPILABELQUANTIFICATIONITKIMPL_H
#define CCPILABELQUANTIFICATIONITKIMPL_H

#include "CCPiImageData.h"
#include "CCPiUserApplicationInterface.h"
#include "CCPiLabelQuantificationResult.h"
#include "Quan3D.hpp"
#include "QuanWorker.hpp"

template<class T>
class CCPiLabelQuantificationITKImpl
{
public:
	CCPiLabelQuantificationITKImpl(CCPiImageData<T>* data, CCPiUserApplicationInterface *userAppInterface, float origin[3], long volumeDims[3], float voxelSize[3], float min,float max,float minFeatureSize, int vtkDataType);
	~CCPiLabelQuantificationITKImpl();
	void Compute();
	CCPiLabelQuantificationResult* GetOutput();
private:
	CCPiImageData<T>* ImageData;
	int   VolumeDims[3];
	float VoxelSize[3];
	float Origin[3];
	float MinDataValue;
	float MaxDataValue;
	float MinFeatureSize;
	int  VtkDataType;
	CCPiLabelQuantificationResult* Result;
	CCPiUserApplicationInterface *UserAppInterface;
    // Create the controller class for calculations
    CCPiQuantification3D<T> quan3D;
};

template <class T>
CCPiLabelQuantificationITKImpl<T>::CCPiLabelQuantificationITKImpl(CCPiImageData<T>* data, CCPiUserApplicationInterface *userAppInterface, float origin[3], long dims[3], float vsize[3], float min,float max, float minFeatureSize, int vtkDataType)
{
	ImageData = data;
	VolumeDims[0] = dims[0];  VolumeDims[1]=dims[1];  VolumeDims[2] = dims[2];
	VoxelSize[0]  = vsize[0]; VoxelSize[1] = vsize[1]; VoxelSize[2] = vsize[2];
	Origin[0] = origin[0];    Origin[1]   = origin[1]; Origin[2] = origin[2];
	MinDataValue = min;
	MaxDataValue = max;
	MinFeatureSize = minFeatureSize;
	VtkDataType = vtkDataType;
	UserAppInterface = userAppInterface;
}

template <class T>
CCPiLabelQuantificationITKImpl<T>::~CCPiLabelQuantificationITKImpl()
{
}

template <class T>
void CCPiLabelQuantificationITKImpl<T>::Compute()
{
	UserAppInterface->SetStatusMessage("Initialising...");
    UserAppInterface->SetProgressValue( 0.01 );

    // Initialise the controller
	quan3D.Initialise(ImageData->GetImage(), VolumeDims, 1,
      VtkDataType, Origin, VoxelSize, 
      MinDataValue, MaxDataValue);

    quan3D.SetMinFeatureSize(MinFeatureSize);
    quan3D.CreateVoxelIndexList();

    quan3D.PrepareForQuantification();

    quan3D.PrintSummaryData();  

    int totalVoxels = quan3D.GetNumVoxelValues(), n = 0;

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < totalVoxels; i++) {

        CCPiQuantificationWorker<T> *worker = NULL;

        // Do the real work
        #pragma omp critical(nextworker)
        {
            worker = quan3D.GetNextWorker();
        }
        if (worker != NULL) {

            if (0 == worker->Run()) {
                #pragma omp critical(writefile)
                {
					quan3D.SetQuantificationResultByWorker(worker->GetId(),worker->GetQuantificationResult());
                }
            }
            delete worker;
        }
        #pragma omp atomic
        n++;
        if (omp_get_thread_num() == 0) {
			UserAppInterface->SetStatusMessage("Quantification underway");
			UserAppInterface->SetProgressValue( (float)n/(float)totalVoxels );
        }
    }
	UserAppInterface->SetStatusMessage("Quantification complete");
	Result = quan3D.GetQuantificationResult();
}

template <class T>
CCPiLabelQuantificationResult* CCPiLabelQuantificationITKImpl<T>::GetOutput()
{
	return Result;
}

#endif