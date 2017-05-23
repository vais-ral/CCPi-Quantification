
//#include <boost/python.hpp>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include "CCPiAccessibleVolumeInputImages.h"
#include "CCPiAccessibleVolumeITKImpl.h"
#include "CCPiConsoleUserInterface.h"

#include "CCPiLabelQuantificationITKImpl.h"
#include "vtkType.h"

namespace bp = boost::python;
namespace np = boost::python::numpy;


void create_accessible_volume();

//Utilities
template <class T>
int getVtkType()
{
	if(typeid(T) == typeid(unsigned char))
		return VTK_SIGNED_CHAR;
	else if(typeid(T) == typeid(char))
		return VTK_CHAR;
	else if(typeid(T) == typeid(short))
		return VTK_SHORT;
	else if(typeid(T) == typeid(unsigned short))
		return VTK_UNSIGNED_SHORT;
	else if(typeid(T) == typeid(int))
		return VTK_INT;
	else if(typeid(T) == typeid(unsigned int))
		return VTK_UNSIGNED_INT;
	else if(typeid(T) == typeid(long))
		return VTK_LONG;
	else if(typeid(T) == typeid(unsigned long))
		return VTK_UNSIGNED_LONG;
	else if(typeid(T) == typeid(float))
		return VTK_FLOAT;
	else if(typeid(T) == typeid(double))
		return VTK_DOUBLE;
	return VTK_ID_TYPE;
}
//Accessible volume
bp::list AccessibleVolumeInputGetDimensions(CCPiAccessibleVolumeInputImages *avii)
{
	const unsigned int* dims = avii->getDimensions();
	bp::list result;
	result.append<int>(dims[0]);
	result.append<int>(dims[1]);
	result.append<int>(dims[2]);
	return result;
}

bp::list AccessibleVolumeInputGetVoxelSize(CCPiAccessibleVolumeInputImages *avii)
{
	const float* voxels = avii->getVoxelSize();
	bp::list result;
	result.append<float>(voxels[0]);
	result.append<float>(voxels[1]);
	result.append<float>(voxels[2]);
	avii->getScafoldVolume();
	return result;	
}

bp::list AccessibleVolumeInputGetOrigin(CCPiAccessibleVolumeInputImages *avii)
{
	const float* origin = avii->getOrigin();
	bp::list result;
	result.append<float>(origin[0]);
	result.append<float>(origin[1]);
	result.append<float>(origin[2]);
	return result;	
}


CCPiAccessibleVolumeInputImages* AccessibleVolumeInputConstructor(np::ndarray npVoxelSize, np::ndarray npOrigin, np::ndarray npVolumeData, np::ndarray npMaskedVolumeData)
{
	int volumeDims[3];
	float voxelSize[3];
	float origin[3];
	long inputDataDims[3];
	origin[0] = bp::extract<int>(npOrigin[0]);
	origin[1] = bp::extract<int>(npOrigin[1]);
	origin[2] = bp::extract<int>(npOrigin[2]);
	voxelSize[0] = bp::extract<float>(npVoxelSize[0]);
	voxelSize[1] = bp::extract<float>(npVoxelSize[1]);
	voxelSize[2] = bp::extract<float>(npVoxelSize[2]);
	volumeDims[0] = npVolumeData.shape(0);
	volumeDims[1] = npVolumeData.shape(1);
	volumeDims[2] = npVolumeData.shape(2);
	inputDataDims[0] = npVolumeData.shape(0);
	inputDataDims[1] = npVolumeData.shape(1);
	inputDataDims[2] = npVolumeData.shape(2);
	
	// Making a deepcopy of the image data. passing a pointer is crashing in itk.
	// TODO:: check why passing a pointer is causing it to crash.
	CCPiImageDataUnsignedChar* inputData = new CCPiImageDataUnsignedChar((unsigned char *)(npVolumeData.get_data()), inputDataDims, true);
	CCPiImageDataUnsignedChar* inputMaskData = new CCPiImageDataUnsignedChar((unsigned char *)(npMaskedVolumeData.get_data()), inputDataDims, true);
	return new CCPiAccessibleVolumeInputImages(volumeDims, voxelSize, origin, inputData, inputMaskData);
}

CCPiAccessibleVolumeITKImpl* AccessibleVolumeITKImplConstructor(CCPiAccessibleVolumeInputImages *input,	float sphereDiameterRangeMin, float sphereDiameterRangeMax, int numberOfSpheres, float imageResolution)
{
		CCPiConsoleUserInterface *ui = new CCPiConsoleUserInterface();
		return new CCPiAccessibleVolumeITKImpl(input, ui, NULL, sphereDiameterRangeMin, sphereDiameterRangeMax, numberOfSpheres, imageResolution);
}

bp::list AccessibleVolumeITKImplGetAccessibleVolume(CCPiAccessibleVolumeITKImpl *avii)
{
	std::map<double,double> outputAV=avii->GetAccessibleVolume();
	bp::list result;
	for(std::map<double,double>::iterator iterator = outputAV.begin(); iterator != outputAV.end(); iterator++) {
		result.append<bp::tuple>(bp::make_tuple(iterator->first, iterator->second));
	}
	return result;
}

void create_accessible_volume()
{
	bp::class_<CCPiAccessibleVolumeInputImages>("AccessibleVolumeInput", bp::no_init)
		.def("__init__", bp::make_constructor(AccessibleVolumeInputConstructor))
		.def("getScafoldVolume", &CCPiAccessibleVolumeInputImages::getScafoldVolume)
		.def("getScafoldPorosity", &CCPiAccessibleVolumeInputImages::getScafoldPorosity)
		.def("getDimensions", AccessibleVolumeInputGetDimensions)
		.def("getVoxelSize", AccessibleVolumeInputGetVoxelSize)
		.def("getOrigin", AccessibleVolumeInputGetOrigin);

	bp::class_<CCPiAccessibleVolumeITKImpl>("AccessibleVolume", bp::no_init)
		.def("__init__", bp::make_constructor(AccessibleVolumeITKImplConstructor))
		.def("compute", &CCPiAccessibleVolumeITKImpl::Compute)
		.def("getAccessibleVolume",AccessibleVolumeITKImplGetAccessibleVolume);
		
}


template<class T>
CCPiLabelQuantificationITKImpl<T>* LabelQuantificationConstructor(np::ndarray input, np::ndarray npOrigin,np::ndarray npVoxelSize,float min, float max, float minFeatureSize)
{
	float voxelSize[3];
	float origin[3];
	long dims[3];
	dims[0] = input.shape(0);
	dims[1] = input.shape(1);
	dims[2] = input.shape(2);	
	CCPiImageData<T>* image = new CCPiImageData<T>((T*) (input.get_data()), dims, true);
	CCPiConsoleUserInterface *ui = new CCPiConsoleUserInterface();	
	origin[0] = bp::extract<int>(npOrigin[0]);
	origin[1] = bp::extract<int>(npOrigin[1]);
	origin[2] = bp::extract<int>(npOrigin[2]);
	voxelSize[0] = bp::extract<float>(npVoxelSize[0]);
	voxelSize[1] = bp::extract<float>(npVoxelSize[1]);
	voxelSize[2] = bp::extract<float>(npVoxelSize[2]);
	return new CCPiLabelQuantificationITKImpl<T>(image, ui, origin, dims, voxelSize, min, max, minFeatureSize, getVtkType<T>());
}

template<class T>
bp::list LabelQuantificationGetOutput(CCPiLabelQuantificationITKImpl<T> *input)
{
	bp::list result;
	CCPiLabelQuantificationResult* output = input->GetOutput();
	std::vector<std::string> columnNames = output->GetQuantityNames();
	bp::list names;
	for(std::vector<std::string>::iterator i= columnNames.begin();i!=columnNames.end();++i)
	{
		names.append<std::string>(*i);
	}
	result.append<bp::list>(names);
	std::list<int> labelIndexes = output->GetLabelIndexes();
	for(std::list<int>::iterator row_itr=labelIndexes.begin();row_itr!=labelIndexes.end();row_itr++)
	{
		bp::list valuesList;
		for(std::vector<std::string>::iterator column_itr = columnNames.begin();column_itr!=columnNames.end();column_itr++)
		{
			double value = output->GetValue(*column_itr, *row_itr);
			valuesList.append<double>(value);
		}
		result.append<bp::list>(valuesList);
	}	
	return result;
}

template<class T>
void export_label_quantification(std::string name)
{
	bp::class_<CCPiLabelQuantificationITKImpl<T> >(name.c_str(), bp::no_init)
		.def("__init__", bp::make_constructor(LabelQuantificationConstructor<T>))
		.def("compute", &CCPiLabelQuantificationITKImpl<T>::Compute)
		.def("getOutput", LabelQuantificationGetOutput<T>);
} 

void create_label_quantification()
{
	export_label_quantification<unsigned int>("LabelQuantificationUInt");
	export_label_quantification<unsigned short>("LabelQuantificationUShort");
	export_label_quantification<unsigned char>("LabelQuantificationUChar");	
	export_label_quantification<int>("LabelQuantificationInt");
	export_label_quantification<short>("LabelQuantificationShort");
	export_label_quantification<char>("LabelQuantificationChar");		
}

BOOST_PYTHON_MODULE(quantification)
{
	np::initialize();
	create_accessible_volume();	
	create_label_quantification();
}