/**
 * This is class to hold Image Data. It will hold image data pointers or it can allocate 
 * memory for image data and copy it from source. Manages the memory for the Algorithms
 *
 * Author: Mr. Srikanth Nagella
 * Date: 08.07.2014
 */

#ifndef CCPIIMAGEDATA_H
#define CCPIIMAGEDATA_H

template <class T>
class CCPiImageData
{

public:
	CCPiImageData(T* data, long dims[3], bool deepCopy);
	~CCPiImageData();
	T* GetImage(){return Image;}
	long* GetDimensions(){return Dimensions;}
private:
	T* Image;
	long Dimensions[3];
	bool isOwner;
};

template<class T>
CCPiImageData<T>::CCPiImageData(T* data, long dims[3], bool deepCopy)
{
	Dimensions[0] = dims[0];
	Dimensions[1] = dims[1];
	Dimensions[2] = dims[2];

	if(deepCopy)
	{
		Image = new T[ Dimensions[0] * Dimensions[1] * Dimensions[2] ];
		for(long index = 0; index < Dimensions[0]*Dimensions[1]*Dimensions[2]; index++)
		{
			Image[index] = data[index];
		}
		isOwner = true;
	}else{
		Image = data;
		isOwner = false;
	}
}

template<class T>
CCPiImageData<T>::~CCPiImageData()
{
	if(isOwner)
		delete[] Image;
}

typedef CCPiImageData<int> CCPiImageDataInt;
typedef CCPiImageData<unsigned char> CCPiImageDataUnsignedChar;
typedef CCPiImageData<char> CCPiImageDataChar;
typedef CCPiImageData<short> CCPiImageDataShort;
typedef CCPiImageData<unsigned short> CCPiImageDataUnsignedShort;
typedef CCPiImageData<unsigned int> CCPiImageDataUnsignedInt;
typedef CCPiImageData<long> CCPiImageDataLong;
typedef CCPiImageData<unsigned long> CCPiImageDataUnsignedLong;
typedef CCPiImageData<float> CCPiImageDataFloat;
typedef CCPiImageData<double> CCPiImageDataDouble;

#endif