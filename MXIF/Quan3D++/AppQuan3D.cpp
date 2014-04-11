/**
\file AppQuan3D.cpp
Source file for command-line application that runs the 3D quantification plug-in

\author David Worth STFC
\date May 2012
*/
#include <iostream>
#include <string>
#include <stdlib.h>

#include "vtkDataArray.h"
#include "vtkImageData.h"
#include "vtkImageReader2.h"
#include "vtkImageReader2Factory.h"
#include "vtkPointData.h"

#include "Quan3D.hpp"
#include "QuanWorker.hpp"

#define vtkTemplateMacroCase(typeN, type, call)     \
  case typeN: { typedef type VTK_TT; call; }; break
#ifndef vtkTemplateMacro
#define vtkTemplateMacro(call)                                    \
  vtkTemplateMacroCase(VTK_DOUBLE, double, call);                 \
  vtkTemplateMacroCase(VTK_FLOAT, float, call);                   \
  vtkTemplateMacroCase(VTK_LONG, long, call);                     \
  vtkTemplateMacroCase(VTK_UNSIGNED_LONG, unsigned long, call);   \
  vtkTemplateMacroCase(VTK_INT, int, call);                       \
  vtkTemplateMacroCase(VTK_UNSIGNED_INT, unsigned int, call);     \
  vtkTemplateMacroCase(VTK_SHORT, short, call);                   \
  vtkTemplateMacroCase(VTK_UNSIGNED_SHORT, unsigned short, call); \
  vtkTemplateMacroCase(VTK_CHAR, char, call);                     \
  vtkTemplateMacroCase(VTK_UNSIGNED_CHAR, unsigned char, call)
#endif

template <class IT>
void Run(vtkImageData *imageData, IT *)
{

  CCPiQuantification3D<IT> quan3D;

  double range[2];
  imageData->GetScalarRange(range);
  
  quan3D.SetMinFeatureSize(1000);

  quan3D.Initialise((IT*)(imageData->GetScalarPointer()),
                    imageData->GetDimensions(),
                    1,
                    imageData->GetScalarType(),
                    imageData->GetOrigin(),
                    imageData->GetSpacing(),
                    range[0],
                    range[1]);

  quan3D.CreateVoxelIndexList();

  quan3D.PrepareForQuantification();

  quan3D.PrintSummaryData();  

  quan3D.WriteCSVData("3D_data.csv");

  int totalVoxels = quan3D.GetNumVoxelValues(), n = 0;

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
          worker->WriteCSVData("3D_data.csv");
        }
      }
      delete worker;
    }
    #pragma omp atomic
    n++;
    #pragma omp critical(progress)
    {
      int count = 40*n/totalVoxels;
      std::cout << '\r';
      for (int j = 1; j < count+1; j++) {
        std::cout << '*';
      }
      for (int j = count+1; j < 41; j++) {
        std::cout << '.';
      }
      std::cout << "|   " << 100*n/totalVoxels << "%";      
    }
  }


}

/**
Main function for running the 3D quantification plug-in.

Arguments are
filename xmin xmax ymin ymax zmin zmax x_spacing y_spacing z_spacing \
  x_origin y_origin z_origin scalar_type is_little_endian
*/
int main (int argc, char **argv) {

  std::string filename;
  int xmin, xmax, ymin, ymax, zmin, zmax;
  double x_spacing, y_spacing, z_spacing;
  double x_origin, y_origin, z_origin;
  int scalar_type;
  bool is_little_endian;

  if (argc == 1) {
    // Interactive mode
    //cout << "What's the name of the input file? "
    //cin >> filename;
  }
  else if (argc != 16) {
    cout << "Usage:" << endl;
    cout << argv[0] << " filename xmin xmax ymin ymax zmin zmax x_spacing " <<
      "y_spacing z_spacing x_origin y_origin z_origin scalar_type " <<
      "is_little_endian" << endl;
    cout << endl << "Note that min/max values are 0 based" << endl;
    return 1;
  }
  else {
    filename = argv[1];
    xmin = atoi(argv[2]);
    xmax = atoi(argv[3]);
    ymin = atoi(argv[4]);
    ymax = atoi(argv[5]);
    zmin = atoi(argv[6]);
    zmax = atoi(argv[7]);
    x_spacing = atof(argv[8]);
    y_spacing = atof(argv[9]);
    z_spacing = atof(argv[10]);
    x_origin = atof(argv[11]);
    y_origin = atof(argv[12]);
    z_origin = atof(argv[13]);
    scalar_type = atoi(argv[14]);
    is_little_endian = atoi(argv[15]);
  }

  vtkImageReader2 *reader = vtkImageReader2::New();

  reader->SetFileDimensionality(3);
  reader->SetFileName(filename.c_str());
  reader->SetDataExtent(xmin, xmax, ymin, ymax, zmin, zmax);
  reader->SetDataSpacing(x_spacing, y_spacing, z_spacing);
  reader->SetDataOrigin(x_origin, y_origin, z_origin);
  reader->SetDataScalarType(scalar_type);
  is_little_endian ? reader->SetDataByteOrderToLittleEndian() :
    reader->SetDataByteOrderToBigEndian();

  reader->UpdateWholeExtent();

  switch (scalar_type) {
    vtkTemplateMacro( Run(reader->GetOutput(), static_cast<VTK_TT*>(0)) );
  }
  reader->Delete();



  return 0;
}
