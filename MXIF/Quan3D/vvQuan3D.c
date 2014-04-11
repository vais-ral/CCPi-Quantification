/*=========================================================================

  Copyright (c) Yue Sheng, University of Manchester and David Worth, STFC
  All rights reserved.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/** \file vvQuan3D.c
\brief This is a VolView plugin to calculate several characteristics from a 
labelled image.

The following characteristics are calculated:
\li Volume by voxel counts
\li Equivalent sphere diameter by voxel counts
\li Bounding box diagonal
\li Principal Component Analysis
\li Ellipsoid fitting by PCA
\li Equivalent circle diameter by PCA
\li Isosurface by marching cube
\li Surface area
\li Surface volume
\li Equivalent sphere diameter from surface volume
\li Sphercity
\li Normalised surface area to volume ratio (Radius*Sa/Vol)
*/

#include "vvHelper.h"
#include "vtkVVPluginAPI.h"
#include <dlfcn.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/**
Struct used when calculating the number of voxels with each voxel value
*/
typedef struct {
  /** The voxel value */
  int value;
  /** The number of voxels with this value */
  int count;
} vvQuan3DVoxelCount;

/**
Comparison function for sorting the image data.
\param a Integer value
\param b Integer value
\return 1 if \c a>b; -1 if \c a<b; 0 if \c a==b.
*/
int cmpfun(void *a, void *b)
{
  return ( *(int*)a - *(int*)b );
}

/**
\brief Calculate sphere diameter from the sphere volume.

Use the formula \f$ 2*(3v/(4\pi))^{\frac{1}{3}} \f$ where \f$ v \f$ is the 
sphere volume.

\param volume The sphere volume
\return Diameter of the sphere
*/
double vol2dia(double volume)
{
  return 2.0 * pow((3.0*volume)/(4.0*3.14159265),(1.0/3.0));
}

/**
\brief Find the (i,j,k) coordinates for each voxel id and also calculate
maximum, minimum and mean coordinates.

The positions pointer must be allocated by the calling function as follows
\code
positions = malloc(3*sizeof(int)*voxelCount);
\endcode

\param voxelSubArray Array of voxel ids
\param voxelCount Number of voxels in voxelSubArray <b>NB This is not the 
length of the array</b>
\param dim The dimsnesions of the image
\param positions Pointer to matrix which will hold the (i,j,k) coordinates
\param min  The minimum coordinate in each dimension
\param max  The maximum coordinate in each dimension
\param mean The mean of the coordinates in each dimension
*/
void getCoords(int *voxelSubArray, int voxelCount, int dim[], 
               gsl_matrix *positions, int min[], int max[], int mean[])
{
  int i;
  int voxelId;
  int icoord, jcoord, kcoord, rem;
  min[0] = min[1] = min[2] = INT_MAX;
  max[0] = max[1] = max[2] = INT_MIN;
  mean[0] = mean[1] = mean[2] = 0;
  
  /* Loop over the voxels and calulate i,j,k coords for each one */
  for (i = 0; i < voxelCount; i++) {
    voxelId = voxelSubArray[i];
    kcoord = floor(voxelId / (dim[0]*dim[1]));
    rem = voxelId % (dim[0]*dim[1]);
    jcoord = floor(rem / dim[0]);
    icoord = rem % dim[0];
    //printf("voxel %d (%d,%d,%d)\n",voxelId,icoord,jcoord,kcoord);
    gsl_matrix_set(positions,i,0,icoord);
    gsl_matrix_set(positions,i,1,jcoord);
    gsl_matrix_set(positions,i,2,kcoord);
    if (icoord > max[0]) max[0] = icoord;
    if (jcoord > max[1]) max[1] = jcoord;
    if (kcoord > max[2]) max[2] = kcoord;
    if (icoord < min[0]) min[0] = icoord;
    if (jcoord < min[1]) min[1] = jcoord;
    if (kcoord < min[2]) min[2] = kcoord;
    mean[0] += icoord;
    mean[1] += jcoord;
    mean[2] += kcoord;
  }
  mean[0] /= voxelCount;
  mean[1] /= voxelCount;
  mean[2] /= voxelCount;
}

/** 
\brief Calculate the data from the image

\author David Worth, STFC
\date May 2012
\todo This only works for unsigned short image type. Make it do other integer
datatypes.

\param inf Pointer to object that holds data for this plugin. The voxel size 
can be obtained from it.
\param pds Pointer to object that carries information on data set to be 
processed. Includes actual buffer of voxel data, number of voxels along each
dimension, voxel spacing and voxel type.
*/
static int ProcessData(void *inf, vtkVVProcessDataStruct *pds)
{
  /* Information about the plugin and the image */
  vtkVVPluginInfo *info = (vtkVVPluginInfo *)inf;
  /* Number of voxels in each direction of image */
  int *dim = info->InputVolumeDimensions;
  /* Size of datatype in the image */
  unsigned int rsize = info->InputVolumeScalarSize;
  /* Voxel size - can be set by user in VolView */
  double voxelSize = atof(info->GetGUIProperty(info, 0, VVP_GUI_VALUE));
  /* Minimum and maximum raw values from image. */
  double imgMin = info->InputVolumeScalarRange[0], 
         imgMax = info->InputVolumeScalarRange[1];
  /* Array of voxel counts for each value. Array index is the value. */
  int maxVoxelVals = (int)ceil(imgMax)-(int)floor(imgMin)+1;
  int *voxelCounts = malloc(sizeof(int)*maxVoxelVals);
  /* Array of arrays - each sub-array holding the ids of voxels with that value */
  int* *voxelArray = malloc(sizeof(int*)*maxVoxelVals);
  /* Pointer to sub array used when creating voxelArray */
  int *voxelSubArray = NULL;
  /* Index array used to write into voxelSubArray */
  int *indexVoxelSubArray = malloc(sizeof(int)*maxVoxelVals);
  /* Number of the non-background voxels */
  int sumVoxels = 0;
  /* Number of different voxel values */
  int totalParticles = 0;
  /* Array of (i,j,k) coordinates for voxels with a particular value */
  gsl_matrix *positions = NULL;
  /* Stats for voxel data for each value */
  int max[3], min[3], mean[3], boundingBoxSize[3];
  double boundingBoxDiam;
  /* Multiplication factors based on physical voxel size */
  double areaFactor = pow(voxelSize,2), volumeFactor = pow(voxelSize,3);
  /* Output data */
  int voxelLabel, voxelCount;
  double voxelCountVolume, voxelCountDia;
  /* Data used when covariance matrix is calculated */
  double* *matrix = NULL;
  double *row = NULL;
  double* *covMatrix = NULL;
  /* Data for eigenanalysis of covariance matrix */
  gsl_matrix *gslCovMatrix = NULL;
  gsl_vector *eval = NULL;
  gsl_matrix *evecs = NULL;
  gsl_eigen_symmv_workspace *workspace = NULL;
  gsl_vector vector;
  int result;
  
  void *gslcblas_handle, *gsl_handle;
  
  /* Buffer for output to user */
  char buffer[1000];
  int i, j = 0;

  /* The raw image data */
  unsigned short *ptr = (unsigned short *)pds->inData;  
  
  /* Check the data type */
  switch (info->InputVolumeScalarType) {
    case VTK_CHAR:                                               
    case VTK_UNSIGNED_CHAR:                                      
    case VTK_FLOAT:                                              
    case VTK_DOUBLE: {                                         
      sprintf(buffer, "Image datatype must be integer or boolean. This \
      algorithm cannot continue");
      info->SetProperty(info, VVP_ERROR, buffer);
      return -1;
      break;  
    }
  }
    
  printf("Image has type %d size %d\n", info->InputVolumeScalarType,rsize);
  printf("The dimensions of the image are: %d, %d, %d \n",
          dim[0], dim[1], dim[2]);
  sprintf(buffer, "The dimensions of the image are: %d, %d, %d \n",
          dim[0], dim[1], dim[2]);
  printf("There are %d components\n",info->InputVolumeNumberOfComponents);
  printf("The minimum value in the image is %d\n", (int)imgMin);
  printf("and the maximum value is %d\n", (int)imgMax);
           
  /* Sort the image values */
  /*qsort(ptr, (size_t)(dim[0]*dim[1]*dim[2]), (size_t)info->InputVolumeScalarSize,
        (__compar_fn_t)cmpfun);*/
 
  /* 
   * Create array of voxel counts for each voxel value except the min value 
   * which we take to be the background
   */
  for (i = 0; i < maxVoxelVals; i++) voxelCounts[i]=0;
  for (i = 0; i < dim[0]*dim[1]*dim[2]; i++) {
    if (ptr[i] != imgMin) { /* Not the minimum value */
      if (ptr[i] > imgMax) printf("%d %d\n",i,ptr[i]); /* Something wrong */
      /* Not seen this value before so it's a new particle */
      if (voxelCounts[ptr[i]] == 0) totalParticles++;
      voxelCounts[ptr[i]]++;
      sumVoxels++;
    }
  }
  
  /* Create the array of arrays with voxel ids for each voxel value */
  for (i = 0; i < maxVoxelVals; i++) {
    if (voxelCounts[i] > 0) {
      /*printf("Non zero voxel count for i = %d - count is %d\n",i,voxelCounts[i]);*/
      voxelArray[i] = malloc(voxelCounts[i]*sizeof(int));
    }
    else voxelArray[i] = NULL;
    indexVoxelSubArray[i] = 0;
  }
  /* Now populate that array with values */
  for (i = 0; i < dim[0]*dim[1]*dim[2]; i++) {
    if (ptr[i] != imgMin) {
      //printf("add entry to voxel array %d - entry id %d\n",ptr[i],i);
      voxelSubArray = voxelArray[ptr[i]];
      voxelSubArray[ indexVoxelSubArray[ptr[i]] ] = i;
      indexVoxelSubArray[ptr[i]]++;
    }
  }
  
  printf("Porosity is %f\n",(double)sumVoxels/(double)(dim[0]*dim[1]*dim[2]));
  sprintf(buffer+strlen(buffer), "Porosity is %f\n",
          (double)sumVoxels/(double)(dim[0]*dim[1]*dim[2]));
  printf("Number of particles is %d\n",totalParticles);

  /* Load the GSL BLAS and GSL itself. RTLD_GLOBAL means we can use the
   * function names from GSL rather than having to get function pointers with
   * dlsym. */
  gslcblas_handle = dlopen("/usr/lib/libgslcblas.so", RTLD_LAZY|RTLD_GLOBAL);
  if (!gslcblas_handle) {
    printf("Error opening libgslcblas.so - %s\n",dlerror());
  }
  gsl_handle = dlopen("/usr/lib/libgsl.so", RTLD_LAZY|RTLD_GLOBAL);
  if (!gsl_handle) {
    printf("Error opening libgsl.so - %s\n",dlerror());
  }
  
  /* Allocate the GSL matrices and vectors for eigenanalysis */ 
  gslCovMatrix = gsl_matrix_alloc(3,3);
  eval = gsl_vector_alloc(3);
  evecs = gsl_matrix_alloc(3,3);
  workspace = gsl_eigen_symmv_alloc(3);

  /* Loop over the particles */
  for (i = 0; i < 1/*totalParticles*/; i++) {
    /* Find the next particle with a non-zero voxel count */
    while ((j < maxVoxelVals) && (voxelCounts[j] == 0)) j++;
    /* Do some simple stuff */
    voxelLabel = j;
    voxelCount = voxelCounts[j];
    voxelCountVolume = voxelCount * volumeFactor;
    voxelCountDia = voxelSize*vol2dia(voxelCount);
    printf("%d, %d, %f, %f",voxelLabel,voxelCount,voxelCountVolume,
      voxelCountDia);
        
    /* Find the (i,j,k) coordinates in the image of the voxels */
    positions = gsl_matrix_alloc(voxelCount,3);
    voxelSubArray = voxelArray[voxelLabel];
    getCoords(voxelSubArray, voxelCount, dim, positions, min, max, mean);
    /*printf(" (%d, %d, %d)",min[0],min[1],min[2]);
    printf(" (%d, %d, %d)",max[0],max[1],max[2]);
    printf(" (%d, %d, %d)",mean[0],mean[1],mean[2]);*/
    boundingBoxSize[0] = max[0]-min[0]+1;
    boundingBoxSize[1] = max[1]-min[1]+1;
    boundingBoxSize[2] = max[2]-min[2]+1;
    boundingBoxDiam = sqrt(boundingBoxSize[0]*boundingBoxSize[0] +
                           boundingBoxSize[1]*boundingBoxSize[1] +
                           boundingBoxSize[2]*boundingBoxSize[2])*voxelSize;
    printf(" (%d, %d, %d)",boundingBoxSize[0],boundingBoxSize[1],
                           boundingBoxSize[2]);
    printf(" %f", boundingBoxDiam);
      
    /* Calculate positions-mean */
    matrix = malloc(voxelCount*sizeof(double*));
    for (j = 0; j < voxelCount; j++) {
      row = malloc(3*sizeof(double));
      row[0] = gsl_matrix_get(positions,j,0)-mean[0];
      row[1] = gsl_matrix_get(positions,j,1)-mean[1];
      row[2] = gsl_matrix_get(positions,j,2)-mean[2];
      matrix[j] = row;
    }
    covMatrix = CovMat(matrix, voxelCount, 3);
    
    for (j = 0; j < voxelCount; j++) {
      free(matrix[j]);
    }
    free(matrix);

    /*printf("\nCovariance matrix\n");
    for (j = 0; j < 3; j++) {
      printf("%f %f %f\n", covMatrix[j][0], covMatrix[j][1], covMatrix[j][2]);
    }*/
    
    /* To do the eigenanalysis use GNU Scientific library so create matrix
     * in that format */
    for (j = 0; j < 3; j++) {
      gsl_matrix_set(gslCovMatrix,0,j,covMatrix[0][j]);
      gsl_matrix_set(gslCovMatrix,1,j,covMatrix[1][j]);
      gsl_matrix_set(gslCovMatrix,2,j,covMatrix[2][j]);
    }
    
    /* Calculate eigenvalues and vectors */
    result = gsl_eigen_symmv(gslCovMatrix, eval, evecs, workspace);
    
    printf("\nEigenvalues\n");
    gsl_vector_fprintf(stdout, eval, "%g");
    printf("\nEigenvectors\n");
    vector = gsl_matrix_const_column(evecs,0).vector;
    gsl_vector_fprintf(stdout, &vector, "%g");
    printf("\n");
    vector = gsl_matrix_const_column(evecs,1).vector;
    gsl_vector_fprintf(stdout, &vector, "%g");
    printf("\n");
    vector = gsl_matrix_const_column(evecs,2).vector;
    gsl_vector_fprintf(stdout, &vector, "%g");

    
    for (j = 0; j < 3; j++) {
      free(covMatrix[j]);
    }
    free(covMatrix);
    
    gsl_matrix_free(positions);
    
    printf("\n");
    /* Move on to next particle */
    j++;
  }
  
  dlclose(gsl_handle);
  dlclose(gslcblas_handle);
  
  info->SetProperty(info, VVP_REPORT_TEXT, buffer);

  free(voxelCounts);
  for (i = 0; i < maxVoxelVals; i++) {
    if (voxelArray[i] != NULL) free(voxelArray[i]);
  }
  
  gsl_matrix_free(gslCovMatrix);
  gsl_vector_free(eval);
  gsl_matrix_free(evecs);
  gsl_eigen_symmv_free(workspace);

  free(voxelArray);
  free(indexVoxelSubArray);
  return 0;
}

/** 
\brief Update the VolView GUI to display user parameters.

Sets one GUI parameter - the physical size of the voxels in the image.

\param inf Pointer to object that should be modified to set up GUI elements for
this plugin. It also contains details of the input and output images.
*/
static int UpdateGUI(void *inf)
{
  int i;
  vtkVVPluginInfo *info = (vtkVVPluginInfo *)inf;

  info->SetGUIProperty(info, 0, VVP_GUI_LABEL, "Voxel Size");
  info->SetGUIProperty(info, 0, VVP_GUI_TYPE, VVP_GUI_SCALE);
  info->SetGUIProperty(info, 0, VVP_GUI_DEFAULT , "1");
  info->SetGUIProperty(info, 0, VVP_GUI_HELP,
               "The physical size of a voxel in the image");
    
  /* what range should we show for possible output values */
  vvPluginSetGUIScaleRange(0);

  info->OutputVolumeScalarType = info->InputVolumeScalarType;
  info->OutputVolumeNumberOfComponents = info->InputVolumeNumberOfComponents;
  for (i = 0; i < 3; i++)
    {
    info->OutputVolumeDimensions[i] = info->InputVolumeDimensions[i];
    info->OutputVolumeSpacing[i] = info->InputVolumeSpacing[i];
    info->OutputVolumeOrigin[i] = info->InputVolumeOrigin[i];
    }

  return 1;
}

/** 
\brief Initialise this plugin.

Sets the name and group for this plugin in the plugin list shown to the VolView
user and gives some documentation. Also defines properties so VolView can judge
the memory requirements and potential for undoing this plugin.

\param info Pointer to object that should be modified to give details about this
plugin.
*/
void VV_PLUGIN_EXPORT vvQuan3DCInit(vtkVVPluginInfo *info)
{
  /* always check the version */
  vvPluginVersionCheck();
  
  /* setup information that never changes */
  info->ProcessData = ProcessData;
  info->UpdateGUI = UpdateGUI;
  info->SetProperty(info, VVP_NAME, "Quantify 3D in C");
  info->SetProperty(info, VVP_GROUP, "Quantification");
  info->SetProperty(info, VVP_TERSE_DOCUMENTATION,
     "Quantify several characteristics from a labelled image");
  info->SetProperty(info, VVP_FULL_DOCUMENTATION,
    "Quantify the following characteristics from a labelled image: Volume by \
voxel counts, Equivalent sphere diameter by voxel counts, Bounding box \
diagonal, Principal Component Analysis, Ellipsoid fitting by PCA, Equivalent \
circle diameter by PCA, Isosurface by marching cube, Surface area, Surface \
volume, Equivalent sphere diameter from surface volume, Sphercity, Normalised \
surface area to volume ratio");
  
  info->SetProperty(info, VVP_SUPPORTS_IN_PLACE_PROCESSING, "1");
  info->SetProperty(info, VVP_SUPPORTS_PROCESSING_PIECES,   "1");
  info->SetProperty(info, VVP_REQUIRED_Z_OVERLAP,           "0");
  info->SetProperty(info, VVP_NUMBER_OF_GUI_ITEMS,          "1");
  info->SetProperty(info, VVP_REQUIRES_SERIES_INPUT,        "0");
  info->SetProperty(info, VVP_SUPPORTS_PROCESSING_SERIES_BY_VOLUMES, "0");
  info->SetProperty(info, VVP_PRODUCES_OUTPUT_SERIES, "0");
  info->SetProperty(info, VVP_PRODUCES_PLOTTING_OUTPUT, "0");
}




