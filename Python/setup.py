#!/usr/bin/env python

import setuptools
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import os
import sys
import numpy
import platform	
import vtk

itk_version="4.11"
vtk_version="%d.%d"%(vtk.vtkVersion.GetVTKMajorVersion(), vtk.vtkVersion.GetVTKMinorVersion())
library_include_path = ""
try:
    library_include_path = os.environ['LIBRARY_INC']
except:
    library_include_path = os.environ['PREFIX']+'/include'
    pass
extra_include_dirs = [numpy.get_include(), library_include_path]
extra_library_dirs = []
extra_compile_args = ['-fopenmp','-O2', '-funsigned-char', '-Wall']
extra_libraries = ['ITKCommon-'+itk_version,'itkvcl-'+itk_version,'itkvnl-'+itk_version, 'itkvnl_algo-'+itk_version, 'itkv3p_netlib-'+itk_version, 'itksys-'+itk_version, 'ITKStatistics-'+itk_version, 'vtkCommonCore-'+vtk_version,  'vtkImagingCore-'+vtk_version, 'vtkImagingHybrid-'+vtk_version, 'vtkFiltersCore-'+vtk_version, 'vtkCommonComputationalGeometry-'+vtk_version,'vtkCommonDataModel-'+vtk_version,'vtkIOImage-'+vtk_version,'vtkCommonExecutionModel-'+vtk_version, 'vtkCommonMisc-'+vtk_version]
if platform.system() == 'Windows':
    extra_compile_args[0:] = ['/DWIN32','/EHsc','/DCCPiCore_EXPORTS','/Ob2','/O2','/DBOOST_ALL_NO_LIB']
    extra_include_dirs += ["..\\Core\\", ".", library_include_path+"\\ITK-"+itk_version, library_include_path+"\\vtk-"+vtk_version]
    if sys.version_info.major == 3 :   
        extra_libraries += ['boost_python3-vc140-mt-1_64', 'boost_numpy3-vc140-mt-1_64']
    else:
        extra_libraries += ['boost_python-vc90-mt-1_64', 'boost_numpy-vc90-mt-1_64']
else:
    extra_include_dirs += ["../Core/", ".", library_include_path+"/ITK-"+itk_version, library_include_path+"/vtk-"+vtk_version]
    extra_include_dirs += [library_include_path]
    if sys.version_info.major == 3 :
        extra_libraries += ['boost_python3', 'boost_numpy3','gomp']
    else:
        extra_libraries += ['boost_python', 'boost_numpy','gomp']
        
setup(
    name='ccpi',
	description='This is a CCPi Core Imaging Library package for Quantification codes',
	version='0.9',
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("ccpi.quantification",
                             sources=["src/ccpi_quantification.cpp",  
                                        "../Core/QuanWorker.cpp",
                                        "../Core/Quan3D.cpp",
                                        "../Core/CCPiAccessibleVolumeInputImages.cpp",
                                        "../Core/CCPiAccessibleVolumeITKImpl.cpp",
                                        "../Core/CCPiConsoleUserInterface.cpp",
                                        "../Core/CCPiLabelQuantificationResult.cpp",
                                        "../Core/CCPiParticle.cpp",
                                        "../Core/CCPiParticleTracker.cpp",
                                        "../Core/CCPiTrack.cpp"],
                             include_dirs=extra_include_dirs, library_dirs=extra_library_dirs, libraries=extra_libraries, extra_compile_args=extra_compile_args),],
    zip_safe = False,    
    packages = {'ccpi'}
)
