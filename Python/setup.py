#!/usr/bin/env python

import setuptools
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import os
import numpy
import platform	

extra_include_dirs = [numpy.get_include()]
extra_library_dirs = []
extra_compile_args = ['-fopenmp','-O2', '-funsigned-char', '-Wall', '-Werror']
extra_libraries = ['ITKCommon-4.11','itkvcl-4.11','itkvnl-4.11', 'itkvnl_algo-4.11', 'itkv3p_netlib-4.11', 'itksys-4.11', 'ITKStatistics-4.11', 'vtkCommonCore-7.1',  'vtkImagingCore-7.1', 'vtkImagingHybrid-7.1', 'vtkFiltersCore-7.1', 'vtkCommonComputationalGeometry-7.1','vtkCommonDataModel-7.1','vtkIOImage-7.1','vtkCommonExecutionModel-7.1']
if platform.system() == 'Windows':
   extra_compile_args[0:] = ['/DWIN32','/EHsc','/DCCPiCore_EXPORTS','/Ob2','/O2']
   extra_include_dirs += ['C:\\Apps\\libs-vs2015\\boost\\1.64\\include\\boost-1_65',"..\\Core\\", ".", "C:\\Apps\\Anaconda3\\envs\\python3.5\\Library\\include\\ITK-4.11","C:\\Apps\\Anaconda3\\envs\\Python3.5\\Library\\include\\vtk-7.1"]
   extra_library_dirs += ['C:\\Apps\\libs-vs2015\\boost\\1.64\\lib',"C:\\Apps\\Anaconda3\\envs\\python3.5\\Library\\lib"]   
setup(
    name='ccpi',
	description='This is a CCPi package for Quantification codes',
	version='0.1',
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
	packages = {'ccpi'}
)