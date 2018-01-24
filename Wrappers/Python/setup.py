#!/usr/bin/env python

import setuptools
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

import os
import sys
import numpy
import platform	
import vtk

cil_version=os.environ['CIL_VERSION']
if  cil_version == '':
    print("Please set the environmental variable CIL_VERSION")
    sys.exit(1)

itk_version="4.11"
vtk_version="%d.%d"%(vtk.vtkVersion.GetVTKMajorVersion(), vtk.vtkVersion.GetVTKMinorVersion())
library_include_path = []
library_lib_path = []
try:
    library_include_path = [ os.environ['LIBRARY_INC'] ]
    library_lib_path = [ os.environ['LIBRARY_LIB'] ]
except:
    if platform.system() == 'Windows':
        pass
    else:
        try:
           library_include_path = [ os.environ['PREFIX']+'/include' ]
           library_lib_path = [ os.environ['PREFiX']+'/lib' ]
        except:
           pass
    pass
extra_include_dirs = [numpy.get_include()]
extra_library_dirs = []
extra_compile_args = []
extra_link_args = []
extra_libraries = ['ITKCommon-'+itk_version,'itkvcl-'+itk_version,'itkvnl-'+itk_version, 'itkvnl_algo-'+itk_version, 'itkv3p_netlib-'+itk_version, 'itksys-'+itk_version, 'ITKStatistics-'+itk_version, 'vtkCommonCore-'+vtk_version,  'vtkImagingCore-'+vtk_version, 'vtkImagingHybrid-'+vtk_version, 'vtkFiltersCore-'+vtk_version, 'vtkCommonComputationalGeometry-'+vtk_version,'vtkCommonDataModel-'+vtk_version,'vtkIOImage-'+vtk_version,'vtkCommonExecutionModel-'+vtk_version, 'vtkCommonMisc-'+vtk_version]
lip = ''
try:
    lip = library_include_path[0]
except:
    pass
if platform.system() == 'Windows':
    extra_compile_args += ['/TP','/DWIN32','/EHsc','/DCCPiCore_EXPORTS','/Ob2','/O2','/DBOOST_ALL_NO_LIB']
    extra_include_dirs += ["..\\..\\Core\\", ".", lip+"\\ITK-"+itk_version, lip+"\\vtk-"+vtk_version]
    if sys.version_info.major == 3 :   
        extra_libraries += ['boost_python3-vc140-mt-1_64', 'boost_numpy3-vc140-mt-1_64']
    else:
        extra_libraries += ['boost_python-vc90-mt-1_64', 'boost_numpy-vc90-mt-1_64']
    extra_include_dirs += library_include_path
    extra_library_dirs += library_lib_path    
else:
    extra_include_dirs += ["../../Core/", ".", lip+"/ITK-"+itk_version, lip+"/vtk-"+vtk_version]
    extra_include_dirs += [library_include_path]
    if sys.version_info.major == 3 :
        extra_libraries += ['boost_python3', 'boost_numpy3','gomp']
    else:
        extra_libraries += ['boost_python', 'boost_numpy','gomp']
    extra_compile_args = ['-fopenmp','-O2', '-funsigned-char', '-Wall']
        
        
setup(
    name='ccpi',
	description='This is a CCPi Core Imaging Library package for Quantification codes',
	version=cil_version,
    cmdclass = {'build_ext': build_ext},
    ext_modules = cythonize([Extension("ccpi.quantification.LabelQuantification",
                             sources=["src/LabelQuantification.pyx", 
                                      "src/LabelQuantification_c.cpp",
                                        "../../Core/QuanWorker.cpp",
                                        "../../Core/Quan3D.cpp",
                                        "../../Core/CCPiConsoleUserInterface.cpp",
                                        "../../Core/CCPiLabelQuantificationResult.cpp",],
                             include_dirs=extra_include_dirs, library_dirs=extra_library_dirs, extra_compile_args=extra_compile_args, libraries=extra_libraries, extra_link_args=extra_link_args),
                             Extension("ccpi.quantification.AccessibleVolume",
                             sources=["src/AccessibleVolume.pyx", 
                                      "src/AccessibleVolume_c.cpp",
                                        "../../Core/CCPiAccessibleVolumeInputImages.cpp",
                                        "../../Core/CCPiAccessibleVolumeITKImpl.cpp",
                                        "../../Core/CCPiConsoleUserInterface.cpp",
                                        "../../Core/CCPiLabelQuantificationResult.cpp",],
                             include_dirs=extra_include_dirs, library_dirs=extra_library_dirs, extra_compile_args=extra_compile_args, libraries=extra_libraries, extra_link_args=extra_link_args)],
    language='c++'),
    zip_safe = False,    
    packages = {'ccpi'}
)
