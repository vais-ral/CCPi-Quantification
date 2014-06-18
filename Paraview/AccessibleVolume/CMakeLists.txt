cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
PROJECT (CCPiParaviewPlugin)

if (NOT ParaView_BINARY_DIR)
  find_package(ParaView REQUIRED)
  include(${PARAVIEW_USE_FILE})
endif()

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})
												
include(ParaViewPlugins)

# create a paraview plugin containing server manager xml and the server
# manager classes to build
# this plugin can be loaded on the server side

ADD_PARAVIEW_PLUGIN(CCPiAccessibleVolume(CCPi) "1.0"
  SERVER_MANAGER_XML CCPiAccessibleVolume.xml 
  SERVER_MANAGER_SOURCES CCPiAccessibleVolumeParaviewImpl.cxx
  SERVER_SOURCES CCPiParaviewUserInterface.cpp ../../Core/CCPiAccessibleVolumeInputImages.cpp ../../Core/CCPiAccessibleVolumeITKImpl.cpp)
ADD_PARAVIEW_PLUGIN(CCPiSimpleHistogramThresholding(CCPi) "1.0"
  SERVER_MANAGER_XML CCPiAccessibleVolume.xml 
  SERVER_MANAGER_SOURCES CCPiSimpleHistogramThresholdingParaviewImpl.cxx
  SERVER_SOURCES CCPiParaviewUserInterface.cpp )
ADD_PARAVIEW_PLUGIN(CCPiLabelQuantification(CCPi) "1.0"
  SERVER_MANAGER_XML CCPiAccessibleVolume.xml
  SERVER_MANAGER_SOURCES CCPiLabelQuantificationParaviewImpl.cxx 
  SERVER_SOURCES CCPiParaviewUserInterface.cpp ../../Core/Quan3D.cpp ../../Core/QuanWorker.cpp )

include_directories(../../Core)  
TARGET_LINK_LIBRARIES(CCPiAccessibleVolume ITKCommon 
							itkvnl itkvnl_algo itkvnl_inst itkv3p_netlib itksys)							
TARGET_LINK_LIBRARIES(CCPiLabelQuantification ITKCommon 
							itkvnl itkvnl_algo itkvnl_inst itkv3p_netlib itksys)
TARGET_LINK_LIBRARIES(CCPiSimpleHistogramThresholding ITKCommon 
							itkvnl itkvnl_algo itkvnl_inst itkv3p_netlib itksys ITKStatistics)  