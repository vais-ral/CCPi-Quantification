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
  SERVER_SOURCES CCPiParaviewUserInterface.cpp)
ADD_PARAVIEW_PLUGIN(CCPiSimpleHistogramThresholding(CCPi) "1.0"
  SERVER_MANAGER_XML CCPiAccessibleVolume.xml 
  SERVER_MANAGER_SOURCES CCPiSimpleHistogramThresholdingParaviewImpl.cxx
  SERVER_SOURCES CCPiParaviewUserInterface.cpp )
ADD_PARAVIEW_PLUGIN(CCPiLabelQuantification(CCPi) "1.0"
  SERVER_MANAGER_XML CCPiAccessibleVolume.xml
  SERVER_MANAGER_SOURCES CCPiLabelQuantificationParaviewImpl.cxx 
  SERVER_SOURCES CCPiParaviewUserInterface.cpp)

ADD_SUBDIRECTORY(../../Core "${CMAKE_CURRENT_BINARY_DIR}/CCPiCore")
 
include_directories(../../Core)  
TARGET_LINK_LIBRARIES(CCPiAccessibleVolume CCPiCore ITKCommon 
							itkvnl itkvnl_algo itkvnl_inst itkv3p_netlib itksys)							
TARGET_LINK_LIBRARIES(CCPiLabelQuantification CCPiCore ITKCommon 
							itkvnl itkvnl_algo itkvnl_inst itkv3p_netlib itksys)
TARGET_LINK_LIBRARIES(CCPiSimpleHistogramThresholding CCPiCore ITKCommon 
							itkvnl itkvnl_algo itkvnl_inst itkv3p_netlib itksys ITKStatistics)  

SET(CPACK_GENERATOR "RPM")
SET(CPACK_RPM_PACKAGE_SUMMARY "CCPi Quantification Code Paraview Plugin")
SET(CPACK_RPM_PACKAGE_NAME "CCPiParaviewPlugin")
SET(CPACK_RPM_PACKAGE_VERSION "4.1")
SET(CPACK_RPM_PACKAGE_RELEASE "0.1")
SET(CPACK_RPM_PACKAGE_VENDOR  "CCPi")
SET(CPACK_RPM_PACKAGE_DESCRIPTION "This package provides CCPi Quantification
code algorithms in Paraview.")
 
INCLUDE(CPack)
