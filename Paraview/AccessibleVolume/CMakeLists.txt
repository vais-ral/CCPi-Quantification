cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)

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

ADD_PARAVIEW_PLUGIN(AccessibleVolume(CCPi) "1.0"
  SERVER_MANAGER_XML CCPiAccessibleVolume.xml 
  SERVER_MANAGER_SOURCES CCPiAccessibleVolumeParaviewImpl.cxx)

