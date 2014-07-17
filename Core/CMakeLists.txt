cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
PROJECT (CCPiCore)

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

# Allow the developer to select if Dynamic or Static libraries are built
OPTION (BUILD_SHARED_LIBS "Build Shared Libraries" OFF)
# Set the LIB_TYPE variable to STATIC
SET (LIB_TYPE STATIC)
IF (BUILD_SHARED_LIBS)
  # User wants to build Dynamic Libraries, so change the LIB_TYPE variable to CMake keyword 'SHARED'
  SET (LIB_TYPE SHARED)
ENDIF (BUILD_SHARED_LIBS)

ADD_LIBRARY(CCPiCore ${LIB_TYPE} QuanWorker.cpp
							Quan3D.cpp
							CCPiAccessibleVolumeInputImages.cpp
							CCPiAccessibleVolumeITKImpl.cpp
							CCPiConsoleUserInterface.cpp
							CCPiLabelQuantificationResult.cpp
							CCPiParticle.cpp
							CCPiParticleTracker.cpp
							CCPiTrack.cpp
			)
						 
TARGET_LINK_LIBRARIES(CCPiCore ITKCommon ITKStatistics
							itkvnl itkvnl_algo itkvnl_inst itkv3p_netlib itksys vtkHybrid)							

OPTION(BUILD_TESTING "Build the testing tree." ON)
IF(BUILD_TESTING)
   ENABLE_TESTING()   
ENDIF(BUILD_TESTING)

SUBDIRS(test)					