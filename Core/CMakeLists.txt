cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
PROJECT (CCPiCore)

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

IF(ParaView_DIR)
FIND_PACKAGE(VTK COMPONENTS vtkCommonCore vtkImagingCore vtkImagingHybrid vtkFiltersCore vtkCommonComputationalGeometry NO_MODULE)
INCLUDE(${VTK_USE_FILE})
ELSE()
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})
ENDIF()


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
			
IF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
  SET_TARGET_PROPERTIES(CCPiCore PROPERTIES COMPILE_FLAGS "-fPIC")
ENDIF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
			
SET (CCPI_VTK_LIBS)
IF (NOT(VTK_MAJOR_VERSION LESS 6))			
SET (CCPI_VTK_LIBS vtkCommonCore vtkImagingCore vtkImagingHybrid vtkFiltersCore vtkCommonComputationalGeometry)
ELSE()
SET (CCPI_VTK_LIBS vtkHybrid)
ENDIF()

 if( ${ITK_VERSION_MAJOR} LESS 4 )
   SET(CCPI_ITK_LIBS ITKCommon itkvnl itkvnl_algo itkvnl_inst itkv3p_netlib itksys ITKStatistics)
 else( "${ITK_VERSION_MAJOR}" )
   SET(CCPI_ITK_LIBS ${ITK_LIBRARIES} )
 endif()


TARGET_LINK_LIBRARIES(CCPiCore  ${CCPI_ITK_LIBS} ${CCPI_VTK_LIBS})							
							

OPTION(BUILD_TESTING "Build the testing tree." ON)
IF(BUILD_TESTING)
   ENABLE_TESTING()   
ENDIF(BUILD_TESTING)

SUBDIRS(test)					

#Allow the developer to select Nexus Widget to be built with the core
OPTION(BUILD_NEXUS_WIDGET "Build Nexus Widget" OFF)
IF(BUILD_NEXUS_WIDGET)
	SUBDIRS(NexusWidget)
ENDIF()
INSTALL(TARGETS CCPiCore
        DESTINATION "lib"
       )
INSTALL(FILES CCPiAccessibleVolumeInputImages.h
              CCPiConsoleUserInterface.h
              CCPiImageData.h
              CCPiLabelQuantificationResult.h
              CCPiParticleTracker.h
              CCPiTrack.h
              CCPiAccessibleVolumeITKImpl.h
              CCPiDefines.h
              CCPiLabelQuantificationITKImpl.h
              CCPiParticle.h
              CCPiSimpleHistogramThresholdingITKImpl.h
              CCPiUserApplicationInterface.h

        DESTINATION "include"
        )
