cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)
PROJECT (cilquant)

FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

IF(ParaView_DIR)
FIND_PACKAGE(VTK COMPONENTS vtkIOCore vtkIOImage vtkIOImport vtkCommonCore vtkImagingCore vtkImagingHybrid vtkFiltersCore vtkCommonComputationalGeometry NO_MODULE)
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

ADD_LIBRARY(cilquant ${LIB_TYPE} QuanWorker.cpp
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
  SET_TARGET_PROPERTIES(cilquant PROPERTIES COMPILE_FLAGS "-fPIC")
ENDIF( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
			
SET (CCPI_VTK_LIBS)
IF (NOT(VTK_MAJOR_VERSION LESS 6))			
SET (CCPI_VTK_LIBS vtkIOImage vtkCommonCore vtkImagingCore vtkImagingHybrid vtkFiltersCore vtkCommonComputationalGeometry)
ELSE()
SET (CCPI_VTK_LIBS vtkHybrid)
ENDIF()

 if( ${ITK_VERSION_MAJOR} LESS 4 )
   SET(CCPI_ITK_LIBS ITKCommon itkvnl itkvnl_algo itkvnl_inst itkv3p_netlib itksys ITKStatistics)
 else( "${ITK_VERSION_MAJOR}" )
   SET(CCPI_ITK_LIBS ${ITK_LIBRARIES} )
 endif()


TARGET_LINK_LIBRARIES(cilquant  ${CCPI_ITK_LIBS} ${CCPI_VTK_LIBS})							
SET_TARGET_PROPERTIES(cilquant  PROPERTIES PUBLIC_HEADER "CCPiAccessibleVolumeInputImages.h;CCPiConsoleUserInterface.h;CCPiImageData.h;CCPiLabelQuantificationResult.h;CCPiParticleTracker.h;CCPiTrack.h;CCPiAccessibleVolumeITKImpl.h;CCPiDefines.h;CCPiLabelQuantificationITKImpl.h;CCPiParticle.h;CCPiSimpleHistogramThresholdingITKImpl.h;CCPiUserApplicationInterface.h;")							

OPTION(BUILD_TESTING "Build the testing tree." ON)
IF(BUILD_TESTING)
   ENABLE_TESTING()   
ENDIF(BUILD_TESTING)

add_subdirectory(test)					


INSTALL(TARGETS cilquant EXPORT CCPiQuantificationTargets
        LIBRARY  DESTINATION "${INSTALL_LIB_DIR}" COMPONENT lib
	PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}" COMPONENT dev
	RUNTIME  DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
	ARCHIVE  DESTINATION "${INSTALL_LIB_DIR}" COMPONENT lib
       )


############Expose include files
set(${PROJECT_NAME}_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}
                                 ${C_INCLUDE_DIRS}
    CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)


