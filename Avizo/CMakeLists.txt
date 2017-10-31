cmake_minimum_required(VERSION 2.8)
project(CCPiQuantificationAvizo)

# Create the source files
configure_file(${CMAKE_SOURCE_DIR}/CMake/version.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/version.h)
configure_file(${CMAKE_SOURCE_DIR}/CMake/api.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/CCPi/api.h)
	#List of function to wrap in Avizo
	set(Avizo_Function_List CCPiAccessibleVolume
				CCPiAvizoNexusCmdReader
				CCPiKMeansFilterITK
				CCPiLabelQuantification
				CCPiParticleTracking
				CCPiSimpleHistogramThresholding
	)

	set(Avizo_Wrapper_Exit_Function_String "")
        set(Avizo_Wrapper_Extern_Function_String "")
        set(Avizo_Wrapper_Init_Function_String "")
	foreach(fname ${Avizo_Function_List})
		set(Avizo_Wrapper_Init_Function_String "${Avizo_Wrapper_Init_Function_String}    mcInitClass_${fname}();")
		set(Avizo_Wrapper_Exit_Function_String "${Avizo_Wrapper_Exit_Function_String}    mcExitClass_${fname}();")
		set(Avizo_Wrapper_Extern_Function_String "${Avizo_Wrapper_Extern_Function_String}void mcExitClass_${fname}();")
		if(WIN32)
			set(Avizo_Wrapper_Init_Function_String "${Avizo_Wrapper_Init_Function_String}\r\n")
			set(Avizo_Wrapper_Exit_Function_String "${Avizo_Wrapper_Exit_Function_String}\r\n")
			set(Avizo_Wrapper_Extern_Function_String "${Avizo_Wrapper_Extern_Function_String}\r\n")
		else()
			set(Avizo_Wrapper_Init_Function_String "${Avizo_Wrapper_Init_Function_String}\n")
			set(Avizo_Wrapper_Exit_Function_String "${Avizo_Wrapper_Exit_Function_String}\n")
			set(Avizo_Wrapper_Extern_Function_String "${Avizo_Wrapper_Extern_Function_String}\n")
		endif()	
		set(Avizo_Wrapper_Extern_Function_String "${Avizo_Wrapper_Extern_Function_String}void mcInitClass_${fname}();")
		if(WIN32)
			set(Avizo_Wrapper_Extern_Function_String "${Avizo_Wrapper_Extern_Function_String}\r\n")
		else()
			set(Avizo_Wrapper_Extern_Function_String "${Avizo_Wrapper_Extern_Function_String}\n")
		endif()			
	endforeach()
configure_file(${CMAKE_SOURCE_DIR}/CMake/init.cpp.in ${CMAKE_CURRENT_SOURCE_DIR}/src/CCPi/init.cpp)

# Set all the source files in the src directory
FILE(GLOB CCPiQuantificationAvizo_SRC ${PROJECT_SOURCE_DIR}/src/CCPi/*.cpp ${PROJECT_SOURCE_DIR}/src/CCPi/*.h)


#MESSAGE( "Avizo Source Files:" ${CCPiAvizo_SRC})
add_library(CCPiQuantificationAvizo SHARED ${CCPiQuantificationAvizo_SRC})

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)



find_package (Avizo)
IF (Avizo_FOUND)
  include_directories( ${Avizo_INCLUDE_DIRS} )
  target_link_libraries (CCPiQuantificationAvizo  ${Avizo_LIBRARIES})
  add_definitions(${Avizo_DEFINITIONS})
ENDIF (Avizo_FOUND)

#Custom definitions
add_definitions(-DPACKAGEID=amirapackage_CCPiQuantificationAvizo)
set(CMAKE_CXX_FLAGS "-shared -fPIC -fopenmp -Wall -Wno-unknown-pragmas -mmmx -msse -msse2 -mfpmath=sse -fmessage-length=0 -fexceptions -frounding-math")


# Find ITK Library
FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

# Find VTK Library
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})
SET(CCPI_VTK_LIBS vtkIOImage vtkCommonCore vtkImagingCore vtkImagingHybrid vtkFiltersCore vtkCommonComputationalGeometry)

# Find HDF5 Library
set(_SAVED_HDF5_DIR ${HDF5_DIR})
FIND_PACKAGE(HDF5 COMPONENTS C NO_MODULE)
include_directories(${HDF5_INCLUDE_DIR})
include_directories(${CCPiCore_INCLUDE_DIRS})
include_directories(${CCPiNexusWidget_INCLUDE_DIRS})
target_link_libraries (CCPiQuantificationAvizo CCPiCore CCPiNexusWidget ${ITK_LIBRARIES} ${CCPI_VTK_LIBS} ${HDF5_LIBRARIES})

