#############################################################################
# CMakeLists.txt file for Building CCPiQuanitifcation Core and Nexus Widget
############################################################################
cmake_minimum_required(VERSION 2.8)
project(CCPiQuantification)

set(CCPIQUANTIFICATION_MAJOR_VERSION 0)
set(CCPIQUANTIFICATION_MINOR_VERSION 9)
set(CCPIQUANTIFCIATION_PATCH_VERSION 1)
set(CCPIQUANTIFICATION_VERSION 
	${CCPIQUANTIFICATION_MAJOR_VERSION}.${CCPIQUANTIFICATION_MINOR_VERSION}.${CCPIQUANTIFCIATION_PATCH_VERSION})

set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")

if(WIN32 AND NOT CYGWIN)
	SET(DEF_INSTALL_CMAKE_DIR cmake)
else()
	SET(DEF_INSTALL_CMAKE_DIR lib/cmake/CCPiQuantification)
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
	"Installation directory for CMake files")

#Make relative paths absolute
foreach(p LIB BIN INCLUDE CMAKE)
	set(var INSTALL_{p}_DIR)
	if(NOT IS_ABSOLUTE "${${var}}")
		set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
	endif()
endforeach()


#setup include directories
include_directories(
	"${PROJECT_SOURCE_DIR}"
	"${PROJECT_BINARY_DIR}")

# Add sub-directories
add_subdirectory(Core)
#Allow the developer to select Nexus Widget to be built with the core
OPTION(BUILD_NEXUS_WIDGET "Build Nexus Widget" OFF)
IF(BUILD_NEXUS_WIDGET)
	add_subdirectory(Core/NexusWidget)
ENDIF()
OPTION(BUILD_AVIZO_PLUGIN "Build Avizo plugin" OFF)
IF(BUILD_AVIZO_PLUGIN)
	add_subdirectory(Wrappers/Avizo)
ENDIF()
OPTION(BUILD_PARAVIEW_PLUGIN "Build Paraview plugin" OFF)
IF(BUILD_PARAVIEW_PLUGIN)
	add_subdirectory(Wrappers/Paraview)
ENDIF()
OPTION(BUILD_IMAGEJ_PLUGIN "Build ImageJ plugin" OFF)
IF(BUILD_IMAGEJ_PLUGIN)
	add_subdirectory(Wrappers/ImageJ)
ENDIF()

IF(BUILD_NEXUS_WIDGET)
SET(CCPI_QUANTIFICATION_BUILD_LIBRARIES "CCPiCore CCPiNexusWidget")
export(TARGETS CCPiCore CCPiNexusWidget
	FILE "${PROJECT_BINARY_DIR}/CCPiQuantificationTargets.cmake")
ELSE()
SET(CCPI_QUANTIFICATION_BUILD_LIBRARIES "CCPiCore")
export(TARGETS CCPiCore
	FILE "${PROJECT_BINARY_DIR}/CCPiQuantificationTargets.cmake")
ENDIF()

export(PACKAGE CCPiQuantification)
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
	"${INSTALL_INCLUDE_DIR}")
# for build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
configure_file(${CMAKE_SOURCE_DIR}/CMake/CCPiQuantificationConfig.cmake.in  
		"${PROJECT_BINARY_DIR}/CCPiQuantificationConfig.cmake" @ONLY)
# for install tree
set(CONF_INCLUDE_DIRS "\${CCPIQUANTIFICATION_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(${CMAKE_SOURCE_DIR}/CMake/CCPiQuantificationConfig.cmake.in
		"${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CCPiQuantificationConfig.cmake" @ONLY)
# for both
configure_file(${CMAKE_SOURCE_DIR}/CMake/CCPiQuantificationConfigVersion.cmake.in
		"${PROJECT_BINARY_DIR}/CCPiQuantificationConfigVersion.cmake" @ONLY)

#Install
install(FILES
	"${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CCPiQuantificationConfig.cmake"
	"${PROJECT_BINARY_DIR}/CCPiQuantificationConfigVersion.cmake"
	DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)

# Install the export set
install(EXPORT CCPiQuantificationTargets DESTINATION
	"${INSTALL_CMAKE_DIR}" COMPONENT dev)		

