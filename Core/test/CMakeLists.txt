IF(BUILD_TESTING)
   ENABLE_TESTING()   
ENDIF(BUILD_TESTING)

INCLUDE_DIRECTORIES (../)
ADD_EXECUTABLE(AccessibleVolumeTest AccessibleVolumeTest.cpp) 
TARGET_LINK_LIBRARIES(AccessibleVolumeTest CCPiCore ITKCommon ITKIO ${CCPI_VTK_LIBS})
	
ADD_TEST(NAME AccessibleVolumeTest 	COMMAND AccessibleVolumeTest "Data/Data128.tif" "Data/DataMask128.tif")	

SUBDIRS(Data)

IF(WIN32)
set(PATH_STRING "$ENV{PATH};${ITK_LIBRARY_DIRS}/Release;${ITK_LIBRARY_DIRS}") 
STRING(REPLACE "\\;" ";" PATH_STRING "${PATH_STRING}")
STRING(REPLACE ";" "\\;" PATH_STRING "${PATH_STRING}")
set_tests_properties(AccessibleVolumeTest PROPERTIES ENVIRONMENT "PATH=${PATH_STRING}") 
ENDIF()