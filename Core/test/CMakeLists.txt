IF(BUILD_TESTING)
   ENABLE_TESTING()   
ENDIF(BUILD_TESTING)

 if( "${ITK_VERSION_MAJOR}" LESS 4 )
   SET(CCPI_ITK_LIBS ITKCommon itkvnl itkvnl_algo itkvnl_inst itkv3p_netlib itksys ITKIO ITKStatistics)
else( "${ITK_VERSION_MAJOR}" )
   SET(CCPI_ITK_LIBS ${ITK_LIBRARIES} )
 endif()
 
INCLUDE_DIRECTORIES (../)
ADD_EXECUTABLE(AccessibleVolumeTest AccessibleVolumeTest.cpp) 
TARGET_LINK_LIBRARIES(AccessibleVolumeTest cilquant ${CCPI_ITK_LIBS} ${CCPI_VTK_LIBS} )
	
ADD_TEST(NAME AccessibleVolumeTest 	COMMAND AccessibleVolumeTest "Data/Data128.tif" "Data/DataMask128.tif")	

SUBDIRS(Data)

IF(WIN32)
set(PATH_STRING "$ENV{PATH};${ITK_LIBRARY_DIRS}/Release;${ITK_LIBRARY_DIRS}") 
STRING(REPLACE "\\;" ";" PATH_STRING "${PATH_STRING}")
STRING(REPLACE ";" "\\;" PATH_STRING "${PATH_STRING}")
set_tests_properties(AccessibleVolumeTest PROPERTIES ENVIRONMENT "PATH=${PATH_STRING}") 
ENDIF()