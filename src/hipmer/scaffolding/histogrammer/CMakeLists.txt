MESSAGE("Building merAlignerAnalyzer ${CMAKE_BUILD_TYPE} UPC code using '${CMAKE_UPC_COMPILER} ${CMAKE_UPC_FLAGS}' to compile UPC code")

ADD_DEFINITIONS(-DUSE_UPC_FOR_COMMON)

SET_SOURCE_FILES_PROPERTIES(merAlignerAnalyzer.c PROPERTIES LANGUAGE "UPC" 
                                                            OBJECT_DEPENDS "${CMAKE_SOURCE_DIR}/src/hipmer/common/common.h")
ADD_EXECUTABLE(merAlignerAnalyzer merAlignerAnalyzer.c $<TARGET_OBJECTS:upc_common> )
ADD_DEPENDENCIES(merAlignerAnalyzer REPLACE_VERSION_H)
INSTALL(TARGETS merAlignerAnalyzer DESTINATION ${CMAKE_INSTALL_PREFIX}/bin/ )
