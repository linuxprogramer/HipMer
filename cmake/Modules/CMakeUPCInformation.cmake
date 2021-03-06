# https://github.com/Kitware/CMake/blob/master/Modules/CMakeAddNewLanguage.txt
#
# set up rule variables for LANG :
#   CMAKE_(LANG)_CREATE_SHARED_LIBRARY
#   CMAKE_(LANG)_CREATE_SHARED_MODULE
#   CMAKE_(LANG)_CREATE_STATIC_LIBRARY
#   CMAKE_(LANG)_COMPILE_OBJECT
#   CMAKE_(LANG)_LINK_EXECUTABLE
#
# variables supplied by the generator at use time
# <TARGET>
# <TARGET_BASE> the target without the suffix
# <OBJECTS>
# <OBJECT>
# <LINK_LIBRARIES>
# <FLAGS>
# <LINK_FLAGS>


message("Building CMake UPC Information and Rules: ${CMAKE_UPC_COMPILER} ${CMAKE_UPC_FLAGS}")
set(CMAKE_UPC_OUTPUT_EXTENSION .o)

set(_INCLUDED_FILE 0)
set(CMAKE_UPC_BASE_NAME)
get_filename_component(CMAKE_UPC_BASE_NAME ${CMAKE_UPC_COMPILER} NAME_WE)

# include specific environment for compiler id / platform
message("Looking for Compiler/${CMAKE_UPC_COMPILER_ID}-UPC or Platform/${CMAKE_SYSTEM_NAME}-${CMAKE_UPC_COMPILER_ID}-UPC or Platform/${CMAKE_SYSTEM_NAME}-${CMAKE_UPC_BASE_NAME} or Platform/${CMAKE_SYSTEM_NAME}")
if(CMAKE_UPC_COMPILER_ID)
  include(Compiler/${CMAKE_UPC_COMPILER_ID}-UPC OPTIONAL)
  include(Platform/${CMAKE_SYSTEM_NAME}-${CMAKE_UPC_COMPILER_ID}-UPC OPTIONAL RESULT_VARIABLE _INCLUDED_FILE)
endif()
if(NOT _INCLUDED_FILE)
  include(Platform/${CMAKE_SYSTEM_NAME}-${CMAKE_UPC_BASE_NAME} OPTIONAL RESULT_VARIABLE _INCLUDED_FILE)
endif()
if(NOT _INCLUDED_FILE)
  include(Platform/${CMAKE_SYSTEM_NAME} OPTIONAL)
endif()

# set UPC languange parameters
set(CMAKE_INCLUDE_FLAG_UPC "-I")
set(CMAKE_UPC_LINKER_PREFERENCE "UPC")

#if(CMAKE_BUILD_TYPE STREQUAL "Release")
#  set(CMAKE_UPC_FLAGS "${CMAKE_UPC_FLAGS} ${CMAKE_UPC_FLAGS_RELEASE}")
#elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
#  set(CMAKE_UPC_FLAGS "${CMAKE_UPC_FLAGS} ${CMAKE_UPC_FLAGS_DEBUG}")
#endif()

set(CMAKE_UPC_FLAGS "${CMAKE_UPC_FLAGS_INIT} ${CMAKE_UPC_FLAGS}" CACHE STRING "Flags used by the UPC compiler during all build types.")

#   CMAKE_(LANG)_CREATE_SHARED_LIBRARY
#   CMAKE_(LANG)_CREATE_SHARED_MODULE
#   CMAKE_(LANG)_CREATE_STATIC_LIBRARY
#   CMAKE_(LANG)_COMPILE_OBJECT
#   CMAKE_(LANG)_LINK_EXECUTABLE

#set(CMAKE_UPC_CREATE_SHARED_LIBRARY
#      "<CMAKE_UPC_COMPILER> <CMAKE_UPC_COMPILER_ARG1> <CMAKE_SHARED_LIBRARY_UPC_FLAGS> <LANGUAGE_COMPILE_FLAGS> <LINK_FLAGS> <CMAKE_SHARED_LIBRARY_CREATE_UPC_FLAGS> <SONAME_FLAG><TARGET_SONAME> -o <TARGET> <OBJECTS> <LINK_LIBRARIES>")
#
#set(CMAKE_UPC_CREATE_SHARED_MODULE ${CMAKE_UPC_CREATE_SHARED_LIBRARY})

set(CMAKE_UPC_COMPILE_OBJECT 
    "<CMAKE_UPC_COMPILER> <DEFINES> <FLAGS> -o <OBJECT> -c <SOURCE>")

set(CMAKE_UPC_LINK_EXECUTABLE
    "<CMAKE_UPC_COMPILER> <FLAGS> <CMAKE_UPC_LINK_FLAGS> <LINK_FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES>")


