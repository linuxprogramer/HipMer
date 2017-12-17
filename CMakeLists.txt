# The name of our project is "HipMer". CMakeLists files in this project can
# refer to the root source directory of the project as ${HIPMER_SOURCE_DIR} and
# to the root binary directory of the project as ${HIPMER_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8)
project (HipMer)
set(HIPMER_SOURCE_DIR ${CMAKE_SOURCE_DIR})

enable_testing()

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
	    FORCE)
endif()
message("Building for ${CMAKE_BUILD_TYPE}")
if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    add_definitions(-DDEBUG)
endif()

if(NOT HIPMER_READ_BUFFER)
  set(HIPMER_READ_BUFFER 1048576)
endif()
add_definitions(-DKSEQ_BUFFER_SIZE=${HIPMER_READ_BUFFER})

# get the version from git
set(HIPMER_VERSION_FILE_TEMPLATE "${CMAKE_SOURCE_DIR}/src/hipmer/common/version.h.in")
set(HIPMER_VERSION_FILE "${CMAKE_SOURCE_DIR}/src/hipmer/common/version.h")
set(HIPMER_GIT_VERSION_FILE "${CMAKE_SOURCE_DIR}/HIPMER_VERSION")
find_package(Git QUIET)
if (GIT_FOUND AND IS_DIRECTORY ${CMAKE_SOURCE_DIR}/.git)
	add_custom_target(GET_GIT_VERSION ALL
			COMMAND ${GIT_EXECUTABLE} describe --tags > ${HIPMER_GIT_VERSION_FILE}
			WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
			COMMENT "Retrieving git version")
else()
	if (NOT EXISTS ${HIPMER_GIT_VERSION_FILE})
		add_custom_target(GET_GIT_VERSION
				COMMAND echo "UNKNOWN-VERSION" > ${HIPMER_GIT_VERSION_FILE})
	else()
		add_custom_target(GET_GIT_VERSION
				COMMAND /bin/true)
	endif()
endif()
set(GIT_PREFIX HIPMER)
configure_file(${CMAKE_SOURCE_DIR}/cmake/Modules/makeVersionFile.cmake.in ${CMAKE_BINARY_DIR}/makeVersionFile.cmake @ONLY)
add_custom_target(REPLACE_VERSION_H ALL
		COMMAND ${CMAKE_COMMAND}
			-DHIPMER_GIT_VERSION_FILE=${HIPMER_GIT_VERSION_FILE}
			-DHIPMER_VERSION_FILE=${HIPMER_VERSION_FILE}
			-DHIPMER_VERSION_FILE_TEMPLATE=${HIPMER_VERSION_FILE_TEMPLATE}
			-P ${CMAKE_BINARY_DIR}/makeVersionFile.cmake
		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
		DEPENDS GET_GIT_VERSION
		COMMENT "Building ${HIPMER_VERSION_FILE}")
add_dependencies(REPLACE_VERSION_H GET_GIT_VERSION)

if(NOT CMAKE_INSTALL_PREFIX)
    message(FATAL_ERROR "Please specify a CMAKE_INSTALL_PREFIX")
endif()
if(NOT IS_ABSOLUTE ${CMAKE_INSTALL_PREFIX})
    get_filename_component(fullpath ${CMAKE_INSTALL_PREFIX} ABSOLUTE)
    set (CMAKE_INSTALL_PREFIX "${fullpath}")
endif()
message("Installing into ${CMAKE_INSTALL_PREFIX}")

# if no prefix was passed, install into the build directory
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}" CACHE PATH "default install path" FORCE )
endif()

# check compiler and set appropriate linker flags
include( CheckCXXCompilerFlag )
check_cxx_compiler_flag("-static-libgcc" CXX_HAS_STATIC_LIBGCC)
check_cxx_compiler_flag("-static-libstdc++" CXX_HAS_STATIC_LIBSTDCXX)
if (CXX_HAS_STATIC_LIBSTDCXX)
    set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++" )
else ()
    set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static" )
endif ()

set(BUILD_STATIC_LIBS ON)
set(BUILD_SHARED_LIBS OFF)

check_cxx_compiler_flag("-std=c++0x" CXX_SUPPORTS_0X)
check_cxx_compiler_flag("-std=c++11" CXX_SUPPORTS_11)

if(CXX_SUPPORTS_11)
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 " )
elseif(CXX_SUPPORTS_0X)
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x " )
else()
    message( FATAL_ERROR "Incompatible c++ compiler.  c++0x or c++11 support is required. CMake will exit." )
endif()

set( CMAKE_VERBOSE_MAKEFILE off)


message(STATUS "Testing the environment..")

# SeqDB with MPI support
find_package(MPI)

# include modules for UPC
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/Modules)
include(CheckLanguage2)
check_language2(UPC)
if (CMAKE_UPC_COMPILER)
    enable_language(UPC OPTIONAL)
    if(CMAKE_UPC_COMPILER_LOADED AND CMAKE_UPC_COMPILER AND EXISTS ${CMAKE_UPC_COMPILER})
        add_subdirectory(src/hipmer)
    else() 
        message(WARNING "UPC compiler was not found to be working, not building HipMer.")
    endif()
else()
    message(WARNING "UPC compiler was not found, not building HipMer")
endif()
if (NOT CMAKE_UPC_COMPILER)
   message(warning "UPC could not be found.  Please make sure it is properly installed and in the PATH. (see http://upc.lbl.gov)
If cmake can not automatically find it, try specifying:
    export HIPMER_BUILD_OPTS=\"-DCMAKE_UPC_COMPILER_ENV_VAR=upcc\"
")
endif()


#unit test scripts
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/test_data.tmp)

set(TEST_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test)

install(DIRECTORY src/perl/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
        USE_SOURCE_PERMISSIONS
        PATTERN "canonical_assembly.pl")
install(DIRECTORY src/perl/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
        USE_SOURCE_PERMISSIONS
        PATTERN "findDMin2.pl")

install(DIRECTORY test/ DESTINATION ${CMAKE_INSTALL_PREFIX}/etc/meraculous PATTERN ".*" EXCLUDE PATTERN "merBlast" EXCLUDE)
install(DIRECTORY doc/ DESTINATION ${CMAKE_INSTALL_PREFIX}/share/doc/meraculous PATTERN ".*" EXCLUDE)
install(FILES ${HIPMER_GIT_VERSION_FILE} README.md LEGAL.txt LICENSE.txt DESTINATION ${CMAKE_INSTALL_PREFIX})

# get the version from git

#	message("Determining version from git...")
#	set(GIT_PREFIX "HIPMER")
#	configure_file(${CMAKE_SOURCE_DIR}/cmake/Modules/GetGitInfo.cmake.in ${CMAKE_BINARY_DIR}/GetGitInfo.cmake @ONLY)
#	set
	#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_BINARY_DIR})
	#include(GetGitInfo)
#	GetGitInfo()
#	add_custom_target(getgit ALL
#			COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}/GetGitInfo.cmake
#			)
#	set_source_files_properties(${version_file} PROPERTIES OBJECT_DEPENDS getgit)
#else()
#	if(IS_FILE ${CMAKE_SOURCE_DIR}/HIPMER_VERSION)
#		message("Determining version from HIPMER_VERSION files...")
#		file(READ "${CMAKE_SOURCE_DIR}/HIPMER_VERSION" GIT_TAG)
#		file(READ "${CMAKE_SOURCE_DIR}/HIPMER_VERSION.txt" GIT_STRING)
#	else()
#		message("No git nor HIPMER_VERSION files are present...")
#	endif()
#endif()

#if (GIT_TAG)
#	# if not in a git repo and HIPMER_VERSION and HIPMER_VERSION.txt do not exist, lets hope that common.h exists!
#	set(HIPMER_VERSION "${GIT_TAG}")
#	message("Configured HipMer version: ${GIT_STRING}")
#	configure_file("${version_file_in}" "${version_file}")
#else()
#	if (NOT EXISTS ${version_file})
#		message(FATAL_ERROR "Could not determine the HIPMER_VERSION through git or the HIPMER_VERSION file(s).  please supply a GIT_TAG")
#	endif()
#endif()
#configure_file("${version_file_in}" "${version_file}" @ONLY)


