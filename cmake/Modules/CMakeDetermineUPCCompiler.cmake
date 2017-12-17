# https://github.com/Kitware/CMake/blob/master/Modules/CMakeAddNewLanguage.txt
#
# this should find the compiler for LANG and configure CMake(LANG)Compiler.cmake.in

# Sets the following variables:
#   CMAKE_UPC_COMPILER
#   CMAKE_UPC_COMPILER_ID
#
# Allows predefined CMAKE_UPC_COMPILER
# or CMAKE_UPC_COMPILER_INIT for multi-argument compiler (such as "cc -h upc")
#

message("CMakeDetermineUPCCompiler starting")
include(${CMAKE_ROOT}/Modules/CMakeDetermineCompiler.cmake)

# set CMAKE_UPC_COMPILER_INIT with potential names
 
if(NOT CMAKE_UPC_COMPILER)
  set(CMAKE_UPC_FLAGS_INIT "")
  if(CMAKE_UPC_COMPILER_INIT)
    message("Using CMAKE_UPC_COMPILER_INIT '${CMAKE_UPC_COMPILER_INIT}' for upc compiler discovery")
    get_filename_component(CMAKE_UPC_COMPILER "${CMAKE_UPC_COMPILER_INIT}" PROGRAM PROGRAM_ARGS CMAKE_UPC_FLAGS_ENV_INIT)

    if(CMAKE_UPC_FLAGS_ENV_INIT)
      set(CMAKE_UPC_FLAGS_INIT "${CMAKE_UPC_FLAGS_INIT} ${CMAKE_UPC_FLAGS_ENV_INIT}" CACHE STRING "First argument to UPC compiler")
    endif()
  else()
    # Load system-specific compiler preferences for this language.
    include(Platform/${CMAKE_SYSTEM_NAME}-UPC OPTIONAL)
    if(NOT CMAKE_UPC_COMPILER_NAMES)
      # guess upcc or the environment variable CC
      get_filename_component(CMAKE_UPC_COMPILER_NAMES upcc PROGRAM)
    endif()

    message("Finding UPC compiler: ${CMAKE_UPC_COMPILER_NAMES}")
    _cmake_find_compiler(UPC)
  endif()

else()
  message("Using user-supplied CMAKE_UPC_COMPILER: ${CMAKE_UPC_COMPILER}")
endif()

if (NOT CMAKE_UPC_COMPILER OR NOT EXISTS ${CMAKE_UPC_COMPILER})
  message("Could not find a UPC compiler. Aborting UPC detection")
  return()
endif()

# set the CMAKE_UPC_COMPILER_ID
if(CMAKE_UPC_COMPILER MATCHES ".*/upcc$")
  set(CMAKE_UPC_COMPILER_ID BUPC)
  message("UPC compiler is Berkeley UPC")
  if(DEFINED BUPC_TRANSLATOR_FLAG)
    set(CMAKE_UPC_FLAGS_INIT "${CMAKE_UPC_FLAGS_INIT} ${BUPC_TRANSLATOR_FLAG}")
  else()
    message("Checking BUPC for -cupc2c translator")
    set(_test_file ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testUPCCompiler.upc)
    file(WRITE ${_test_file}
    "#include <upc.h>\n"
    "#include <stdio.h>\n"
    "#include <stdlib.h>\n"
    "int main() {\n"
    "     printf(\"Hello from %d of %d\\n\", MYTHREAD, THREADS);\n"
    "     return 0;\n"
    "}\n")
    set(test_args ${CMAKE_UPC_FLAGS_INIT} -cupc2c -o ${_test_file}.a.out ${_test_file})
    message("${CMAKE_UPC_COMPILER} ${test_args}")
    execute_process(COMMAND ${CMAKE_UPC_COMPILER} ${test_args}
              RESULT_VARIABLE __CUPC2C_WORKS
              OUTPUT_VARIABLE __CMAKE_UPC_COMPILER_OUTPUT)
    if(__CUPC2C_WORKS EQUAL 0)
      message("Using upc2c Berkeley UPC translator")
      set(CMAKE_UPC_FLAGS_INIT "${CMAKE_UPC_FLAGS_INIT} -cupc2c")
      set(BUPC_TRANSLATOR_FLAG "-cupc2c" CACHE STRING "Berkeley UPC compiler flag")
    else()
      message("Could not use upc2c Berkeley UPC translator: ${__CMAKE_UPC_COMPILER_OUTPUT}")
      set(BUPC_TRANSLATOR_FLAG "" CACHE STRING "Berkeley UPC compiler flag")
    endif()
  endif()
  
elseif(CMAKE_UPC_COMPILER MATCHES ".*/gupc$")
  set(CMAKE_UPC_COMPILER_ID GUPC)
  message("UPC compiler is Gnu UPC")
elseif(CMAKE_UPC_COMPILER MATCHES ".*/cc$")
  set(CMAKE_UPC_COMPILER_ID CrayUPC)
  message("UPC compiler is Cray UPC")
else()
  message("Unknown UPC compiler.  Aborting UPC detection")
  set(CMAKE_UPC_COMPILER)
  set(CMAKE_UPC_COMPILER_ID)
  return()
endif()


if (NOT CMAKE_UPC_FLAGS_INIT)
  set(CMAKE_UPC_FLAGS_INIT "")
endif()

set(CMAKE_UPC_SOURCE_FILE_EXTENSIONS "upc")

message("Discovered UPC Compiler (${CMAKE_UPC_COMPILER_ID}): ${CMAKE_UPC_COMPILER} ${CMAKE_UPC_FLAGS_INIT}")

set(_UPC_CMAKE_IN ${CMAKE_SOURCE_DIR}/cmake/Modules/CMakeUPCCompiler.cmake.in)
if(NOT EXISTS ${_UPC_CMAKE_IN})
  set(_UPC_CMAKE_IN ${HIPMER_SOURCE_DIR}/cmake/Modules/CMakeUPCCompiler.cmake.in)
endif()
configure_file(${_UPC_CMAKE_IN}
               ${CMAKE_PLATFORM_INFO_DIR}/CMakeUPCCompiler.cmake @ONLY)
set(CMAKE_UPC_COMPILER_ENV_VAR "UPC_COMPILER")

