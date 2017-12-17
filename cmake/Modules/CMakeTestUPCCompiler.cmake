# https://github.com/Kitware/CMake/blob/master/Modules/CMakeAddNewLanguage.txt
#
#test the compiler and set:

if(CMAKE_UPC_COMPILER_FORCED)
  message("Not testing UPC compiler: CMAKE_UPC_COMPILER_FORCED is set")
  set(CMAKE_UPC_COMPILER_WORKS TRUE)
  return()
endif()

message("Testing UPC Compiler")
include(CMakeTestCompilerCommon)

# Remove any cached result from an older CMake version.
# We now store this in CMakeCXXCompiler.cmake.
unset(CMAKE_UPC_COMPILER_WORKS CACHE)

if(NOT CMAKE_UPC_COMPILER_WORKS) 
  PrintTestCompilerStatus("UPC" "")

  EXECUTE_PROCESS(COMMAND ${CMAKE_UPC_COMPILER} --help
                          RESULT_VARIABLE CAN_UPCC
                          OUTPUT_VARIABLE LOG_UPCC)
  IF(NOT CAN_UPCC EQUAL 0)
    MESSAGE("Could not execute '${CMAKE_UPC_COMPILER} ${CMAKE_UPC_FLAGS} --help': ${LOG_UPCC}")
    return()
  ENDIF()

  set(_test_file ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testUPCCompiler.upc)
  file(WRITE ${_test_file}
  "#include <upc.h>\n"
  "#include <stdio.h>\n"
  "#include <stdlib.h>\n"
  "int main() {\n"
  "	printf(\"Hello from %d of %d\\n\", MYTHREAD, THREADS);\n"
  "     return 0;\n"
  "}\n")
  try_compile(CMAKE_UPC_COMPILER_WORKS ${CMAKE_BINARY_DIR} ${_test_file}
              OUTPUT_VARIABLE __CMAKE_UPC_COMPILER_OUTPUT)
  set(CMAKE_UPC_COMPILER_WORKS ${CMAKE_UPC_COMPILER_WORKS})
  unset(CMAKE_UPC_COMPILER_WORKS CACHE)
  set(UPC_TEST_WAS_RUN 1)
endif()

if(NOT CMAKE_UPC_COMPILER_WORKS)
  PrintTestCompilerStatus("UPC" " -- broken")
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determining if the UPC compiler works failed with "
    "the following output:\n${__CMAKE_UPC_COMPILER_OUTPUT}\n\n")
  message("The UPC compiler \"${CMAKE_UPC_COMPILER}\" "
    "is not able to compile a simple test program.\nIt fails "
    "with the following output:\n ${__CMAKE_UPC_COMPILER_OUTPUT}\n\n"
    "CMake will not be able to correctly generate this project.")
  return()
else()
  if(UPC_TEST_WAS_RUN)
    PrintTestCompilerStatus("UPC" " -- works")
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
      "Determining if the UPC compiler works passed with "
      "the following output:\n${__CMAKE_UPC_COMPILER_OUTPUT}\n\n")
  endif()
endif()

unset(__CMAKE_UPC_COMPILER_OUTPUT)

