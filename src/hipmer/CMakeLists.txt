message("Building HipMer")

# SET MPI / UPC definitions and variables
if (NOT HIPMER_KMER_LENGTHS)
#    set(HIPMER_KMER_LENGTHS 15 19 21 31 41 51 63 95 127)
    set(HIPMER_KMER_LENGTHS 19 21 31 51 61)
endif ()
if (NOT HIPMER_MAX_FILE_PATH)
  set(HIPMER_MAX_FILE_PATH 384)
endif()
if (NOT HIPMER_MAX_READ_LEN)
  set(HIPMER_MAX_READ_LEN 1408)
endif()
if (NOT HIPMER_MAX_READ_NAME)
#  set(HIPMER_MAX_READ_NAME 160)
  set(HIPMER_MAX_READ_NAME 80)
endif()
if (NOT HIPMER_BLOOM)
  set(HIPMER_BLOOM64 1)
else()
  if (${HIPMER_BLOOM} STREQUAL "64")
    set(HIPMER_BLOOM64 1)
  endif()
endif()

add_definitions(-DBIN_INPUT) 
add_definitions(-DPROFILE)
ADD_DEFINITIONS(-DCONFIG_USE_COLORS)
ADD_DEFINITIONS(-DCONFIG_SANITY_CHECK)
ADD_DEFINITIONS(-DCONFIG_SHOW_PROGRESS)
add_definitions(-DCONFIG_CHECK_SEQS)
add_definitions(-DMAX_READ_LEN=${HIPMER_MAX_READ_LEN})
add_definitions(-DMAX_READ_NAME_LEN=${HIPMER_MAX_READ_NAME})
add_definitions(-DMAX_FILE_PATH=${HIPMER_MAX_FILE_PATH})
add_definitions(-DNO_PAD)

# include modules for UPC
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")
include(CheckLanguage2)
check_language2(UPC)

#macro(get_kmer_opts)
#    math(EXPR KMER_PACKED_LENGTH "(${HIPMER_KMER_LEN} + 3) / 4")
#    math(EXPR MAX_KMER_SIZE "((${HIPMER_KMER_LEN} + 31) / 32) * 32")
#    set(HIPMER_KMER_DEFS -DMAX_KMER_SIZE=${MAX_KMER_SIZE}
#        -DKMER_LENGTH=${HIPMER_KMER_LEN} -DKMER_PACKED_LENGTH=${KMER_PACKED_LENGTH})
#endmacro()

macro(ADD_KMER_DEFS_TO_LIB TARGET)
#    get_kmer_opts(${HIPMER_KMER_LEN})
    math(EXPR KMER_PACKED_LENGTH "(${HIPMER_KMER_LEN} + 3) / 4")
    math(EXPR MAX_KMER_SIZE "((${HIPMER_KMER_LEN} + 31) / 32) * 32")
    set(HIPMER_KMER_DEFS -DMAX_KMER_SIZE=${MAX_KMER_SIZE}
        -DKMER_LENGTH=${HIPMER_KMER_LEN} -DKMER_PACKED_LENGTH=${KMER_PACKED_LENGTH})
    target_compile_definitions(${TARGET} PUBLIC ${HIPMER_KMER_DEFS})
endmacro()

macro(ADD_KMER_DEFS TARGET)
    add_kmer_defs_to_lib(${TARGET})
    get_target_property(_lang ${TARGET} LINKER_LANGUAGE)
    if (${_lang} STREQUAL "UPC")
        target_add_definitions_to_link_flags(${TARGET})
    endif ()
endmacro()

enable_language(UPC)
if ((NOT CMAKE_UPC_COMPILER) OR (NOT CMAKE_UPC_COMPILER_ID) OR (NOT EXISTS ${CMAKE_UPC_COMPILER}))
    message(FATAL_ERROR "No working UPC compiler found.")
endif()
if (CMAKE_UPC_COMPILER_ID STREQUAL "BUPC")
    message("Detected Berkeley UPC")
    add_definitions( -DUSE_BUPC )
elseif (CMAKE_UPC_COMPILER_ID STREQUAL "CrayUPC")
    message("Detected Cray UPC")
    add_definitions( -DUSE_CRAY )
else()
    message("Unknown UPC version, assuming Berkeley UPC")
    add_definitions( -DUSE_BUPC )
endif ()

if (CMAKE_UPC_USE_PTHREADS)
    SET(CMAKE_UPC_FLAGS "${CMAKE_UPC_FLAGS} -pthreads=${CMAKE_UPC_USE_PTHREADS}")
else()
    SET(CMAKE_UPC_FLAGS "${CMAKE_UPC_FLAGS} -pthreads=0")
endif()

message("Build type: ${CMAKE_BUILD_TYPE}")

if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(CMAKE_UPC_FLAGS "${CMAKE_UPC_FLAGS} -g")
endif()

message("Using UPC: ${CMAKE_UPC_COMPILER} ${CMAKE_UPC_FLAGS}")

add_subdirectory(fqreader)
add_subdirectory(kmercount)
add_subdirectory(contigs)
add_subdirectory(meraligner)
add_subdirectory(scaffolding)
add_subdirectory(loadfq)
add_subdirectory(gapclosing)
add_subdirectory(common)

# FIXME: the perl scripts shouldn't need to be installed
install(PROGRAMS 
    run_hipmer.sh 
    rerun_stage.sh
    rerun_single_stage.sh
    get_stage_timings.sh
    get_diagnostics.sh
    ${CMAKE_SOURCE_DIR}/test/hipmer/run_hipmer_validation.sh
    ${CMAKE_SOURCE_DIR}/test/hipmer/check_validation_results.sh
    ${CMAKE_SOURCE_DIR}/test/hipmer/check_ecoli_results.sh
    ${CMAKE_SOURCE_DIR}/test/hipmer/check_results.sh
    ${CMAKE_SOURCE_DIR}/test/hipmer/hipmer_setup_validation_data.sh
    ${CMAKE_SOURCE_DIR}/test/hipmer/run_hipmer_validation.sh
    ${CMAKE_SOURCE_DIR}/test/hipmer/hipmer_setup_ecoli_data.sh 
    ${CMAKE_SOURCE_DIR}/test/hipmer/hipmer_setup_human_data.sh
    ${CMAKE_SOURCE_DIR}/test/hipmer/hipmer_setup_human_subsets_data.sh
    ${CMAKE_SOURCE_DIR}/test/hipmer/hipmer_setup_chr14_data.sh
    ${CMAKE_SOURCE_DIR}/src/hipmer/scaffolding/merger/merge_ono_inputs.sh
    DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
