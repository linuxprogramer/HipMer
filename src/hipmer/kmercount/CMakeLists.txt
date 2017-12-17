
find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
find_package(ZLIB 1.2.3 REQUIRED)

if (NOT MPI_CXX_FOUND)
    message(FATAL_ERROR "MPI C++ not found. Cannot build ufx.")
endif()

include_directories(${MPI_CXX_INCLUDE_PATH})
set(HIPMER_MPI_DEFINITIONS MPIPARALLEL)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wshadow -Wreturn-type -Wno-unused-result")
set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -std=c99")

message("Building kMerCount MPI utilities: ${CMAKE_CXX_COMPILER} ${CMAKE_CXX_FLAGS}")
add_subdirectory(DataStructures/libbloom)

add_library(MPIType OBJECT MPIType.cpp)
set_target_properties(MPIType PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")

# fqreader
#add_definitions(-DCHECK_SEQS)
add_definitions(-DUSE_MPI_FOR_COMMON -DHEAVYHITTERS -DHISTOGRAM )
#set_source_files_properties(${MPI_FQ_FILES} PROPERTIES LANGUAGE "CXX")
#add_library(MPI_FQ_OBJS OBJECT ${MPI_FQ_FILES})

# zlib
include_directories(${ZLIB_INCLUDE_DIRS})
get_filename_component(ZLIB_BASEPATH ${ZLIB_INCLUDE_DIRS} PATH)
link_directories(${ZLIB_BASEPATH}/lib64 ${ZLIB_BASEPATH}/lib)

# unlink
#set(myTargets unlink)
#add_executable(unlink unlinker.cpp)
#set_target_properties(unlink PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}" LINK_FLAGS "${MPI_CXX_LINK_FLAGS}" LINK_LIBRARIES "${MPI_CXX_LIBRARIES}")
#target_link_libraries(unlink -lrt)

FOREACH(HIPMER_KMER_LEN ${HIPMER_KMER_LENGTHS})
    add_library(KmerObjects-${HIPMER_KMER_LEN} OBJECT Kmer.cpp KmerIterator.cpp hash.cpp KmerMiddle.cpp)
    add_dependencies(KmerObjects-${HIPMER_KMER_LEN} REPLACE_VERSION_H)
    add_kmer_defs_to_lib(KmerObjects-${HIPMER_KMER_LEN})

    add_library(readufx-${HIPMER_KMER_LEN} OBJECT readufx.cpp)
    add_dependencies(readufx-${HIPMER_KMER_LEN} REPLACE_VERSION_H)
    set_target_properties(readufx-${HIPMER_KMER_LEN} PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
    add_kmer_defs_to_lib(readufx-${HIPMER_KMER_LEN})

    add_library(readufx-mmap-${HIPMER_KMER_LEN} OBJECT readufx.cpp)
    add_dependencies(readufx-mmap-${HIPMER_KMER_LEN} REPLACE_VERSION_H)
    set_target_properties(readufx-mmap-${HIPMER_KMER_LEN} PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
    target_compile_definitions(readufx-mmap-${HIPMER_KMER_LEN} PUBLIC -DUSE_MMAP)
    add_kmer_defs_to_lib(readufx-mmap-${HIPMER_KMER_LEN})

    add_library(UFXextended-${HIPMER_KMER_LEN} OBJECT UFXextended.cpp)
    add_dependencies(UFXextended-${HIPMER_KMER_LEN} REPLACE_VERSION_H)
    set_target_properties(UFXextended-${HIPMER_KMER_LEN} PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}" )
    add_kmer_defs_to_lib(UFXextended-${HIPMER_KMER_LEN})

    add_library(bin2textUFX-${HIPMER_KMER_LEN} OBJECT bin2textUFX.cpp)
    add_dependencies(bin2textUFX-${HIPMER_KMER_LEN} REPLACE_VERSION_H)
    set_target_properties(bin2textUFX-${HIPMER_KMER_LEN} PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
    add_kmer_defs_to_lib(bin2textUFX-${HIPMER_KMER_LEN})

    add_library(bin2textCount-${HIPMER_KMER_LEN} OBJECT bin2textCount.cpp)
    add_dependencies(bin2textCount-${HIPMER_KMER_LEN} REPLACE_VERSION_H)
    set_target_properties(bin2textCount-${HIPMER_KMER_LEN} PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
    add_kmer_defs_to_lib(bin2textCount-${HIPMER_KMER_LEN})

    # SNPidentify
    add_executable(SNPidentify-${HIPMER_KMER_LEN} SNPidentify.cpp $<TARGET_OBJECTS:KmerObjects-${HIPMER_KMER_LEN}>)
    add_dependencies(SNPidentify-${HIPMER_KMER_LEN} REPLACE_VERSION_H)
    set_target_properties(SNPidentify-${HIPMER_KMER_LEN} PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
    target_link_libraries(SNPidentify-${HIPMER_KMER_LEN} ${MPI_CXX_LIBRARIES} -lz)
    add_kmer_defs(SNPidentify-${HIPMER_KMER_LEN})

    # ufx
    if (HIPMER_BLOOM64)
      add_definitions(-DHIPMER_BLOOM64)
    endif()
    add_executable(ufx-${HIPMER_KMER_LEN}
        $<TARGET_OBJECTS:MURMURHASH> $<TARGET_OBJECTS:LIBBLOOM>
        $<TARGET_OBJECTS:UFXextended-${HIPMER_KMER_LEN}>
        $<TARGET_OBJECTS:KmerObjects-${HIPMER_KMER_LEN}> 
        $<TARGET_OBJECTS:Buffer>
        $<TARGET_OBJECTS:MPI_FQ_OBJS>
	$<TARGET_OBJECTS:MPIType>
        )
    set_target_properties(ufx-${HIPMER_KMER_LEN}
        PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}"
        COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};${HIPMER_MPI_DEFINITIONS}"
        LINK_FLAGS "${CMAKE_LINK_FLAGS} ${OpenMP_C_FLAGS} ${MPI_CXX_LINK_FLAGS}")
    target_link_libraries(ufx-${HIPMER_KMER_LEN} ${MPI_CXX_LIBRARIES} -lz -lrt)
    add_dependencies(ufx-${HIPMER_KMER_LEN} REPLACE_VERSION_H)
    add_kmer_defs(ufx-${HIPMER_KMER_LEN})
    set(myTargets ${myTargets} ufx-${HIPMER_KMER_LEN})

    # b2tufx
    add_executable(b2tufx-${HIPMER_KMER_LEN}
        $<TARGET_OBJECTS:bin2textUFX-${HIPMER_KMER_LEN}>
        $<TARGET_OBJECTS:KmerObjects-${HIPMER_KMER_LEN}>
        $<TARGET_OBJECTS:Buffer>
        $<TARGET_OBJECTS:MPI_FQ_OBJS>)
    set_target_properties(b2tufx-${HIPMER_KMER_LEN}
        PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};${HIPMER_MPI_DEFINITIONS}"
        COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}"
        LINK_FLAGS "${MPI_CXX_LINK_FLAGS}") 
    add_dependencies(b2tufx-${HIPMER_KMER_LEN} REPLACE_VERSION_H)
    target_link_libraries(b2tufx-${HIPMER_KMER_LEN} ${MPI_LIBRARIES} -lz -lrt)
    add_kmer_defs(b2tufx-${HIPMER_KMER_LEN})

    # b2tcnt
    add_executable(b2tcnt-${HIPMER_KMER_LEN} 
        $<TARGET_OBJECTS:bin2textCount-${HIPMER_KMER_LEN}> 
        $<TARGET_OBJECTS:KmerObjects-${HIPMER_KMER_LEN}>)
    set_target_properties(b2tcnt-${HIPMER_KMER_LEN}
        PROPERTIES COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};${HIPMER_MPI_DEFINITIONS}"
        COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}"
        LINK_FLAGS "${MPI_CXX_LINK_FLAGS}") 
    add_dependencies(b2tcnt-${HIPMER_KMER_LEN} REPLACE_VERSION_H)
    target_link_libraries(b2tcnt-${HIPMER_KMER_LEN} ${MPI_CXX_LIBRARIES} -lz -lrt)
    add_kmer_defs(b2tcnt-${HIPMER_KMER_LEN})

    # ufxdrive
    add_executable(ufxdrive-${HIPMER_KMER_LEN} ufxdriver.c $<TARGET_OBJECTS:readufx-${HIPMER_KMER_LEN}> $<TARGET_OBJECTS:KmerObjects-${HIPMER_KMER_LEN}> )
    set_target_properties(ufxdrive-${HIPMER_KMER_LEN}
        PROPERTIES COMPILE_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} ${MPI_CXX_COMPILE_FLAGS}"
        COMPILE_DEFINITIONS "${COMPILE_DEFINITIONS};USEMPI"
        LINK_FLAGS "${CMAKE_LINK_FLAGS} ${OpenMP_C_FLAGS} ${MPI_CXX_LINK_FLAGS}")
    add_dependencies(ufxdrive-${HIPMER_KMER_LEN} REPLACE_VERSION_H)
    target_link_libraries(ufxdrive-${HIPMER_KMER_LEN} ${MPI_CXX_LIBRARIES} -lz -lrt)
    add_kmer_defs(ufxdrive-${HIPMER_KMER_LEN})

    set(myTargets ${myTargets}
        SNPidentify-${HIPMER_KMER_LEN}
        b2tufx-${HIPMER_KMER_LEN}
        b2tcnt-${HIPMER_KMER_LEN}
        ufxdrive-${HIPMER_KMER_LEN})

ENDFOREACH()

FOREACH(targetName ${myTargets})
    # TODO make tests for these...
    # add_test(test_${targetName} ${targetName} )
    install(TARGETS ${targetName} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin/ )
ENDFOREACH(targetName)

