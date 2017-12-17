MESSAGE("Building contigMerDepth ${CMAKE_BUILD_TYPE} UPC code using '${CMAKE_UPC_COMPILER} ${CMAKE_UPC_FLAGS}' to compile UPC code")

# for c++ ufx reading code
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wshadow -Wreturn-type")
SET(CMAKE_EXE_LINKER_FLAGS)

ADD_DEFINITIONS(-DMERDEPTH)
ADD_DEFINITIONS(-DUSE_UPC_FOR_COMMON)

include_directories(${CMAKE_SOURCE_DIR}/src/hipmer/contigs/)

SET_SOURCE_FILES_PROPERTIES(contigMerDepth.c PROPERTIES LANGUAGE "UPC" OBJECT_DEPENDS "${CMAKE_SOURCE_DIR}/src/hipmer/common/common.h")

FOREACH(HIPMER_KMER_LEN ${HIPMER_KMER_LENGTHS})
    ADD_EXECUTABLE(contigMerDepth-${HIPMER_KMER_LEN} contigMerDepth.c
                                                     $<TARGET_OBJECTS:KmerObjects-${HIPMER_KMER_LEN}> 
                                                     $<TARGET_OBJECTS:readufx-${HIPMER_KMER_LEN}>
                                                     $<TARGET_OBJECTS:upc_common>
                                                     $<TARGET_OBJECTS:Buffer>
                  )
    SET_TARGET_PROPERTIES(contigMerDepth-${HIPMER_KMER_LEN} PROPERTIES LINKER_LANGUAGE "UPC" LINK_LIBRARIES "z" )
    ADD_DEPENDENCIES(contigMerDepth-${HIPMER_KMER_LEN} REPLACE_VERSION_H)
    TARGET_LINK_LIBRARIES(contigMerDepth-${HIPMER_KMER_LEN} z)
    add_kmer_defs(contigMerDepth-${HIPMER_KMER_LEN})
    INSTALL(TARGETS contigMerDepth-${HIPMER_KMER_LEN} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin/ )
ENDFOREACH()
