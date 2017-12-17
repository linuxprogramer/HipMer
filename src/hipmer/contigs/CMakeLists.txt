MESSAGE("Building contigs ${CMAKE_BUILD_TYPE} UPC code using '${CMAKE_UPC_COMPILER} ${CMAKE_UPC_FLAGS}' to compile UPC code")

# for c++ ufx reading code
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wshadow -Wreturn-type")

# fix linking when objects are from multiple languages
set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "")
set(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "")
set(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "stdc++")
set(CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "")

SET(CMAKE_EXE_LINKER_FLAGS)
SET_SOURCE_FILES_PROPERTIES(meraculous.c PROPERTIES LANGUAGE "UPC" OBJECT_DEPENDS "${CMAKE_SOURCE_DIR}/src/hipmer/common/common.h")

ADD_DEFINITIONS(-DSTORE_OPT -DDETAILED_BUILD_PROFILE)
ADD_DEFINITIONS(-DUSE_UPC_FOR_COMMON)
ADD_DEFINITIONS(-DMERACULOUS -DSYNC_PROTOCOL)
IF(HIPMER_VERBOSE)
    ADD_DEFINITIONS(-DVERBOSE=${HIPMER_VERBOSE})
ENDIF()

FOREACH(HIPMER_KMER_LEN ${HIPMER_KMER_LENGTHS})
    ADD_EXECUTABLE(meraculous-${HIPMER_KMER_LEN} meraculous.c 
        $<TARGET_OBJECTS:KmerObjects-${HIPMER_KMER_LEN}> 
        $<TARGET_OBJECTS:readufx-${HIPMER_KMER_LEN}>  
        $<TARGET_OBJECTS:upc_common>
        $<TARGET_OBJECTS:Buffer> )
    SET_TARGET_PROPERTIES(meraculous-${HIPMER_KMER_LEN} PROPERTIES LINKER_LANGUAGE "UPC" LINK_LIBRARIES "z") 
    ADD_DEPENDENCIES(meraculous-${HIPMER_KMER_LEN} REPLACE_VERSION_H)
    INSTALL(TARGETS meraculous-${HIPMER_KMER_LEN} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin/ )
    add_kmer_defs(meraculous-${HIPMER_KMER_LEN})
ENDFOREACH()

#SET(TEST_THREADS 4)
#test: meraculous-19
#	echo "Starting Unit Test" \
#	&& rm -f {output,log}_*_validate.fa \
#	&& upcrun -n 1 ./meraculous-19 -i test.ufx -o validate.fa -m 19 \
#	&& ./canonical_assembly.pl output_*_validate.fa | diff - test.ufx.fa \
#	&& rm -f {output,log}_*_validate.fa \
#	&& upcrun -n $(TEST_THREADS) ./meraculous-19 -i test.ufx -o validate.fa -m 19 \
#	&& ./canonical_assembly.pl output_*_validate.fa | diff - test.ufx.fa \
#	&& echo Test Passed
#
#% : %.gz
#	gunzip -c $< > $@
#
#test-ecoli: ecoli_in.19.ufx test_ecoli-19.fa test_ecoli-200.fa test_ecoli-250.fa meraculous-19
#	echo "Starting Ecoli test assembly" \
#	&& rm -f {output,log}_*_test_ecoli-{19,200,250}.fa \
#	&& upcrun -n $(TEST_THREADS) ./meraculous-19 -i ecoli_in.19.ufx -o test_ecoli-19.fa -m 19 \
#	&& ./canonical_assembly.pl output_*_test_ecoli-19.fa | diff -q - test_ecoli-19.fa \
#	&& upcrun -n $(TEST_THREADS) ./meraculous-19 -i ecoli_in.19.ufx -o test_ecoli-200.fa -m 200 \
#	&& ./canonical_assembly.pl output_*_test_ecoli-200.fa | diff -q - test_ecoli-200.fa  \
#	&& upcrun -n $(TEST_THREADS) ./meraculous-19 -i ecoli_in.19.ufx -o test_ecoli-250.fa -m 250 \
#	&& ./canonical_assembly.pl output_*_test_ecoli-250.fa | diff -q - test_ecoli-250.fa \
#	&& echo "Ecoli assemblies passed"
#
#
#clean :
#	rm -f *.o
#	rm -rf mem_reporter-[1-9]*
#	rm -rf meraculous-[1-9]*
#	rm -rf contigMerDepth-[1-9]*
#	rm -rf contigEndAnalyzer-[1-9]*
#	rm -f ecoli_in.19.txt.ufx test_ecoli-19.fa test_ecoli-200.fa test_ecoli-250.fa
