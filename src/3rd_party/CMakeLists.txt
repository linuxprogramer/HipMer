
add_executable (fasta_stats 
	       fasta_stats/fasta_stats.cc
	       fasta_stats/itoa.cc
	       fasta_stats/open_compressed.cc
	       fasta_stats/pretty_print.cc)


set_target_properties(fasta_stats
		PROPERTIES COMPILE_FLAGS "-I${CMAKE_CURRENT_SOURCE_DIR}/fasta_stats/include -g -O2 -W -Wall -Wno-unused-parameter -Wpointer-arith -Wshadow -Wundef" )


#installation
install(TARGETS fasta_stats DESTINATION ${CMAKE_INSTALL_PREFIX}/bin/)
