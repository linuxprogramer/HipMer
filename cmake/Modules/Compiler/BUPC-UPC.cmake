message("Loading Compiler/BUPC-UPC")
if(NOT CMAKE_UPC_FLAGS_INIT)
  set(CMAKE_UPC_FLAGS_INIT "")
endif()
set(CMAKE_UPC_FLAGS_DEBUG "-Wc,-g -Wl,-g")
set(CMAKE_UPC_FLAGS_RELEASE "-DNDEBUG -O")
set(CMAKE_UPC_FLAGS_MINSIZEREL "${CMAKE_UPC_FLAGS_RELEASE}")
set(CMAKE_UPC_FLAGS_RELWITHDEBINFO "${CMAKE_UPC_FLAGS_RELEASE} ${CMAKE_UPC_FLAGS_DEBUG}")

