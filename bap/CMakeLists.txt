CMAKE_MINIMUM_REQUIRED(VERSION 3.5)

set(SRC_PATH ${PROJECT_SOURCE_DIR})
#set(CXX_INCLUDE_PATH ${SRC_PATH}/cxx/include)
#set(CXX_LIB_PATH ${SRC_PATH}/cxx/lib)

#set(CMAKE_CXX_FLAGS "-std=c++14 -O3 -DNDEBUG -W -Wall -pedantic -fopenmp -lpthread -lrt")

# set(CMAKE_CXX_FLAGS "-std=c++17 -O2 -DNDEBUG -W -Wall -pedantic -fopenmp -lpthread -lrt")

#set(CMAKE_CXX_FLAGS "-std=c++17 -pg -O2 -DNDEBUG -W -Wall -pedantic -fopenmp -lpthread -lrt")

set(CMAKE_CXX_FLAGS "-std=c++17 -g -ggdb -gstabs+ -pedantic -fopenmp -lpthread -lrt")
#set(CMAKE_CXX_FLAGS "-std=c++17 -gdwarf-2 -pedantic -fopenmp -lpthread -lrt")

link_directories(${INSTALL_PATH}/lib)

set (src
    # assemble_fragments.cpp
    # annotate_fragments.cpp
    # gz_io.cpp
    # compute_bap_stat_chr.cpp
    utility.cpp
    # determine_barcode_merges.cpp
    # reannotate_fragments_get_summary_stats_chr.cpp
    # reannotate_bam_file.cpp
    # final_qc_se.cpp
    bap.cpp
    samReader.cpp
    bamCat.cpp

    main.cpp
    )

set(EXECUTABLE_OUTPUT_PATH ${INSTALL_PATH}/bin)
add_executable(app_bap ${src})
#add_definitions(-DTIMING)
#add_definitions(-DDEBUG_)

# CLI11 use std::filesystem, but for gcc<9,it needs a seperate library libstdc++fs
# set CLI11_HAS_FILESYSTEM=0 for prevent using std::filesystem
add_definitions(-DCLI11_HAS_FILESYSTEM=0)

target_link_libraries(app_bap
    z
    bz2
    hts
    )

set_target_properties(app_bap PROPERTIES OUTPUT_NAME "bap")
install(FILES 10b_knee_execute.R 00_knee_CL.R DESTINATION bin)