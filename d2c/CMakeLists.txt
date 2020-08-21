CMAKE_MINIMUM_REQUIRED(VERSION 3.5)

set(SRC_PATH ${PROJECT_SOURCE_DIR})
#set(CXX_INCLUDE_PATH ${SRC_PATH}/cxx/include)
#set(CXX_LIB_PATH ${SRC_PATH}/cxx/lib)


set(CMAKE_CXX_FLAGS "-std=c++17 -O2 -DNDEBUG -W -Wall -pedantic -fopenmp -lpthread -lrt")
#set(CMAKE_CXX_FLAGS "-std=c++17 -g -ggdb -gstabs+ -pedantic -fopenmp -lpthread -lrt")


include_directories(../)
link_directories(${INSTALL_PATH}/lib)

set (src
    utility.cpp
    d2c.cpp
    samReader.cpp
    bamCat.cpp
    saturation.cpp

    main.cpp
    )

set(EXECUTABLE_OUTPUT_PATH ${INSTALL_PATH}/bin)
add_executable(app ${src})
#add_definitions(-DTIMING)
#add_definitions(-DDEBUG_)

# CLI11 use std::filesystem, but for gcc<9,it needs a seperate library libstdc++fs
# set CLI11_HAS_FILESYSTEM=0 for prevent using std::filesystem
add_definitions(-DCLI11_HAS_FILESYSTEM=0)

target_link_libraries(app
    density
    z
    bz2
    hts
    )

set_target_properties(app PROPERTIES OUTPUT_NAME "d2c")
install(FILES 19_makeKneePlots.R barcode.list DESTINATION bin)