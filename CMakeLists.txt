cmake_minimum_required(VERSION 2.8)
project(spin_glass_solver)

#--------------------------------------------------
# Setup Output directories
#--------------------------------------------------
SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single Directory for all Executables."
)
SET (CMAKE_LIBRARY_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single Directory for all Libraries"
)
SET (CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single Directory for all static libraries."
)

#--------------------------------------------------
# Flags
#--------------------------------------------------
# needed for all compilers
set(CMAKE_CXX_FLAGS "-std=c++11")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_BUILD_TYPE RELEASE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -pedantic ")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_BUILD_TYPE RELEASE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -pedantic ")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
endif()

#--------------------------------------------------
# WhiplashDB
#--------------------------------------------------
include(FindPkgConfig)
pkg_check_modules(MONGOCXX REQUIRED libmongocxx)
include_directories($ENV{WDB_HOME}/src ${MONGOCXX_INCLUDE_DIRS})
link_directories(${MONGOCXX_LIBRARY_DIRECTORIES} /Users/ethan/lib)
link_libraries(dl)
set(WDB_FLAGS "-w -Wfatal-errors -fPIC -shared")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WDB_FLAGS}")

#--------------------------------------------------
# Exe
#--------------------------------------------------
add_executable(main src/main.cpp)
target_link_libraries(main ${MONGOCXX_LIBRARIES})

#--------------------------------------------------
# Install
#--------------------------------------------------
install(TARGETS main
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        )
SET(CMAKE_INSTALL_PREFIX ${CMAKE_HOME_DIRECTORY})

