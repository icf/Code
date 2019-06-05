cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(bcs C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_testing()

find_package(AFQMCLAB REQUIRED)

set(CMAKE_CXX_COMPILER ${AFQMCLAB_COMPILER})
include_directories(${AFQMCLAB_INCLUDE_DIRS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${AFQMCLAB_FLAGS}")
add_definitions(${AFQMCLAB_DEFINITIONS})

set(SOURCE_HEADERS
        include/bcs.h
        include/bcsMethod.h
        )

set(SOURCE_FILES
        source/bcs.cpp
        source/bcsDetail.cpp
        source/bcsMethod.cpp
        source/bcsTuneMu.cpp
        source/bcsTunePairing.cpp
        source/main.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_HEADERS} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} ${AFQMCLAB_LIBRARIES})

add_custom_target(rmdat COMMAND rm -rf *dat)