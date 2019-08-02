cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(afqmcConstraintPath C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_testing()

find_package(AFQMCLAB REQUIRED)

set(CMAKE_CXX_COMPILER ${AFQMCLAB_COMPILER})
include_directories(${AFQMCLAB_INCLUDE_DIRS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${AFQMCLAB_FLAGS}")
add_definitions(${AFQMCLAB_DEFINITIONS})

set(SOURCE_HEADERS
        include/ghf.h
        include/ghfMethod.h
        include/eigen
        include/LBFGSpp-master/include/LBFGS.h
        include/vmpBCS.h
        include/BCS.h
        include/BCSSDOperation.h
        include/HubbardSOCMeasureCommuteBCSSD.h
        include/HubbardSOCMeasureObserveBCSSD.h
        include/HubbardSOCMeasureFixedBCSSD.h
        include/BCSXiao.h
        include/NiupNidnBCSOperation.h
        include/hopBCSOperation.h
        include/afqmcBackPropagationPop.h
        include/afqmcConstraintPath.h
        include/afqmcConstraintPathDefine.h
        include/afqmcConstraintPathMethod.h
        include/afqmcWalkerPop.h
        )

set(SOURCE_FILES
        source/ghf.cpp
        source/ghfDetail.cpp
        source/ghfMethod.cpp
        source/vmpBCS.cpp
        source/NiupNidnBCSOperation.cpp
        source/hopBCSOperation.cpp
        source/BCS.cpp
        source/HubbardSOCMeasureCommuteBCSSD.cpp
        source/HubbardSOCMeasureObserveBCSSD.cpp
        source/HubbardSOCMeasureFixedBCSSD.cpp
        source/BCSSDOperation.cpp
        source/afqmcBackPropagationPop.cpp
        source/afqmcConstraintPath.cpp
        source/afqmcConstraintPathMeasure.cpp
        source/afqmcConstraintPathMethod.cpp
        source/afqmcConstraintPathProjection.cpp
        source/afqmcConstraintPathSetWalker.cpp
        source/afqmcConstraintPathBackPropagation.cpp
        source/afqmcWalkerPop.cpp
        source/main.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_HEADERS} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} ${AFQMCLAB_LIBRARIES})

add_custom_target(rmdat COMMAND rm -rf *dat)
