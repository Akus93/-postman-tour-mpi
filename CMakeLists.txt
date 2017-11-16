cmake_minimum_required(VERSION 3.5)
project(pr_graphs)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_CXX_STANDARD 11)
set(SOURCE_FILES
        src/main.cpp
        )

include_directories(.)
include_directories(src)

add_executable(pr_graphs ${SOURCE_FILES} ${HEADER_FILES})

target_link_libraries(pr_graphs ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
    set_target_properties(pr_graphs PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(pr_graphs PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()