cmake_minimum_required(VERSION 3.7...3.29)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(PhononXFEM  VERSION 1.0
                    DESCRIPTION "Project to evaluate dispersion curves of phononic crystals using XFEM"
                    LANGUAGES CXX)

set(SOURCE_FILES1
    src/main.cpp
    src/Quad.cpp
    src/Material.cpp
    src/Quad2.cpp
    src/Quad3.cpp
    src/Quad4.cpp
    src/Quad5.cpp
)

set(SOURCE_FILES2
    src/plot.cpp
)

set(INCLUDE_FILES
    include/Quad.h
    include/Material.h
    include/matplotlibcpp.h
    include/Quad2.h
    include/Quad3.h
    include/Quad4.h
    include/Quad5.h
)

add_executable(PhononXFEM ${SOURCE_FILES1} ${INCLUDE_FILES} )

add_executable(PlotDispersion ${SOURCE_FILES2} ${INCLUDE_FILES})

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
target_link_libraries(PhononXFEM ${PYTHON_LIBRARIES})
target_link_libraries(PlotDispersion ${PYTHON_LIBRARIES})