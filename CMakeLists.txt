cmake_minimum_required(VERSION 3.12)
project(Diplom)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)
include_directories(.idea)
include_directories(cmake-build-debug)

add_executable(Diplom
        .idea/codeStyles/Project.xml
        .idea/Diplom.iml
        .idea/misc.xml
        .idea/modules.xml
        .idea/workspace.xml
        cmake-build-debug/cmake_install.cmake
        cmake-build-debug/CMakeCache.txt
        cmake-build-debug/Diplom.cbp
        cmake-build-debug/Diplom.exe
        cmake-build-debug/Makefile
        alanin.psf
        alanin.txt
        alanin.txt.txt
        CMakeLists.txt
        dna.psf
        gradient_method.cpp
        gradient_method.h
        graph.cpp
        graph.h
        main.cpp)
