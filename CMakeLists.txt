cmake_minimum_required(VERSION 3.7)
project(courseProject_sdz3)

set(CMAKE_CXX_STANDARD 14)

SET(CMAKE_CXX_FLAGS_DEBUG )

set(SOURCE_FILES
        AlgorithmTest.cpp
        AlgorithmTest.h
        main.cpp
        TimeMeasurement.cpp
        TimeMeasurement.h
        TravellingSalesmanProblem.cpp
        TravellingSalesmanProblem.h
        TSPLIB_Parser.cpp TSPLIB_Parser.h)

add_executable(courseProject_sdz3 ${SOURCE_FILES})