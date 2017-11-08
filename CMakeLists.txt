cmake_minimum_required(VERSION 3.7)
project(courseProject_pea1)

SET(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_FLAGS_RELEASE "-fno-rtti -O2")
SET(CMAKE_CXX_FLAGS_DEBUG "-fno-rtti -O2")

set(SOURCE_FILES
        AlgorithmTest.cpp
        AlgorithmTest.h
        main.cpp
        TimeMeasurement.cpp
        TimeMeasurement.h
        TravellingSalesmanProblem.cpp
        TravellingSalesmanProblem.h
        TSPLIB_Parser.cpp TSPLIB_Parser.h)

add_executable(courseProject_pea1 ${SOURCE_FILES})