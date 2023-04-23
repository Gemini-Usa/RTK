cmake_minimum_required(VERSION 3.24)
project(RealTimeKinematic)

set(CMAKE_CXX_STANDARD 20)

include_directories(.)

add_executable(RealTimeKinematic
        constants.h
        coordtrans.cpp
        errorprocess.cpp
        function.h
        main.cpp
        matrix.cpp
        novdecode.cpp
        satellite.cpp
        sockets.h
        sppspv.cpp
        timesync.cpp
        serial.h
        singlediff.cpp
        doublediff.cpp
        rtk.cpp
        lambda.cpp)
