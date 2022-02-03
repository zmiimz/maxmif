cmake_minimum_required(VERSION 3.14)
project(MaxMIF_zmi)

set(CMAKE_CXX_STANDARD 17)

include_directories(.)
link_libraries(stdc++fs)
add_executable(MaxMIF_zmi
        MaxMIF.cpp
        MaxMIF.h)
