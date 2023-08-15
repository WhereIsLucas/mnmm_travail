cmake_minimum_required(VERSION 3.15)
project(mnmm)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

include_directories(.)

add_executable(mnmm
        Container.cpp
        Container.h
        Cell.cpp
        Cell.h
        Grain.cpp
        Grain.h
        main.cpp
        Plan.cpp
        Plan.h GrainPrinter.cpp GrainPrinter.h Vector2.cpp Vector2.h Domain.cpp Domain.h collisions.cpp collisions.h CollisionSettings.cpp CollisionSettings.h)
