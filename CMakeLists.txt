cmake_minimum_required(VERSION 3.15)
project(travail2)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

include_directories(.)

add_executable(travail2
        Container.cpp
        Container.h
        Cell.cpp
        Cell.h
        Grain.cpp
        Grain.h
        main.cpp
        Plan.cpp
        Plan.h GrainPrinter.cpp GrainPrinter.h Vector2.cpp Vector2.h Domain.cpp Domain.h Barrel.cpp Barrel.h BarrelPrinter.cpp BarrelPrinter.h)
