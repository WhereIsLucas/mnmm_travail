cmake_minimum_required(VERSION 3.15)
project(mn_boris)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

include_directories(.)

add_executable(mn_boris
        Container.cpp
        Container.h
        Cell.cpp
        Cell.h
        Grain.cpp
        Grain.h
        main.cpp
        Plane.cpp
        Plane.h GrainPrinter.cpp GrainPrinter.h Vector2.cpp Vector2.h Domain.cpp Domain.h Ball.cpp Ball.h BallPrinter.cpp BallPrinter.h CollisionSettings.cpp CollisionSettings.h Collisions.cpp Collisions.h
        PlanePrinter.cpp)
