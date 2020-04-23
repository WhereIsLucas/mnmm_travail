//
// Created by lucas on 22/04/2020.
//


#include <iostream>
#include "Vector2.h"
#include "Grain.h"

int main(){
    Grain grains[2];
    grains[0].initDisk(1,.2,.2,Vector2(1,0),Vector2(0));
    grains[1].initDisk(1,.2,.2,Vector2(-1,0),Vector2(0));
    Vector2 normal = (grains[0].getPosition() - grains[1].getPosition());
    std::cout << normal.normalize().getNorm() << std::endl;
}