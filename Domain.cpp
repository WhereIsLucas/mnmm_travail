//
// Created by hugo on 14/05/2020.
//

#include <string>
#include <iostream>
#include <fstream>
#include "Domain.h"

double Domain::getX() const {
    return x;
}

void Domain::setX(double x) {
    Domain::x = x;
}

double Domain::getY() const {
    return y;
}

void Domain::setY(double y) {
    Domain::y = y;
}

Domain::Domain(double x, double y) : x(x), y(y) {
    Domain::x = x;
    Domain::y = y;
}

void Domain::printDomainInfos(std::string fileName) {
    std::ofstream file;
    file.open(fileName.c_str());
    file.precision(10);
    file << Domain::getX() << "," << Domain::getY() << std::endl;
    file.close();
}
