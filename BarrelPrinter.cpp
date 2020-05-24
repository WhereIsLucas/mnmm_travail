//
// Created by hugo on 18/05/2020.
//

#include <fstream>
#include "BarrelPrinter.h"

const std::string &BarrelPrinter::getPath() const {
    return path;
}

void BarrelPrinter::setPath(const std::string path) {
    BarrelPrinter::path = path;
}

void BarrelPrinter::print(Barrel barrel, int frameNumber) {
    std::string fileName = BarrelPrinter::getPath() + "barrel" + std::to_string(frameNumber) + ".txt";
    std::ofstream file;
    file.open(fileName.c_str(), std::ios::app);
    file.precision(10);
    file << barrel.getX() << "," << barrel.getY() << "," << barrel.getVx() << "," << barrel.getVy()
         << "," << barrel.getTheta() << "," << barrel.getRadius() << std::endl;
    file.close();
}

void BarrelPrinter::clearPrint(int frameNumber) {
    std::string fileName = BarrelPrinter::getPath() + "barrel" + std::to_string(frameNumber) + ".txt";
    remove(fileName.c_str());
}


BarrelPrinter::BarrelPrinter(const std::string &path) : path(path) {
    BarrelPrinter::path = path;
}
