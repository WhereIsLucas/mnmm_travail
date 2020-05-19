#include <fstream>
#include <iostream>
#include "GrainPrinter.h"

void GrainPrinter::print(Grain grain, int frameNumber) {
    std::string fileName = GrainPrinter::getPath()+"grain" + std::to_string(frameNumber) + ".txt";
    std::ofstream file;
    file.open(fileName.c_str(), std::ios::app);
    file.precision(10);
    file << grain.index() << "," << grain.getX() << "," << grain.getY() << "," << grain.getVx() << "," << grain.getVy() << "," << grain.getOmega() << "," << grain.getRadius() << std::endl;
    file.close();
}

void GrainPrinter::print(Grain *grains, int size, int frameNumber) {
    for (int i = 0; i < size; ++i) {
        print(grains[i], frameNumber);
    }
}

const std::string &GrainPrinter::getPath() const {
    return path;
}

void GrainPrinter::setPath(const std::string path) {
    GrainPrinter::path = path;
}

void GrainPrinter::clearPrint(int frameNumber) {
    std::string fileName = GrainPrinter::getPath()+"grain" + std::to_string(frameNumber) + ".txt";
    remove(fileName.c_str());
}

GrainPrinter::GrainPrinter(const std::string &path) : path(path) {
    GrainPrinter::path = path;
}

void GrainPrinter::clearPrints(int frameNumber) {
    for (int l = 0; l <= frameNumber; l++) {
        GrainPrinter::clearPrint(l);
    }
}
