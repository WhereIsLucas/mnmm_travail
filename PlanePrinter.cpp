#include <fstream>
#include <iostream>
#include "PlanePrinter.h"
#include "Plane.h"

void PlanePrinter::print(Plane plane, int frameNumber) {
    std::string fileName = PlanePrinter::getPath()+"plane" + std::to_string(frameNumber) + ".txt";
    std::ofstream file;
    file.open(fileName.c_str(), std::ios::app);
    file.precision(10);
    Vector2 positionVector = plane.getPosition();
    Vector2 normalVector = plane.getNormal();
    double m = -normalVector.getX() / normalVector.getY();
    double p = positionVector.getY() - positionVector.getX() * (m);
    file << m << "," << p << std::endl;
    file.close();
}

const std::string &PlanePrinter::getPath() const {
    return path;
}

void PlanePrinter::setPath(const std::string path) {
    PlanePrinter::path = path;
}

void PlanePrinter::clearPrint(int frameNumber) const {
    std::string fileName = PlanePrinter::getPath() + "plane" + std::to_string(frameNumber) + ".txt";
    remove(fileName.c_str());
}

PlanePrinter::PlanePrinter(const std::string &path) : path(path) {
    PlanePrinter::path = path;
}
