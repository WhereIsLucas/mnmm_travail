//
// Created by hugo on 18/05/2020.
//

#include <fstream>
#include "BallPrinter.h"

const std::string &BallPrinter::getPath() const {
    return path;
}

void BallPrinter::setPath(const std::string path) {
    BallPrinter::path = path;
}

void BallPrinter::print(Ball barrel, int frameNumber) {
    std::string fileName = BallPrinter::getPath() + "barrel" + std::to_string(frameNumber) + ".txt";
    std::ofstream file;
    file.open(fileName.c_str(), std::ios::app);
    file.precision(10);
    file << barrel.getX() << "," << barrel.getY() << "," << barrel.getVx() << "," << barrel.getVy()
         << "," << barrel.getOmega() << "," << barrel.getRadius() << std::endl;
    file.close();
}

void BallPrinter::clearPrint(int frameNumber) {
    std::string fileName = BallPrinter::getPath() + "barrel" + std::to_string(frameNumber) + ".txt";
    remove(fileName.c_str());
}


BallPrinter::BallPrinter(const std::string &path) : path(path) {
    BallPrinter::path = path;
}
