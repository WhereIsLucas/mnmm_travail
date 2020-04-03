#include <fstream>
#include "GrainPrinter.h"

void GrainPrinter::print(Grain grain, int frameNumber) {
    std::string fileName = "grain" + std::to_string(frameNumber) + ".txt";
    std::ofstream file;
    file.open(fileName.c_str(), std::ios::app);
    file.precision(10);
    file << grain.index() << "," << grain.x() << "," << grain.y() << "," << grain.vx() << "," << grain.vy() << "," << grain.theta() << "," << grain.radius() << std::endl;
    file.close();

}

const std::string &GrainPrinter::getPath() const {
    return path;
}

void GrainPrinter::setPath(const std::string &path) {
    GrainPrinter::path = path;
}
