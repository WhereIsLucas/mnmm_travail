//
// Created by lucas on 02/04/2020.
//

#ifndef TRAVAIL2_GRAINPRINTER_H
#define TRAVAIL2_GRAINPRINTER_H


#include <string>
#include "Grain.h"

class GrainPrinter {
private:
    std::string path;
public:
    GrainPrinter(const std::string &path);

    const std::string &getPath() const;
    void setPath(const std::string path);

    void print(Grain grain, int frameNumber);
    void print(Grain *grains, int size, int frameNumber);
    void clearPrint(int frameNumber);
    void clearPrints(int frameNumber);

};


#endif //TRAVAIL2_GRAINPRINTER_H
