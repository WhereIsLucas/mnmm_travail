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
    const std::string &getPath() const;
    void print(Grain grain, int frameNumber);
    void clearPrint(int frameNumber);
    void setPath(const std::string path);
};


#endif //TRAVAIL2_GRAINPRINTER_H
