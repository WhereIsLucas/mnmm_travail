//
// Created by lucas on 02/04/2020.
//

#ifndef MNMM_GRAINPRINTER_H
#define MNMM_GRAINPRINTER_H


#include <string>
#include "Grain.h"

class GrainPrinter {
private:
    std::string path;
public:
    GrainPrinter(const std::string &path);

    const std::string &getPath() const;
    void setPath(const std::string path);

    void print(Grain grain, int frameNumber) const;
    void clearPrint(int frameNumber) const;
};


#endif //MNMM_GRAINPRINTER_H
