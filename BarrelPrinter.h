//
// Created by hugo on 18/05/2020.
//

#ifndef TRAVAIL2_BARRELPRINTER_H
#define TRAVAIL2_BARRELPRINTER_H

#include <string>
#include "Barrel.h"

class BarrelPrinter {
private:
    std::string path;
public:
    BarrelPrinter(const std::string &path);

    const std::string &getPath() const;

    void setPath(const std::string path);

    void print(Barrel barrel, int frameNumber);

    void clearPrint(int frameNumber);
};


#endif //TRAVAIL2_BARRELPRINTER_H