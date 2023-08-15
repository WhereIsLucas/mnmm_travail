
#ifndef MN_BORIS_GRAINPRINTER_H
#define MN_BORIS_GRAINPRINTER_H


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

    void clearPrint(int frameNumber);
};


#endif //MN_BORIS_GRAINPRINTER_H
