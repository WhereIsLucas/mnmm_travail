
#ifndef MN_BORIS_PLANEPRINTER_H
#define MN_BORIS_PLANEPRINTER_H


#include <string>
#include "Grain.h"
#include "Plane.h"

class PlanePrinter {
private:
    std::string path;
public:
    PlanePrinter(const std::string &path);

    const std::string &getPath() const;

    void setPath(const std::string path);

    void print(Plane plane, int frameNumber);

    void clearPrint(int frameNumber) const;
};


#endif //MN_BORIS_PLANEPRINTER_H
