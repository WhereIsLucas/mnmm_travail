//
// Created by hugo on 18/05/2020.
//

#ifndef MN_BORIS_BALLPRINTER_H
#define MN_BORIS_BALLPRINTER_H

#include <string>
#include "Ball.h"

class BallPrinter {
private:
    std::string path;
public:
    BallPrinter(const std::string &path);

    const std::string &getPath() const;

    void setPath(const std::string path);

    void print(Ball barrel, int frameNumber);

    void clearPrint(int frameNumber);
};


#endif //MN_BORIS_BALLPRINTER_H