
#ifndef MN_BORIS_PLANE_H
#define MN_BORIS_PLANE_H

#include <string>
#include "Vector2.h"

class Plane {

private:
    Vector2 position;
public:
    const Vector2 &getPosition() const;

    void setPosition(const Vector2 &position);

    const Vector2 &getNormal() const;

    void setNormal(const Vector2 &normal);

    void printPlanInfos(std::string fileName);

    Vector2 getPointFromX(double xValue);

private:
    Vector2 normal;

public:
    Plane();

    ~Plane();

    void initPlan(Vector2, Vector2);

    void initPlanFromCoordinates(Vector2 a, Vector2 b);
};

#endif //MN_BORIS_PLANE_H
