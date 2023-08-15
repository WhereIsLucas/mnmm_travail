//
// Created by lucas on 28/03/2020.
//

#ifndef MN_BORIS_CONTAINER_H
#define MN_BORIS_CONTAINER_H


#include "Vector2.h"

class Container {
private:
    double radius{};
    Vector2 center;
public:
    const Vector2 &getCenter() const;

    void setCenter(const Vector2 &centerVector);

    Container(double radiusValue, Vector2 centerVector);

    double getRadius() const;

    void setRadius(double radiusValue);
};

#endif //MN_BORIS_CONTAINER_H
