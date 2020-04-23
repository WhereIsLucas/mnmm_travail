//
// Created by lucas on 28/03/2020.
//

#include "Container.h"
#include "Vector2.h"

void Container::initContainer(double radiusValue, Vector2 centerPosition) {
    Container::radius = radiusValue;
    Container::center = centerPosition;
}

double Container::getRadius() const {
    return radius;
}

void Container::setRadius(double radiusValue) {
    Container::radius = radiusValue;
}

const Vector2 &Container::getCenter() const {
    return center;
}

void Container::setCenter(const Vector2 &centerVector) {
    Container::center = centerVector;
}

Container::Container(double radiusValue, Vector2 centerVector) {
    Container::setRadius(radiusValue);
    Container::setCenter(centerVector);
}

