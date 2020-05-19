
#include "Container.h"
#include "Vector2.h"

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

Container::Container(double radiusValue, Vector2 centerVector, double omegaValue) {
    Container::setRadius(radiusValue);
    Container::setOmega(omegaValue);
    Container::setCenter(centerVector);
}

double Container::getOmega() const {
    return omega;
}

void Container::setOmega(double omega) {
    Container::omega = omega;
}

double Container::getRadialSpeed() {
    return Container::omega * Container::getRadius();
}

