#include <math.h>
#include <iostream>
#include <fstream>
#include "Plane.h"

Plane::Plane() {

}

Plane::~Plane() {
    
}

void Plane::initPlan(Vector2 positionVector, Vector2 normalVector) {
    Plane::position = positionVector;
    Plane::normal = normalVector.normalize();
}

void Plane::initPlanFromCoordinates(Vector2 a, Vector2 b) {
    Plane::position = 0.5 * (a + b);
    Plane::normal = Vector2(-(b - a).getY(), (b - a).getX()).normalize();
}

const Vector2 &Plane::getPosition() const {
    return position;
}

void Plane::setPosition(const Vector2 &position) {
    Plane::position = position;
}

const Vector2 &Plane::getNormal() const {
    return normal;
}

void Plane::setNormal(const Vector2 &normal) {
    Plane::normal = normal;
}

void Plane::printPlanInfos(std::string fileName) {
    std::ofstream file;
    file.open(fileName.c_str());
    file.precision(10);
    Vector2 positionVector = Plane::getPosition();
    Vector2 normalVector = Plane::getNormal();
    double m = -normalVector.getX() / normalVector.getY();
    double p = positionVector.getY() - positionVector.getX() * (m);
    file << m << "," << p << std::endl;
    file.close();
}

Vector2 Plane::getPointFromX(double xValue) {
    Vector2 positionVector = Plane::getPosition();
    Vector2 normalVector = Plane::getNormal();
    double m = -normalVector.getX() / normalVector.getY();
    double p = positionVector.getY() - positionVector.getX() * (m);
    return Vector2(xValue, m * xValue + p);
}
