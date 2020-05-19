#include <math.h>
#include <iostream>
#include <fstream>
#include "Plan.h"

Plan::Plan() {

}

Plan::~Plan() {
    
}

void Plan::initPlan(Vector2 positionVector, Vector2 normalVector) {
    Plan::position = positionVector;
    Plan::normal = normalVector.normalize();
}

void Plan::initPlanFromCoordinates(Vector2 a, Vector2 b) {
    Plan::position = 0.5 * (a + b);
    Plan::normal = Vector2(-(b - a).getY(), (b - a).getX()).normalize();
}

const Vector2 &Plan::getPosition() const {
    return position;
}

void Plan::setPosition(const Vector2 &position) {
    Plan::position = position;
}

const Vector2 &Plan::getNormal() const {
    return normal;
}

void Plan::setNormal(const Vector2 &normal) {
    Plan::normal = normal;
}

void Plan::printPlanInfos(std::string fileName) {
    std::ofstream file;
    file.open(fileName.c_str());
    file.precision(10);
    Vector2 positionVector = Plan::getPosition();
    Vector2 normalVector = Plan::getNormal();
    double m = -normalVector.getX() / normalVector.getY();
    double p = positionVector.getY() - positionVector.getX() * (m);
    file << m << "," << p << std::endl;
    file.close();
}

Vector2 Plan::getPointFromX(double xValue) {
    Vector2 positionVector = Plan::getPosition();
    Vector2 normalVector = Plan::getNormal();
    double m = -normalVector.getX() / normalVector.getY();
    double p = positionVector.getY() - positionVector.getX() * (m);
    return Vector2(xValue, m * xValue + p);
}
