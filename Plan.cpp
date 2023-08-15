#include "Plan.h"

Plan::Plan()
{

}

Plan::~Plan()
{
    
}

void Plan::initPlan(Vector2 positionVector, Vector2 normalVector) {
    Plan::position = positionVector;
    Plan::normal = normalVector;
}

void Plan::initPlanFromCoordinates(Vector2 a, Vector2 b) {
    Plan::position = 0.5 * (a + b);
    Plan::normal = Vector2(-(b - a).getY(), (b - a).getX());
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