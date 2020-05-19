//
// Created by hugo on 18/05/2020.
//

#include "Barrel.h"

Barrel::Barrel() = default;

Barrel::~Barrel() {

}

void Barrel::initBarrel(double radius, double mass, Vector2 positionVector, Vector2 velocityVector) {
    Barrel::setRadius(radius);
    Barrel::setMass(mass);
    inertia = 0.5 * mass * radius * radius;

    Barrel::setPosition(positionVector);
    Barrel::setVelocity(velocityVector);
    Barrel::setAcceleration(Vector2(0.));

    Barrel::setTheta(0.);
    Barrel::setW(0.);
    Barrel::setAlpha(0.);

    Barrel::setForce(Vector2(0.));

    M = 0.;
}

const Vector2 &Barrel::getPosition() const {
    return position;
}

void Barrel::setPosition(const Vector2 &positionVector) {
    Barrel::position = positionVector;
}

const Vector2 &Barrel::getVelocity() const {
    return velocity;
}

void Barrel::setVelocity(const Vector2 &velocityVector) {
    Barrel::velocity = velocityVector;
}

const Vector2 &Barrel::getAcceleration() const {
    return acceleration;
}

void Barrel::setAcceleration(const Vector2 &accelerationVector) {
    Barrel::acceleration = accelerationVector;
}

const Vector2 &Barrel::getForce() const {
    return force;
}

void Barrel::setForce(const Vector2 &forceVector) {
    Barrel::force = forceVector;
}

void Barrel::updateVelocity(double dt) {
    Barrel::setAcceleration((1 / mass) * Barrel::force);
    Barrel::setVelocity(Barrel::getVelocity() + dt * Barrel::getAcceleration());
    alpha = M / inertia;
    w += alpha * dt;
}

void Barrel::updatePosition(double dt) {
    Barrel::position = Barrel::getPosition() + dt * Barrel::getVelocity();
    theta += w * dt;
}

void Barrel::resetForce() {
    Barrel::setForce(Vector2(0.));
    M = 0.;
}

void Barrel::addForce(Vector2 addedForce) {
    Barrel::force = Barrel::force + addedForce;
}

void Barrel::addMomentum(double addM) {
    Barrel::M = Barrel::M + addM;
}

void Barrel::addGravityForce(Vector2 gravityDirection) {
    Barrel::addForce(mass * gravityDirection);

}

double Barrel::getRadius() {
    return radius;
}

double Barrel::getX() {
    return Barrel::position.getX();
}

double Barrel::getY() {
    return Barrel::position.getY();
}

double Barrel::getVy() {
    return Barrel::velocity.getY();
}

double Barrel::getOmega() {
    return w;
}

double Barrel::getMass() {
    return mass;
}

double Barrel::getTheta() {
    return theta;
}

double Barrel::getVx() {
    return Barrel::velocity.getX();
}

void Barrel::setRadius(double radius) {
    Barrel::radius = radius;
}

void Barrel::setMass(double mass) {
    Barrel::mass = mass;
}

double Barrel::getInertia() const {
    return inertia;
}

void Barrel::setInertia(double inertia) {
    Barrel::inertia = inertia;
}

void Barrel::setTheta(double theta) {
    Barrel::theta = theta;
}

double Barrel::getW() const {
    return w;
}

void Barrel::setW(double w) {
    Barrel::w = w;
}

double Barrel::getAlpha() const {
    return alpha;
}

void Barrel::setAlpha(double alpha) {
    Barrel::alpha = alpha;
}

double Barrel::getM() const {
    return M;
}

void Barrel::setM(double m) {
    M = m;
}
