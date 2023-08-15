//
// Created by hugo on 18/05/2020.
//

#include "Ball.h"

Ball::Ball() = default;

Ball::~Ball() {

}

void Ball::initBarrel(double radius, double mass, Vector2 positionVector, Vector2 velocityVector) {
    Ball::setRadius(radius);
    Ball::setMass(mass);
    inertia = 0.5 * mass * radius * radius;

    Ball::setPosition(positionVector);
    Ball::setVelocity(velocityVector);
    Ball::setAcceleration(Vector2(0.));

    Ball::setTheta(0.);
    Ball::setW(0.);
    Ball::setAlpha(0.);

    Ball::setForce(Vector2(0.));

    M = 0.;
}

const Vector2 &Ball::getPosition() const {
    return position;
}

void Ball::setPosition(const Vector2 &positionVector) {
    Ball::position = positionVector;
}

const Vector2 &Ball::getVelocity() const {
    return velocity;
}

void Ball::setVelocity(const Vector2 &velocityVector) {
    Ball::velocity = velocityVector;
}

const Vector2 &Ball::getAcceleration() const {
    return acceleration;
}

void Ball::setAcceleration(const Vector2 &accelerationVector) {
    Ball::acceleration = accelerationVector;
}

const Vector2 &Ball::getForce() const {
    return force;
}

void Ball::setForce(const Vector2 &forceVector) {
    Ball::force = forceVector;
}

void Ball::updateVelocity(double dt) {
    Ball::setAcceleration((1 / mass) * Ball::force);
    Ball::setVelocity(Ball::getVelocity() + dt * Ball::getAcceleration());
    alpha = M / inertia;
    w += alpha * dt;
}

void Ball::updatePosition(double dt) {
    Ball::position = Ball::getPosition() + dt * Ball::getVelocity();
    theta += w * dt;
}

void Ball::resetForce() {
    Ball::setForce(Vector2(0.));
    M = 0.;
}

void Ball::addForce(Vector2 addedForce) {
    Ball::force = Ball::force + addedForce;
}

void Ball::addMomentum(double addM) {
    Ball::M = Ball::M + addM;
}

void Ball::addGravityForce(Vector2 gravityDirection) {
    Ball::addForce(mass * gravityDirection);

}

double Ball::getRadius() {
    return radius;
}

double Ball::getX() {
    return Ball::position.getX();
}

double Ball::getY() {
    return Ball::position.getY();
}

double Ball::getVy() {
    return Ball::velocity.getY();
}

double Ball::getOmega() {
    return w;
}

double Ball::getMass() {
    return mass;
}

double Ball::getTheta() {
    return theta;
}

double Ball::getVx() {
    return Ball::velocity.getX();
}

void Ball::setRadius(double radius) {
    Ball::radius = radius;
}

void Ball::setMass(double mass) {
    Ball::mass = mass;
}

double Ball::getInertia() const {
    return inertia;
}

void Ball::setInertia(double inertia) {
    Ball::inertia = inertia;
}

void Ball::setTheta(double theta) {
    Ball::theta = theta;
}

double Ball::getW() const {
    return w;
}

void Ball::setW(double w) {
    Ball::w = w;
}

double Ball::getAlpha() const {
    return alpha;
}

void Ball::setAlpha(double alpha) {
    Ball::alpha = alpha;
}

double Ball::getM() const {
    return M;
}

void Ball::setM(double m) {
    M = m;
}
