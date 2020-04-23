#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "Grain.h"

using namespace std;

Grain::Grain() = default;

void Grain::initDisk(int i_index, double i_radius, double i_mass, Vector2 positionVector, Vector2 velocityVector)
{
    m_index = i_index;
    m_linkedDisk = -9;
    m_linkedCell = -9;
    m_radius = i_radius;
    m_mass = i_mass;
    m_inertia = 0.5*m_mass*m_radius*m_radius;

    Grain::setPosition(positionVector);
    Grain::setVelocity(velocityVector);
    Grain::setAcceleration(Vector2(0.));

    m_theta = 0.;
    m_w = 0.;
    m_alpha = 0.;

    Grain::setForce(Vector2(0.));

    m_M = 0.;
}

void Grain::resetForce()
{
    Grain::setForce(Vector2(0.));
    m_M = 0.;
}

void Grain::addForce(Vector2 addedForce)
{
    Grain::force = Grain::force + addedForce;
}

void Grain::addMomentum(double i_M)
{
    m_M += i_M;
}

void Grain::addGravityForce(Vector2 gravityDirection)
{
    Grain::addForce(m_mass * gravityDirection);
}

void Grain::updatePosition(double dt)
{
    Grain::position = Grain::getPosition() + dt*Grain::getVelocity();
    m_theta += m_w*dt;
}

void Grain::updateVelocity(double dt)
{
    Grain::setAcceleration((1/m_mass)* Grain::force);
    Grain::setVelocity(Grain::getVelocity() + dt*Grain::getAcceleration());
    m_alpha = m_M/m_inertia;
    m_w += m_alpha*dt;
}

void Grain::setLinkedDisk(int i_linkedDisk)
{
    m_linkedDisk = i_linkedDisk;
}

void Grain::setLinkedCell(int i_linkedCell)
{
    m_linkedCell = i_linkedCell;
}

int Grain::linkedDisk()
{
    return m_linkedDisk;
}

int Grain::linkedCell()
{
    return m_linkedCell;
}

double Grain::getRadius()
{
    return m_radius;
}

double Grain::getMass()
{
    return m_mass;
}

double Grain::getX()
{
    return Grain::position.getX();
}

double Grain::getY()
{
    return Grain::position.getY();
}

double Grain::getVx()
{
    return Grain::velocity.getX();
}

double Grain::getVy()
{
    return Grain::velocity.getY();
}

double Grain::w()
{
    return m_w;
}

int Grain::index() {
    return m_index;
}

double Grain::theta() {
    return m_theta;
}

void Grain::setVelocity(const Vector2 &velocityVector) {
    Grain::velocity = velocityVector;
}

const Vector2 &Grain::getAcceleration() const {
    return acceleration;
}

void Grain::setAcceleration(const Vector2 &accelerationVector) {
    Grain::acceleration = accelerationVector;
}

const Vector2 &Grain::getForce() const {
    return force;
}

void Grain::setForce(const Vector2 &forceVector) {
    Grain::force = forceVector;
}

const Vector2 &Grain::getPosition() const {
    return position;
}

void Grain::setPosition(const Vector2 &positionVector) {
    Grain::position = positionVector;
}

const Vector2 &Grain::getVelocity() const {
    return velocity;
}

Grain::~Grain() {

}

double getDistanceBetweenGrains(Grain grain1, Grain grain2) {
    getDistanceBetweenVectors(grain1.getPosition(),grain1.getPosition()) - (grain1.getRadius() + grain2.getRadius());
}

