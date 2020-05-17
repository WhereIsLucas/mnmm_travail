#include "Grain.h"

using namespace std;

Grain::Grain() = default;

void Grain::initDisk(int i_index, double i_radius, double i_mass, Vector2 positionVector, Vector2 velocityVector)
{
    m_index = i_index;
    m_linkedDisk = -9;
    m_linkedCell = -9;
    radius = i_radius;
    mass = i_mass;
    inertia = 0.5 * mass * radius * radius;

    Grain::setPosition(positionVector);
    Grain::setVelocity(velocityVector);
    Grain::setAcceleration(Vector2(0.));
    Grain::setTheta(    0.);
    Grain::setOmega(    0.);
    Grain::setAlpha(    0.);

    Grain::resetForce();
}

void Grain::resetForce()
{
    Grain::setForce(Vector2(0.));
    Grain::setMomentum(0);
}

void Grain::addForce(Vector2 addedForce)
{
    Grain::force = Grain::force + addedForce;
}

void Grain::addMomentum(double addedMomentum)
{
    Grain::momentum = Grain::momentum + addedMomentum;
}

void Grain::addGravityForce(Vector2 gravityDirection)
{
    Grain::addForce(mass * gravityDirection);
}

void Grain::updatePosition(double dt)
{
    Grain::position = Grain::getPosition() + dt*Grain::getVelocity();
    theta += omega * dt;
}

void Grain::updateVelocity(double dt)
{
    Grain::setAcceleration((1 / mass) * Grain::force);
    Grain::setVelocity(Grain::getVelocity() + dt*Grain::getAcceleration());

    Grain::setAlpha((1./Grain::getInertia())*Grain::getMomentum());
    Grain::setOmega(Grain::getOmega() + dt*Grain::getAlpha());
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
    return radius;
}

double Grain::getMass()
{
    return mass;
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

double Grain::getOmega()
{
    return omega;
}

int Grain::index() {
    return m_index;
}

double Grain::getTheta() {
    return theta;
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

void Grain::setTheta(double theta) {
    Grain::theta = theta;
}

void Grain::setOmega(double omega) {
    Grain::omega = omega;
}

double Grain::getAlpha() const {
    return alpha;
}

void Grain::setAlpha(double alpha) {
    Grain::alpha = alpha;
}

void Grain::setRadius(double radius) {
    Grain::radius = radius;
}

void Grain::setMass(double mass) {
    Grain::mass = mass;
}

double Grain::getInertia() const {
    return inertia;
}

void Grain::setInertia(double inertia) {
    Grain::inertia = inertia;
}

double Grain::getMomentum() const {
    return momentum;
}

void Grain::setMomentum(double momentum) {
    Grain::momentum = momentum;
}

double getDistanceBetweenGrains(Grain grain1, Grain grain2) {
    return getDistanceBetweenVectors(grain1.getPosition(),grain2.getPosition()) - (grain1.getRadius() + grain2.getRadius());
}

