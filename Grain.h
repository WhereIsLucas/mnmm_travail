#ifndef MNMM_GRAIN_H
#define MNMM_GRAIN_H

#include "Vector2.h"

class Grain {

private:
    int m_index;
    Vector2 position;
    Vector2 velocity;
    Vector2 acceleration;
    Vector2 force = Vector2(0);
public:
    const Vector2 &getPosition() const;
    void setPosition(const Vector2 &positionVector);
    const Vector2 &getVelocity() const;
    void setVelocity(const Vector2 &velocity);
    const Vector2 &getAcceleration() const;
    void setAcceleration(const Vector2 &accelerationVector);
    const Vector2 &getForce() const;
    void setForce(const Vector2 &forceVector);
private:
    int m_linkedCell;
    int m_linkedDisk;

    double radius, mass, inertia;
public:
    void setRadius(double radius);

    void setMass(double mass);

    double getInertia() const;

    void setInertia(double inertia);

public:
    void setTheta(double theta);

    void setOmega(double omega);

    double getAlpha() const;

    void setAlpha(double alpha);

private:
    double theta;
    double omega;
    double alpha;
    double momentum;
public:
    double getMomentum() const;

    void setMomentum(double momentum);

public:
    Grain();

    ~Grain();


    void updateVelocity(double);

    void updatePosition(double);

    void resetForce();

    void addForce(Vector2 addedForce);

    void addMomentum(double);

    void addGravityForce(Vector2 gravityDirection);

    void setLinkedDisk(int);

    void setLinkedCell(int);

    void print(int);

    int linkedCell();

    int linkedDisk();

    double getRadius();

    double getX();

    double getY();

    double getVx();

    double getVy();


    double getOmega();

    double getMass();

    int index();

    double getTheta();

    void initDisk(int i_index, double i_radius, double i_mass, Vector2 positionVector, Vector2 velocityVector);

};

double getDistanceBetweenGrains(Grain grain1, Grain grain2);

#endif