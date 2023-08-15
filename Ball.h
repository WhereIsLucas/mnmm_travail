#ifndef MN_BORIS_BALL_H
#define MN_BORIS_BALL_H

#include "Vector2.h"

class Ball {

private:
    Vector2 position;
    Vector2 velocity;
    Vector2 acceleration;
    Vector2 force = Vector2(0);

    double radius, mass, inertia;
public:
    void setRadius(double radius);

    void setMass(double mass);

    double getInertia() const;

    void setInertia(double inertia);

    void setTheta(double theta);

    double getW() const;

    void setW(double w);

    double getAlpha() const;

    void setAlpha(double alpha);

    double getM() const;

    void setM(double m);

private:
    double theta;
    double w;
    double alpha;
    double M;

public:
    Ball();

    ~Ball();

    const Vector2 &getPosition() const;

    void setPosition(const Vector2 &positionVector);

    const Vector2 &getVelocity() const;

    void setVelocity(const Vector2 &velocityVector);

    const Vector2 &getAcceleration() const;

    void setAcceleration(const Vector2 &accelerationVector);

    const Vector2 &getForce() const;

    void setForce(const Vector2 &forceVector);


    void updateVelocity(double);

    void updatePosition(double);

    void resetForce();

    void addForce(Vector2 addedForce);

    void addMomentum(double);

    void addGravityForce(Vector2 gravityDirection);


    double getRadius();

    double getX();

    double getY();

    double getVx();

    double getVy();


    double getOmega();


    double getMass();

    double getTheta();

    void initBarrel(double radius, double mass, Vector2 positionVector, Vector2 velocityVector);
};

#endif