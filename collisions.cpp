

#include <iostream>
#include "collisions.h"
#include "Grain.h"
#include "Container.h"
#include "CollisionSettings.h"


void computeCollisionWithGrain(Grain *pGrain1, Grain *pGrain2, CollisionSettings *collisionSettings) {
    double delta = getDistanceBetweenGrains(*pGrain1, *pGrain2);
    if (delta < 0) {
        Vector2 normalVector = (pGrain1->getPosition() - pGrain2->getPosition()).normalize();
        double vx = pGrain1->getVx() - pGrain2->getVx()
                    + pGrain1->getRadius() * pGrain1->getOmega() * normalVector.getY()
                    + pGrain2->getRadius() * pGrain2->getOmega() * normalVector.getY();
        double vy = pGrain1->getVy() - pGrain2->getVy()
                    - pGrain1->getRadius() * pGrain1->getOmega() * normalVector.getX()
                    - pGrain2->getRadius() * pGrain2->getOmega() * normalVector.getX();

        Vector2 velocityAtContactPoint(vx, vy);

        Vector2 normalVelocity = projectOntoVector(velocityAtContactPoint, normalVector);
        Vector2 tangentVelocity = velocityAtContactPoint - normalVelocity;
        Vector2 tangentVector(0.);
        if (tangentVelocity.getNorm() != 0.) {
            tangentVector = tangentVelocity.normalize();
        }

        //contact forces and torque
        double effectiveMass = (pGrain1->getMass() * pGrain2->getMass()) /
                               (pGrain1->getMass() + pGrain2->getMass());
        double eta = collisionSettings->getEta(effectiveMass);
        double normalForceNorm = -1.*(collisionSettings->getKn() * delta + eta * normalVelocity.getNorm());
        double tangentForceNorm = -1.* collisionSettings->getKt() * tangentVelocity.getNorm();
        Vector2 tangentForce(tangentForceNorm * tangentVector);
        // check if normal force is repulsive
        if (normalForceNorm > 0.) {
            pGrain1->addForce(normalForceNorm * normalVector);
            pGrain2->addForce(-1. * normalForceNorm * normalVector);
        } else {
            normalForceNorm = 0.;
        }

        if (tangentForce.getNorm() > collisionSettings->getMu() * normalForceNorm) {
            tangentForce = -1. * collisionSettings->getMu() * normalForceNorm * tangentVector;
        }
        pGrain1->addForce(tangentForce);
        pGrain2->addForce(-1.* tangentForce);

//torque
        double M = pGrain1->getRadius() *
                   (-1.*normalVector.getX() * tangentForce.getY()
                    + (normalVector.getY() * tangentForce.getX()));
        pGrain1->addMomentum(M);
        double M2 = pGrain2->getRadius() *
                   (-1.*normalVector.getX() * tangentForce.getY()
                    + (normalVector.getY() * tangentForce.getX()));
        pGrain2->addMomentum(M2);
        return;
    }
}


void computeCollisionWithContainer(Grain *pGrain1, Container *container, CollisionSettings *collisionSettings) {
    double delta = container->getRadius() - pGrain1->getRadius() -
                   getDistanceBetweenVectors(pGrain1->getPosition(), container->getCenter());
    if (delta < 0) {
        Vector2 normalVector = (container->getCenter() - pGrain1->getPosition()).normalize();
        double radialSpeed = container->getRadialSpeed();
        Vector2 containerSpeedVector(-radialSpeed*normalVector.getY(),radialSpeed*normalVector.getX());
        double vx = pGrain1->getVx() - containerSpeedVector.getX() + pGrain1->getRadius() * pGrain1->getOmega() * normalVector.getY();
        double vy = pGrain1->getVy() - containerSpeedVector.getY() - pGrain1->getRadius() * pGrain1->getOmega() * normalVector.getX();

        Vector2 velocityAtContactPoint(vx, vy);
        Vector2 normalVelocity = projectOntoVector(velocityAtContactPoint,normalVector);
        Vector2 tangentVelocity = velocityAtContactPoint - normalVelocity;
        Vector2 tangentVector(0);
        if (tangentVelocity.getNorm() != 0.) {
            tangentVector = tangentVelocity.normalize();
        }
        //contact forces and torque
        double effectiveMass = pGrain1->getMass();
        double eta = collisionSettings->getEta(effectiveMass);
        double normalForceNorm = -1. * ((collisionSettings->getKn() * delta) + (eta * normalVelocity.getNorm()));
        double tangentForceNorm = -1.* collisionSettings->getKt() * tangentVelocity.getNorm();
        Vector2 tangentForce(tangentForceNorm * tangentVector);
        if (normalForceNorm > 0.) {
            pGrain1->addForce(normalForceNorm * normalVector);
        } else {
            normalForceNorm = 0.;
        }

        if (tangentForce.getNorm() > collisionSettings->getMu() * normalForceNorm) {
            tangentForce = -1. * collisionSettings->getMu() * normalForceNorm * tangentVector;
        }
        pGrain1->addForce(tangentForce);
        //torque
        double M = pGrain1->getRadius() *
                     (-1.*normalVector.getX() * tangentForce.getY()
                      + (normalVector.getY() * tangentForce.getX()));
        pGrain1->addMomentum(M);
        return;
    }


}
