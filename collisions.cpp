//
// Created by lucas on 16/05/2020.
//

#include "collisions.h"
#include "Grain.h"
#include "Container.h"
#include "CollisionSettings.h"


//void computeCollision(Grain *pGrain1, Grain *pGrain2) {
//    Vector2 normalVector = (pGrain2->getPosition() - pGrain1->getPosition()).normalize();
//    double delta = getDistanceBetweenGrains(*pGrain1, *pGrain2);
//    if (delta < 0) {
//        double vy = pGrain1->getVy() - pGrain2->getVy();
////                    - pGrain1->getRadius() * pGrain1->getOmega() * normalVector.getX() -
////                    pGrain2->getRadius() * pGrain2->getOmega() * normalVector.getX();
//        double vx = pGrain1->getVx() - pGrain2->getVx();
////                    + pGrain1->getRadius() * pGrain1->getOmega() * normalVector.getY() +
////                    pGrain2->getRadius() * pGrain2->getOmega() * normalVector.getY();
//        Vector2 velocityAtContactPoint(vx, vy);
//        Vector2 normalVelocity(velocityAtContactPoint.getX() * normalVector.getX(),
//                               velocityAtContactPoint.getY() * normalVector.getY());
//        Vector2 tangentVelocity = velocityAtContactPoint - normalVelocity;
//        Vector2 tangentVector;
//        if (tangentVelocity.getNorm() != 0.) {
//            tangentVector = tangentVelocity.normalize();
//        } else {
//            tangentVector.setComponents(0, 0);
//        }
//
//        //contact forces and torque
//        double effectiveMass = (pGrain1->getMass() * pGrain2->getMass()) /
//                               (pGrain1->getMass() + pGrain2->getMass());
//        double eta = -2. * log(e) * sqrt(effectiveMass * kn / (log(e) * log(e) + M_PI * M_PI));
//        double normalForceNorm = -1.*(kn * delta + eta * normalVelocity.getNorm());
//        double tangentForceNorm = -kt * tangentVelocity.getNorm();
//        Vector2 tangentForce(-kt * tangentVelocity);
//
//        // check if normal force is repulsive
//        if (normalForceNorm > 0) {
//            Vector2 normalForce(tangentForceNorm * normalVector);
//            pGrain1->addForce(normalForce);
//            pGrain2->addForce(-1 * normalForce);
//        } else {
//            double fn = 0.;
//        }
//
//
//        if (tangentForce.getNorm() > mu * normalForceNorm) {
////            ftx = -mu * normalForce.getNorm() * normalizedTangentVelocity();
////            fty = -mu * fn * ty;
//            pGrain1->addForce(Vector2(0));
//            pGrain1->addForce(-1. * Vector2(0));
//        } else {
//            pGrain1->addForce(Vector2(0));
//            pGrain2->addForce(Vector2(0));
//        }
//
////torque
////        double M = -pGrain1->getRadius() * nx * fty + pGrain1->getRadius() * ny * ftx;
////        pGrain1->addMomentum(M);
////        double M2 = -pGrain2->getRadius() * nx * fty + pGrain2->getRadius() * ny * ftx;
////        pGrain2->addMomentum(M2);
//    }
//}


void computeCollisionWithContainer(Grain *pGrain1, Container *container, CollisionSettings *collisionSettings) {
    double delta = container->getRadius() - pGrain1->getRadius() -
                   getDistanceBetweenVectors(pGrain1->getPosition(), container->getCenter());
    if (delta < 0) {
        Vector2 normalVector = (container->getCenter() - pGrain1->getPosition()).normalize();

        double vy = pGrain1->getVy();
        double vx = pGrain1->getVx();

        Vector2 velocityAtContactPoint(vx, vy);
        Vector2 normalVelocity(velocityAtContactPoint.getX() * normalVector.getX(),
                               velocityAtContactPoint.getY() * normalVector.getY());
        Vector2 tangentVelocity = velocityAtContactPoint - normalVelocity;
        Vector2 tangentVector(0);
        if (tangentVelocity.getNorm() != 0.) {
            tangentVector = tangentVelocity.normalize();
        }

        //contact forces and torque
        double effectiveMass = pGrain1->getMass();
        double eta = collisionSettings->getEta(effectiveMass);
        double normalForceNorm = -1. * (collisionSettings->getKn() * delta) + (eta * normalVelocity.getNorm());
        double tangentForceNorm = -collisionSettings->getKt() * tangentVelocity.getNorm();
        Vector2 tangentForce(tangentForceNorm * tangentVector.getX(), (tangentForceNorm * tangentVector.getY()));

        if (normalForceNorm > 0) {
            pGrain1->addForce(normalForceNorm * normalVector);
//            std::cout << (normalForceNorm * normalVector).getX() << "" << (normalForceNorm * normalVector).getY() << " N" << std::endl;
        } else {
            normalForceNorm = 0.;
        }

        if (tangentForce.getNorm() > collisionSettings->getMu() * normalForceNorm) {
            pGrain1->addForce(-1 * collisionSettings->getMu() * normalForceNorm * tangentVector);
        } else {
            pGrain1->addForce(tangentForce);
        }
        Vector2 forceVector = pGrain1->getForce();
//        std::cout << forceVector.getX() << " " << forceVector.getY() << std::endl;

        //torque
        double M = (-1 * pGrain1->getRadius() * normalVector.getX() * tangentForce.getY())
                   + (pGrain1->getRadius() * normalVector.getY() * tangentForce.getX());
//        std::cout << M << std::endl;

        pGrain1->addMomentum(M);
    }


}
