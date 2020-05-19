//
// Created by lucas on 16/05/2020.
//

#include <iostream>
#include <tgmath.h>
#include "Collisions.h"
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

        Vector2 normalVelocity = projectOntoVector(velocityAtContactPoint, normalVector).getNorm() * normalVector;
        Vector2 tangentVelocity = velocityAtContactPoint - normalVelocity;
        Vector2 tangentVector(0);
        if (tangentVelocity.getNorm() != 0.) {
            tangentVector = tangentVelocity.normalize();
        }

        //contact forces and torque
        double effectiveMass = (pGrain1->getMass() * pGrain2->getMass()) /
                               (pGrain1->getMass() + pGrain2->getMass());
        double eta = collisionSettings->getEta(effectiveMass);
        double normalForceNorm = -1. * (collisionSettings->getKn() * delta + eta * normalVelocity.getNorm());
        double tangentForceNorm = -1. * collisionSettings->getKt() * tangentVelocity.getNorm();
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
        pGrain1->addForce(-1. * tangentForce);

//torque
        double M = pGrain1->getRadius() *
                   (-1. * normalVector.getX() * tangentForce.getY()
                    + (normalVector.getY() * tangentForce.getX()));
        pGrain1->addMomentum(M);
        double M2 = pGrain2->getRadius() *
                    (-1. * normalVector.getX() * tangentForce.getY()
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
        double vx = pGrain1->getVx() + pGrain1->getRadius() * pGrain1->getOmega() * normalVector.getY();
        double vy = pGrain1->getVy() - pGrain1->getRadius() * pGrain1->getOmega() * normalVector.getX();

        Vector2 velocityAtContactPoint(vx, vy);
        Vector2 normalVelocity = projectOntoVector(velocityAtContactPoint, normalVector).getNorm() * normalVector;
        Vector2 tangentVelocity = velocityAtContactPoint - normalVelocity;
        Vector2 tangentVector(0);
        if (tangentVelocity.getNorm() != 0.) {
            tangentVector = tangentVelocity.normalize();
        }
        //contact forces and torque
        double effectiveMass = pGrain1->getMass();
        double eta = collisionSettings->getEta(effectiveMass);
        double normalForceNorm = -1. * ((collisionSettings->getKn() * delta) + (eta * normalVelocity.getNorm()));
        double tangentForceNorm = -1. * collisionSettings->getKt() * tangentVelocity.getNorm();
        Vector2 tangentForce(tangentForceNorm * tangentVector);
        if (normalForceNorm > 0) {
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
                   (-1. * normalVector.getX() * tangentForce.getY()
                    + (normalVector.getY() * tangentForce.getX()));
        pGrain1->addMomentum(M);
        return;
    }


}

void computeCollisionBetweenBarrelAndPlan(Barrel *pBarrel, Plan *plan, CollisionSettings *collisionSettings) {
    Vector2 normalVecteur = plan->getNormal();
    Vector2 vecteur = plan->getPosition() - pBarrel->getPosition();
    double delta = vecteur.getX() * normalVecteur.getX() + vecteur.getY() * normalVecteur.getY();
    normalVecteur.display();
//    std::cout << delta << std::endl;

    if (delta < 0) {

        Vector2 normalVector = normalVecteur;

        double vy = pBarrel->getVy();
        double vx = pBarrel->getVx();

        Vector2 velocityAtContactPoint(vx, vy);
        Vector2 normalVelocity(velocityAtContactPoint.getX() * normalVector.getX(),
                               velocityAtContactPoint.getY() * normalVector.getY());
        Vector2 tangentVelocity = velocityAtContactPoint - normalVelocity;
        Vector2 tangentVector(0);
        if (tangentVelocity.getNorm() != 0.) {
            tangentVector = tangentVelocity.normalize();
        }

        //contact forces and torque
        double effectiveMass = pBarrel->getMass();
        double eta = -2. * log(collisionSettings->getE()) * sqrt(effectiveMass * collisionSettings->getKn() /
                                                                 (pow(log(collisionSettings->getE()), 2) +
                                                                  pow(M_PI, 2)));
        double normalForceNorm = -1. * (collisionSettings->getKn() * delta) + (eta * normalVelocity.getNorm());
        double tangentForceNorm = -collisionSettings->getKt() * tangentVelocity.getNorm();
        Vector2 tangentForce(tangentForceNorm * tangentVector.getX(), (tangentForceNorm * tangentVector.getY()));

        if (normalForceNorm > 0) {
            pBarrel->addForce(normalForceNorm * normalVector);
//         std::cout << (normalForceNorm * normalVector).getX() << "" << (normalForceNorm * normalVector).getY() << " N" << std::endl;
        } else {
            normalForceNorm = 0.;
        }

        if (tangentForce.getNorm() > collisionSettings->getMu() * normalForceNorm) {
            pBarrel->addForce(-1 * collisionSettings->getMu() * normalForceNorm * tangentVector);
            std::cout << "adding force" << std::endl;
        } else {
            pBarrel->addForce(tangentForce);
        }
        Vector2 forceVector = pBarrel->getForce();
//        std::cout << forceVector.getX() << " " << forceVector.getY() << std::endl;


//CE QUE J AI AJOUTE
        // check if normal force is repulsive
        if (normalForceNorm > 0) {
            Vector2 normalForce(tangentForceNorm * normalVector);
            pBarrel->addForce(normalForce);
        } else {
            double fn = 0.;
        }


        if (tangentForce.getNorm() > collisionSettings->getMu() * normalForceNorm) {
//            ftx = -mu * normalForce.getNorm() * normalizedTangentVelocity();
//            fty = -mu * fn * ty;
            pBarrel->addForce(Vector2(0));
            pBarrel->addForce(-1. * Vector2(0));
        } else {
            pBarrel->addForce(Vector2(0));
        }
//JUSQUICI

        //torque
        double M = (-1 * pBarrel->getRadius() * normalVector.getX() * tangentForce.getY())
                   + (pBarrel->getRadius() * normalVector.getY() * tangentForce.getX());
//        std::cout << M << std::endl;

        pBarrel->addMomentum(M);
    }


}