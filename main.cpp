#include <iostream>
#include <random>
#include <map>
#include <omp.h>
#include "Grain.h"
#include "Barrel.h"
#include "BarrelPrinter.h"
#include "Plan.h"
#include "Cell.h"
#include "GrainPrinter.h"
#include "Container.h"
#include "Vector2.h"
#include "Domain.h"
#include <chrono>

void computeCollision(Grain *pGrain1, Grain *pGrain2);

void computeCollisionWithContainer(Grain *pGrain1, Container *pContainer);

void computeCollisionWithPlan(Grain *pGrain1, Plan *pplan);

void computeCollisionBetweenBarrelAndPlan(Barrel *pBarrel, Plan *pplan);

// CONTACT PARAMETERS
double e = 0.99;
double mu = 0.1;
double kn = 1000.;
double kt = 100000.;
double dt = 10. * 0.000001;


int main(int argc, char **argv) {
    auto t1 = omp_get_wtime();

    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<double> uniformRealDistribution(0, 1);

    int i, j, k;


    // VIDEO OPTIONS
    int fps = 50;
    double tStartCapture = 0.;
    double totalTime = 2;
    int totalFrames = (int) ((totalTime - tStartCapture) * fps);
    double recTime;
    //TODO make nice constructor
    GrainPrinter grainPrinter;
    grainPrinter.setPath("datas/");
    BarrelPrinter barrelPrinter;
    barrelPrinter.setPath("datas/");

// GRAINS
    int numberOfGrains = 0;
    double radius = 0.05;
    double mass;
    double rho = 2000.;
    auto *grains = new Grain[numberOfGrains];
    int numberOfPlacedGrains = 0;
    int numberOfOverlaps;

    //Barrel
    double barrelRadius = 0.05;
    double barrelMass = 1.;
    auto *barrel = new Barrel;

//this setups the getRadius distribution
    double radiusMean = radius;
    std::uniform_real_distribution<double> radiusDistribution(radiusMean - (radiusMean / 2),
                                                              radiusMean + (radiusMean / 2));


//container
    double containerRadius = .3;

    int numberOfRevolution = 3;

    //plan et domaine
    double alphaDegree = 30;
    double alpha = alphaDegree * M_PI / 180;

    double xDomain = 2. * containerRadius * numberOfRevolution;

    Vector2 a(0, xDomain * tan(alpha));
    Vector2 b(xDomain, 0);
    Plan plan;
    plan.initPlanFromCoordinates(a, b);
    plan.printPlanInfos("datas/plan.txt");

    double yDomain = 2. * containerRadius * numberOfRevolution + plan.getPointFromX(0).getY();
    Domain domain(xDomain, yDomain);
    domain.printDomainInfos("datas/domain.txt");



//CONTAINER
    Vector2 containerCenter = plan.getPointFromX(containerRadius) + Vector2(0, containerRadius);
    Container container(containerRadius, containerCenter);

//ON place les grains et le barrel
    barrel->initBarrel(barrelRadius, barrelMass, plan.getPointFromX(barrelRadius) + Vector2(0, barrelRadius),
                       Vector2(0.));

    while (numberOfPlacedGrains < numberOfGrains) {
        radius = fabs(radiusDistribution(gen));
        numberOfOverlaps = 0;
        double direction = (double) (uniformRealDistribution(gen)) * 2 * M_PI;
        double randomRadius =
                (double) sqrt(uniformRealDistribution(gen)) * (barrelRadius) * .95;//sqrt pour que ce soit uniforme
        Vector2 randomPosition(randomRadius * cos(direction), randomRadius * sin(direction));
        randomPosition = randomPosition + barrel->getPosition();
        randomPosition.display();


//        randomPosition.setComponents(0.3-randomRadius, 0);
        for (i = 0; i < numberOfPlacedGrains; i++) { //Regarder si il n'y a pas d'overlap
            if (getDistanceBetweenVectors(randomPosition, grains[i].getPosition()) <
                radius + grains[i].getRadius()) {
                //BREAK
                numberOfOverlaps++;
                i = numberOfPlacedGrains;
            }
        }
        if (numberOfOverlaps == 0) {
            mass = 4. / 3. * M_PI * pow(radius, 3) * rho;
            grains[numberOfPlacedGrains].initDisk(numberOfPlacedGrains, radius, mass, randomPosition, Vector2(0.));
            numberOfPlacedGrains++;
        }
    }

    std::cout << "Grains and Barrel are placed" << std::endl;

//CLEAR FILES
    for (int l = 0; l <= totalFrames; l++) {
        grainPrinter.clearPrint(l);
    }

//record initial state
    for (i = 0; i < numberOfGrains; i++) {
        grainPrinter.print(grains[i], 0);
    }
//    barrelPrinter.print(reinterpret_cast<const Barrel &>(barrel), 0);
    std::cout << "Initial state in saved" << std::endl;


//linked cells
//TODO get max size
    double cellSize = domain.getX() / 5.;
    int nCellX = (int) ((domain.getX()) / cellSize);
    int nCellY = (int) ((domain.getY()) / cellSize);
    int nCell = nCellX * nCellY;
    Cell *cells = new Cell[nCell];
//    double dx = (containerRadius * 2.) / nCellX;
//    double dy = (containerRadius * 2.) / nCellY;
    std::cout << "Cells values are initialized" << std::endl;
    std::cout << "nCellX " << nCellX << std::endl;
    std::cout << "nCellY " << nCellY << std::endl;

    int ix, iy, jx, jy;
    for (i = 0; i < nCell; i++) {
        cells[i].initCell(i);
        iy = i / nCellX;
        ix = i % nCellX;
        if (ix > 0) { //if not first on the column
            cells[i].addNeighbor(ix + (iy * nCellX) - 1);
        }
        if (ix < nCellX - 1) { //if not last on the column
            cells[i].addNeighbor(ix + (iy * nCellX) + 1);
        }
        if (iy > 0) { //if not on first line
            cells[i].addNeighbor(ix + (iy - 1) * nCellX);
            if (ix > 0) { //if not first on the line
                cells[i].addNeighbor(ix + (iy - 1) * nCellX - 1);
            }
            if (ix < nCellX - 1) { //if not last on the line
                cells[i].addNeighbor(ix + (iy - 1) * nCellX + 1);
            }
        }
        if (iy < nCellY - 1) { //if not on last line
            cells[i].addNeighbor(ix + ((iy + 1) * nCellX));
            if (ix > 0) { //if not first on the line
                cells[i].addNeighbor(ix + ((iy + 1) * nCellX) - 1);
            }
            if (ix < nCellX - 1) { //if not last on the line
                cells[i].addNeighbor(ix + ((iy + 1) * nCellX) + 1);
            }
        }
    }

    std::cout << "Cells are positioned" << std::endl;
    //variables
    int cellIndex, hol;
    int neighborCellIndex, nNeighbors;
    double t;

    for (t = 0.; t < totalTime; t += dt) {
        /*** refresh and update position***/

        // reset linked cells
        for (i = 0; i < nCell; i++) {
            cells[i].setHeadOfList(-9);
        }


        //loop on grains

        for (i = 0; i < numberOfGrains; i++) {

            grains[i].updatePosition(dt / 2.);
            cellIndex = (int) (grains[i].getX() / cellSize) +
                        (int) ((grains[i].getY() / cellSize) * nCellX);
            grains[i].setLinkedCell(cellIndex);
            hol = cells[cellIndex].getHeadOfList();
            grains[i].setLinkedDisk(hol);
            cells[cellIndex].setHeadOfList(i);
            grains[i].resetForce();
            grains[i].addGravityForce(Vector2(0, -9.81));
        }


        /*** contact detection and forces ***/
        for (i = 0; i < numberOfGrains; i++) {
            // In cell
            cellIndex = grains[i].linkedCell();
            j = cells[cellIndex].getHeadOfList();
            while (j != -9) {
                if (i < j) {
                    computeCollision(&grains[i], &grains[j]);
                }
                j = grains[j].linkedDisk();
            }

            // In neighbor cells
            nNeighbors = cells[cellIndex].numberOfNeighbors();
            for (k = 0; k < nNeighbors; k++) {
                neighborCellIndex = cells[cellIndex].neighbor(k);
                j = cells[neighborCellIndex].getHeadOfList();
                while (j != -9) {
                    computeCollision(&grains[i], &grains[j]);
                    j = grains[j].linkedDisk();
                }
            }

            //Collisions with the container
            //computeCollisionWithContainer(&grains[i], &container);
            computeCollisionWithPlan(&grains[i], &plan);
            computeCollisionBetweenBarrelAndPlan(&barrel, &plan);
        }

        //update velocity and position
        for (i = 0; i < numberOfGrains; i++) {
            grains[i].updateVelocity(dt);
            grains[i].updatePosition(dt / 2.);
        }

        //record
        recTime = t - tStartCapture;
        if (recTime >= 0.) {
            if ((int) ((recTime + dt) * fps) > (int) (recTime * fps)) {
                for (i = 0; i < numberOfGrains; i++) {
                    grainPrinter.print(grains[i], (int) ((recTime + dt) * fps));
                }
                std::cout << "PRINTED IMAGE : " << (int) ((recTime + dt) * fps) << std::endl;
            }
        }
    }

    delete[] cells;
    delete[] grains;

    // Recording end time.
    auto t2 = omp_get_wtime();

    printf("Work took %f seconds", t2 - t1);
    return 0;

}

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
        pGrain2->addForce(-1. * tangentForce);

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

void computeCollisionWithPlan(Grain *pGrain1, Plan *plan) {
    Vector2 normalVecteur = -1 * plan->getNormal();
    Vector2 vecteur = plan->getPosition() - pGrain1->getPosition();
    double delta = vecteur.getX() * normalVecteur.getX() + vecteur.getY() * normalVecteur.getY();


    if (delta < 0) {

        Vector2 normalVector = normalVecteur;

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
        double eta = -2. * log(e) * sqrt(effectiveMass * kn / (pow(log(e), 2) + pow(M_PI, 2)));
        double normalForceNorm = -1. * (kn * delta) + (eta * normalVelocity.getNorm());
        double tangentForceNorm = -kt * tangentVelocity.getNorm();
        Vector2 tangentForce(tangentForceNorm * tangentVector.getX(), (tangentForceNorm * tangentVector.getY()));

        if (normalForceNorm > 0) {
            pGrain1->addForce(normalForceNorm * normalVector);
//         std::cout << (normalForceNorm * normalVector).getX() << "" << (normalForceNorm * normalVector).getY() << " N" << std::endl;
        } else {
            normalForceNorm = 0.;
        }

        if (tangentForce.getNorm() > mu * normalForceNorm) {
            pGrain1->addForce(-1 * mu * normalForceNorm * tangentVector);
            std::cout << "adding force" << std::endl;
        } else {
            pGrain1->addForce(tangentForce);
        }
        Vector2 forceVector = pGrain1->getForce();
//        std::cout << forceVector.getX() << " " << forceVector.getY() << std::endl;


//CE QUE J AI AJOUTE
        // check if normal force is repulsive
        if (normalForceNorm > 0) {
            Vector2 normalForce(tangentForceNorm * normalVector);
            pGrain1->addForce(normalForce);
        } else {
            double fn = 0.;
        }


        if (tangentForce.getNorm() > mu * normalForceNorm) {
//            ftx = -mu * normalForce.getNorm() * normalizedTangentVelocity();
//            fty = -mu * fn * ty;
            pGrain1->addForce(Vector2(0));
            pGrain1->addForce(-1. * Vector2(0));
        } else {
            pGrain1->addForce(Vector2(0));
        }
//JUSQUICI

        //torque
        double M = (-1 * pGrain1->getRadius() * normalVector.getX() * tangentForce.getY())
                   + (pGrain1->getRadius() * normalVector.getY() * tangentForce.getX());
//        std::cout << M << std::endl;

        pGrain1->addMomentum(M);
    }


}

void computeCollisionBetweenBarrelAndPlan(Barrel *pBarrel, Plan *plan) {
    Vector2 normalVecteur = -1 * plan->getNormal();
    Vector2 vecteur = plan->getPosition() - pBarrel->getPosition();
    double delta = vecteur.getX() * normalVecteur.getX() + vecteur.getY() * normalVecteur.getY();


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
        double eta = -2. * log(e) * sqrt(effectiveMass * kn / (pow(log(e), 2) + pow(M_PI, 2)));
        double normalForceNorm = -1. * (kn * delta) + (eta * normalVelocity.getNorm());
        double tangentForceNorm = -kt * tangentVelocity.getNorm();
        Vector2 tangentForce(tangentForceNorm * tangentVector.getX(), (tangentForceNorm * tangentVector.getY()));

        if (normalForceNorm > 0) {
            pBarrel->addForce(normalForceNorm * normalVector);
//         std::cout << (normalForceNorm * normalVector).getX() << "" << (normalForceNorm * normalVector).getY() << " N" << std::endl;
        } else {
            normalForceNorm = 0.;
        }

        if (tangentForce.getNorm() > mu * normalForceNorm) {
            pBarrel->addForce(-1 * mu * normalForceNorm * tangentVector);
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


        if (tangentForce.getNorm() > mu * normalForceNorm) {
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

