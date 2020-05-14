#include <iostream>
#include <random>
#include <map>
#include <omp.h>
#include "Grain.h"
#include "Plan.h"
#include "Cell.h"
#include "GrainPrinter.h"
#include "Container.h"
#include "Vector2.h"
#include <chrono>



void computeCollision(Grain *pGrain1, Grain *pGrain2);

void computeCollisionWithContainer(Grain *pGrain1, Container *pContainer);

// CONTACT PARAMETERS
double e = 0.99;
double mu = 0.1;
double kn = 200.;
double kt = 100000.;
double dt = 4. * 0.000001;


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

// GRAINS
    int numberOfGrains = 200;
    double radius = 0.0005;
    double mass;
    double rho = 2000.;
    auto *grains = new Grain[numberOfGrains];
    int numberOfPlacedGrains = 0;
    int numberOfOverlaps;

//this setups the getRadius distribution
    double radiusMean = radius;
    std::uniform_real_distribution<double> radiusDistribution(radiusMean - (radiusMean / 3),
                                                              radiusMean + (radiusMean / 3));


//container
    double containerRadius = .3;
    Vector2 containerCenter(0, 0);
    Container container(containerRadius, containerCenter);

    while (numberOfPlacedGrains < numberOfGrains) {
        radius = fabs(radiusDistribution(gen));
        numberOfOverlaps = 0;
        double direction = (double) (uniformRealDistribution(gen)) * 2 * M_PI;
        double randomRadius = (double) sqrt(uniformRealDistribution(gen)) * (containerRadius - 2. * radius) * .95;
        Vector2 randomPosition(randomRadius * cos(direction), randomRadius * sin(direction));

//        randomPosition.setComponents(0.3-randomRadius, 0);
        for (i = 0; i < numberOfPlacedGrains; i++) {
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

    std::cout << "Grains are placed" << std::endl;

//CLEAR FILES
    for (int l = 0; l <= totalFrames; l++) {
        grainPrinter.clearPrint(l);
    }

//record initial state
    for (i = 0; i < numberOfGrains; i++) {
        grainPrinter.print(grains[i], 0);
    }
    std::cout << "Initial state in saved" << std::endl;


//linked cells
//TODO get max size
    double cellSize = 2.2 * radiusMean;
    int nCellX = (int) ((containerRadius) / cellSize);
    int nCellY = (int) ((containerRadius) / cellSize);
    int nCell = nCellX * nCellY;
    Cell *cells = new Cell[nCell];
    double dx = (containerRadius * 2.) / nCellX;
    double dy = (containerRadius * 2.) / nCellY;
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
    double rx, ry, nx, ny, tx, ty;
    double vn, vnx, vny, vt, vtx, vty;
    double rij, delta;
    double effectiveMass;
    double fn, fnx, fny;
    double ft, ftx, fty;
    double M, t;

    for (t = 0.; t < totalTime; t += dt) {
        /*** refresh and update position***/

        // reset linked cells
        for (i = 0; i < nCell; i++) {
            cells[i].setHeadOfList(-9);
        }

        //loop on grains
        for (i = 0; i < numberOfGrains; i++) {

            grains[i].updatePosition(dt / 2.);
            cellIndex = (int) ((grains[i].getX() + container.getRadius()) / dx) +
                        (int) ((grains[i].getY() + container.getRadius()) / dy) * nCellX;
            grains[i].setLinkedCell(cellIndex);
            hol = cells[cellIndex].getHeadOfList();
            grains[i].setLinkedDisk(hol);
            cells[cellIndex].setHeadOfList(i);
            grains[i].resetForce();
            grains[i].addGravityForce(Vector2(0, -9.81));
        }


        /*** contact detection and forces ***/
        for (i = 0; i < numberOfGrains; i++) {
//            In cell
            cellIndex = grains[i].linkedCell();
            j = cells[cellIndex].getHeadOfList();
            while (j != -9) {
                if (i < j) {
                    computeCollision(&grains[i], &grains[j]);
                }
                j = grains[j].linkedDisk();
            }
//in neighbor cells
            nNeighbors = cells[cellIndex].numberOfNeighbors();
            for (k = 0; k < nNeighbors; k++) {
                neighborCellIndex = cells[cellIndex].neighbor(k);
                j = cells[neighborCellIndex].getHeadOfList();
                while (j != -9) {
                    computeCollision(&grains[i], &grains[j]);
                    j = grains[j].linkedDisk();
                }
            }

//            Collisions with the container
            computeCollisionWithContainer(&grains[i], &container);
//
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

void computeCollision(Grain *pGrain1, Grain *pGrain2) {
    Vector2 normalVector = (pGrain2->getPosition() - pGrain1->getPosition()).normalize();
    double delta = getDistanceBetweenGrains(*pGrain1, *pGrain2);
    if (delta < 0) {
        double vy = pGrain1->getVy() - pGrain2->getVy();
//                    - pGrain1->getRadius() * pGrain1->getOmega() * normalVector.getX() -
//                    pGrain2->getRadius() * pGrain2->getOmega() * normalVector.getX();
        double vx = pGrain1->getVx() - pGrain2->getVx();
//                    + pGrain1->getRadius() * pGrain1->getOmega() * normalVector.getY() +
//                    pGrain2->getRadius() * pGrain2->getOmega() * normalVector.getY();
        Vector2 velocityAtContactPoint(vx, vy);
        Vector2 normalVelocity(velocityAtContactPoint.getX() * normalVector.getX(),
                               velocityAtContactPoint.getY() * normalVector.getY());
        Vector2 tangentVelocity = velocityAtContactPoint - normalVelocity;
        Vector2 tangentVector;
        if (tangentVelocity.getNorm() != 0.) {
            tangentVector = tangentVelocity.normalize();
        } else {
            tangentVector.setComponents(0, 0);
        }

        //contact forces and torque
        double effectiveMass = (pGrain1->getMass() * pGrain2->getMass()) /
                               (pGrain1->getMass() + pGrain2->getMass());
        double eta = -2. * log(e) * sqrt(effectiveMass * kn / (log(e) * log(e) + M_PI * M_PI));
        double normalForceNorm = -1.*(kn * delta + eta * normalVelocity.getNorm());
        double tangentForceNorm = -kt * tangentVelocity.getNorm();
        Vector2 tangentForce(-kt * tangentVelocity);

        // check if normal force is repulsive
        if (normalForceNorm > 0) {
            Vector2 normalForce(tangentForceNorm * normalVector);
            pGrain1->addForce(normalForce);
            pGrain2->addForce(-1 * normalForce);
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
            pGrain2->addForce(Vector2(0));
        }

//torque
//        double M = -pGrain1->getRadius() * nx * fty + pGrain1->getRadius() * ny * ftx;
//        pGrain1->addMomentum(M);
//        double M2 = -pGrain2->getRadius() * nx * fty + pGrain2->getRadius() * ny * ftx;
//        pGrain2->addMomentum(M2);
    }
}

void computeCollisionWithContainer(Grain *pGrain1, Container *container) {
    double delta = container->getRadius() - pGrain1->getRadius() -
                   getDistanceBetweenVectors(pGrain1->getPosition(), container->getCenter());
    if (delta < 0) {
        Vector2 normalVector = (container->getCenter() - pGrain1->getPosition()).normalize();
//        std::cout << normalVector.getX() << " " << normalVector.getY() << std::endl;
//        std::cout << pGrain1->getX() << " " << pGrain1->getY() << std::endl;
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
//            std::cout << (normalForceNorm * normalVector).getX() << "" << (normalForceNorm * normalVector).getY() << " N" << std::endl;
        } else {
            normalForceNorm = 0.;
        }

        if (tangentForce.getNorm() > mu * normalForceNorm) {
            pGrain1->addForce(-1 * mu * normalForceNorm * tangentVector);
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


