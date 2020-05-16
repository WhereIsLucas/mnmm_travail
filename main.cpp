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
#include "Domain.h"
#include "collisions.h"
#include <chrono>

// CONTACT PARAMETERS

double dt = 20. * 0.000001;


int main(int argc, char **argv) {
    auto t1 = omp_get_wtime();

    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<double> uniformRealDistribution(0, 1);

    int i, j, k;
    // COLLISIONS SETTINGS
    auto containerCollisionSettings = new CollisionSettings(.99,.6,200.,100000.);
    auto particlesCollisionSettings = new CollisionSettings(.99,.6,200.,100000.);

    // VIDEO OPTIONS
    int fps = 50;
    double tStartCapture = 0.;
    double totalTime = 2;
    int totalFrames = (int) ((totalTime - tStartCapture) * fps);
    double recTime;
    GrainPrinter grainPrinter("datas/");

    // GRAINS
    int numberOfGrains = 20;
    double radius = 0.0005;
    double mass;
    double rho = 2000.;
    auto *grains = new Grain[numberOfGrains];
    int numberOfPlacedGrains = 0;
    int numberOfOverlaps;

    //this setups the getRadius distribution
    double radiusMean = radius;
    std::uniform_real_distribution<double> radiusDistribution(radiusMean - (radiusMean / 2),
                                                              radiusMean + (radiusMean / 2));


    //container
    double containerRadius = .3;
    Vector2 containerCenter(containerRadius);
    Container container(containerRadius, containerCenter);

    //DOMAIN
    double xDomain = 2. * containerRadius * 1.2;
    double yDomain = 2. * containerRadius * 1.2;
    Domain domain(xDomain, yDomain);
    domain.printDomainInfos("datas/domain.txt");

    //Placing the grains
    while (numberOfPlacedGrains < numberOfGrains) {
        radius = fabs(radiusDistribution(gen));
        numberOfOverlaps = 0;
        //We choose a random direction and an random radius
        double direction = (double) (uniformRealDistribution(gen)) * 2 * M_PI;
        double randomRadius = (double) sqrt(uniformRealDistribution(gen)) * (containerRadius - 2. * radius) *
                              .95;
        Vector2 randomPosition(randomRadius * cos(direction), randomRadius * sin(direction));
        randomPosition = randomPosition + containerCenter;
//        randomPosition.setComponents(0.3-randomRadius, 0);
        // We check for an overlap
        for (i = 0; i < numberOfPlacedGrains; i++) { //Regarder si il n'y a pas d'overlap
            if (getDistanceBetweenVectors(randomPosition, grains[i].getPosition()) <
                radius + grains[i].getRadius()) {
                numberOfOverlaps++;
                // exiting the loop
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
    double cellSize = 2.2 * radiusMean;
    int nCellX = (int) ((domain.getX()) / cellSize);
    int nCellY = (int) ((domain.getY()) / cellSize);
    int nCell = nCellX * nCellY;
    Cell *cells = new Cell[nCell];

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
//            cellIndex = grains[i].linkedCell();
//            j = cells[cellIndex].getHeadOfList();
//            while (j != -9) {
//                if (i < j) {
//                    computeCollision(&grains[i], &grains[j]);
//                }
//                j = grains[j].linkedDisk();
//            }

//            // In neighbor cells
//            nNeighbors = cells[cellIndex].numberOfNeighbors();
//            for (k = 0; k < nNeighbors; k++) {
//                neighborCellIndex = cells[cellIndex].neighbor(k);
//                j = cells[neighborCellIndex].getHeadOfList();
//                while (j != -9) {
//                    computeCollision(&grains[i], &grains[j]);
//                    j = grains[j].linkedDisk();
//                }
//            }

            //Collisions with the container
            computeCollisionWithContainer(&grains[i], &container, containerCollisionSettings);

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


