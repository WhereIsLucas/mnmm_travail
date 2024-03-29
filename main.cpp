#include <iostream>
#include <random>
#include <omp.h>
#include "Grain.h"
#include "Plan.h"
#include "Cell.h"
#include "GrainPrinter.h"
#include "Container.h"
#include "Vector2.h"
#include "Domain.h"
#include "collisions.h"

// CONTACT PARAMETERS



int main(int argc, char **argv) {

    auto t1 = omp_get_wtime();

    //initialize the random librairies
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<double> uniformRealDistribution(0, 1);

    //iterators
    int i, j, k;

    // COLLISIONS SETTINGS
    auto containerCollisionSettings = new CollisionSettings(.9,.6,1000.,1000000.);
    auto grainCollisionSettings = new CollisionSettings(.9,.6,8000.,1000000.);

    // VIDEO SETTINGS for output data
    int fps = 25;
    double tStartCapture = 0.;
    double totalTime = 10;
    int totalFrames = (int) ((totalTime - tStartCapture) * fps);
    double recTime;
    //dt is the delta used for calculate the interframes
    double dt = 1. * 0.00001;

    GrainPrinter grainPrinter("data/");

    // GRAINS
    int numberOfGrains = 20;
    double radius = 0.015;
    double mass;
    double rho = 2000.;
    auto *grains = new Grain[numberOfGrains];
    int numberOfPlacedGrains = 0;
    int numberOfOverlaps;

    // Container is the ball
    double containerRadius = .3;
    Vector2 containerCenter(containerRadius);
    Container container(containerRadius, containerCenter);

    //DOMAIN
    double xDomain = 2. * containerRadius * 1.2;
    double yDomain = 4. * containerRadius * 1.2;
    Domain domain(xDomain, yDomain);
    domain.printDomainInfos("data/domain.txt");

    //Placing the grains
    while (numberOfPlacedGrains < numberOfGrains) {
        numberOfOverlaps = 0;
        //We choose a random direction and a random radius
        double direction = (double) (uniformRealDistribution(gen)) * 2 * M_PI;
        double randomRadius = (double) sqrt(uniformRealDistribution(gen)) * (containerRadius - 2. * radius) *
                              .95;
        Vector2 randomPosition(randomRadius * cos(direction), randomRadius * sin(direction));
//        randomPosition.setComponents( (-.2)+(.2*numberOfPlacedGrains),0);
//        randomPosition.setComponents(-.1,0.);
        randomPosition = randomPosition + containerCenter;
        // We check for an overlap
        for (i = 0; i < numberOfPlacedGrains; i++) {
            if (getDistanceBetweenVectors(randomPosition, grains[i].getPosition()) <
                radius + grains[i].getRadius()) {
                numberOfOverlaps++;
                // exiting the loop
                i = numberOfPlacedGrains;
                std::cout <<  "Overlap" << std::endl;
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
    double cellSize = 2.2 * radius;
    cellSize = domain.getX()/5.;
    int nCellX = (int) ((domain.getX()) / cellSize);
    int nCellY = (int) ((domain.getY()) / cellSize);
    int nCell = nCellX * nCellY;
    Cell *cells = new Cell[nCell];

    std::cout << "Cells values are initialized" << std::endl;
    std::cout << "nCellX " << nCellX << std::endl;
    std::cout << "nCellY " << nCellY << std::endl;
    std::cout << "nCells " << nCell << std::endl;

    int ix, iy;
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
    double t = 0.;
    auto t1_image = omp_get_wtime();
    int iteration = 0;
    while (true) {
        iteration++;
        t += dt;
        if(t > totalTime){
            break;
        }

        // reset linked cells
        for (i = 0; i < nCell; i++) {
            cells[i].setHeadOfList(-9);
        }

        //loop on grains
        for (i = 0; i < numberOfGrains; i++) {

            /* Leap frog step 1 */
            grains[i].updatePosition(dt / 2.);
            cellIndex = (int) (grains[i].getX() / cellSize) +
                        (int) ((grains[i].getY() / cellSize) * nCellX);
            if(abs(cellIndex) > nCell-1){
                grains[i].getPosition().display();
                std::cout << "Out of domain limits" << std::endl;
                exit(1);
            }

//            std::cout << grains[i].getY() << std::endl;
            grains[i].setLinkedCell(cellIndex);
            hol = cells[cellIndex].getHeadOfList();
            grains[i].setLinkedDisk(hol);
            cells[cellIndex].setHeadOfList(i);
            grains[i].resetForce();
            grains[i].addGravityForce(Vector2(0, -9.81));
        }


        /*** contact detection and forces ***/
        for (i = 0; i < numberOfGrains; i++) {

//             In cell
            cellIndex = grains[i].linkedCell();
            j = cells[cellIndex].getHeadOfList();
            while (j != -9) {
                if (i < j) {
                    computeCollisionWithGrain(&grains[i], &grains[j], grainCollisionSettings);
                }
                j = grains[j].linkedDisk();
            }

            // In neighbor cells
            nNeighbors = cells[cellIndex].numberOfNeighbors();
            for (k = 0; k < nNeighbors; k++) {

                neighborCellIndex = cells[cellIndex].neighbor(k);
                j = cells[neighborCellIndex].getHeadOfList();
                while (j != -9) {
                    computeCollisionWithGrain(&grains[i], &grains[j], grainCollisionSettings);
                    j = grains[j].linkedDisk();
                }
            }

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
                auto t2_image = omp_get_wtime();
                std::cout << "printing image " << (int) ((recTime + dt) * fps) << " | " <<  t2_image - t1_image << " seconds" << std::endl;
                t1_image = omp_get_wtime();
            }
        }
    }

    delete[] cells;
    delete[] grains;
    std::cout << "PRINTED " << (int) ((recTime + dt) * fps) << " images" << std::endl;

    // Recording end time.
    auto t2 = omp_get_wtime();

    printf("Work took %f seconds", t2 - t1);
    return 0;

}


