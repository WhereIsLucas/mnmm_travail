#include <iostream>
#include <random>
#include <map>
#include <omp.h>
#include "Grain.h"
#include "Plan.h"
#include "Cell.h"
#include "GrainPrinter.h"
#include <chrono>

int main(int argc, char **argv) {
    auto t1 = omp_get_wtime();

    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<double> uniformRealDistribution(0, 1);

    int i, j, k;

    // CONTACT PARAMETERS
    double e = 0.9;
    double mu = 0.6;
    double kn = 200.;
    double kt = 100000.;
    double dt = 0.000001;

    // VIDEO OPTIONS
    int fps = 50;
    double tStartCapture = 0.;
    double totalTime = .2;
    int totalFrames = (int) ((totalTime - tStartCapture) * fps);
    double recTime;
    GrainPrinter grainPrinter;
    grainPrinter.setPath("datas/");

    // GRAINS
    int numberOfGrains = 100;
    double radius = 0.0005;
    double mass;
    double x, y, vx, vy;
    double eta;
    double rho = 2000.;
    auto *grains = new Grain[numberOfGrains];
    int numberOfPlacedGrains = 0;
    int numberOfOverlaps;

    //this setups the radius distribution
    double radiusMean = radius;
    std::uniform_real_distribution<double> radiusDistribution(radiusMean - (radiusMean / 3),
                                                              radiusMean + (radiusMean / 3));


    //container
//    double containerRadius = .3;
    double lx = 0.03;
    double ly = 0.06;
    double amplx = 0.;
    double amply = 0.;
    double freqx = 0.;
    double freqy = 0.;
    Plan *wall = new Plan[4];
    wall[0].initPlan(0., -ly / 2., 0., 1., amplx, amply, freqx, freqy);
    wall[1].initPlan(0., ly / 2., 0., -1., amplx, amply, freqx, freqy);
    wall[2].initPlan(-lx / 2., 0., 1., 0., amplx, amply, freqx, freqy);
    wall[3].initPlan(lx / 2., 0., -1., 0., amplx, amply, freqx, freqy);

    while (numberOfPlacedGrains < numberOfGrains) {
        radius = fabs(radiusDistribution(gen));
        numberOfOverlaps = 0;
        x = -lx / 2. + radius + ((double) (uniformRealDistribution(gen)) * (lx - 2. * radius));
        y = -ly / 2. + radius + ((double) (uniformRealDistribution(gen)) * (ly - 2. * radius));
        for (i = 0; i < numberOfPlacedGrains; i++) {
            if ((x - grains[i].x()) * (x - grains[i].x()) + (y - grains[i].y()) * (y - grains[i].y()) <
                (radius + grains[i].radius()) * (radius + grains[i].radius())) {
                numberOfOverlaps++;
                i = numberOfPlacedGrains;
            }
        }
        if (numberOfOverlaps == 0) {
            mass = 4. / 3. * M_PI * radius * radius * radius * rho;
            vx = 0.;
            vy = 0.;
            grains[numberOfPlacedGrains].initDisk(numberOfPlacedGrains, radius, mass, x, y, vx, vy);
            numberOfPlacedGrains++;
        }
    }

    //CLEAR FILES
    for (int l = 0; l <= totalFrames; l++) {
        grainPrinter.clearPrint(l);
    }

    //record initial state
    for (i = 0; i < numberOfGrains; i++) {
        grainPrinter.print(grains[i], 0);
    }

    //linked cells
    double cellSize = 2.2 * radius;
    int nCellx = (int) ((lx + 2. * amplx) / cellSize);
    int nCelly = (int) ((ly + 2. * amply) / cellSize);
    int nCell = nCellx * nCelly;
    Cell *cellule = new Cell[nCell];
    double dx = (lx + 2. * amplx) / nCellx;
    double dy = (ly + 2. * amply) / nCelly;

    int ix, iy, jx, jy;
    for (i = 0; i < nCell; i++) {
        iy = i / nCellx;
        ix = i % nCellx;
        cellule[i].initCell(i);
        for (j = i + 1; j < nCell; j++) {
            jy = j / nCellx;
            jx = j % nCellx;
            if (abs(ix - jx) < 2 && abs(iy - jy) < 2) {
                cellule[i].addNeighbor(j);
            }
        }
    }

    /********************/

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

        // linked cells
        for (i = 0; i < nCell; i++) {
            cellule[i].setHeadOfList(-9);
        }

        //grains
        for (i = 0; i < numberOfGrains; i++) {
            //update positions
            grains[i].updatePosition(dt / 2.);

            //class in linked cells
            cellIndex = (int) ((grains[i].x() + lx / 2. + amplx) / dx) +
                        (int) ((grains[i].y() + ly / 2. + amply) / dy) * nCellx;
            grains[i].setLinkedCell(cellIndex);
            hol = cellule[cellIndex].headOfList();
            grains[i].setLinkedDisk(hol);
            cellule[cellIndex].setHeadOfList(i);

            //reset force
            grains[i].resetForce();

            //set gravity force
            grains[i].addGravityForce(0., -9.81);
        }

        //container
        for (i = 0; i < 4; i++) {
            wall[i].updatePosition(t + dt / 2.);
        }

        /*** contact detection and forces ***/
        for (i = 0; i < numberOfGrains; i++) {
            //in same Cell
            cellIndex = grains[i].linkedCell();
            j = cellule[cellIndex].headOfList();
            while (j != -9) {
                if (i < j) {
                    nx = grains[i].x() - grains[j].x();
                    ny = grains[i].y() - grains[j].y();
                    rij = sqrt(nx * nx + ny * ny);
                    nx /= rij;
                    ny /= rij;

                    delta = rij - (grains[i].radius() + grains[j].radius());
                    if (delta < 0) {
                        //contact base
                        vx = grains[i].vx() - grains[j].vx() + grains[i].radius() * grains[i].w() * ny +
                             grains[j].radius() * grains[j].w() * ny;
                        vy = grains[i].vy() - grains[j].vy() - grains[i].radius() * grains[i].w() * nx -
                             grains[j].radius() * grains[j].w() * nx;
                        vn = vx * nx + vy * ny;
                        vnx = vn * nx;
                        vny = vn * ny;
                        vtx = vx - vnx;
                        vty = vy - vny;
                        vt = sqrt(vtx * vtx + vty * vty);
                        if (vt != 0.) {
                            tx = vtx / vt;
                            ty = vty / vt;
                        } else {
                            tx = 0.;
                            ty = 0.;
                        }

                        //contact forces and torque
                        effectiveMass = (grains[i].mass() * grains[j].mass()) / (grains[i].mass() + grains[j].mass());
                        eta = -2. * log(e) * sqrt(effectiveMass * kn / (log(e) * log(e) + M_PI * M_PI));
                        fn = -kn * delta - eta * vn;
                        ft = -kt * vt;

                        //forces
                        if (fn > 0) {
                            fnx = fn * nx;
                            fny = fn * ny;
                            grains[i].addForce(fnx, fny);
                            grains[j].addForce(-fnx, -fny);
                        } else {
                            fn = 0.;
                        }

                        if (fabs(ft) > mu * fn) {
                            ftx = -mu * fn * tx;
                            fty = -mu * fn * ty;
                            grains[i].addForce(ftx, fty);
                            grains[j].addForce(-ftx, -fty);
                        } else {
                            ftx = ft * tx;
                            fty = ft * ty;
                            grains[i].addForce(ftx, fty);
                            grains[j].addForce(-ftx, -fty);
                        }

                        //torque
                        M = -grains[i].radius() * nx * fty + grains[i].radius() * ny * ftx;
                        grains[i].addMomentum(M);
                        M = -grains[j].radius() * nx * fty + grains[j].radius() * ny * ftx;
                        grains[j].addMomentum(M);
                    }
                }
                j = grains[j].linkedDisk();
            }

            //in neighbor cells
            nNeighbors = cellule[cellIndex].numberOfNeighbors();
            for (k = 0; k < nNeighbors; k++) {
                neighborCellIndex = cellule[cellIndex].neighbor(k);
                j = cellule[neighborCellIndex].headOfList();
                while (j != -9) {
                    nx = grains[i].x() - grains[j].x();
                    ny = grains[i].y() - grains[j].y();
                    rij = sqrt(nx * nx + ny * ny);
                    nx /= rij;
                    ny /= rij;

                    delta = rij - (grains[i].radius() + grains[j].radius());
                    if (delta < 0) {
                        //contact base
                        vx = grains[i].vx() - grains[j].vx() + grains[i].radius() * grains[i].w() * ny +
                             grains[j].radius() * grains[j].w() * ny;
                        vy = grains[i].vy() - grains[j].vy() - grains[i].radius() * grains[i].w() * nx -
                             grains[j].radius() * grains[j].w() * nx;
                        vn = vx * nx + vy * ny;
                        vnx = vn * nx;
                        vny = vn * ny;
                        vtx = vx - vnx;
                        vty = vy - vny;
                        vt = sqrt(vtx * vtx + vty * vty);
                        if (vt != 0.) {
                            tx = vtx / vt;
                            ty = vty / vt;
                        } else {
                            tx = 0.;
                            ty = 0.;
                        }

                        //contact forces and torque
                        effectiveMass = (grains[i].mass() * grains[j].mass()) / (grains[i].mass() + grains[j].mass());
                        eta = -2. * log(e) * sqrt(effectiveMass * kn / (log(e) * log(e) + M_PI * M_PI));
                        fn = -kn * delta - eta * vn;
                        ft = -kt * vt;

                        //forces
                        if (fn > 0) {
                            fnx = fn * nx;
                            fny = fn * ny;
                            grains[i].addForce(fnx, fny);
                            grains[j].addForce(-fnx, -fny);
                        } else {
                            fn = 0.;
                        }

                        if (fabs(ft) > mu * fn) {
                            ftx = -mu * fn * tx;
                            fty = -mu * fn * ty;
                            grains[i].addForce(ftx, fty);
                            grains[j].addForce(-ftx, -fty);
                        } else {
                            ftx = ft * tx;
                            fty = ft * ty;
                            grains[i].addForce(ftx, fty);
                            grains[j].addForce(-ftx, -fty);
                        }

                        //torque
                        M = -grains[i].radius() * nx * fty + grains[i].radius() * ny * ftx;
                        grains[i].addMomentum(M);
                        M = -grains[j].radius() * nx * fty + grains[j].radius() * ny * ftx;
                        grains[j].addMomentum(M);
                    }
                    j = grains[j].linkedDisk();
                }
            }

            //with walls
            for (j = 0; j < 4; j++) {
                nx = wall[j].nx();
                ny = wall[j].ny();
                rx = grains[i].x() - wall[j].x();
                ry = grains[i].y() - wall[j].y();

                delta = (rx * nx + ry * ny) - grains[i].radius();
                if (delta < 0) {
                    vx = grains[i].vx() - wall[j].vx() + grains[i].radius() * grains[i].w() * ny;
                    vy = grains[i].vy() - wall[j].vy() - grains[i].radius() * grains[i].w() * nx;
                    vn = vx * nx + vy * ny;
                    vnx = vn * nx;
                    vny = vn * ny;
                    vtx = vx - vnx;
                    vty = vy - vny;
                    vt = sqrt(vtx * vtx + vty * vty);
                    if (vt != 0.) {
                        tx = vtx / vt;
                        ty = vty / vt;
                    } else {
                        tx = 0.;
                        ty = 0.;
                    }

                    //contact forces and torque
                    effectiveMass = grains[i].mass();
                    eta = -2. * log(e) * sqrt(effectiveMass * kn / (log(e) * log(e) + M_PI * M_PI));
                    fn = -kn * delta - eta * vn;
                    ft = -kt * vt;

                    //forces
                    if (fn > 0) {
                        fnx = fn * nx;
                        fny = fn * ny;
                        grains[i].addForce(fnx, fny);
                    } else {
                        fn = 0.;
                    }


                    if (fabs(ft) > mu * fn) {
                        ftx = -mu * fn * tx;
                        fty = -mu * fn * ty;
                        grains[i].addForce(ftx, fty);
                    } else {
                        ftx = ft * tx;
                        fty = ft * ty;
                        grains[i].addForce(ftx, fty);
                    }

                    //torque
                    M = -grains[i].radius() * nx * fty + grains[i].radius() * ny * ftx;
                    grains[i].addMomentum(M);
                }
            }
        }

        //update velocity and position
        for (i = 0; i < numberOfGrains; i++) {
            grains[i].updateVelocity(dt);
            grains[i].updatePosition(dt / 2.);
        }
        for (i = 0; i < 4; i++) {
            wall[i].updateVelocity(t + dt);
            wall[i].updatePosition(t + dt);
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

    delete[] cellule;
    delete[] wall;
    delete[] grains;

    // Recording end time.
    auto t2 = omp_get_wtime();

    printf("Work took %f seconds", t2 - t1);


    return 0;

}