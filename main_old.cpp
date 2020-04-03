#include <cmath>
#include <cstdlib>
#include <string>
#include <ctime>
#include <unistd.h>
#include <sstream>
#include <iostream>

#include "Cell.h"
#include "Grain.h"
#include "Plan.h"

int main(int argc, char **argv) {
    srand(time(NULL));
    int i, j, k;

    //contact
    double e = 0.9;
    double mu = 0.6;
    double kn = 200.;
    double kt = 100000.;
    double dt = 0.000001;

    //video
    int fps = 100;
    double tStartCapture = 0.;
    double totalTime = 1;
    int totalFrames = (int) (totalTime - tStartCapture) * fps;
    double recTime;

    //container
    double containerRadius = .3;
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

    //grains
    int numberOfGrains = 100;
    double radius;
    double mass;
    double x, y, vx, vy;
    double eta;
    double rho = 2000.;
    Disk *grain = new Disk[numberOfGrains];
    int numberOfPlacedGrains = 0;
    int numberOfOverlaps;

    //place grains
    while (numberOfPlacedGrains < numberOfGrains) {
        numberOfOverlaps = 0;
        radius = 0.0005;
        mass = 4. / 3. * M_PI * radius * radius * radius * rho;
        x = -lx / 2. + radius + ((double) (rand()) / RAND_MAX) * (lx - 2. * radius);
        y = -ly / 2. + radius + ((double) (rand()) / RAND_MAX) * (ly - 2. * radius);
        vx = 0.;
        vy = 0.;
        for (i = 0; i < numberOfPlacedGrains; i++) {
            if ((x - grain[i].x()) * (x - grain[i].x()) + (y - grain[i].y()) * (y - grain[i].y()) <
                (radius + grain[i].radius()) * (radius + grain[i].radius())) {
                numberOfOverlaps++;
                i = numberOfPlacedGrains;
            }
        }
        if (numberOfOverlaps == 0) {
            grain[numberOfPlacedGrains].initDisk(numberOfPlacedGrains, radius, mass, x, y, vx, vy);
            numberOfPlacedGrains++;
        }
    }

    //CLEAR FILES


    for (int l = 0; l <= totalFrames; l++) {
        std::string fileName = "grain" + std::to_string(l) + ".txt";
        remove(fileName.c_str());
    }

    //record initial state
    for (i = 0; i < numberOfGrains; i++) {

        grain[i].print(0);
    }

    //linked cells
    double cellSize = 2.2 * radius;
    int nCellx = (lx + 2. * amplx) / cellSize;
    int nCelly = (ly + 2. * amply) / cellSize;
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
    double Fn, Fnx, Fny;
    double Ft, Ftx, Fty;
    double M, t;

    for (t = 0.; t < totalTime; t += dt) {
        /*** refresh and update position***/

        //linked cells
        for (i = 0; i < nCell; i++) {
            cellule[i].setHeadOfList(-9);
        }

        //grains
        for (i = 0; i < numberOfGrains; i++) {
            //update positions
            grain[i].updatePosition(dt / 2.);

            //class in linked cells
            cellIndex = (int) ((grain[i].x() + lx / 2. + amplx) / dx) +
                        (int) ((grain[i].y() + ly / 2. + amply) / dy) * nCellx;
            grain[i].setLinkedCell(cellIndex);
            hol = cellule[cellIndex].headOfList();
            grain[i].setLinkedDisk(hol);
            cellule[cellIndex].setHeadOfList(i);

            //reset force
            grain[i].resetForce();

            //set gravity force
            grain[i].addGravityForce(0., -9.81);
        }

        //container
        for (i = 0; i < 4; i++) {
            wall[i].updatePosition(t + dt / 2.);
        }

        /*** contact detection and forces ***/

        for (i = 0; i < numberOfGrains; i++) {
            //in same Cell
            cellIndex = grain[i].linkedCell();
            j = cellule[cellIndex].headOfList();
            while (j != -9) {
                if (i < j) {
                    nx = grain[i].x() - grain[j].x();
                    ny = grain[i].y() - grain[j].y();
                    rij = sqrt(nx * nx + ny * ny);
                    nx /= rij;
                    ny /= rij;

                    delta = rij - (grain[i].radius() + grain[j].radius());
                    if (delta < 0) {
                        //contact base
                        vx = grain[i].vx() - grain[j].vx() + grain[i].radius() * grain[i].w() * ny +
                             grain[j].radius() * grain[j].w() * ny;
                        vy = grain[i].vy() - grain[j].vy() - grain[i].radius() * grain[i].w() * nx -
                             grain[j].radius() * grain[j].w() * nx;
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
                        effectiveMass = (grain[i].mass() * grain[j].mass()) / (grain[i].mass() + grain[j].mass());
                        eta = -2. * log(e) * sqrt(effectiveMass * kn / (log(e) * log(e) + M_PI * M_PI));
                        Fn = -kn * delta - eta * vn;
                        Ft = -kt * vt;

                        //forces
                        if (Fn > 0) {
                            Fnx = Fn * nx;
                            Fny = Fn * ny;
                            grain[i].addForce(Fnx, Fny);
                            grain[j].addForce(-Fnx, -Fny);
                        } else {
                            Fn = 0.;
                        }

                        if (fabs(Ft) > mu * Fn) {
                            Ftx = -mu * Fn * tx;
                            Fty = -mu * Fn * ty;
                            grain[i].addForce(Ftx, Fty);
                            grain[j].addForce(-Ftx, -Fty);
                        } else {
                            Ftx = Ft * tx;
                            Fty = Ft * ty;
                            grain[i].addForce(Ftx, Fty);
                            grain[j].addForce(-Ftx, -Fty);
                        }

                        //torque
                        M = -grain[i].radius() * nx * Fty + grain[i].radius() * ny * Ftx;
                        grain[i].addMomentum(M);
                        M = -grain[j].radius() * nx * Fty + grain[j].radius() * ny * Ftx;
                        grain[j].addMomentum(M);
                    }
                }
                j = grain[j].linkedDisk();
            }

            //in neighbor cells
            nNeighbors = cellule[cellIndex].numberOfNeighbors();
            for (k = 0; k < nNeighbors; k++) {
                neighborCellIndex = cellule[cellIndex].neighbor(k);
                j = cellule[neighborCellIndex].headOfList();
                while (j != -9) {
                    nx = grain[i].x() - grain[j].x();
                    ny = grain[i].y() - grain[j].y();
                    rij = sqrt(nx * nx + ny * ny);
                    nx /= rij;
                    ny /= rij;

                    delta = rij - (grain[i].radius() + grain[j].radius());
                    if (delta < 0) {
                        //contact base
                        vx = grain[i].vx() - grain[j].vx() + grain[i].radius() * grain[i].w() * ny +
                             grain[j].radius() * grain[j].w() * ny;
                        vy = grain[i].vy() - grain[j].vy() - grain[i].radius() * grain[i].w() * nx -
                             grain[j].radius() * grain[j].w() * nx;
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
                        effectiveMass = (grain[i].mass() * grain[j].mass()) / (grain[i].mass() + grain[j].mass());
                        eta = -2. * log(e) * sqrt(effectiveMass * kn / (log(e) * log(e) + M_PI * M_PI));
                        Fn = -kn * delta - eta * vn;
                        Ft = -kt * vt;

                        //forces
                        if (Fn > 0) {
                            Fnx = Fn * nx;
                            Fny = Fn * ny;
                            grain[i].addForce(Fnx, Fny);
                            grain[j].addForce(-Fnx, -Fny);
                        } else {
                            Fn = 0.;
                        }

                        if (fabs(Ft) > mu * Fn) {
                            Ftx = -mu * Fn * tx;
                            Fty = -mu * Fn * ty;
                            grain[i].addForce(Ftx, Fty);
                            grain[j].addForce(-Ftx, -Fty);
                        } else {
                            Ftx = Ft * tx;
                            Fty = Ft * ty;
                            grain[i].addForce(Ftx, Fty);
                            grain[j].addForce(-Ftx, -Fty);
                        }

                        //torque
                        M = -grain[i].radius() * nx * Fty + grain[i].radius() * ny * Ftx;
                        grain[i].addMomentum(M);
                        M = -grain[j].radius() * nx * Fty + grain[j].radius() * ny * Ftx;
                        grain[j].addMomentum(M);
                    }
                    j = grain[j].linkedDisk();
                }
            }

            //with walls
            for (j = 0; j < 4; j++) {
                nx = wall[j].nx();
                ny = wall[j].ny();
                rx = grain[i].x() - wall[j].x();
                ry = grain[i].y() - wall[j].y();

                delta = (rx * nx + ry * ny) - grain[i].radius();
                if (delta < 0) {
                    vx = grain[i].vx() - wall[j].vx() + grain[i].radius() * grain[i].w() * ny;
                    vy = grain[i].vy() - wall[j].vy() - grain[i].radius() * grain[i].w() * nx;
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
                    effectiveMass = grain[i].mass();
                    eta = -2. * log(e) * sqrt(effectiveMass * kn / (log(e) * log(e) + M_PI * M_PI));
                    Fn = -kn * delta - eta * vn;
                    Ft = -kt * vt;

                    //forces
                    if (Fn > 0) {
                        Fnx = Fn * nx;
                        Fny = Fn * ny;
                        grain[i].addForce(Fnx, Fny);
                    } else {
                        Fn = 0.;
                    }


                    if (fabs(Ft) > mu * Fn) {
                        Ftx = -mu * Fn * tx;
                        Fty = -mu * Fn * ty;
                        grain[i].addForce(Ftx, Fty);
                    } else {
                        Ftx = Ft * tx;
                        Fty = Ft * ty;
                        grain[i].addForce(Ftx, Fty);
                    }

                    //torque
                    M = -grain[i].radius() * nx * Fty + grain[i].radius() * ny * Ftx;
                    grain[i].addMomentum(M);
                }
            }
        }

        //update velocity and position
        for (i = 0; i < numberOfGrains; i++) {
            grain[i].updateVelocity(dt);
            grain[i].updatePosition(dt / 2.);
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
                    grain[i].print((int) ((recTime + dt) * fps));
                }
            }
        }
    }

    delete[] cellule;
    delete[] wall;
    delete[] grain;
    return 0;
}
