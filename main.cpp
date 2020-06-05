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
#include "Collisions.h"
#include <chrono>
#include <fstream>

int main(int argc, char **argv) {
    //RANDOM
    auto t1 = omp_get_wtime();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> uniformRealDistribution(0, 1);
    //Variables
    double dt = 10. * 0.000001;
    int i;
    double m_tot = 0;

//Collisions settings
    auto grainCollisionsSettings = new CollisionSettings(.9, .6, 100000., 1000.);
    auto barrelCollisionsSettings = new CollisionSettings(.9, .6, 100. * 100000., 1000.);
    auto barrelGrainCollisionsSettings = new CollisionSettings(.9, .6, 1000000., 10000.);
    // VIDEO OPTIONS
    int fps = 25;
    double tStartCapture = 0.;
    double totalTime = 10;
    int totalFrames = (int) ((totalTime - tStartCapture) * fps);
    double recTime;

    GrainPrinter grainPrinter("datas/grain/");
    BarrelPrinter barrelPrinter("datas/barrel/");

// GRAINS
    double prop = 0.7; //Proportion de remplissage voulue (entre 0 et 1)
    double radiusRandom = 0.5; //Distribution des radius (entre 0 et 1)
    double radius = 0.075; //Radius moyen des grains
    int numberOfGrains = (int) (prop / pow(radius, 2));
    std::ofstream myfile;
    myfile.open("datas/infos.txt");
    myfile << prop << "," << radius << ',' << radiusRandom << ',' << numberOfGrains << std::endl;;
    myfile.close();
    std::cout << numberOfGrains << std::endl;

    double mass; //Define later
    double rho = 1600.; //masse volumique du sable
    auto *grains = new Grain[numberOfGrains];
    int numberOfPlacedGrains = 0;
    int numberOfOverlaps;

// BARREL
    double barrelRadius = 1.;
    double barrelRho = 1400.; //masse volumique du PVC
    double barrelMass = 2. * M_PI + barrelRadius * 0.01 * barrelRho;
    Barrel barrel;

//this setups the getRadius distribution
    double radiusMean = radius;
    std::uniform_real_distribution<double> radiusDistribution(radiusMean - (radiusRandom * radiusMean),
                                                              radiusMean + (radiusRandom
                                                                            * radiusMean));

//plan et domaine
    int numberOfRevolution = 4;  //Nbr de révoluton du barrel sur le plan

    double alphaDegree = 15; //Angle d'inclinaison du plan
    double alpha = alphaDegree * M_PI / 180;

    double xDomain = 4. * barrelRadius * numberOfRevolution;

    Vector2 a(0, xDomain * tan(alpha));
    Vector2 b(xDomain, 0);

    Plan plan;
    plan.initPlanFromCoordinates(a, b);
    plan.printPlanInfos("datas/plan.txt");

    double yDomain = plan.getPointFromX(0).getY() + 2. * barrelRadius;

    //Plan 2 pour stopper le tonneau
    Vector2 c(xDomain, 0);
    Vector2 d(xDomain + 1, yDomain);
    Plan plan2;
    plan2.initPlanFromCoordinates(c, d);
    plan.printPlanInfos("datas/plan2.txt");

    //Définition du domaine et impression pour le code python
    Domain domain(xDomain, yDomain);
    std::cout << xDomain << " " << yDomain << std::endl;
    domain.printDomainInfos("datas/domain.txt");


//ON place les grains et le barrel
    double xBarrel = barrelRadius - barrelRadius * sin(alpha);
    Vector2 barrelPosition(plan.getPointFromX(xBarrel));
    double deltaBarrel = barrelPosition.getY() - plan.getPointFromX(barrelRadius).getY();
    barrel.initBarrel(barrelRadius, barrelMass,
                      plan.getPointFromX(barrelRadius) + Vector2(0, barrelRadius) + Vector2(0, deltaBarrel),
                      Vector2(0.));

    while (numberOfPlacedGrains < numberOfGrains) {
        radius = fabs(radiusDistribution(gen));
        numberOfOverlaps = 0;
        double direction = (double) (uniformRealDistribution(gen)) * 2 * M_PI;
        double randomRadius =
                (double) sqrt(uniformRealDistribution(gen)) * (barrelRadius) * .90;//sqrt pour que ce soit uniforme
        Vector2 grainPosition(randomRadius * cos(direction), randomRadius * sin(direction));
        grainPosition = grainPosition + barrel.getPosition();
        for (i = 0; i < numberOfPlacedGrains; i++) { //Regarder si il n'y a pas d'overlap
            if (getDistanceBetweenVectors(grainPosition, grains[i].getPosition()) <
                radius + grains[i].getRadius()) {
                numberOfOverlaps++;
                i = numberOfPlacedGrains;
            }
        }
        if (numberOfOverlaps == 0) {
            mass = 4. / 3. * M_PI * pow(radius, 3) * rho;
            m_tot += mass;
            grains[numberOfPlacedGrains].initDisk(numberOfPlacedGrains, radius, mass, grainPosition, Vector2(0.));
            numberOfPlacedGrains++;
            std::cout << numberOfPlacedGrains << std::endl;
        }
    }
    //Collisions Settings
    std::cout << "Grains and Barrel are placed" << std::endl;


//CLEAR FILES
    for (int l = 0; l <= totalFrames; l++) {
        grainPrinter.clearPrint(l);
        barrelPrinter.clearPrint(l);
    }

//record initial state
    for (i = 0; i < numberOfGrains; i++) {
        grainPrinter.print(grains[i], 0);
    }
    barrelPrinter.print(barrel, 0);
    std::cout << "Initial state in saved" << std::endl;

    double t = 0.;
    while (true) {
        t += dt;
        if (t > totalTime) {
            break;
        }
        /*** refresh and update position***/

        //loop on grains

        for (i = 0; i < numberOfGrains; i++) {
            grains[i].updatePosition(dt / 2.);
            grains[i].resetForce();
            grains[i].addGravityForce(Vector2(0, -9.81));
        }
        barrel.updatePosition(dt / 2.);
        barrel.resetForce();
        barrel.addGravityForce(Vector2(0, -9.81));


        /*** contact detection and forces ***/
        for (i = 0; i < numberOfGrains; i++) {
            for (int j = i + 1; j < numberOfGrains; ++j) {
                computeCollisionWithGrain(&grains[i], &grains[j], grainCollisionsSettings);
            }
            //Collisions with the barrel
            computeCollusionBetweenGrainAndBarrel(&grains[i], &barrel, barrelGrainCollisionsSettings);
        }
        computeCollisionBetweenBarrelAndPlan(&barrel, &plan, barrelCollisionsSettings);
        computeCollisionBetweenBarrelAndPlan(&barrel, &plan2, barrelCollisionsSettings);

        //update velocity and position
        for (i = 0; i < numberOfGrains; i++) {
            grains[i].updateVelocity(dt);
            grains[i].updatePosition(dt / 2.);
        }
        barrel.updateVelocity(dt);
        barrel.updatePosition(dt / 2.);
        //record
        recTime = t - tStartCapture;
        if (recTime >= 0.) {
            if ((int) ((recTime + dt) * fps) > (int) (recTime * fps)) {
                for (i = 0; i < numberOfGrains; i++) {
                    grainPrinter.print(grains[i], (int) ((recTime + dt) * fps));
                }
                barrelPrinter.print(barrel, (int) ((recTime + dt) * fps));
                std::cout << "PRINTED IMAGE : " << (int) ((recTime + dt) * fps) << std::endl;
            }
        }
    }

    delete[] grains;

    // Recording end time.
    auto t2 = omp_get_wtime();

    printf("Work took %f seconds", t2 - t1);
    return 0;

}

