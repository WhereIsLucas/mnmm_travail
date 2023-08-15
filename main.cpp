#include <iostream>
#include <random>
#include <map>
#include <omp.h>
#include "Grain.h"
#include "Ball.h"
#include "BallPrinter.h"
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
    //RANDOM INIT
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
    auto ballCollisionsSettings = new CollisionSettings(.9, .6, 100. * 100000., 1000.);
    auto ballGrainCollisionsSettings = new CollisionSettings(.9, .6, 1000000., 10000.);

    // VIDEO OPTIONS
    int fps = 25;
    double tStartCapture = 0.;
    double totalTime = 10;
    int totalFrames = (int) ((totalTime - tStartCapture) * fps);
    double recTime;

    GrainPrinter grainPrinter("data/grain/");
    BallPrinter ballPrinter("data/ball/");

// GRAINS
    double prop = 0.4; //Proportion de remplissage voulue (entre 0 et 1)
    double radius = 0.075; //Radius moyen des grains
    int numberOfGrains = (int) (prop / pow(radius, 2));
    std::ofstream infoFile;
    infoFile.open("data/infos.txt");
    infoFile << prop << "," << radius <<  ',' << numberOfGrains << std::endl;;
    infoFile.close();

    double mass; //Define later
    double rho = 1600.; //masse volumique du sable
    auto *grains = new Grain[numberOfGrains];
    int numberOfPlacedGrains = 0;
    int numberOfOverlaps;

// BARREL
    double ballRadius = 1.;
    double ballRho = 1400.; //masse volumique du PVC
    double ballMass = 2. * M_PI + ballRadius * 0.01 * ballRho;
    Ball ball;

    //this setups the getRadius distribution
    double radiusMean = radius;

    // setup the domain and the base plane

    double xDomain = 1.2 * ballRadius;
    Vector2 a(0, xDomain);
    Vector2 b(xDomain, 0);
    Plan vibratingPlane;

    vibratingPlane.initPlanFromCoordinates(a, b);
    vibratingPlane.printPlanInfos("data/vibratingPlane.txt");

    double yDomain = 5.* ballRadius;

    //Définition du domain et impression pour le code python
    Domain domain(xDomain, yDomain);
    domain.printDomainInfos("data/domain.txt");


    //ON place les grains et le ball
    double xBall = ballRadius;
    Vector2 ballPosition(vibratingPlane.getPointFromX(xBall));
    ball.initBarrel(ballRadius, ballMass,Vector2(0, ballRadius),Vector2(0.));

    while (numberOfPlacedGrains < numberOfGrains) {
        numberOfOverlaps = 0;
        double direction = (double) (uniformRealDistribution(gen)) * 2 * M_PI;
        double randomPositionRadius = (double) sqrt(uniformRealDistribution(gen)) * (ballRadius) * .90;//sqrt pour que ce soit uniforme
        Vector2 grainPosition(randomPositionRadius * cos(direction), randomPositionRadius * sin(direction));
        grainPosition = grainPosition + ball.getPosition();
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
        }
    }
    //Collisions Settings
    std::cout << "Grains and Ball are placed" << std::endl;


    //CLEAR FILES
    for (int l = 0; l <= totalFrames; l++) {
        grainPrinter.clearPrint(l);
        ballPrinter.clearPrint(l);
    }

    //record initial state
    for (i = 0; i < numberOfGrains; i++) {
        grainPrinter.print(grains[i], 0);
    }
    ballPrinter.print(ball, 0);
    std::cout << "Initial state in saved" << std::endl;

    double t = 0.;
    while (true) {
        t += dt;
        if (t > totalTime) {
            break;
        }

        //loop on grains

        for (i = 0; i < numberOfGrains; i++) {
            grains[i].updatePosition(dt / 2.);
            grains[i].resetForce();
            grains[i].addGravityForce(Vector2(0, -9.81));
        }
        ball.updatePosition(dt / 2.);
        ball.resetForce();
        ball.addGravityForce(Vector2(0, -9.81));


        /*** contact detection and forces ***/
        for (i = 0; i < numberOfGrains; i++) {
            for (int j = i + 1; j < numberOfGrains; ++j) {
                computeCollisionWithGrain(&grains[i], &grains[j], grainCollisionsSettings);
            }
            //Collisions with the ball
            computeCollusionBetweenGrainAndBarrel(&grains[i], &ball, ballGrainCollisionsSettings);
        }
        computeCollisionBetweenBarrelAndPlan(&ball, &vibratingPlane, ballCollisionsSettings);

        //update velocity and position
        for (i = 0; i < numberOfGrains; i++) {
            grains[i].updateVelocity(dt);
            grains[i].updatePosition(dt / 2.);
        }
        ball.updateVelocity(dt);
        ball.updatePosition(dt / 2.);
        //record
        recTime = t - tStartCapture;
        if (recTime >= 0.) {
            if ((int) ((recTime + dt) * fps) > (int) (recTime * fps)) {
                for (i = 0; i < numberOfGrains; i++) {
                    grainPrinter.print(grains[i], (int) ((recTime + dt) * fps));
                }
                ballPrinter.print(ball, (int) ((recTime + dt) * fps));
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

