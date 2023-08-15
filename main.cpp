#include <iostream>
#include <random>
#include <map>
#include <omp.h>
#include "Grain.h"
#include "Ball.h"
#include "BallPrinter.h"
#include "Plane.h"
#include "Cell.h"
#include "GrainPrinter.h"
#include "Vector2.h"
#include "Domain.h"
#include "Collisions.h"
#include "PlanePrinter.h"
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

    GrainPrinter grainPrinter("data/grains/");
    BallPrinter ballPrinter("data/ball/");
    PlanePrinter planePrinter("data/plane/");

// GRAINS
    double prop = 0.; //Proportion de remplissage voulue (entre 0 et 1)
    double radius = 0.002; //Radius moyen des grains
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

    // Ball
    double ballRadius = 0.0375;
    double ballMass = 0.0019;
    Ball ball;

    //this setups the getRadius distribution
    double radiusMean = radius;

    // setup the domain and the vibrating plane
    double xDomain = 1.2 * ballRadius;
    Plane vibratingPlane;
    vibratingPlane.initPlanFromCoordinates(Vector2(-2., 0), Vector2(2., 0));
    double frequency = 6.;
    double amplitude = 0.0095;

    double yDomain = 5.* ballRadius;

    //Définition du domain et impression pour le code python
    Domain domain(xDomain, yDomain);
    domain.printDomainInfos("data/domain.txt");


    //ON place les grains et le ball
    ball.initBall(ballRadius, ballMass, Vector2(0, ballRadius), Vector2(0.));

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
        planePrinter.clearPrint(l);
    }

    //record initial state
    for (i = 0; i < numberOfGrains; i++) {
        grainPrinter.print(grains[i], 0);
    }
    ballPrinter.print(ball, 0);
    planePrinter.print(vibratingPlane, 0);
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

        //we update the position of the vibrating plane. It is a sin wave based on the paremeters above
        vibratingPlane.setPosition(Vector2(0., amplitude * sin(2. * M_PI * frequency * t)));


        /*** contact detection and forces ***/
        for (i = 0; i < numberOfGrains; i++) {
            for (int j = i + 1; j < numberOfGrains; ++j) {
                computeCollisionWithGrain(&grains[i], &grains[j], grainCollisionsSettings);
            }
            //Collisions with the ball
            computeCollusionBetweenGrainAndBall(&grains[i], &ball, ballGrainCollisionsSettings);
        }
        computeCollisionBetweenBallAndPlane(&ball, &vibratingPlane, ballCollisionsSettings);

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
                planePrinter.print(vibratingPlane, (int) ((recTime + dt) * fps));
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

