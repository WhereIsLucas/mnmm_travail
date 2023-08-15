//
// Created by lucas on 16/05/2020.
//

#ifndef MN_BORIS_COLLISIONS_H
#define MN_BORIS_COLLISIONS_H

#include "Grain.h"
#include "Container.h"
#include "CollisionSettings.h"
#include "Plan.h"
#include "Ball.h"

void computeCollisionWithGrain(Grain *pGrain1, Grain *pGrain2, CollisionSettings *collisionSettings);

void computeCollisionWithContainer(Grain *pGrain1, Container *container, CollisionSettings *collisionSettings);

void computeCollusionBetweenGrainAndBarrel(Grain *pGrain1, Ball *pBarrel, CollisionSettings *collisionSettings);

void computeCollisionBetweenBarrelAndPlan(Ball *pBarrel, Plan *plan, CollisionSettings *collisionSettings);


#endif //MN_BORIS_COLLISIONS_H