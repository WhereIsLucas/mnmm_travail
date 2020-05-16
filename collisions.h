//
// Created by lucas on 16/05/2020.
//

#ifndef TRAVAIL2_COLLISIONS_H
#define TRAVAIL2_COLLISIONS_H

#include "Grain.h"
#include "Container.h"
#include "CollisionSettings.h"

void computeCollision(Grain *pGrain1, Grain *pGrain2);

void computeCollisionWithContainer(Grain *pGrain1, Container *container, CollisionSettings *collisionSettings);


#endif //TRAVAIL2_COLLISIONS_H
