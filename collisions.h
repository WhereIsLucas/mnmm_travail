//
// Created by lucas on 16/05/2020.
//

#ifndef MNMM_COLLISIONS_H
#define MNMM_COLLISIONS_H

#include "Grain.h"
#include "Container.h"
#include "CollisionSettings.h"

void computeCollisionWithGrain(Grain *pGrain1, Grain *pGrain2, CollisionSettings *collisionSettings);

void computeCollisionWithContainer(Grain *pGrain1, Container *container, CollisionSettings *collisionSettings);


#endif //MNMM_COLLISIONS_H
