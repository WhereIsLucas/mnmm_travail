//
// Created by lucas on 01/06/2020.
//

#include "cohesion.h"
#include "CohesionSettings.h"

void  addCohesionForce(Grain *pGrain1, Grain *pGrain2, CohesionSettings cohesionSettings){
    double delta = getDistanceBetweenGrains(*pGrain1, *pGrain2);
    if (delta < cohesionSettings.getDistanceThreshold()) {
        Vector2 normalVector = (pGrain1->getPosition() - pGrain2->getPosition()).normalize();
        double forceValue = cohesionSettings.getForceValue() - cohesionSettings.getForceValue()*delta/cohesionSettings.getDistanceThreshold();
        pGrain1->addForce(-1.*forceValue * normalVector);
        pGrain2->addForce(forceValue * normalVector);
    }
}
