//
// Created by lucas on 01/06/2020.
//

#include "CohesionSettings.h"

double CohesionSettings::getDistanceThreshold() const {
    return distanceThreshold;
}

void CohesionSettings::setDistanceThreshold(double distanceThreshold) {
    CohesionSettings::distanceThreshold = distanceThreshold;
}

double CohesionSettings::getCohesionConstant() const {
    return cohesionConstant;
}

void CohesionSettings::setCohesionConstant(double forceValue) {
    CohesionSettings::cohesionConstant = forceValue;
}

CohesionSettings::CohesionSettings(double distanceThreshold, double forceValue) : distanceThreshold(distanceThreshold),
                                                                                  cohesionConstant(forceValue) {
    CohesionSettings::setDistanceThreshold(distanceThreshold);
    CohesionSettings::setCohesionConstant(forceValue);
}
