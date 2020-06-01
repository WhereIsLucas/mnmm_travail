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

double CohesionSettings::getForceValue() const {
    return forceValue;
}

void CohesionSettings::setForceValue(double forceValue) {
    CohesionSettings::forceValue = forceValue;
}

CohesionSettings::CohesionSettings(double distanceThreshold, double forceValue) : distanceThreshold(distanceThreshold),
                                                                                  forceValue(forceValue) {
    CohesionSettings::setDistanceThreshold(distanceThreshold);
    CohesionSettings::setForceValue(forceValue);
}
