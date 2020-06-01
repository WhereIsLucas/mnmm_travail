//
// Created by lucas on 01/06/2020.
//

#ifndef TRAVAIL2_COHESIONSETTINGS_H
#define TRAVAIL2_COHESIONSETTINGS_H


class CohesionSettings {
private:
    double distanceThreshold;
    double forceValue;
public:
    CohesionSettings(double distanceThreshold, double forceValue);

public:
    double getDistanceThreshold() const;

    void setDistanceThreshold(double distanceThreshold);

    double getForceValue() const;

    void setForceValue(double forceValue);
};


#endif //TRAVAIL2_COHESIONSETTINGS_H
