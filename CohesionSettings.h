//
// Created by lucas on 01/06/2020.
//

#ifndef TRAVAIL2_COHESIONSETTINGS_H
#define TRAVAIL2_COHESIONSETTINGS_H


class CohesionSettings {
private:
    double distanceThreshold;
    double cohesionConstant;
public:
    CohesionSettings(double distanceThreshold, double forceValue);

public:
    double getDistanceThreshold() const;

    void setDistanceThreshold(double distanceThreshold);

    double getCohesionConstant() const;

    void setCohesionConstant(double forceValue);
};


#endif //TRAVAIL2_COHESIONSETTINGS_H
