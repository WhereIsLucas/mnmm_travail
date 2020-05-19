#include <cmath>
#include <iostream>
#include "CollisionSettings.h"

double CollisionSettings::getE() const {
    return e;
}

void CollisionSettings::setE(double e) {
    CollisionSettings::e = e;
}

double CollisionSettings::getMu() const {
    return mu;
}

void CollisionSettings::setMu(double mu) {
    CollisionSettings::mu = mu;
}

double CollisionSettings::getKn() const {
    return kn;
}

void CollisionSettings::setKn(double kn) {
    CollisionSettings::kn = kn;
}

double CollisionSettings::getKt() const {
    return kt;
}

void CollisionSettings::setKt(double kt) {
    CollisionSettings::kt = kt;
}

CollisionSettings::CollisionSettings(double e, double mu, double kn, double kt) : e(e), mu(mu), kn(kn), kt(kt) {
    CollisionSettings::setE(e);
    CollisionSettings::setMu(mu);
    CollisionSettings::setKn(kn);
    CollisionSettings::setKt(kt);
}

double CollisionSettings::getEta(double effectiveMass) {
    return -2. * log(e) * sqrt(effectiveMass * kn / (pow(log(e), 2) + pow(M_PI, 2)));
}
