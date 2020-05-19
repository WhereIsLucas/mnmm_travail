#ifndef TRAVAIL2_COLLISIONSETTINGS_H
#define TRAVAIL2_COLLISIONSETTINGS_H


class CollisionSettings {
private:
    double e = 0.99;
    double mu = 0.1;
    double kn = 200.;
    double kt = 100000.;
public:
    CollisionSettings(double e, double mu, double kn, double kt);

    double getE() const;

    void setE(double e);

    double getMu() const;

    void setMu(double mu);

    double getKn() const;

    void setKn(double kn);

    double getKt() const;

    void setKt(double kt);

    double getEta(double effectiveMass);
};


#endif //TRAVAIL2_COLLISIONSETTINGS_H
