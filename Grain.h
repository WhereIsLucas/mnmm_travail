#ifndef TRAVAIL2_GRAIN_H
#define TRAVAIL2_GRAIN_H

class Grain {

private:
    int m_index;

    int m_linkedCell;
    int m_linkedDisk;

    double m_radius, m_mass, m_inertia;
    double m_x, m_y;
    double m_vx, m_vy, m_v;
    double m_ax, m_ay;
    double m_theta;
    double m_w;
    double m_alpha;
    double m_Fx, m_Fy;
    double m_M;

public:
    Grain();

    ~Grain();

    void initDisk(int, double, double, double, double, double, double);

    void updateVelocity(double);

    void updatePosition(double);

    void resetForce();

    void addForce(double, double);

    void addMomentum(double);

    void addGravityForce(double, double);

    void setLinkedDisk(int);

    void setLinkedCell(int);

    void print(int);

    int linkedCell();

    int linkedDisk();

    double radius();

    double x();

    double y();

    double vx();

    double vy();

    double v();

    double w();


    double mass();

    int index();

    double theta();

};

#endif