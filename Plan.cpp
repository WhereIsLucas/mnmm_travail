#include <math.h>
#include "Plan.h"

Plan::Plan()
{

}

Plan::~Plan()
{
    
}

void Plan::initPlan(double i_x, double i_y, double i_nx, double i_ny, double i_amplx, double i_amply, double i_freqx, double i_freqy)
{
    double omegax=2.*M_PI*i_freqx;
    double omegay=2.*M_PI*i_freqy;
    double norm=sqrt(i_nx*i_nx+i_ny*i_ny);

    m_x0 = i_x;
    m_y0 = i_y;
    m_x = m_x0;
    m_y = m_y0;
    m_nx = i_nx/norm;
    m_ny = i_ny/norm;
    m_amplx = i_amplx;
    m_amply = i_amply;
    m_freqx = i_freqx;
    m_freqy = i_freqy;
    m_vx = m_amplx*omegax;
    m_vy = m_amply*omegay;
}

void Plan::updatePosition(double i_t)
{
    m_x = m_x0+m_amplx*sin(2.*M_PI*m_freqx*i_t);
    m_y = m_y0+m_amply*sin(2.*M_PI*m_freqy*i_t);
}

void Plan::updateVelocity(double i_t)
{
    double omegax=2.*M_PI*m_freqx;
    double omegay=2.*M_PI*m_freqy;
    m_vx = m_amplx*omegax*cos(omegax*i_t);
    m_vy = m_amply*omegay*cos(omegay*i_t);
}

double Plan::x()
{
    return m_x;
}

double Plan::y()
{
    return m_y;
}

double Plan::vx()
{
    return m_vx;
}

double Plan::vy()
{
    return m_vy;
}

double Plan::nx()
{
    return m_nx;
}

double Plan::ny()
{
    return m_ny;
}
