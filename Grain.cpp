#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "Grain.h"

using namespace std;

Grain::Grain()
{

}

Grain::~Grain()
{
    
}

void Grain::initDisk(int i_index, double i_radius, double i_mass, double i_x, double i_y, double i_vx, double i_vy)
{
    m_index = i_index;
    m_linkedDisk = -9;
    m_linkedCell = -9;
    m_radius = i_radius;
    m_mass = i_mass;
    m_inertia = 0.5*m_mass*m_radius*m_radius;
    
    m_x = i_x;
    m_y = i_y;
    m_vx = i_vx;
    m_vy = i_vy;
    m_v = sqrt(m_vx*m_vx+m_vy*m_vy);
    m_ax = 0.;
    m_ay = 0.;
    
    m_theta = 0.;
    m_w = 0.;
    m_alpha = 0.;
    
    m_Fx = 0.;
    m_Fy = 0.;
    m_M = 0.;
}

void Grain::resetForce()
{
    m_Fx = 0.;
    m_Fy = 0.;
    m_M = 0.;
}

void Grain::addForce(double i_Fx, double i_Fy)
{
    m_Fx += i_Fx;
    m_Fy += i_Fy;
}

void Grain::addMomentum(double i_M)
{
    m_M += i_M;
}

void Grain::addGravityForce(double i_gx, double i_gy)
{
    m_Fx += m_mass*i_gx;
    m_Fy += m_mass*i_gy;
}

void Grain::updatePosition(double dt)
{
    m_x += m_vx*dt;
    m_y += m_vy*dt;
    
    m_theta += m_w*dt;
}

void Grain::updateVelocity(double dt)
{
    m_ax = m_Fx/m_mass;
    m_ay = m_Fy/m_mass;
    m_vx += m_ax*dt;
    m_vy += m_ay*dt;
    m_v = sqrt(m_vx*m_vx+m_vy*m_vy);
    
    m_alpha = m_M/m_inertia;
    m_w += m_alpha*dt;
}

void Grain::setLinkedDisk(int i_linkedDisk)
{
    m_linkedDisk = i_linkedDisk;
}

void Grain::setLinkedCell(int i_linkedCell)
{
    m_linkedCell = i_linkedCell;
}

void Grain::print(int i_num)
{
    std::string fileName = "grain" + std::to_string(i_num) + ".txt";
    ofstream myfile;
    myfile.open(fileName.c_str(),ios::app);
    myfile.precision(10);
    myfile<<m_index<<","<<m_x<<","<<m_y<<","<<m_vx<<","<<m_vy<<","<<m_theta<<","<<m_radius<<endl;
    myfile.close();
}

int Grain::linkedDisk()
{
    return m_linkedDisk;
}

int Grain::linkedCell()
{
    return m_linkedCell;
}

double Grain::radius()
{
    return m_radius;
}

double Grain::mass()
{
    return m_mass;
}

double Grain::x()
{
    return m_x;
}

double Grain::y()
{
    return m_y;
}

double Grain::vx()
{
    return m_vx;
}

double Grain::vy()
{
    return m_vy;
}

double Grain::v()
{
    return m_v;
}

double Grain::w()
{
    return m_w;
}

int Grain::index() {
    return m_index;
}

double Grain::theta() {
    return m_theta;
}
