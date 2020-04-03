class Plan
{
    
private:
    
    double m_t;
    double m_x0, m_y0;
    double m_x, m_y;
    double m_vx, m_vy;
    double m_nx, m_ny;
    double m_amplx, m_amply;
    double m_freqx, m_freqy;
    
public:
    Plan();
    ~Plan();
    
    void initPlan(double,double,double,double,double,double,double,double);
    void updateVelocity(double);
    void updatePosition(double);
    double x();
    double y();
    double vx();
    double vy();
    double nx();
    double ny();
};
