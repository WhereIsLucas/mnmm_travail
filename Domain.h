//
// Created by hugo on 14/05/2020.
//

#ifndef MNMM_DOMAIN_H
#define MNMM_DOMAIN_H


class Domain {

private:
    double x;
    double y;

public:
    double getX() const;

    void setX(double x);

    double getY() const;

    void setY(double y);

    Domain(double x, double y);

    void printDomainInfos(std::string fileName);
};


#endif //MNMM_DOMAIN_H
