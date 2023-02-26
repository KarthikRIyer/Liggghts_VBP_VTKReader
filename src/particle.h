//
// Created by karthik on 26/10/20.
//

#ifndef LIGGGHTS_VBP_VTKREADER_PARTICLE_H
#define LIGGGHTS_VBP_VTKREADER_PARTICLE_H

#define four_by_three_pi 4.18879020479

#include "Eigen/Dense"

class Particle {
public:

    Particle() {}

    Particle(double mass, double x, double y, double z, double vx=0, double vy=0, double vz=0) : mass(mass), x(x), y(y), z(z),
                                                                                           vx(vx), vy(vy), vz(vz) {}

    virtual double getTopZ();

    virtual double getBottomZ();

    virtual double volume();

    virtual double getMaxDimLength();

    virtual Eigen::MatrixXd getQMatrix();

    double mass = 0;
    double vol = 0;
    double x = 0;
    double y = 0;
    double z = 0;
    double vx = 0;
    double vy = 0;
    double vz = 0;
    int type = -1;
};

#endif //LIGGGHTS_VBP_VTKREADER_PARTICLE_H
