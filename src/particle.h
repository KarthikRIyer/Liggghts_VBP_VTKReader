//
// Created by karthik on 26/10/20.
//

#ifndef LIGGGHTS_VBP_VTKREADER_PARTICLE_H
#define LIGGGHTS_VBP_VTKREADER_PARTICLE_H

#define four_by_three_pi 4.18879020479

class Particle {
public:

    Particle() {}

    Particle(double mass, double x, double y, double z) : mass(mass), x(x), y(y), z(z) {}

    virtual double getTopZ();

    virtual double getBottomZ();

    virtual double volume();

    double mass = 0;
    double vol = 0;
    double x = 0;
    double y = 0;
    double z = 0;
};

#endif //LIGGGHTS_VBP_VTKREADER_PARTICLE_H