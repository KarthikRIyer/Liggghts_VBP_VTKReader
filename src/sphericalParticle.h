//
// Created by karthik on 26/10/20.
//
#include "particle.h"

#ifndef LIGGGHTS_VBP_VTKREADER_SPHERICALPARTICLE_H
#define LIGGGHTS_VBP_VTKREADER_SPHERICALPARTICLE_H

class SphericalParticle : public Particle {
public:
    SphericalParticle() : Particle() {}

    SphericalParticle(double mass, double radius, double x, double y, double z, double vx = 0, double vy = 0,
                      double vz = 0) : radius(radius),
                                       Particle(mass, x, y, z, vx, vy, vz) {
        vol = radius * radius * radius * four_by_three_pi;
    }

    virtual double volume() const;

    double radius = 0;

    double getTopZ() override;

    double getBottomZ() override;

    double getMaxDimLength() override;

    Eigen::MatrixXd getQMatrix() override;
};

#endif //LIGGGHTS_VBP_VTKREADER_SPHERICALPARTICLE_H
