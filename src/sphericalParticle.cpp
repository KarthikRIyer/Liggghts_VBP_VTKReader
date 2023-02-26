//
// Created by karthik on 26/10/20.
//
#include "sphericalParticle.h"


double SphericalParticle::volume() const {
    return vol;
}

double SphericalParticle::getTopZ() {
    return z + radius;
}

double SphericalParticle::getBottomZ() {
    return z - radius;
}

double SphericalParticle::getMaxDimLength() {
    return radius * 2.0;
}

Eigen::MatrixXd SphericalParticle::getQMatrix() {
    Eigen::MatrixXd m(3,3);
    m.setZero();
    return m;
}
