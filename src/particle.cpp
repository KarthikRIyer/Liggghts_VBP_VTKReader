//
// Created by karthik on 26/10/20.
//

#include "particle.h"

double Particle::getTopZ() {
    return z;
}

double Particle::getBottomZ() {
    return z;
}

double Particle::volume() {
    return vol;
}

double Particle::getMaxDimLength() {
    return 0;
}

Eigen::MatrixXd Particle::getQMatrix() {
    Eigen::MatrixXd m(3,3);
    m.setZero();
    return m;
}
