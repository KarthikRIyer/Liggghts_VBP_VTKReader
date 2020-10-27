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