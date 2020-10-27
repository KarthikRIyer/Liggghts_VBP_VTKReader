//
// Created by karthik on 26/10/20.
//
#include "superquadricParticle.h"

double SuperquadricParticle::volume() {
    return vol;
}

double SuperquadricParticle::getTopZ() {
    return z + zDiff;
}

double SuperquadricParticle::getBottomZ() {
    return z - zDiff;
}
