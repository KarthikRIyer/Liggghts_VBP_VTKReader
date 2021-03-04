//
// Created by karthik on 27/02/21.
//
#include "sphSpikyKernel.h"

SphSpikyKernel::SphSpikyKernel() : h(0), h2(0), h3(0), h4(0), h5(0) {}

SphSpikyKernel::SphSpikyKernel(double kernelRadius) : h(kernelRadius), h2(h * h), h3(h2 * h), h4(h2 * h2), h5(h3 * h2) {
}

double SphSpikyKernel::operator()(double distance) const {
    if (distance >= h)
        return 0;
    else {
        double x = 1.0 - distance / h;
        return (15.0 / (kPiD * h3)) * x * x * x;
    }
}

double SphSpikyKernel::firstDerivative(double distance) const {
    if (distance >= h)
        return 0;
    else {
        double x = 1.0 - distance / h;
        return -45.0 / (kPiD * h4) * x * x;
    }
}

glm::dvec3 SphSpikyKernel::gradient(double distance, const glm::dvec3 &direction) const {
    return -firstDerivative(distance) * direction;
}
