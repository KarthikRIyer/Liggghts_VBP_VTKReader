//
// Created by karthik on 27/02/21.
//
#include "sphStdKernel.h"

SphStdKernel::SphStdKernel() : h(0), h2(0), h3(0), h4(0), h5(0) {

}

SphStdKernel::SphStdKernel(double kernelRadius) : h(kernelRadius), h2(h * h), h3(h2 * h), h4(h2 * h2),
                                                         h5(h3 * h2) {

}

double SphStdKernel::operator()(double distance) const {
    if (distance * distance >= h2)
        return 0;
    else {
        double x = 1.0 - distance * distance / h2;
        return 315.0 / (64 * kPiD * h3) * x * x * x;
    }
}

double SphStdKernel::firstDerivative(double distance) const {
    if (distance * distance >= h2)
        return 0;
    else {
        double x = 1.0 - distance * distance / h2;
        return -945.0 / (32 * kPiD * h5) * distance * x * x * x;
    }
}

glm::dvec3 SphStdKernel::gradient(double distance, const glm::dvec3 &direction) const {
    return -firstDerivative(distance) * direction;
}
