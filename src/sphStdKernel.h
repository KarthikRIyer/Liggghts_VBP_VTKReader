//
// Created by karthik on 27/02/21.
//
#pragma once

#include <glm/glm.hpp>

struct SphStdKernel {
    static constexpr double kPiD = 3.14159265358979323846264338327950288;

    double h, h2, h3, h4, h5;

    SphStdKernel();

    explicit SphStdKernel(double kernelRadius);

    double operator()(double distance) const;

    double firstDerivative(double distance) const;

    glm::dvec3 gradient(double distance, const glm::dvec3 &direction) const;
};
