//
// Created by karthik on 26/10/20.
//
#include <cmath>
#include "particle.h"

#ifndef LIGGGHTS_VBP_VTKREADER_SUPERQUADRICPARTICLE_H
#define LIGGGHTS_VBP_VTKREADER_SUPERQUADRICPARTICLE_H

class SuperquadricParticle : public Particle {
public:
    SuperquadricParticle() : Particle() {}

    SuperquadricParticle(double mass, double shapex, double shapey, double shapez, double x, double y, double z,
                         double quat1, double quat2, double quat3, double quat4, double vx = 0, double vy = 0,
                         double vz = 0)
            : Particle(mass, x, y, z, vx, vy, vz), shapex(shapex),
              shapey(shapey), shapez(shapez), quat1(quat1), quat2(quat2), quat3(quat3), quat4(quat4) {
        vol = shapex * shapey * shapez * four_by_three_pi;

        /*Quaternion to rotation matrix from: https://link.springer.com/article/10.1007/s40571-016-0131-6*/
        double a1 = 1 - 2.0 * (quat3 * quat3 + quat4 * quat4);
        double b1 = 2.0 * (quat2 * quat3 - quat1 * quat4);
        double c1 = 2.0 * (quat2 * quat4 - quat1 * quat3);

        double a2 = 2.0 * (quat2 * quat3 + quat1 * quat4);
        double b2 = 1 - 2.0 * (quat2 * quat2 + quat4 * quat4);
        double c2 = 2.0 * (quat3 * quat4 - quat1 * quat2);

        double a3 = 2.0 * (quat2 * quat4 - quat1 * quat3);
        double b3 = 2.0 * (quat3 * quat4 + quat1 * quat2);
        double c3 = 1 - 2.0 * (quat2 * quat2 + quat3 * quat3);

        this->q00 = (c1 * c1 - 1.0 / 3.0) * (3.0 / 2.0);
        this->q01 = (c1 * c2) * (3.0 / 2.0);
        this->q02 = (c1 * c3) * (3.0 / 2.0);

        this->q10 = (c2 * c1) * (3.0 / 2.0);
        this->q11 = (c2 * c2 - 1.0 / 3.0) * (3.0 / 2.0);
        this->q12 = (c2 * c3) * (3.0 / 2.0);

        this->q20 = (c3 * c1) * (3.0 / 2.0);
        this->q21 = (c3 * c2) * (3.0 / 2.0);
        this->q22 = (c3 * c3 - 1.0 / 3.0) * (3.0 / 2.0);

        double D = a1 * (b2 * c3 - c2 * b3) - b1 * (a2 * c3 - a3 * c2) + c1 * (a2 * b2 - a3 * b2);
        double one_by_D_squared = 1.0 / (D * D);
        double one_by_a_squared = 1.0 / (shapex * shapex);
        double one_by_b_squared = 1.0 / (shapey * shapey);
        double one_by_c_squared = 1.0 / (shapez * shapez);

        double p = one_by_D_squared * (((b2 * c3 - c2 * b3) * (b2 * c3 - c2 * b3) * one_by_a_squared) +
                                       ((a3 * c2 - a2 * c3) * (a3 * c2 - a2 * c3) * one_by_b_squared) +
                                       ((a2 * b3 - a3 * b2) * (a2 * b3 - a3 * b2) * one_by_c_squared));
        double q = one_by_D_squared * (((b3 * c1 - b2 * c3) * (b3 * c1 - b2 * c3) * one_by_a_squared) +
                                       ((a1 * c3 - a3 * c1) * (a1 * c3 - a3 * c1) * one_by_b_squared) +
                                       ((a3 * b1 - a1 * b3) * (a3 * b1 - a1 * b3) * one_by_c_squared));
        double r = one_by_D_squared * (((b1 * c2 - b2 * c1) * (b1 * c2 - b2 * c1) * one_by_a_squared) +
                                       ((a2 * c1 - a1 * c2) * (a2 * c1 - a1 * c2) * one_by_b_squared) +
                                       ((a1 * b2 - a2 * b1) * (a1 * b2 - a2 * b1) * one_by_c_squared));
        double h = one_by_D_squared * (((b2 * c3 - b3 * c2) * (b3 * c1 - b1 * c3) * one_by_a_squared) +
                                       ((a3 * c2 - a2 * c3) * (a1 * c3 - a3 * c1) * one_by_b_squared) +
                                       ((a2 * b3 - a3 * b2) * (a3 * b1 - a1 * b3) * one_by_c_squared));
        double g = one_by_D_squared * (((b1 * c2 - b2 * c1) * (b2 * c3 - b3 * c2) * one_by_a_squared) +
                                       ((a2 * c1 - a1 * c2) * (a3 * c2 - a2 * c3) * one_by_b_squared) +
                                       ((a1 * b2 - a2 * b1) * (a2 * b3 - a3 * b2) * one_by_c_squared));
        double f = one_by_D_squared * (((b3 * c1 - b1 * c3) * (b1 * c2 - b2 * c1) * one_by_a_squared) +
                                       ((a1 * c3 - a3 * c1) * (a2 * c1 - a1 * c2) * one_by_b_squared) +
                                       ((a3 * b1 - a1 * b3) * (a1 * b2 - a2 * b1) * one_by_c_squared));
        double A = (p * (q * r - f * f)) - (h * (h * r - f * g)) + (g * (h * f - g * q));
        /*Max Z from : https://math.stackexchange.com/questions/3610652/find-x-y-z-maximum-and-minimum-points-of-ellipsoid*/
        zDiff = std::sqrt((p * q - h * h) / A);
    }

    double volume() override;

    double getTopZ() override;

    double getBottomZ() override;

    double getMaxDimLength() override;

    double shapex = 0;
    double shapey = 0;
    double shapez = 0;
    double quat1 = 0;
    double quat2 = 0;
    double quat3 = 0;
    double quat4 = 1;
    double blockiness1 = 2;
    double blockiness2 = 2;
    double zDiff = 0;

    double q00 =0; //Qij matrix to calculate order parameter S
    double q01 =0;
    double q02 =0;
    double q10 =0;
    double q11 =0;
    double q12 =0;
    double q20 =0;
    double q21 =0;
    double q22 =0;
};

#endif //LIGGGHTS_VBP_VTKREADER_SUPERQUADRICPARTICLE_H
