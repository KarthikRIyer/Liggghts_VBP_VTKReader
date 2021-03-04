//
// Created by karthik on 26/02/21.
//
#pragma once

#include <vector>
#include <memory>
#include "particle.h"
#include <glm/glm.hpp>
#include <functional>

class ParticleSearcher {
public:
    ParticleSearcher(std::vector<std::shared_ptr<Particle>> &particles, double rightX, double leftX, double frontY,
                     double backY, double topZ, double bottomZ);

    void forEachNearbyPoint(glm::dvec3 &origin, double radius,
                            const std::function<void(size_t, std::shared_ptr<Particle>)> &callback);

private:
    double rightX;
    double leftX;
    double frontY;
    double backY;
    double topZ;
    double bottomZ;
    double gridSpacing;
    std::vector<std::shared_ptr<Particle>> _particles;
    glm::ivec3 _resolution;
    std::vector<std::vector<size_t>> _buckets;

    size_t getHashKeyFromPosition(glm::dvec3 position);

    void getNearbyKeys(glm::dvec3 &position, size_t *nearbyKeys);

    inline int getBin(double lo, double val, double increment) {
        return std::floor((val - lo) / increment);
    }

    template<typename T>
    inline T clamp(const T &n, const T &lower, const T &upper) {
        return n <= lower ? lower : n >= upper ? upper : n;
    }
};