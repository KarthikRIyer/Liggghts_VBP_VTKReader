//
// Created by karthik on 26/02/21.
//
#include <iostream>
#include "particleSearcher.h"

ParticleSearcher::ParticleSearcher(std::vector<std::shared_ptr<Particle>> &particles, double rightX, double leftX,
                                   double frontY, double backY, double topZ, double bottomZ) : rightX(rightX),
                                                                                               leftX(leftX),
                                                                                               frontY(frontY),
                                                                                               backY(backY), topZ(topZ),
                                                                                               bottomZ(bottomZ) {
    _buckets.clear();
    _particles.clear();

    gridSpacing = 4 * particles[0]->getMaxDimLength();
    _resolution.x = std::floor(std::abs(rightX - leftX) / gridSpacing);
    _resolution.y = std::floor(std::abs(frontY - backY) / gridSpacing);
    _resolution.z = std::floor(std::abs(topZ - bottomZ) / gridSpacing);

    //allocate memory
    _buckets.resize(_resolution.x * _resolution.y * _resolution.z);
    _particles.resize(particles.size());

    //put points into buckets
    for (size_t i = 0; i < particles.size(); ++i) {
        _particles[i] = particles[i];
        size_t key = getHashKeyFromPosition(glm::dvec3(particles[i]->x, particles[i]->y, particles[i]->z));
        _buckets[key].push_back(i);
    }
}

void ParticleSearcher::forEachNearbyPoint(glm::dvec3 &origin, double radius,
                                          const std::function<void(size_t, std::shared_ptr<Particle>)> &callback) {

    if (_buckets.empty())return;

    size_t nearbyKeys[8];
    getNearbyKeys(origin, nearbyKeys);

    const double queryRadiusSquared = radius * radius;

    for (int i = 0; i < 8; ++i) {
        const auto &bucket = _buckets[nearbyKeys[i]];
        size_t numberOfPointsInBucket = bucket.size();
        for (int j = 0; j < numberOfPointsInBucket; ++j) {
            size_t pointIndex = bucket[j];
            double rSquared = glm::length(
                    glm::dvec3(_particles[pointIndex]->x, _particles[pointIndex]->y, _particles[pointIndex]->z) -
                    origin);
            rSquared = rSquared * rSquared;
            if (rSquared <= queryRadiusSquared) {
                callback(pointIndex, _particles[pointIndex]);
            }
        }
    }
}

size_t ParticleSearcher::getHashKeyFromPosition(glm::dvec3 position) {

//    int binI = getBin(bottomZ, particle->z, zIncrement);
//    int bin = clamp(binI, 0, sliceCount - 1);

    glm::ivec3 bucketIndex;
    bucketIndex.x = clamp(getBin(std::min(leftX, rightX), position.x, gridSpacing), 0, _resolution.x - 1);
    bucketIndex.y = clamp(getBin(std::min(frontY, backY), position.y, gridSpacing), 0, _resolution.y - 1);
    bucketIndex.z = clamp(getBin(std::min(topZ, bottomZ), position.z, gridSpacing), 0, _resolution.z - 1);

    return (bucketIndex.z * _resolution.y + bucketIndex.y) * _resolution.x + bucketIndex.x;
}

void ParticleSearcher::getNearbyKeys(glm::dvec3 &position, size_t *nearbyKeys) {
    glm::ivec3 originIndex;
    originIndex.x = clamp(getBin(std::min(leftX, rightX), position.x, gridSpacing), 0, _resolution.x - 1);
    originIndex.y = clamp(getBin(std::min(frontY, backY), position.y, gridSpacing), 0, _resolution.y - 1);
    originIndex.z = clamp(getBin(std::min(topZ, bottomZ), position.z, gridSpacing), 0, _resolution.z - 1);

    glm::ivec3 nearbyBucketIndices[8];

    for (int i = 0; i < 8; ++i) {
        nearbyBucketIndices[i] = originIndex;
    }

    if ((originIndex.x + 0.5f) * gridSpacing <= position.x) {
        nearbyBucketIndices[4].x += 1;
        nearbyBucketIndices[5].x += 1;
        nearbyBucketIndices[6].x += 1;
        nearbyBucketIndices[7].x += 1;
    } else {
        nearbyBucketIndices[4].x -= 1;
        nearbyBucketIndices[5].x -= 1;
        nearbyBucketIndices[6].x -= 1;
        nearbyBucketIndices[7].x -= 1;
    }

    if ((originIndex.y + 0.5f) * gridSpacing <= position.x) {
        nearbyBucketIndices[2].y += 1;
        nearbyBucketIndices[3].y += 1;
        nearbyBucketIndices[6].y += 1;
        nearbyBucketIndices[7].y += 1;
    } else {
        nearbyBucketIndices[2].y -= 1;
        nearbyBucketIndices[3].y -= 1;
        nearbyBucketIndices[6].y -= 1;
        nearbyBucketIndices[7].y -= 1;
    }

    if ((originIndex.z + 0.5f) * gridSpacing <= position.z) {
        nearbyBucketIndices[1].z += 1;
        nearbyBucketIndices[3].z += 1;
        nearbyBucketIndices[5].z += 1;
        nearbyBucketIndices[7].z += 1;
    } else {
        nearbyBucketIndices[1].z -= 1;
        nearbyBucketIndices[3].z -= 1;
        nearbyBucketIndices[5].z -= 1;
        nearbyBucketIndices[7].z -= 1;
    }

    for (int i = 0; i < 8; ++i) {
        nearbyBucketIndices[i].x = clamp(nearbyBucketIndices[i].x, 0, _resolution.x - 1);
        nearbyBucketIndices[i].y = clamp(nearbyBucketIndices[i].y, 0, _resolution.y - 1);
        nearbyBucketIndices[i].z = clamp(nearbyBucketIndices[i].z, 0, _resolution.z - 1);
        nearbyKeys[i] = (nearbyBucketIndices[i].z * _resolution.y + nearbyBucketIndices[i].y) * _resolution.x +
                        nearbyBucketIndices[i].x;
    }
}
