//
// Created by karthik on 23/11/20.
//

#ifndef LIGGGHTS_VBP_VTKREADER_BEDLIMS_H
#define LIGGGHTS_VBP_VTKREADER_BEDLIMS_H

#include <vtkGenericDataObjectReader.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <memory>

#include "particle.h"
#include "superquadricParticle.h"
#include "sphericalParticle.h"

std::vector<std::pair<double, double>> getBedLims(std::string filePath, bool isSuperQ) {
    vtkSmartPointer<vtkGenericDataObjectReader> reader =
            vtkSmartPointer<vtkGenericDataObjectReader>::New();
    reader->SetFileName(filePath.c_str());
    reader->Update();
    if (reader->IsFilePolyData()) {
        vtkPolyData *output = reader->GetPolyDataOutput();
        vtkPointData *pd = output->GetPointData();

        vtkDoubleArray *radiusData;
        vtkDoubleArray *shapeXData;
        vtkDoubleArray *shapeYData;
        vtkDoubleArray *shapeZData;
        vtkDoubleArray *quat1Data;
        vtkDoubleArray *quat2Data;
        vtkDoubleArray *quat3Data;
        vtkDoubleArray *quat4Data;
        vtkDoubleArray *rotationMatrixData;
        if (isSuperQ) {
            shapeXData = vtkDoubleArray::SafeDownCast(pd->GetArray("shapex"));
            shapeYData = vtkDoubleArray::SafeDownCast(pd->GetArray("shapey"));
            shapeZData = vtkDoubleArray::SafeDownCast(pd->GetArray("shapez"));
            quat1Data = vtkDoubleArray::SafeDownCast(pd->GetArray("quat1"));
            quat2Data = vtkDoubleArray::SafeDownCast(pd->GetArray("quat2"));
            quat3Data = vtkDoubleArray::SafeDownCast(pd->GetArray("quat3"));
            quat4Data = vtkDoubleArray::SafeDownCast(pd->GetArray("quat4"));
            rotationMatrixData = vtkDoubleArray::SafeDownCast(pd->GetArray("TENSOR"));
        } else {
            radiusData = vtkDoubleArray::SafeDownCast(pd->GetArray("radius"));
        }

        double radius;
        double shapex, shapey, shapez;
        double quat1, quat2, quat3, quat4;
        double rotationMatrix[9];

        double p[3];
        double leftX = VTK_DOUBLE_MAX;
        double rightX = VTK_DOUBLE_MIN;
        double backY = VTK_DOUBLE_MAX;
        double frontY = VTK_DOUBLE_MIN;
        double topZ = VTK_DOUBLE_MIN;
        double bottomZ = VTK_DOUBLE_MAX;

        for (vtkIdType i = 0; i < output->GetNumberOfPoints(); ++i) {
            output->GetPoint(i, p);
            std::shared_ptr<Particle> particlePtr = nullptr;
            if (isSuperQ) {
                shapex = shapeXData->GetValue(i);
                shapey = shapeYData->GetValue(i);
                shapez = shapeZData->GetValue(i);
                rotationMatrixData->GetTupleValue(int(i), rotationMatrix);
                if (quat1Data != nullptr) quat1 = quat1Data->GetValue(i);
                else quat1 = 1;
                if (quat2Data != nullptr) quat2 = quat2Data->GetValue(i);
                else quat2 = 0;
                if (quat3Data != nullptr) quat3 = quat3Data->GetValue(i);
                else quat3 = 0;
                if (quat4Data != nullptr) quat4 = quat4Data->GetValue(i);
                else quat4 = 0;
                Eigen::MatrixXd rotMat(3,3);
                rotMat(0,0) = rotationMatrix[0];
                rotMat(0,1) = rotationMatrix[1];
                rotMat(0,2) = rotationMatrix[2];
                rotMat(1,0) = rotationMatrix[3];
                rotMat(1,1) = rotationMatrix[4];
                rotMat(1,2) = rotationMatrix[5];
                rotMat(2,0) = rotationMatrix[6];
                rotMat(2,1) = rotationMatrix[7];
                rotMat(2,2) = rotationMatrix[8];
                particlePtr = std::make_shared<SuperquadricParticle>(0, shapex, shapey, shapez, p[0], p[1], p[2],
                                                                     quat1, quat2, quat3, quat4, rotMat, 0,
                                                                     0, 0);
            } else {
                radius = radiusData->GetValue(i);
                particlePtr = std::make_shared<SphericalParticle>(0, radius, p[0], p[1], p[2], 0,
                                                                  0, 0);
            }
            topZ = std::max(topZ, particlePtr->getTopZ());
            bottomZ = std::min(bottomZ, particlePtr->getBottomZ());
            leftX = std::min(leftX, particlePtr->x);
            rightX = std::max(rightX, particlePtr->x);
            backY = std::min(backY, particlePtr->y);
            frontY = std::max(frontY, particlePtr->y);
        }
        return std::vector<std::pair<double, double>>{{leftX,   rightX},
                                                      {backY,   frontY},
                                                      {bottomZ, topZ}};
    }
    return std::vector<std::pair<double, double>>{{0, 0},
                                                  {0, 0},
                                                  {0, 0}};
}

#endif //LIGGGHTS_VBP_VTKREADER_BEDLIMS_H
