//
// Created by karthik on 23/11/20.
//

#ifndef LIGGGHTS_VBP_VTKREADER_BEDHEIGHT_H
#define LIGGGHTS_VBP_VTKREADER_BEDHEIGHT_H

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

std::pair<double, double> getBedHeight(std::string filePath, bool isSuperQ) {
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
        if (isSuperQ) {
            shapeXData = vtkDoubleArray::SafeDownCast(pd->GetArray("shapex"));
            shapeYData = vtkDoubleArray::SafeDownCast(pd->GetArray("shapey"));
            shapeZData = vtkDoubleArray::SafeDownCast(pd->GetArray("shapez"));
            quat1Data = vtkDoubleArray::SafeDownCast(pd->GetArray("quat1"));
            quat2Data = vtkDoubleArray::SafeDownCast(pd->GetArray("quat2"));
            quat3Data = vtkDoubleArray::SafeDownCast(pd->GetArray("quat3"));
            quat4Data = vtkDoubleArray::SafeDownCast(pd->GetArray("quat4"));
        } else {
            radiusData = vtkDoubleArray::SafeDownCast(pd->GetArray("radius"));
        }

        double radius;
        double shapex, shapey, shapez;
        double quat1, quat2, quat3, quat4;

        double p[3];
        double topZ = VTK_DOUBLE_MIN;
        double bottomZ = VTK_DOUBLE_MAX;

        for (vtkIdType i = 0; i < output->GetNumberOfPoints(); ++i) {
            output->GetPoint(i, p);
            std::shared_ptr<Particle> particlePtr = nullptr;
            if (isSuperQ) {
                shapex = shapeXData->GetValue(i);
                shapey = shapeYData->GetValue(i);
                shapez = shapeZData->GetValue(i);
                if (quat1Data != nullptr) quat1 = quat1Data->GetValue(i);
                else quat1 = 1;
                if (quat2Data != nullptr) quat2 = quat2Data->GetValue(i);
                else quat2 = 0;
                if (quat3Data != nullptr) quat3 = quat3Data->GetValue(i);
                else quat3 = 0;
                if (quat4Data != nullptr) quat4 = quat4Data->GetValue(i);
                else quat4 = 0;
                particlePtr = std::make_shared<SuperquadricParticle>(0, shapex, shapey, shapez, p[0], p[1], p[2],
                                                                     quat1, quat2, quat3, quat4, 0,
                                                                     0, 0);
            } else {
                radius = radiusData->GetValue(i);
                particlePtr = std::make_shared<SphericalParticle>(0, radius, p[0], p[1], p[2], 0,
                                                                  0, 0);
            }
            topZ = std::max(topZ, particlePtr->getTopZ());
            bottomZ = std::min(bottomZ, particlePtr->getBottomZ());
        }
        return {bottomZ, topZ};
    }
    return {0, 0};
}

#endif //LIGGGHTS_VBP_VTKREADER_BEDHEIGHT_H
