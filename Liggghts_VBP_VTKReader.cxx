#include <vtkGenericDataObjectReader.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <string>
#include <set>
#include <map>

inline int getBin(double lo, double val, double increment) {
    return std::floor((val - lo) / increment);
}

int main(int argc, char *argv[]) {
    // Ensure a filename was specified
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " InputFilename" << endl;
        return EXIT_FAILURE;
    }

    // Get the filename from the command line
    std::string inputFilename = argv[1];

    // Get all data from the file
    vtkSmartPointer<vtkGenericDataObjectReader> reader =
            vtkSmartPointer<vtkGenericDataObjectReader>::New();
    reader->SetFileName(inputFilename.c_str());
    reader->Update();

    // All of the standard data types can be checked and obtained like this:
    if (reader->IsFilePolyData()) {
        std::cout << "output is a polydata" << std::endl;
        vtkPolyData *output = reader->GetPolyDataOutput();
        vtkPointData *pd = output->GetPointData();
        vtkDoubleArray *velData = vtkDoubleArray::SafeDownCast(pd->GetArray("v"));
        vtkDoubleArray *forceData = vtkDoubleArray::SafeDownCast(pd->GetArray("f"));
        vtkDoubleArray *omegaData = vtkDoubleArray::SafeDownCast(pd->GetArray("omega"));
        vtkDoubleArray *radiusData = vtkDoubleArray::SafeDownCast(pd->GetArray("radius"));
        vtkDoubleArray *massData = vtkDoubleArray::SafeDownCast(pd->GetArray("mass"));
        double velocity[3];
        double omega[3];
        double force[3];
        double radius;
        double mass;
        double p[3];
        double topZ = VTK_DOUBLE_MIN;
        double topZRadius;
        double bottomZ = VTK_DOUBLE_MAX;
        double bottomZRadius;
        std::set<double> radii;
        std::map<double, int> radiiIndexMap;

        std::cout << "output has " << output->GetNumberOfPoints() << " points." << std::endl;
        for (vtkIdType i = 0; i < output->GetNumberOfPoints(); ++i) {
            //velData->GetTupleValue(int(i), velocity);
            //forceData->GetTupleValue(int(i), force);
            //omegaData->GetTupleValue(int(i), omega);
            radius = radiusData->GetValue(i);
            radii.insert(radius);
            //mass = massData->GetValue(i);
            output->GetPoint(i, p);
            if (p[2] > topZ) {
                topZ = p[2];
                topZRadius = radius;
            }
            if (p[2] < bottomZ) {
                bottomZ = p[2];
                bottomZRadius = radius;
            }
            /*
            std::cout << "Point " << i << " : " << p[0] << " " << p[1] << " " << p[2] << "\n";
            std::cout << "Velocity " << i << " : " << velocity[0] << " " << velocity[1] << " " << velocity[2] << "\n";
            std::cout << "Omega " << i << " : " << omega[0] << " " << omega[1] << " " << omega[2] << "\n";
            std::cout << "Force " << i << " : " << force[0] << " " << force[1] << " " << force[2] << "\n";
            std::cout << "Radius " << i << " : " << radius << "\n";
            std::cout << "Mass " << i << " : " << mass << "\n";
            std::cout << "\n";
             */
        }
        bottomZ -= bottomZRadius;
        topZ += topZRadius;
        double zIncrement = (topZ - bottomZ) / 10.0;

        std::vector<long long> row(radii.size(), 0);
        std::vector<std::vector<long long >> particlesBin(10, row);
        int radiiIndex = 0;
        for (double it : radii) {
            radiiIndexMap[it] = radiiIndex++;
        }

        for (vtkIdType i = 0; i < output->GetNumberOfPoints(); ++i) {
            output->GetPoint(i, p);
            radius = radiusData->GetValue(i);
            int bin = getBin(bottomZ, p[2], zIncrement);
            int sizeType = radiiIndexMap[radius];
            particlesBin[bin][sizeType]++;
        }

        std::cout << "\nPARTICLE TYPES:\n";

        for (auto &it : radiiIndexMap) {
            std::cout << "Type " << it.second << " --> Radius = " << it.first << "\n";
        }

        std::cout << "\nPARTICLE COUNT: \n";
        for (int i = 9; i >= 0; --i) {
            std::cout << "Bin " << i << " : ";
            for (int j = 0; j < particlesBin[i].size(); ++j) {
                std::cout << "Type " << j << " = " << particlesBin[i][j] << " ; ";
            }
            std::cout << "\n";
        }
    }

    return EXIT_SUCCESS;
}
