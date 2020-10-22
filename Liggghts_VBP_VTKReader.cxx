#include <vtkGenericDataObjectReader.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <string>
#include <set>
#include <map>

//#include "matplotlibcpp.h"

//namespace plt = matplotlibcpp;

inline int getBin(double lo, double val, double increment) {
    return std::floor((val - lo) / increment);
}

struct Particle {
    bool operator()(const Particle &lhs, const Particle &rhs) const {
        return lhs.radius < rhs.radius;
    }

    bool operator<(const Particle &p) const {
        return this->radius < p.radius;
    }

    Particle(double radius, double mass) : radius(radius), mass(mass) {}

    Particle() {};
    double radius = 0;
    double mass = 0;
};

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
        std::set<Particle, Particle> particleSet;
        std::map<Particle, int, Particle> particleIndexMap;
        std::map<int, Particle> particleMap;
        int sliceCount = 10;

        std::cout << "output has " << output->GetNumberOfPoints() << " points." << std::endl;
        for (vtkIdType i = 0; i < output->GetNumberOfPoints(); ++i) {
            //velData->GetTupleValue(int(i), velocity);
            //forceData->GetTupleValue(int(i), force);
            //omegaData->GetTupleValue(int(i), omega);
            radius = radiusData->GetValue(i);
            mass = massData->GetValue(i);
            particleSet.insert(Particle(radius, mass));
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
        double zIncrement = (topZ - bottomZ) / sliceCount;

        std::vector<std::vector<long long >> particlesBin(sliceCount, std::vector<long long>(particleSet.size(), 0));
        std::vector<std::vector<double >> particlesMassBin(sliceCount, std::vector<double>(particleSet.size(), 0));
        std::vector<std::vector<double >> particlesMassFractionBin(sliceCount,
                                                                   std::vector<double>(particleSet.size(), 0));
        int particleIndex = 0;
        for (Particle particle : particleSet) {
            particleMap[particleIndex] = particle;
            particleIndexMap[particle] = particleIndex;
            particleIndex++;
        }

        for (vtkIdType i = 0; i < output->GetNumberOfPoints(); ++i) {
            output->GetPoint(i, p);
            radius = radiusData->GetValue(i);
            mass = massData->GetValue(i);
            int bin = getBin(bottomZ, p[2], zIncrement);
            int sizeType = particleIndexMap[Particle(radius, mass)];
            particlesBin[bin][sizeType]++;
        }

        double segregationIndex = 0;

        for (int i = 0; i < particlesBin.size(); ++i) {
            double totalMass = 0;
            for (int j = 0; j < particlesMassBin[i].size(); ++j) {
                particlesMassBin[i][j] = particlesBin[i][j] * particleMap[j].mass;
                totalMass += particlesMassBin[i][j];
            }
            for (int j = 0; j < particlesMassFractionBin[i].size(); ++j) {
                particlesMassFractionBin[i][j] = particlesMassBin[i][j] / totalMass;
            }
            segregationIndex += (particlesMassFractionBin[i][0] - 1.0) * (particlesMassFractionBin[i][0] - 1.0);
        }
        segregationIndex /= particlesBin.size();
        segregationIndex = std::sqrt(segregationIndex);

        std::cout << "\nPARTICLE TYPES:\n";

        for (auto &it : particleIndexMap) {
            std::cout << "Type " << it.second << " --> Radius = " << it.first.radius << "; Mass = " << it.first.mass
                      << "\n";
        }

        std::cout << "\nPARTICLE COUNT: \n";
        for (int i = sliceCount - 1; i >= 0; --i) {
            std::cout << "Bin " << i << " : \n";
            for (int j = 0; j < particlesBin[i].size(); ++j) {
                std::cout << "Type " << j << "; Count = " << particlesBin[i][j] << " ; Slice Mass Fraction = "
                          << particlesMassFractionBin[i][j] << " ;\n";
            }
            std::cout << "\n";
        }
        std::cout << "\nOVERALL SEGREGATION INDEX = " << segregationIndex << "\n";

        // Plot fine fraction with distance
        std::vector<double> dist(particlesMassFractionBin.size(), 0);
        std::vector<double> fineMassFraction(particlesMassFractionBin.size(), 0);
        for (int i = 0; i < particlesMassFractionBin.size(); ++i) {
            dist[i] = (i + 0.5) * zIncrement;
            fineMassFraction[i] = particlesMassFractionBin[i][0];
        }
        /*
        plt::figure_size(1200, 800);
        plt::plot(dist, fineMassFraction);
        plt::title("DISTANCE FROM BOTTOM vs FINE MASS FRACTION");
        plt::xlabel("Distance from bottom of packed bed (m)");
        plt::ylabel("Fine mass fraction");
        plt::annotate("Overall segregation index = " + std::to_string(segregationIndex), 0.01, 0.9);
        plt::show();
         */
    }

    return EXIT_SUCCESS;
}
