#include <vtkGenericDataObjectReader.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <string>
#include <set>
#include <map>
#include <sys/stat.h>
#include <algorithm>
#include <Eigen/Dense>

#ifdef _WIN32
#include "dirent.h"
#endif
#ifdef __unix__

#include <dirent.h>
#include <memory>
#include <cfloat>
#include "particle.h"
#include "superquadricParticle.h"
#include "sphericalParticle.h"
#include "bedLims.h"
#include "particleSearcher.h"
#include "sphStdKernel.h"
#include "sphSpikyKernel.h"


#endif
//#include "matplotlibcpp.h"
//
//namespace plt = matplotlibcpp;

inline int getBin(double lo, double val, double increment) {
    return std::floor((val - lo) / increment);
}

template<typename T>
inline T clamp(const T &n, const T &lower, const T &upper) {
    return n <= lower ? lower : n >= upper ? upper : n;
}

//struct Particle {
//    bool operator()(const Particle &lhs, const Particle &rhs) const {
//        return lhs.radius < rhs.radius;
//    }
//
//    bool operator<(const Particle &p) const {
//        return this->radius < p.radius;
//    }
//
//    Particle(double radius, double mass) : radius(radius), mass(mass) {}
//
//    Particle() {};
//    double radius = 0;
//    double mass = 0;
//};

int isFile(const std::string &path) {
    DIR *directory = opendir(path.c_str());
    if (directory != NULL) {
        closedir(directory);
        return 0;
    }
    if (errno == ENOTDIR) {
        return 1;
    }
    return -1;
}

struct fileNameComparator {
    inline bool operator()(const std::string &str1, const std::string &str2) {
        int sIndex = 0;
        for (; sIndex < str1.length(); sIndex++) { if (isdigit(str1[sIndex])) break; }
        std::string str1_cp = str1.substr(sIndex, str1.length() - 1);
        long long fn1 = atoll(str1_cp.c_str());
        sIndex = 0;
        for (; sIndex < str1.length(); sIndex++) { if (isdigit(str2[sIndex])) break; }
        std::string str2_cp = str2.substr(sIndex, str2.length() - 1);
        long long fn2 = atoll(str2_cp.c_str());
        return fn1 < fn2;
    }
};

int main(int argc, char *argv[]) {
    // Ensure a filename was specified
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " InputFilename" << endl;
        return EXIT_FAILURE;
    }

    std::vector<std::string> filesVector;

    // Get the filename from the command line
    std::string inputFilename = argv[1];

    std::map<std::string, unsigned int> flagDispatchTable{
            {"-notall",    1u << 1u},
            {"-onlylast",  1u << 2u},
            {"-superq",    1u << 3u},
            {"-initframe", 1u << 4u},
            {"-velgrad",   1u << 5u}
    };
    unsigned int flags = 0;
    long long init_frame_count = LONG_LONG_MIN;
    std::string initFilePath = "";
    double bedBottom = -DBL_MAX;
    double bedTop = -DBL_MAX;
    double bedLeft = -DBL_MAX;
    double bedRight = -DBL_MAX;
    double bedBack = -DBL_MAX;
    double bedFront = -DBL_MAX;

    if (argc > 2) {
        for (int i = 2; i < argc; ++i) {
            std::string flg = argv[i];
            unsigned int flgVal = flagDispatchTable[flg];
            if (flgVal == flagDispatchTable["-initframe"]) {
                try {
                    if (i >= argc) {
                        std::cout << "Enter init frame number after -initframe flag!!\n";
                        exit(0);
                    }
                    std::string frame_count = argv[i + 1];
                    init_frame_count = std::stoll(frame_count);
                    i++;
                } catch (...) {
                    std::cout << "Enter init frame number after -initframe flag!!\n";
                    exit(0);
                }
            }
            flags = flags | flgVal;
        }
    }
    std::cout << "init_frame_count " << init_frame_count << "\n";
    std::string newDirName;
    setlocale(LC_ALL, "");
    if (isFile(inputFilename) == 1) {
        filesVector.push_back(inputFilename);
    } else {
#ifdef _WIN32
        if (inputFilename[inputFilename.size()-1] != '\\')
            inputFilename.push_back('\\');
        newDirName = inputFilename + "_post_processed\\";
        _mkdir(dir.data());
#endif
#ifdef __unix__
        if (inputFilename[inputFilename.size() - 1] != '/')
            inputFilename.push_back('/');
        newDirName = inputFilename + "_post_processed/";
        if (mkdir(newDirName.data(), 0777) == -1) {
            std::cout << "Unable to create directory: " << newDirName << "\n";
        } else {
            std::cout << "Created directory: " << newDirName << "\n";
        }
#endif

        DIR *dir = opendir(inputFilename.c_str());

        struct dirent *ent = nullptr;
        while ((ent = readdir(dir)) != nullptr) {
            if (ent->d_type == DT_REG) {
                std::string prospectiveFile = ent->d_name;
                if (prospectiveFile.substr(prospectiveFile.rfind('.') + 1) == "vtk" &&
                    prospectiveFile.find("boundingBox") == std::string::npos) {
                    filesVector.push_back(inputFilename + prospectiveFile);
                    if (initFilePath.empty() && init_frame_count != LONG_LONG_MIN) {
//                        std::cout<<"here\n";
                        std::string initFile = filesVector.back();
//                        std::cout<<initFile<<"\n";
                        int sIndex = 0;
                        for (; sIndex < initFile.length(); sIndex++) {
                            if (isdigit(initFile[sIndex]) && initFile[sIndex - 1] == 'p')break;
                        }
                        initFile = initFile.substr(sIndex, initFile.length() - 1);
                        long long fn = atoll(initFile.c_str());
                        if (fn == init_frame_count)
                            initFilePath = filesVector.back();
                    }
                }
            }
        }
        if (!initFilePath.empty() && init_frame_count != LONG_LONG_MIN) {
            std::vector<std::pair<double, double>> bedLims = getBedLims(initFilePath,
                                                                        flags & flagDispatchTable["-superq"]);
            bedLeft = bedLims[0].first;
            bedRight = bedLims[0].second;
            bedBack = bedLims[1].first;
            bedFront = bedLims[1].second;
            bedBottom = bedLims[2].first;
            bedTop = bedLims[2].second;
        } else if (initFilePath.empty() && init_frame_count != LONG_LONG_MIN) {
            std::cout << "Enter correct init frame number!!\n";
            exit(0);
        }
    }

    std::sort(filesVector.begin(), filesVector.end(), fileNameComparator());
    if (flags & flagDispatchTable["-notall"] && filesVector.size() >= 3) {
        std::string firstFile = filesVector[0];
        std::string lastFile = filesVector[filesVector.size() - 1];
        filesVector.clear();
        filesVector.push_back(firstFile);
        filesVector.push_back(lastFile);
    }
    if (flags & flagDispatchTable["-onlylast"] && filesVector.size() > 1) {
        std::string lastFile = filesVector[filesVector.size() - 1];
        filesVector.clear();
        filesVector.push_back(lastFile);
    }
    std::cout << "File count = " << filesVector.size() << "\n";
    int fc = 0;
    for (auto &filePath : filesVector) {
        std::cout << ++fc << "\n";
        std::string postProcessedFileName;
        std::string postProcessedVelFileName;
        if (!newDirName.empty()) {
#if __unix__
            postProcessedFileName = newDirName + filePath.substr(filePath.rfind('/') + 1) + ".postprocessed";
            postProcessedVelFileName = newDirName + filePath.substr(filePath.rfind('/') + 1) + ".vel";
//            std::cout << postProcessedFileName << "\nHola\n";
#endif
#ifdef _WIN32
            postProcessedFileName = newDirName + filePath.substr(filePath.rfind('\\') + 1) + ".postprocessed";
            postProcessedVelFileName = newDirName + filePath.substr(filePath.rfind('\\') + 1) + ".vel";
#endif
        } else {
            postProcessedFileName = filePath + ".postprocessed";
            postProcessedVelFileName = filePath + ".vel";
        }
        ofstream postProcessedFile;
        if (!postProcessedFileName.empty()) {
            postProcessedFile.open(postProcessedFileName.c_str());
        }
        ofstream postProcessedVelFile;
        if (!postProcessedVelFileName.empty()) {
            postProcessedVelFile.open(postProcessedVelFileName.c_str());
        }
        // Get all data from the file
        vtkSmartPointer<vtkGenericDataObjectReader> reader =
                vtkSmartPointer<vtkGenericDataObjectReader>::New();
        reader->SetFileName(filePath.c_str());
        reader->Update();

        // All of the standard data types can be checked and obtained like this:
        if (reader->IsFilePolyData()) {
            postProcessedFile << "output is a polydata\n";
            vtkPolyData *output = reader->GetPolyDataOutput();
            vtkPointData *pd = output->GetPointData();

            bool isSuperQ = flags & flagDispatchTable["-superq"];


            vtkDoubleArray *velData = vtkDoubleArray::SafeDownCast(pd->GetArray("v"));
            /*
            vtkDoubleArray *forceData = vtkDoubleArray::SafeDownCast(pd->GetArray("f"));
            vtkDoubleArray *omegaData = vtkDoubleArray::SafeDownCast(pd->GetArray("omega"));
             */

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
            vtkDoubleArray *massData = vtkDoubleArray::SafeDownCast(pd->GetArray("mass"));


            double velocity[3];
            double rotationMatrix[9];
            /*
            double omega[3];
            double force[3];
             */
            double radius;
            double shapex, shapey, shapez;
            double quat1, quat2, quat3, quat4;
            double mass;
            double p[3];
            double leftX = VTK_DOUBLE_MAX;
            double rightX = VTK_DOUBLE_MIN;
            double backY = VTK_DOUBLE_MAX;
            double frontY = VTK_DOUBLE_MIN;
            double topZ = VTK_DOUBLE_MIN;
            double bottomZ = VTK_DOUBLE_MAX;
            std::set<std::pair<double, double>> particleVolumeMassSet;
            std::map<std::pair<double, double>, int> particleVolumeMassIndexMap;
            std::map<int, std::pair<double, double>> particleIndexVolumeMassMap;
            int sliceCount = 10;
            std::vector<std::shared_ptr<Particle>> particleVector(output->GetNumberOfPoints(), nullptr);
            std::map<std::pair<double, double>, std::vector<std::shared_ptr<Particle>>> particleVolumeMassParticleMap;

            postProcessedFile << "output has " << output->GetNumberOfPoints() << " points.\n";
            for (vtkIdType i = 0; i < output->GetNumberOfPoints(); ++i) {
                velData->GetTupleValue(int(i), velocity);
                //forceData->GetTupleValue(int(i), force);
                //omegaData->GetTupleValue(int(i), omega);
                mass = massData->GetValue(i);
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
                    particlePtr = std::make_shared<SuperquadricParticle>(mass, shapex, shapey, shapez, p[0], p[1], p[2],
                                                                         quat1, quat2, quat3, quat4, rotMat, velocity[0],
                                                                         velocity[1], velocity[2]);
                } else {
                    radius = radiusData->GetValue(i);
                    particlePtr = std::make_shared<SphericalParticle>(mass, radius, p[0], p[1], p[2], velocity[0],
                                                                      velocity[1], velocity[2]);
                }
                particleVolumeMassSet.insert({particlePtr->volume(), particlePtr->mass});
                particleVolumeMassParticleMap[{particlePtr->volume(), particlePtr->mass}].push_back(particlePtr);

//                std::cout<<"bedTop = "<<bedTop<<"\n";
//                std::cout<<"bedBottom = "<<bedBottom<<"\n";
//                std::cout<<"bedFront = "<<bedFront<<"\n";
//                std::cout<<"bedBack = "<<bedBack<<"\n";
//                std::cout<<"bedLeft = "<<bedLeft<<"\n";
//                std::cout<<"bedRight = "<<bedRight<<"\n";

                if (bedTop == -DBL_MAX && bedBottom == -DBL_MAX &&
                    bedFront == -DBL_MAX && bedBack == -DBL_MAX &&
                    bedLeft == -DBL_MAX && bedRight == -DBL_MAX) {
                    topZ = std::max(topZ, particlePtr->getTopZ());
                    bottomZ = std::min(bottomZ, particlePtr->getBottomZ());
                    leftX = std::min(leftX, particlePtr->x);
                    rightX = std::max(rightX, particlePtr->x);
                    backY = std::min(backY, particlePtr->y);
                    frontY = std::max(frontY, particlePtr->y);
                }
                particleVector[int(i)] = particlePtr;
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

            double xIncrement = 0;
            double yIncrement = 0;
            double zIncrement = 0;
            if (bedTop == -DBL_MAX && bedBottom == -DBL_MAX &&
                bedFront == -DBL_MAX && bedBack == -DBL_MAX &&
                bedLeft == -DBL_MAX && bedRight == -DBL_MAX) {
                xIncrement = (rightX - leftX) / sliceCount;
                yIncrement = (frontY - backY) / sliceCount;
                zIncrement = (topZ - bottomZ) / sliceCount;
            } else {
                rightX = bedRight;
                leftX = bedLeft;
                frontY = bedFront;
                backY = bedBack;
                topZ = bedTop;
                bottomZ = bedBottom;
                xIncrement = (rightX - leftX) / sliceCount;
                yIncrement = (frontY - backY) / sliceCount;
                zIncrement = (topZ - bottomZ) / sliceCount;
            }

            std::vector<std::vector<long long >> particlesBin(sliceCount,
                                                              std::vector<long long>(particleVolumeMassSet.size(), 0));
            std::vector<std::vector<double >> particlesMassBin(sliceCount,
                                                               std::vector<double>(particleVolumeMassSet.size(), 0));
            std::vector<std::vector<double >> particlesMassFractionBin(sliceCount,
                                                                       std::vector<double>(particleVolumeMassSet.size(),
                                                                                           0));
            int particleIndex = 0;
            for (std::pair<double, double> particleVolMass : particleVolumeMassSet) {
                particleVolumeMassIndexMap[particleVolMass] = particleIndex;
                particleIndexVolumeMassMap[particleIndex] = particleVolMass;
                particleIndex++;
            }
            std::cout << "Particle count: " << particleVector.size() << "\n";
            int pc = 0;

            for (auto &particle : particleVector) {
                int binI = getBin(bottomZ, particle->z, zIncrement);
                int bin = clamp(binI, 0, sliceCount - 1);
                int sizeType = particleVolumeMassIndexMap[{particle->volume(), particle->mass}];
                particle->type = sizeType;
                particlesBin[bin][sizeType]++;
            }
            double segregationIndex = 0;

            for (int i = 0; i < particlesBin.size(); ++i) {
                double totalMass = 0;
                for (int j = 0; j < particlesMassBin[i].size(); ++j) {
                    particlesMassBin[i][j] = particlesBin[i][j] * particleIndexVolumeMassMap[j].second;
                    totalMass += particlesMassBin[i][j];
                }
                for (int j = 0; j < particlesMassFractionBin[i].size(); ++j) {
                    particlesMassFractionBin[i][j] = totalMass == 0 ? 0 : particlesMassBin[i][j] / totalMass;
//                    particlesMassFractionBin[i][j] = particlesMassBin[i][j] / totalMass;
                }
                segregationIndex +=
                        (particlesMassFractionBin[i][0] * 2.0 - 1.0) * (particlesMassFractionBin[i][0] * 2.0 - 1.0);
            }
            segregationIndex /= particlesBin.size();
            segregationIndex = std::sqrt(segregationIndex);


            postProcessedFile << "\nPARTICLE TYPES:\n";

            for (auto &it : particleIndexVolumeMassMap) {
                postProcessedFile << "Type " << it.first << " --> Volume = " << it.second.first << "; Mass = "
                                  << it.second.second
                                  << "\n";
            }

            postProcessedFile << "\nPARTICLE COUNT: \n";
            for (int i = sliceCount - 1; i >= 0; --i) {
                postProcessedFile << "Bin " << i << " : \n";
                for (int j = 0; j < particlesBin[i].size(); ++j) {
                    postProcessedFile << "Type " << j << "; Count = " << particlesBin[i][j]
                                      << " ; Slice Mass Fraction = "
                                      << particlesMassFractionBin[i][j] << " ;\n";
                }
                postProcessedFile << "\n";
            }

            postProcessedFile << "\nOVERALL SEGREGATION INDEX = " << segregationIndex << "\n";

            /*SMI CALCULATION*/

            int pTypeCount = particlesBin[0].size();
            std::vector<double> f(pTypeCount, 0);
            std::vector<long long> totalTypeCount(pTypeCount, 0);
            long long maxTypeCount = LONG_LONG_MIN;
            long long totalParticles = 0;
            for (int i = 0; i < totalTypeCount.size(); ++i) {
                for (int j = 0; j < particlesBin.size(); ++j) {
                    totalTypeCount[i] += particlesBin[j][i];
                }
                maxTypeCount = std::max(totalTypeCount[i], maxTypeCount);
                totalParticles += totalTypeCount[i];
            }
            for (int i = 0; i < f.size(); ++i) {
                f[i] = double(maxTypeCount) / double(totalTypeCount[i]);
            }
            std::vector<std::vector<double>> P(particlesBin.size(), std::vector<double>(pTypeCount, 0));
            for (int i = 0; i < P.size(); ++i) {
                double max_nf = LONG_LONG_MIN;
                for (int j = 0; j < P[i].size(); ++j) {
                    max_nf = std::max(max_nf, particlesBin[i][j] * f[j]);
                }
                for (int j = 0; j < P[i].size(); ++j) {
                    P[i][j] = (double(particlesBin[i][j]) * f[j]) / max_nf;
                }
            }
            std::vector<double> SMI_bin(particlesBin.size(), 0);
            double one_by_type_minus_one = 1.0 / (double(pTypeCount) - 1.0);
            for (int i = 0; i < SMI_bin.size(); ++i) {
                double Psum = 0;
                for (int j = 0; j < P[i].size(); ++j) {
                    Psum += P[i][j];
                }
                SMI_bin[i] = one_by_type_minus_one * (Psum - 1.0);
            }
            double SMI = 0;
            for (int i = 0; i < particlesBin.size(); ++i) {
                long long binParticleCount = 0;
                for (int j = 0; j < particlesBin[i].size(); ++j) {
                    binParticleCount += particlesBin[i][j];
                }
                SMI += SMI_bin[i] * double(binParticleCount);
            }
            SMI /= double(totalParticles);
            postProcessedFile << "\nSMI = " << SMI << "\n";
            //////////////////////////////////////////////////////////////////////////
            // S (order parameter calculation)
//            std::map<std::pair<double, double>, int> particleVolumeMassIndexMap;
//            std::map<int, std::pair<double, double>> particleIndexVolumeMassMap;
//            int sliceCount = 10;
//            std::vector<std::shared_ptr<Particle>> particleVector(output->GetNumberOfPoints(), nullptr);
//            std::map<std::pair<double, double>, std::vector<std::shared_ptr<Particle>>> particleVolumeMassParticleMap

            for (auto it = particleIndexVolumeMassMap.begin(); it != particleIndexVolumeMassMap.end(); it++) {
                std::pair<double, double> volumeMassPair = it->second;
                Eigen::MatrixXd Q(3, 3);
                Q.setZero();
                for (auto &particle: particleVolumeMassParticleMap[volumeMassPair]) {
                    Q += particle->getQMatrix();
                }
                Q /= particleVolumeMassParticleMap[volumeMassPair].size();
                Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(Q, false);
                double maxEigenVal = std::numeric_limits<double>::min();
                for (double ev: eigenSolver.eigenvalues().real()) {
                    maxEigenVal = std::max(maxEigenVal, ev);
                }
                postProcessedFile <<"\nType "<<it->first<< "; S (order parameter)= " << maxEigenVal << "\n";
            }

            Eigen::MatrixXd Q(3,3);
            Q.setZero();
            for (auto &particle : particleVector) {
                Q += particle->getQMatrix();
            }
            Q /= (particleVector.size());

            Eigen::EigenSolver<Eigen::MatrixXd> eigenSolver(Q, false);
//            postProcessedFile << "\n Q = " << Q << "\n";

            double maxEigenVal = std::numeric_limits<double>::min();
            for (double ev: eigenSolver.eigenvalues().real()) {
                maxEigenVal = std::max(maxEigenVal, ev);
            }
//            postProcessedFile << "\neigen vals = " << eigenSolver.eigenvalues() << "\n";
            postProcessedFile << "\nS (order parameter)= " << maxEigenVal << "\n";

            //////////////////////////////////////////////////////////////////////////

            // Plot fine fraction with distance
            std::vector<double>
                    dist(particlesMassFractionBin.size(), 0);
            std::vector<double> fineMassFraction(particlesMassFractionBin.size(), 0);
            for (int i = 0; i < particlesMassFractionBin.size(); ++i) {
                dist[i] = (i + 0.5) * zIncrement;
                fineMassFraction[i] = particlesMassFractionBin[i][0];
            }
            std::vector<double> x, y, z, u, v, w;
            double sum_vx, sum_vy, sum_vz;
            sum_vx = sum_vy = sum_vz = 0;
            std::vector<std::vector<double>> velXByType(particleIndexVolumeMassMap.size(), std::vector<double>(2, 0));
            std::vector<std::vector<double>> velYByType(particleIndexVolumeMassMap.size(), std::vector<double>(2, 0));
            std::vector<std::vector<double>> velZByType(particleIndexVolumeMassMap.size(), std::vector<double>(2, 0));
            std::vector<std::vector<double>> velXByTypeAbs(particleIndexVolumeMassMap.size(), std::vector<double>(2, 0));
            std::vector<std::vector<double>> velYByTypeAbs(particleIndexVolumeMassMap.size(), std::vector<double>(2, 0));
            std::vector<std::vector<double>> velZByTypeAbs(particleIndexVolumeMassMap.size(), std::vector<double>(2, 0));
            std::vector<std::vector<double>> velBinZ(particlesBin.size(), std::vector<double>(4, 0));
            std::vector<std::vector<double>> velBinX(particlesBin.size(), std::vector<double>(4, 0));
            std::vector<std::vector<double>> velBinY(particlesBin.size(), std::vector<double>(4, 0));
//            std::vector<std::vector<double>> velBinZAbs(particlesBin.size(), std::vector<double>(4, 0));
//            std::vector<std::vector<double>> velBinXAbs(particlesBin.size(), std::vector<double>(4, 0));
//            std::vector<std::vector<double>> velBinYAbs(particlesBin.size(), std::vector<double>(4, 0));
            std::vector<std::vector<std::vector<double>>> velBinZByType(particleIndexVolumeMassMap.size(),
                                                                        std::vector<std::vector<double>>(
                                                                                particlesBin.size(),
                                                                                std::vector<double>(4, 0)));
            std::vector<std::vector<std::vector<double>>> velBinXByType(particleIndexVolumeMassMap.size(),
                                                                        std::vector<std::vector<double>>(
                                                                                particlesBin.size(),
                                                                                std::vector<double>(4, 0)));
            std::vector<std::vector<std::vector<double>>> velBinYByType(particleIndexVolumeMassMap.size(),
                                                                        std::vector<std::vector<double>>(
                                                                                particlesBin.size(),
                                                                                std::vector<double>(4, 0)));
            std::vector<std::vector<std::vector<double>>> velBinZAbsByType(particleIndexVolumeMassMap.size(),
                                                                        std::vector<std::vector<double>>(
                                                                                particlesBin.size(),
                                                                                std::vector<double>(4, 0)));
            std::vector<std::vector<std::vector<double>>> velBinXAbsByType(particleIndexVolumeMassMap.size(),
                                                                        std::vector<std::vector<double>>(
                                                                                particlesBin.size(),
                                                                                std::vector<double>(4, 0)));
            std::vector<std::vector<std::vector<double>>> velBinYAbsByType(particleIndexVolumeMassMap.size(),
                                                                        std::vector<std::vector<double>>(
                                                                                particlesBin.size(),
                                                                                std::vector<double>(4, 0)));
            for (auto &i : particleVector) {
                sum_vx += i->vx;
                sum_vy += i->vy;
                sum_vz += i->vz;

                velXByType[i->type][0] += i->vx;
                velYByType[i->type][0] += i->vy;
                velZByType[i->type][0] += i->vz;
                velXByType[i->type][1] += 1;
                velYByType[i->type][1] += 1;
                velZByType[i->type][1] += 1;
                velXByTypeAbs[i->type][0] += std::abs(i->vx);
                velYByTypeAbs[i->type][0] += std::abs(i->vy);
                velZByTypeAbs[i->type][0] += std::abs(i->vz);
                velXByTypeAbs[i->type][1] += 1;
                velYByTypeAbs[i->type][1] += 1;
                velZByTypeAbs[i->type][1] += 1;

                int vBin = clamp(getBin(bottomZ, i->z, zIncrement), 0, sliceCount - 1);
                velBinZ[vBin][0] += i->vx;
                velBinZ[vBin][1] += i->vy;
                velBinZ[vBin][2] += i->vz;
                velBinZ[vBin][3] += 1;
                velBinZByType[i->type][vBin][0] += i->vx;
                velBinZByType[i->type][vBin][1] += i->vy;
                velBinZByType[i->type][vBin][2] += i->vz;
                velBinZByType[i->type][vBin][3] += 1;
//                velBinZAbs[vBin][0] += std::abs(i->vx);
//                velBinZAbs[vBin][1] += std::abs(i->vy);
//                velBinZAbs[vBin][2] += std::abs(i->vz);
//                velBinZAbs[vBin][3] += 1;
                velBinZAbsByType[i->type][vBin][0] += std::abs(i->vx);
                velBinZAbsByType[i->type][vBin][1] += std::abs(i->vy);
                velBinZAbsByType[i->type][vBin][2] += std::abs(i->vz);
                velBinZAbsByType[i->type][vBin][3] += 1;

                vBin = clamp(getBin(leftX, i->x, xIncrement), 0, sliceCount - 1);
                velBinX[vBin][0] += i->vx;
                velBinX[vBin][1] += i->vy;
                velBinX[vBin][2] += i->vz;
                velBinX[vBin][3] += 1;
                velBinXByType[i->type][vBin][0] += i->vx;
                velBinXByType[i->type][vBin][1] += i->vy;
                velBinXByType[i->type][vBin][2] += i->vz;
                velBinXByType[i->type][vBin][3] += 1;
//                velBinXAbs[vBin][0] += std::abs(i->vx);
//                velBinXAbs[vBin][1] += std::abs(i->vy);
//                velBinXAbs[vBin][2] += std::abs(i->vz);
//                velBinXAbs[vBin][3] += 1;
                velBinXAbsByType[i->type][vBin][0] += std::abs(i->vx);
                velBinXAbsByType[i->type][vBin][1] += std::abs(i->vy);
                velBinXAbsByType[i->type][vBin][2] += std::abs(i->vz);
                velBinXAbsByType[i->type][vBin][3] += 1;

                vBin = clamp(getBin(backY, i->y, yIncrement), 0, sliceCount - 1);
                velBinY[vBin][0] += i->vx;
                velBinY[vBin][1] += i->vy;
                velBinY[vBin][2] += i->vz;
                velBinY[vBin][3] += 1;
                velBinYByType[i->type][vBin][0] += i->vx;
                velBinYByType[i->type][vBin][1] += i->vy;
                velBinYByType[i->type][vBin][2] += i->vz;
                velBinYByType[i->type][vBin][3] += 1;
//                velBinYAbs[vBin][0] += std::abs(i->vx);
//                velBinYAbs[vBin][1] += std::abs(i->vy);
//                velBinYAbs[vBin][2] += std::abs(i->vz);
//                velBinYAbs[vBin][3] += 1;
                velBinYAbsByType[i->type][vBin][0] += std::abs(i->vx);
                velBinYAbsByType[i->type][vBin][1] += std::abs(i->vy);
                velBinYAbsByType[i->type][vBin][2] += std::abs(i->vz);
                velBinYAbsByType[i->type][vBin][3] += 1;
//                if (i->y <= 6.0 / 1000.0 && i->y >= -6.0 / 1000.0) {
                x.push_back(i->x);
                y.push_back(i->y);
                z.push_back(i->z);
                u.push_back(i->vx);
                v.push_back(i->vy);
                w.push_back(i->vz);
                postProcessedVelFile << i->x << " " << i->y << " " << i->z << " " << i->vx << " " << i->vy << " "
                                     << i->vz << "\n";
//                }
            }
            postProcessedFile << "\n";
            for (int i = velBinX.size() - 1; i >= 0; i--) {
                double v_avg =
                        sqrt(velBinX[i][0] * velBinX[i][0] + velBinX[i][1] * velBinX[i][1] +
                             velBinX[i][2] * velBinX[i][2]) /
                        velBinX[i][3];
                postProcessedVelFile << "AVG VEL | BIN-X = "<< i <<" | VALUE = " << v_avg << "\n";
                postProcessedFile << "AVG VEL | BIN-X = "<< i <<" | VALUE = " << v_avg << "\n";
            }
            postProcessedFile << "\n";
            for (int i = velBinY.size() - 1; i >= 0; i--) {
                double v_avg =
                        sqrt(velBinY[i][0] * velBinY[i][0] + velBinY[i][1] * velBinY[i][1] +
                             velBinY[i][2] * velBinY[i][2]) /
                        velBinY[i][3];
                postProcessedVelFile << "AVG VEL | BIN-Y = "<< i <<" | VALUE = " << v_avg << "\n";
                postProcessedFile << "AVG VEL | BIN-Y = "<< i <<" | VALUE = " << v_avg << "\n";
            }
            postProcessedFile << "\n";
            for (int i = velBinZ.size() - 1; i >= 0; i--) {
                double v_avg =
                        sqrt(velBinZ[i][0] * velBinZ[i][0] + velBinZ[i][1] * velBinZ[i][1] +
                             velBinZ[i][2] * velBinZ[i][2]) /
                        velBinZ[i][3];
                postProcessedVelFile << "AVG VEL | BIN-Z = "<< i <<" | VALUE = " << v_avg << "\n";
                postProcessedFile << "AVG VEL | BIN-Z = "<< i <<" | VALUE = " << v_avg << "\n";
            }

            postProcessedFile << "\n";
            for (int i = 0; i < velBinXByType.size(); i++) {
                for (int j = velBinXByType[i].size() - 1; j >= 0; j--) {
                    double v_avg = velBinXByType[i][j][3] == 0 ? 0 :
                                   sqrt(velBinXByType[i][j][0] * velBinXByType[i][j][0] +
                                        velBinXByType[i][j][1] * velBinXByType[i][j][1] +
                                        velBinXByType[i][j][2] * velBinXByType[i][j][2]) /
                                   velBinXByType[i][j][3];
                    postProcessedVelFile << "AVG VEL | TYPE = "<< i <<" | BIN-X = "<< j <<" | VALUE = " << v_avg << "\n";
                    postProcessedFile << "AVG VEL | TYPE = "<< i <<" | BIN-X = "<< j <<" | VALUE = " << v_avg << "\n";

                }
            }
            postProcessedFile << "\n";
            for (int i = 0; i < velBinYByType.size(); i++) {
                for (int j = velBinYByType[i].size() - 1; j >= 0; j--) {
                    double v_avg = velBinYByType[i][j][3] == 0 ? 0 :
                                   sqrt(velBinYByType[i][j][0] * velBinYByType[i][j][0] +
                                        velBinYByType[i][j][1] * velBinYByType[i][j][1] +
                                        velBinYByType[i][j][2] * velBinYByType[i][j][2]) /
                                   velBinYByType[i][j][3];
                    postProcessedVelFile << "AVG VEL | TYPE = "<< i <<" | BIN-Y = "<< j <<" | VALUE = " << v_avg << "\n";
                    postProcessedFile << "AVG VEL | TYPE = "<< i <<" | BIN-Y = "<< j <<" | VALUE = " << v_avg << "\n";
                }
            }
            postProcessedFile << "\n";
            for (int i = 0; i < velBinZByType.size(); i++) {
                for (int j = velBinZByType[i].size() - 1; j >= 0; j--) {
                    double v_avg = velBinZByType[i][j][3] == 0 ? 0 :
                                   sqrt(velBinZByType[i][j][0] * velBinZByType[i][j][0] +
                                        velBinZByType[i][j][1] * velBinZByType[i][j][1] +
                                        velBinZByType[i][j][2] * velBinZByType[i][j][2]) /
                                   velBinZByType[i][j][3];
                    postProcessedVelFile << "AVG VEL | TYPE = "<< i <<" | BIN-Z = "<< j <<" | VALUE = " << v_avg << "\n";
                    postProcessedFile << "AVG VEL | TYPE = "<< i <<" | BIN-Z = "<< j <<" | VALUE = " << v_avg << "\n";
                }
            }

            // CONVECTION VELOCITY CALCULATION
            postProcessedFile << "=================================================================\n";
            std::vector<std::vector<std::vector<double>>> convectionVelocityBinXByType(velBinXByType.size(), std::vector<std::vector<double>>(velBinXByType[0].size(), std::vector<double>(3, 0)));
            std::vector<std::vector<std::vector<double>>> convectionVelocityBinYByType(velBinXByType.size(), std::vector<std::vector<double>>(velBinYByType[0].size(), std::vector<double>(3, 0)));
            std::vector<std::vector<std::vector<double>>> convectionVelocityBinZByType(velBinXByType.size(), std::vector<std::vector<double>>(velBinZByType[0].size(), std::vector<double>(3, 0)));
            // bins along x
            for (int i = velBinXByType.size()-1; i >=0 ; i--) {
                for (int j = velBinXByType[i].size()-1; j >=0 ; j--) {
                    convectionVelocityBinXByType[i][j][0] = (velBinXAbsByType[i][j][0] - std::abs(velBinXByType[i][j][0]))/velBinXByType[i][j][3];
                    convectionVelocityBinXByType[i][j][1] = (velBinXAbsByType[i][j][1] - std::abs(velBinXByType[i][j][1]))/velBinXByType[i][j][3];
                    convectionVelocityBinXByType[i][j][2] = (velBinXAbsByType[i][j][2] - std::abs(velBinXByType[i][j][2]))/velBinXByType[i][j][3];

                    postProcessedFile << "CONVECTION VEL | TYPE = " << i << " | BIN-X = "<< j <<" | X-COMPONENT = " << convectionVelocityBinXByType[i][j][0] << "\n";
                    postProcessedFile << "CONVECTION VEL | TYPE = " << i << " | BIN-X = "<< j <<" | Y-COMPONENT = " << convectionVelocityBinXByType[i][j][1] << "\n";
                    postProcessedFile << "CONVECTION VEL | TYPE = " << i << " | BIN-X = "<< j <<" | Z-COMPONENT = " << convectionVelocityBinXByType[i][j][2] << "\n";
                }
                postProcessedFile << "\n";
            }
            postProcessedFile << "=================================================================\n";
            // bins along y
            for (int i = velBinYByType.size()-1; i >=0 ; i--) {
                for (int j = velBinYByType[i].size()-1; j >=0 ; j--) {
                    convectionVelocityBinYByType[i][j][0] = (velBinYAbsByType[i][j][0] - std::abs(velBinYByType[i][j][0]))/velBinYByType[i][j][3];
                    convectionVelocityBinYByType[i][j][1] = (velBinYAbsByType[i][j][1] - std::abs(velBinYByType[i][j][1]))/velBinYByType[i][j][3];
                    convectionVelocityBinYByType[i][j][2] = (velBinYAbsByType[i][j][2] - std::abs(velBinYByType[i][j][2]))/velBinYByType[i][j][3];

                    postProcessedFile << "CONVECTION VEL | TYPE = " << i << " | BIN-Y = "<< j <<" | X-COMPONENT = " << convectionVelocityBinYByType[i][j][0] << "\n";
                    postProcessedFile << "CONVECTION VEL | TYPE = " << i << " | BIN-Y = "<< j <<" | Y-COMPONENT = " << convectionVelocityBinYByType[i][j][1] << "\n";
                    postProcessedFile << "CONVECTION VEL | TYPE = " << i << " | BIN-Y = "<< j <<" | Z-COMPONENT = " << convectionVelocityBinYByType[i][j][2] << "\n";
                }
                postProcessedFile << "\n";
            }
            postProcessedFile << "=================================================================\n";
            // bins along z
            for (int i = velBinZByType.size()-1; i >=0 ; i--) {
                for (int j = velBinZByType[i].size()-1; j >=0 ; j--) {
                    convectionVelocityBinZByType[i][j][0] = (velBinZAbsByType[i][j][0] - std::abs(velBinZByType[i][j][0]))/velBinZByType[i][j][3];
                    convectionVelocityBinZByType[i][j][1] = (velBinZAbsByType[i][j][1] - std::abs(velBinZByType[i][j][1]))/velBinZByType[i][j][3];
                    convectionVelocityBinZByType[i][j][2] = (velBinZAbsByType[i][j][2] - std::abs(velBinZByType[i][j][2]))/velBinZByType[i][j][3];

                    postProcessedFile << "CONVECTION VEL | TYPE = " << i << " | BIN-Z = "<< j <<" | X-COMPONENT = " << convectionVelocityBinZByType[i][j][0] << "\n";
                    postProcessedFile << "CONVECTION VEL | TYPE = " << i << " | BIN-Z = "<< j <<" | Y-COMPONENT = " << convectionVelocityBinZByType[i][j][1] << "\n";
                    postProcessedFile << "CONVECTION VEL | TYPE = " << i << " | BIN-Z = "<< j <<" | Z-COMPONENT = " << convectionVelocityBinZByType[i][j][2] << "\n";
                }
                postProcessedFile << "\n";
            }
            postProcessedFile << "=================================================================\n\n";

            std::vector<double> bulkXConvVelByType(velXByType.size(), 0);
            std::vector<double> bulkYConvVelByType(velYByType.size(), 0);
            std::vector<double> bulkZConvVelByType(velZByType.size(), 0);
            for (int i = bulkXConvVelByType.size() - 1; i >= 0 ; i--) {
                bulkXConvVelByType[i] = (velXByTypeAbs[i][0] - std::abs(velXByType[i][0])) / velXByTypeAbs[i][1];
//                postProcessedFile << "BULK CONVECTION VEL | TYPE = "<< i <<" | X-COMPONENT = " << bulkXConvVelByType[i] << "\n";
            }
            for (int i = bulkYConvVelByType.size() - 1; i >= 0 ; i--) {
                bulkYConvVelByType[i] = (velYByTypeAbs[i][0] - std::abs(velYByType[i][0])) / velYByTypeAbs[i][1];
//                postProcessedFile << "BULK CONVECTION VEL | TYPE = "<< i <<" | Y-COMPONENT = " << bulkYConvVelByType[i] << "\n";
            }
            for (int i = bulkZConvVelByType.size() - 1; i >= 0 ; i--) {
                bulkZConvVelByType[i] = (velZByTypeAbs[i][0] - std::abs(velZByType[i][0])) / velZByTypeAbs[i][1];
//                postProcessedFile << "BULK CONVECTION VEL | TYPE = "<< i <<" | Z-COMPONENT = " << bulkZConvVelByType[i] << "\n";
            }
            for (int i = bulkXConvVelByType.size() - 1; i >= 0 ; i--) {
                postProcessedFile << "BULK CONVECTION VEL | TYPE = "<< i <<" | X-COMPONENT = " << bulkXConvVelByType[i] << "\n";
                postProcessedFile << "BULK CONVECTION VEL | TYPE = "<< i <<" | Y-COMPONENT = " << bulkYConvVelByType[i] << "\n";
                postProcessedFile << "BULK CONVECTION VEL | TYPE = "<< i <<" | Z-COMPONENT = " << bulkZConvVelByType[i] << "\n";
            }
            postProcessedFile << "\n=================================================================\n\n";

            double v_avg = sqrt(sum_vx * sum_vx + sum_vy * sum_vy + sum_vz * sum_vz) / particleVector.size();
            postProcessedFile << "\nAVERAGE VELOCITY = " << v_avg << "\n";
            postProcessedVelFile << "\nAVERAGE VELOCITY = " << v_avg << "\n";
            postProcessedVelFile << "\nSMI = " << SMI << "\n";
            postProcessedVelFile.close();

            if (flags & flagDispatchTable["-velgrad"]) {
                int particleType = -1;
//                std::cout<<"here\n";
//                std::cout<<"size "<<particleVolumeMassParticleMap.size()<<"\n";
                for (auto volMassParticle: particleVolumeMassParticleMap) {
                    particleType++;
                    std::vector<std::vector<size_t>> neighborLists;
                    ParticleSearcher particleSearcher(volMassParticle.second, rightX, leftX, frontY, backY, topZ,
                                                      bottomZ);

                    neighborLists.resize(volMassParticle.second.size());
//                    std::cout<<"mapsize "<<volMassParticle.second.size()<<"\n";
                    double kernelRadius = 4 * volMassParticle.second.at(0)->getMaxDimLength();
//                    std::cout<<"kernelRadius "<<kernelRadius<<"\n";
                    double _mass = volMassParticle.second.at(0)->mass;
//                    std::cout<<"_mass "<<_mass<<"\n";

                    for (size_t i = 0; i < volMassParticle.second.size(); ++i) {
                        glm::dvec3 origin(volMassParticle.second.at(i)->x, volMassParticle.second.at(i)->y,
                                          volMassParticle.second.at(i)->z);
                        neighborLists[i].clear();
                        particleSearcher.forEachNearbyPoint(origin, kernelRadius,
                                                            [&](size_t j, const std::shared_ptr<Particle> &particle) {
                                                                if (i != j)neighborLists[i].push_back(j);
//                                                                std::cout<<"foreach "<<i<<" "<<j<<"\n";
                                                            });
//                        if (neighborLists[i].size() != 0) {
//                            std::cout << "NL " << i << "\n";
//                        }
//                        std::cout<<"neigborcoount "<<neighborLists[i].size()<<"\n";
                    }

                    //calculate densities TODO PARALLELIZE
                    std::vector<double> d(volMassParticle.second.size(), 0);
                    for (int i = 0; i < d.size(); ++i) {
                        double sum = 0;
                        SphStdKernel kernel(kernelRadius);
                        glm::dvec3 origin(volMassParticle.second.at(i)->x, volMassParticle.second.at(i)->y,
                                          volMassParticle.second.at(i)->z);
                        particleSearcher.forEachNearbyPoint(origin, kernelRadius,
                                                            [&](size_t j, const std::shared_ptr<Particle> &particle) {
                                                                glm::dvec3 neighborPos(particle->x, particle->y,
                                                                                       particle->z);
                                                                double dist = glm::distance(origin, neighborPos);
                                                                sum += kernel(dist);
                                                            });
                        d[i] = _mass * sum;
                    }

                    std::vector<double> duz_by_dz;
//                    std::cout<<"volMassParticle.second.size() "<<volMassParticle.second.size()<<"\n";
                    for (size_t i = 0; i < volMassParticle.second.size(); ++i) {
                        auto particle = volMassParticle.second.at(i);

                        const auto &neighbors = neighborLists[i];
                        if (neighbors.empty()) continue;

                        glm::dvec3 origin(particle->x, particle->y, particle->z);
                        SphSpikyKernel kernel(kernelRadius);
                        glm::dvec3 zVelGradSum(0, 0, 0);

                        for (size_t j = 0; j < neighbors.size(); ++j) {
                            glm::dvec3 neighborPos(volMassParticle.second.at(neighbors.at(j))->x,
                                                   volMassParticle.second.at(neighbors.at(j))->y,
                                                   volMassParticle.second.at(neighbors.at(j))->z);
                            double distance = glm::distance(origin, neighborPos);
                            if (distance > 0.0) {
                                glm::dvec3 dir = (neighborPos - origin) / distance;
//                                zVelGradSum +=
//                                        volMassParticle.second.at(neighbors.at(j))->vz * _mass / d[neighbors.at(j)] *
//                                        kernel.gradient(distance, dir);
                                zVelGradSum += d[i] * _mass *
                                               ((volMassParticle.second.at(i)->vz / (d[i] * d[i])) +
                                                (volMassParticle.second.at(neighbors.at(j))->vz /
                                                 (d[neighbors.at(j)] * d[neighbors.at(j)]))) *
                                               kernel.gradient(distance, dir);
                            }
                        }

                        duz_by_dz.push_back(zVelGradSum.z);
//                        std::cout << "velgrad "<<zVelGradSum.z << "\n";
                    }

                    //calculate mean gradient
                    int n = 0;
                    double meanGrad = 0;
                    for (auto val: duz_by_dz) {
                        double delta = val - meanGrad;
                        meanGrad += delta / ++n;
                    }
                    postProcessedFile << "\nduz/dz TYPE " << particleType << " = " << meanGrad << "\n";
//                    std::cout<<"\nduz/dz TYPE " << particleType << " = " << meanGrad << "\n";
                }

            }

//            plt::figure_size(1200, 800);
//            plt::plot(dist, fineMassFraction);
//            plt::title("DISTANCE FROM BOTTOM vs FINE MASS FRACTION");
//            plt::xlabel("Distance from bottom of packed bed (m)");
//            plt::ylabel("Fine mass fraction");
//            plt::annotate("Overall segregation index = " + std::to_string(segregationIndex), 0.01, 0.9);
//            plt::quiver(x,z,u,w);
//            plt::show();
        }
        if (postProcessedFile.is_open())
            postProcessedFile.close();
    }

    return EXIT_SUCCESS;
}
