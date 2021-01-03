#ifndef FLUID_NODE_H_
#define FLUID_NODE_H_

#include "gloo/SceneNode.hpp"
#include "ParticleState.hpp"
#include "gloo/debug/PrimitiveFactory.hpp"
#include "gloo/components/RenderingComponent.hpp"
#include "gloo/components/ShadingComponent.hpp"
#include "gloo/components/MaterialComponent.hpp"
#include "gloo/shaders/PhongShader.hpp"
#include "gloo/shaders/SimpleShader.hpp"
#include "gloo/InputManager.hpp"
#include "FluidGrid.hpp"
#include "FluidSystem.hpp"
#include "WaterCollectionNode.hpp"
#include <map>
#include <set>
#include <unordered_set>

struct PairTriples {
    GLOO::Triple t1;
    GLOO::Triple t2;

    bool operator== (const PairTriples &other) const {
        return (t1 == other.t1 && t2 == other.t2) || (t2 == other.t1 && t1 == other.t2);
    }
};

struct TriplePairTriples {
    PairTriples p1;
    PairTriples p2;
    PairTriples p3;

    bool operator== (const TriplePairTriples &other) const {
        return (p1 == other.p1 && p2 == other.p2 && p3 == other.p3) ||
               (p1 == other.p1 && p3 == other.p2 && p2 == other.p3) ||
               (p2 == other.p1 && p1 == other.p2 && p3 == other.p3) ||
               (p3 == other.p1 && p2 == other.p2 && p1 == other.p3) ||
               (p2 == other.p1 && p3 == other.p2 && p1 == other.p3) ||
               (p3 == other.p1 && p1 == other.p2 && p2 == other.p3);
    }
};

namespace std {

  template <>
  struct hash<PairTriples>
  {
    std::size_t operator()(const PairTriples& pair) const
    {
      using std::hash;

      return ((hash<int>()(pair.t1.i + pair.t2.i)
               ^ (hash<int>()(pair.t1.j + pair.t2.j) << 1)) >> 1)
               ^ (hash<int>()(pair.t1.k + pair.t2.k) << 1);
    }
  };

  template <>
  struct hash<TriplePairTriples>
  {
    std::size_t operator()(const TriplePairTriples& pair) const
    {
      using std::hash;

      return ((hash<PairTriples>()(pair.p1)
               ^ (hash<PairTriples>()(pair.p2) << 1)) >> 1)
               ^ (hash<PairTriples>()(pair.p3) << 1);
    }
  };

}

namespace GLOO {

static const std::vector<int> edgesToCheck = {0, 1, 0, 2, 0, 3, 1, 4, 1, 6, 2, 4, 2, 5, 3, 5, 5, 7, 6, 7, 4, 7, 3, 6};

class FluidNode: public SceneNode {

    public:

        FluidNode(double integrator_step) {
            cellWidth = 0.1f;
            gridSize = {20, 300, 20};
            int radius = 6;
            std::vector<Triple> solidCells, fluidCells;
            keyPressedE = false;
            keyPressedR = false;
            int i_c = gridSize.i / 2;
            int k_c = gridSize.k / 2;
            // int j_c = radius + 18;
            int j_c = 280;
            for (int i = 0; i < gridSize.i; i++) {
                for (int j = 0; j < gridSize.j; j++) {
                    for (int k = 0; k < gridSize.k; k++) {
                        if (i == 0 || j == 0 || k == 0 || i == gridSize.i - 1 || j == gridSize.j - 1 || k == gridSize.k - 1) {
                            solidCells.push_back({i,j,k});
                        // } else if (abs(k_c - k) < 5 - j && abs(i_c - i) < 5 - j && j > 1) {
                            // fluidCells.push_back({i,j,k});
                        // } else if (j == 7 && abs(i_c - i) < 3 ) {
                        //     solidCells.push_back({i,j,k});
                        } else if (j <= 40) {
                            fluidCells.push_back({i,j,k});  
                        } else if ( (i_c - i)*(i_c - i) + (j_c - j)*(j_c - j) + (k_c - k)*(k_c - k) <= radius * radius) {
                            fluidCells.push_back({i,j,k});  
                        // } else if ( (i_c - i)*(i_c - i) + (22 - j)*(22 - j) + (k_c - k)*(k_c - k) <= 2 * 2) {
                        //     fluidCells.push_back({i,j,k});  
                        }
                    }
                }
            }
            // fluidCells.push_back({i_c, 21, k_c});
            // for (int i = 0; i < gridSize.i; i++) {
            //     for (int j = 0; j < gridSize.j; j++) {
            //         for (int k = 0; k < gridSize.k; k++) {
            //             if (i == 0 || j == 0 || k == 0 || i == gridSize.i - 1 || j == gridSize.j - 1 || k == gridSize.k - 1) {
            //                 solidCells.push_back({i,j,k});
            //             // } else if (abs(k_c - k) < 5 - j && abs(i_c - i) < 5 - j && j > 1) {
            //             //     fluidCells.push_back({i,j,k});
            //             } else if (j <= 1) {
            //                 fluidCells.push_back({i,j,k});  
            //             } else if ( abs(i_c - i) < radius && abs(j_c - j) < radius && abs(k_c - k) < radius ) {
            //                 fluidCells.push_back({i,j,k});  
            //             }
            //         }
            //     }
            // }


            FluidGrid grid = FluidGrid(solidCells, fluidCells, gridSize, cellWidth, integrator_step);
            fluidSystem = make_unique<FluidSystem>(grid);

            for (int i = 0; i < gridSize.i; i++) {
                for (int j = 0; j < gridSize.j; j++) {
                    for (int k = 0; k < gridSize.k; k++) {
                        float phi = sqrt((i - i_c)*(i - i_c) + (j - j_c)*(j - j_c) + (k - k_c)*(k - k_c) + 0.00001f) - (float)radius;
                        float level_phi = j - 40.5f;
                        float best_phi = abs(phi) > abs(level_phi) ? level_phi : phi;

                        // float other_phi = sqrt((i - i_c)*(i - i_c) + (j - 22)*(j - 22) + (k - k_c)*(k - k_c) + 0.00001f) - (float)2;
                        // best_phi = best_phi > other_phi ? best_phi : other_phi;

                        fluidSystem->grid_ptr->SetPhi({i,j,k}, best_phi);
                        if (abs(best_phi) < 15.f * cellWidth && best_phi > 0.f) {
                            fluidSystem->grid_ptr->airCellsAtSurface->push_back({i,j,k});
                        } else if (abs(best_phi) < 15.f * cellWidth && best_phi < 0.f) {
                            fluidSystem->grid_ptr->liquidCellsAtSurface->push_back({i,j,k});
                        }
                    }
                }
            }
            fluidSystem->grid_ptr->FlipPhiStorage();

            vtx_obj_ptr = std::make_shared<VertexObject>();
            auto positions = make_unique<PositionArray>();
            auto normals = make_unique<NormalArray>();
            auto indices = make_unique<IndexArray>();

            // glm::vec3 pos1(0.f);
            // glm::vec3 pos2(0.f, 1.f, 0.f);
            // glm::vec3 pos3(1.f, 0.f, 0.f);
            // glm::vec3 pos4(1.f, 1.f, 0.f);
            // positions->push_back(pos1);
            // positions->push_back(pos2);
            // positions->push_back(pos3);
            // positions->push_back(pos4);

            // positions->push_back(pos2);
            // positions->push_back(pos3);
            // // positions->push_back(pos4);

            // normals->push_back(glm::vec3(0.f, 0.f, 1.f));
            // normals->push_back(glm::vec3(0.f, 0.f, 1.f));
            // std::vector<int> vertices = {0, 0, 0, 0, 0};
            // int a, b=1, c=2;
            // for (a = 0; a < b; a++) {
            //     for (b = a+1; b < c; b++) {
            //         for (c = b+1; c < vertices.size(); c++) {
            //             std::cout << "A: " << a << ", B: " << b << ", C: " << c << std::endl;
            //         }
            //     }
            // }
            // indices->push_back(0);
            // indices->push_back(1);
            // indices->push_back(2);
            // indices->push_back(3);
            // indices->push_back(4);
            // indices->push_back(5);
            // indices->push_back(3);
            // indices->push_back(1);
            // indices->push_back(2);
            // indices->push_back(0);
            // indices->push_back(1);
            // indices->push_back(2);
            // const Triple t1 = {1,2,3};
            // const Triple t2 = {4,5,6};
            // const PairTriples p1 = {t1,t2};
            // const PairTriples p2 = {t2, t1};
            // std::cout << "Pairtriples are equivalent: " << (p1 == p2) << std::endl;

            vtx_obj_ptr->UpdatePositions(std::move(positions));
            vtx_obj_ptr->UpdateNormals(std::move(normals));
            vtx_obj_ptr->UpdateIndices(std::move(indices));

            auto& rc = CreateComponent<RenderingComponent>(vtx_obj_ptr);
            rc.SetDrawMode(DrawMode::Triangles);
            CreateComponent<ShadingComponent>(std::make_shared<PhongShader>());
            glm::vec3 color(0.7f, 0.7f, 1.f);
            CreateComponent<MaterialComponent>(std::make_shared<Material>(color, color, color, 0));
            // WaterCollectionNode* water = new WaterCollectionNode(fluidSystem->grid_ptr->cellTypes_ptr, cellWidth, gridSize);
            RenderFluid();


            
            // auto waterNode = std::unique_ptr<WaterCollectionNode>(new WaterCollectionNode(fluidSystem->grid_ptr->cellTypes_ptr, cellWidth, gridSize));
            // AddChild(std::move(waterNode));  
            
            
            // std::cout << "Initialization finished" << std::endl;
            // Update(0.01f);
        }

        void Update(double delta_time) override {
            current_time += delta_time;
            // RenderFluid();

            if (InputManager::GetInstance().IsKeyPressed(GLFW_KEY_E)) {
                if (keyPressedE) {
                    return;
                } else {
                    keyPressedE = true;
                }
            } else {
                keyPressedE = false;
                return;
            }
            // if (!keyPressedE) {
            //     return;
            // }
            float step = 0.03f;
            std::cout << "Max time step is: " << fluidSystem->grid_ptr->max_time_step << std::endl;
            if (fluidSystem->grid_ptr->max_time_step < 0.03f) {
                std::cout << "Set new Step..." << std::endl;
                step = fluidSystem->grid_ptr->max_time_step;
            }
            fluidSystem->UpdateBodyForces(step);
            // std::cout << "Body forces finished" << std::endl;

            // for (int j = 0; j < gridSize.j; j++) {
            //     int i = gridSize.i / 2;
            //     int k = gridSize.k /2;
            //     std::cout << "PHI AT " << i << ", " << j << ", " << k << ": " << fluidSystem->grid_ptr->GetPhi({i,j,k}) << std::endl;
            //     std::cout << "  V: " << fluidSystem->grid_ptr->GetVMinus({i,j,k}) << std::endl;

            // }
            fluidSystem->AdvectPhi(step);
            // std::cout << "Phi Advection finished" << std::endl;

            // for (int j = 0; j < gridSize.j; j++) {
            //     int i = gridSize.i / 2;
            //     int k = gridSize.k /2;
            //     std::cout << "PHI AT " << i << ", " << j << ", " << k << ": " << fluidSystem->grid_ptr->GetPhi({i,j,k}) << std::endl;
            //     std::cout << "  V: " << fluidSystem->grid_ptr->GetVMinus({i,j,k}) << std::endl;

            // }
            fluidSystem->AdvectVelocity(step);
            // std::cout << "Advection finished" << std::endl;
            // for (int j = 0; j < gridSize.j; j++) {
            //     int i = gridSize.i / 2;
            //     int k = gridSize.k /2;
            //     std::cout << "PHI AT " << i << ", " << j << ", " << k << ": " << fluidSystem->grid_ptr->GetPhi({i,j,k}) << std::endl;
            //     std::cout << "  V: " << fluidSystem->grid_ptr->GetVMinus({i,j,k}) << std::endl;

            // }
            // fluidSystem->UpdateBodyForces(step);
            // std::cout << "Body Forces finished" << std::endl;

            fluidSystem->UpdatePressureSOE(step);
            // std::cout << "SOE finished" << std::endl;

            fluidSystem->UpdatePressure(100);
            // std::cout << "Pressure Solver finished" << std::endl;

            fluidSystem->ProjectVelocity(step);
            // std::cout << "Projection finished" << std::endl;

            // for (int j = 0; j < gridSize.j; j++) {
            //     int i = gridSize.i / 2;
            //     int k = gridSize.k /2;
            //     std::cout << "PHI pre AT " << i << ", " << j << ", " << k << ": " << fluidSystem->grid_ptr->GetPhi({i,j,k}) << std::endl;

            // }
            // std::cout << "Projection finished" << std::endl;

            // fluidSystem->UpdatePhi();
            // fluidSystem->UpdatePhi();
            fluidSystem->UpdatePhi();
            // for (int j = 0; j < gridSize.j; j++) {
            //     int i = gridSize.i / 2;
            //     int k = gridSize.k /2;
            //     std::cout << "PHI Post update AT " << i << ", " << j << ", " << k << ": " << fluidSystem->grid_ptr->GetPhi({i,j,k}) << std::endl;
            // }
            // std::cout << "Update Phi Finished" << std::endl;

            RenderFluid();


            // fluidSystem->UpdateGridMarkerCells(step);
            // std::cout << "Grid Marker Projection finished" << std::endl;
            // if (fluidSystem->grid_ptr->max_time_step < max_time_step) {
            //     max_time_step = fluidSystem->grid_ptr->max_time_step;
            // }
            // std::cout << "MAX_TIME_STEP: " << max_time_step << std::endl;


            // if (InputManager::GetInstance().IsKeyPressed(GLFW_KEY_R)) {
            //     if (keyPressedR) {
            //         return;
            //     } else {
            //         keyPressedR = true;
            //     }
            // } else {
            //     keyPressedR = false;
            //     return;
            // }


            // for (int i = 0; i < GetChildrenCount(); i++) {
            //     if (GetChild(i).IsActive()) {
            //         GetChild(i).SetActive(false);
            //         for (int j = 0; j < GetChild(i).GetChildrenCount(); j++) {
            //             GetChild(i).GetChild(j).SetActive(false);
            //         }
            //     }
                
            //     // delete &GetChild(i); //.SetActive(false)
            // }
            // auto newWaterNode = std::unique_ptr<WaterCollectionNode>(new WaterCollectionNode(fluidSystem->grid_ptr->cellTypes_ptr, cellWidth, gridSize));
            // AddChild(std::move(newWaterNode));
            // std::cout << "------------Render Update finished---------------" << std::endl;
            
            
            
            // for (int i = 0; i < fluidSystem->grid_ptr->gridWidth; i++) {
            //     for (int j = 0; j < fluidSystem->grid_ptr->gridHeight; j++) {
            //         for (int k = 0; k < fluidSystem->grid_ptr->gridDepth; k++) {
            //             if (fluidSystem->grid_ptr->GetCellType({i,j,k}) == CellType::LIQUID) {
            //                 std::cout << i << ", " << j << ", " << k << std::endl;
            //             }
            //         }
            //     }
            // }
            // update the position of the one node
            // GetTransform().SetPosition(state.positions[0]);
        }
    private:

        void RenderFluid() {
            // if (InputManager::GetInstance().IsKeyPressed(GLFW_KEY_R)) {
            //     if (keyPressedR) {
            //         return;
            //     } else {
            //         keyPressedR = true;
            //     }
            // } else {
            //     keyPressedR = false;
            //     return;
            // }
            // std::vector<int> edgesToCheck = {0, 1, 0, 2, 0, 3, 1, 4, 1, 6, 2, 4, 2, 5, 3, 5, 5, 7, 6, 7, 4, 7, 3, 6};

            auto edge_to_index_ptr = make_unique<std::unordered_map<PairTriples, int>>();
            auto set_surfaces = make_unique<std::unordered_set<TriplePairTriples>>();
            float cellWidth = fluidSystem->grid_ptr->cellWidth;
            //find all grid cells within 3 cellWidths of surface
            auto cellsNearSurface = make_unique<std::vector<Triple>>();
            for (int i = 0; i < fluidSystem->grid_ptr->gridWidth; i++) {
                for (int j = 0; j < fluidSystem->grid_ptr->gridHeight; j++) {
                    for (int k = 0; k < fluidSystem->grid_ptr->gridDepth; k++) {
                        if (abs(fluidSystem->grid_ptr->GetPhi({i,j,k})) <  3 * cellWidth) {
                            cellsNearSurface->push_back({i,j,k});
                        }
                    }
                }
            }
            // std::cout << "Cells Near surface: " << cellsNearSurface->size() << std::endl;
            //for every cell, look at 7 neighbors up and back and to the right, generate Vertex Object
            auto positions = make_unique<PositionArray>();
            auto normals = make_unique<NormalArray>();
            auto indices = make_unique<IndexArray>();
            std::cout << "Cells Near Surface: " << cellsNearSurface->size() << std::endl;
            for (int m = 0; m < cellsNearSurface->size(); m++) {
                // // if (m % 100 == 0) {
                // std::cout << m << " cells complete" << std::endl;
                // // }
                Triple triple = cellsNearSurface->at(m);
                int i = triple.i, j = triple.j, k = triple.k;
                std::vector<Triple> neighbors = {{i,j,k}, {i+1,j,k}, {i,j+1,k}, {i,j,k+1},
                                        {i+1,j+1,k}, {i+1,j,k+1}, {i,j+1,k+1}, {i+1,j+1,k+1}};
                std::vector<std::vector<PairTriples>> triangles/* = ComputeTriangles(triple, neighbors)*/;
                //check if on boundary with solid, if so return triangles of face no in solid.
                std::vector<PairTriples> edges = {{neighbors[0], neighbors[1]}, {neighbors[0], neighbors[2]}, {neighbors[0], neighbors[3]},
                                                {neighbors[1], neighbors[4]}, {neighbors[1], neighbors[5]}, {neighbors[2], neighbors[4]},
                                                {neighbors[2], neighbors[6]}, {neighbors[3], neighbors[5]}, {neighbors[3], neighbors[6]},
                                                {neighbors[4], neighbors[7]}, {neighbors[5], neighbors[7]}, {neighbors[6], neighbors[7]}};
                //face above
                if (fluidSystem->grid_ptr->GetCellType({i,j+1,k}) == CellType::SOLID) {
                    continue;
                } /* face behind */ else if (fluidSystem->grid_ptr->GetCellType({i,j,k+1}) == CellType::SOLID) {
                    continue;
                } /* face to right */ else if (fluidSystem->grid_ptr->GetCellType({i+1,j,k}) == CellType::SOLID) {
                    continue;
                }
                int gridNumber = 0;
                int powerOf2 = 128;
                for (int z = 0; z < neighbors.size(); z++) {
                    if (fluidSystem->grid_ptr->GetPhi(neighbors[z]) < 0.f) {
                        gridNumber += powerOf2;
                    }
                    powerOf2 /= 2;
                }

                // handle no triangles case (all fluid or all air)
                

                // handle case with alone triangle  

                if (gridNumber == 128 || 255 - gridNumber == 128 || gridNumber == 137 || 255 - gridNumber == 137 
                || gridNumber == 133 || 255 - gridNumber == 133 || gridNumber == 131 || 255 - gridNumber == 131
                || gridNumber == 140 || 255 - gridNumber == 140 || gridNumber == 134 || 255 - gridNumber == 134 
                || gridNumber == 138 || 255 - gridNumber == 138 || gridNumber == 136 || 255 - gridNumber == 136
                || gridNumber == 132 || 255 - gridNumber == 132 || gridNumber == 130 || 255 - gridNumber == 130 || gridNumber == 142 
                || gridNumber == 129 || 255 - gridNumber == 129 || gridNumber == 141 || gridNumber == 139 || gridNumber == 135) {
                    triangles.push_back({edges[0], edges[1], edges[2]});
                } 
                if (gridNumber == 64 || 255 - gridNumber == 64 || gridNumber == 98 || 255 - gridNumber == 98 
                || gridNumber == 82 || 255 - gridNumber == 82 || gridNumber == 67 || 255 - gridNumber == 67 
                || gridNumber == 112 || 255 - gridNumber == 112 || gridNumber == 81 || 255 - gridNumber == 81 
                || gridNumber == 97 || 255 - gridNumber == 97 || gridNumber == 96 || 255 - gridNumber == 96 
                || gridNumber == 80 || 255 - gridNumber == 80 || gridNumber == 66 || 255 - gridNumber == 66 || gridNumber == 113 
                || gridNumber == 65 || 255 - gridNumber == 65 || gridNumber == 83 || gridNumber == 99 || gridNumber == 114) {
                    triangles.push_back({edges[0], edges[3], edges[4]}); 
                } 
                if (gridNumber == 32 || 255 - gridNumber == 32 || gridNumber == 100 || 255 - gridNumber == 100 
                || gridNumber == 52 || 255 - gridNumber == 52 || gridNumber == 37 || 255 - gridNumber == 37 
                || gridNumber == 112 || 255 - gridNumber == 112 || gridNumber == 49 || 255 - gridNumber == 49 
                || gridNumber == 97 || 255 - gridNumber == 97 || gridNumber == 96 || 255 - gridNumber == 96 
                || gridNumber == 48 || 255 - gridNumber == 48 || gridNumber == 36 || 255 - gridNumber == 36 || gridNumber == 113 
                || gridNumber == 33 || 255 - gridNumber == 33 || gridNumber == 116 || gridNumber == 101 || gridNumber == 53) {
                    triangles.push_back({edges[1], edges[5], edges[6]}); 
                } 
                if (gridNumber == 16 || 255 - gridNumber == 16 || gridNumber == 88 || 255 - gridNumber == 88 
                || gridNumber == 56 || 255 - gridNumber == 56 || gridNumber == 25 || 255 - gridNumber == 25 
                || gridNumber == 112 || 255 - gridNumber == 112 || gridNumber == 81 || 255 - gridNumber == 81 
                || gridNumber == 49 || 255 - gridNumber == 49 || gridNumber == 80 || 255 - gridNumber == 80 
                || gridNumber == 48 || 255 - gridNumber == 48 || gridNumber == 24 || 255 - gridNumber == 24 || gridNumber == 113 
                || gridNumber == 17 || 255 - gridNumber == 17 || gridNumber == 120 || gridNumber == 89 || gridNumber == 57) {
                    triangles.push_back({edges[2], edges[7], edges[8]}); 
                } 
                if (gridNumber == 8 || 255 - gridNumber == 8 || gridNumber == 152 || 255 - gridNumber == 152 
                || gridNumber == 28 || 255 - gridNumber == 28 || gridNumber == 26 || 255 - gridNumber == 26 
                || gridNumber == 140 || 255 - gridNumber == 140 || gridNumber == 138 || 255 - gridNumber == 138 
                || gridNumber == 14 || 255 - gridNumber == 14 || gridNumber == 136 || 255 - gridNumber == 136 
                || gridNumber == 24 || 255 - gridNumber == 24 || gridNumber == 12 || 255 - gridNumber == 12 || gridNumber == 142
                || gridNumber == 10 || 255 - gridNumber == 10 || gridNumber == 30 || gridNumber == 156 || gridNumber == 154) {
                    triangles.push_back({edges[3], edges[5], edges[9]}); 
                } 
                if (gridNumber == 4 || 255 - gridNumber == 4 || gridNumber == 164 || 255 - gridNumber == 164 
                || gridNumber == 44 || 255 - gridNumber == 44 || gridNumber == 38 || 255 - gridNumber == 38 
                || gridNumber == 133 || 255 - gridNumber == 133 || gridNumber == 131 || 255 - gridNumber == 131 
                || gridNumber == 133 || 255 - gridNumber == 133 || gridNumber == 132 || 255 - gridNumber == 132
                || gridNumber == 36 || 255 - gridNumber == 36 || gridNumber == 12 || 255 - gridNumber == 12 || gridNumber == 142 
                || gridNumber == 6 || 255 - gridNumber == 6 || gridNumber == 166 || gridNumber == 172 || gridNumber == 46) {
                    triangles.push_back({edges[4], edges[7], edges[10]}); 
                } 
                if (gridNumber == 2 || 255 - gridNumber == 2 || gridNumber == 194 || 255 - gridNumber == 194 
                || gridNumber == 74 || 255 - gridNumber == 74 || gridNumber == 70 || 255 - gridNumber == 70 
                || gridNumber == 134 || 255 - gridNumber == 134 || gridNumber == 138 || 255 - gridNumber == 138 
                || gridNumber == 14 || 255 - gridNumber == 14 || gridNumber == 130 || 255 - gridNumber == 130 
                || gridNumber == 66 || 255 - gridNumber == 66 || gridNumber == 10 || 255 - gridNumber == 10 || gridNumber == 142 
                || gridNumber == 6 || 255 - gridNumber == 6 || gridNumber == 202 || gridNumber == 198 || gridNumber == 78) {
                    triangles.push_back({edges[6], edges[8], edges[11]}); 
                } 
                if (gridNumber == 1 || 255 - gridNumber == 1 || gridNumber == 193 || 255 - gridNumber == 193 
                || gridNumber == 161 || 255 - gridNumber == 161 || gridNumber == 145 || 255 - gridNumber == 145 
                || gridNumber == 81 || 255 - gridNumber == 81 || gridNumber == 49 || 255 - gridNumber == 49 
                || gridNumber == 97 || 255 - gridNumber == 97 || gridNumber == 129 || 255 - gridNumber == 129 
                || gridNumber == 65 || 255 - gridNumber == 65 || gridNumber == 33 || 255 - gridNumber == 33 || gridNumber == 113 
                || gridNumber == 17 || 255 - gridNumber == 17 || gridNumber == 225 || gridNumber == 209 || gridNumber == 177) {
                    triangles.push_back({edges[9], edges[10], edges[11]}); 
                } 

                // handle case of adjacent two nodes
                if (gridNumber == 192 || 255 - gridNumber == 192 || gridNumber == 194 || 255 - gridNumber == 194
                    || gridNumber == 193 || 255 - gridNumber == 193 || gridNumber == 195) {
                    triangles.push_back({edges[1], edges[2], edges[4]});
                    triangles.push_back({edges[1], edges[3], edges[4]});
                } 
                if (gridNumber == 160 || 255 - gridNumber == 160 || gridNumber == 164 || 255 - gridNumber == 164
                    || gridNumber == 161 || 255 - gridNumber == 161 || gridNumber == 165) {
                    triangles.push_back({edges[0], edges[2], edges[5]});
                    triangles.push_back({edges[2], edges[5], edges[6]});
                } 
                if (gridNumber == 144 || 255 - gridNumber == 144 || gridNumber == 152 || 255 - gridNumber == 152
                    || gridNumber == 145 || 255 - gridNumber == 145 || gridNumber == 153) {
                    triangles.push_back({edges[0], edges[1], edges[8]});
                    triangles.push_back({edges[0], edges[7], edges[8]});
                } 
                if (gridNumber == 72 || 255 - gridNumber == 72 || gridNumber == 88 || 255 - gridNumber == 88
                    || gridNumber == 74 || 255 - gridNumber == 74 || gridNumber == 90) {
                    triangles.push_back({edges[0], edges[4], edges[5]});
                    triangles.push_back({edges[4], edges[5], edges[9]});
                } 
                if (gridNumber == 68 || 255 - gridNumber == 68 || gridNumber == 100 || 255 - gridNumber == 100
                    || gridNumber == 70 || 255 - gridNumber == 70 || gridNumber == 102) {
                    triangles.push_back({edges[0], edges[3], edges[7]});
                    triangles.push_back({edges[3], edges[7], edges[10]});
                } 
                if (gridNumber == 40 || 255 - gridNumber == 40 || gridNumber == 56 || 255 - gridNumber == 56
                    || gridNumber == 44 || 255 - gridNumber == 44 || gridNumber == 60) {
                    triangles.push_back({edges[1], edges[3], edges[6]});
                    triangles.push_back({edges[3], edges[6], edges[9]});
                } 
                if (gridNumber == 34 || 255 - gridNumber == 34 || gridNumber == 98 || 255 - gridNumber == 98
                    || gridNumber == 38 || 255 - gridNumber == 38 || gridNumber == 102) {
                    triangles.push_back({edges[1], edges[5], edges[8]});
                    triangles.push_back({edges[5], edges[8], edges[11]});
                }
                if (gridNumber == 20 || 255 - gridNumber == 20 || gridNumber == 52 || 255 - gridNumber == 52
                    || gridNumber == 28 || 255 - gridNumber == 28 || gridNumber == 60) {
                    triangles.push_back({edges[2], edges[4], edges[8]});
                    triangles.push_back({edges[4], edges[8], edges[10]});
                } 
                if (gridNumber == 18 || 255 - gridNumber == 18 || gridNumber == 82 || 255 - gridNumber == 82
                    || gridNumber == 26 || 255 - gridNumber == 26 || gridNumber == 90) {
                    triangles.push_back({edges[2], edges[6], edges[7]});
                    triangles.push_back({edges[6], edges[7], edges[11]});
                } 
                if (gridNumber == 9 || 255 - gridNumber == 9 || gridNumber == 137 || 255 - gridNumber == 137
                    || gridNumber == 25 || 255 - gridNumber == 25 || gridNumber == 153) {
                    triangles.push_back({edges[3], edges[5], edges[10]});
                    triangles.push_back({edges[5], edges[10], edges[11]});
                }
                if (gridNumber == 5 || 255 - gridNumber == 5 || gridNumber == 133 || 255 - gridNumber == 133
                    || gridNumber == 37 || 255 - gridNumber == 37 || gridNumber == 165) {
                    triangles.push_back({edges[4], edges[7], edges[9]});
                    triangles.push_back({edges[7], edges[9], edges[11]});
                } else if (gridNumber == 3 || 255 - gridNumber == 3 || gridNumber == 131 || 255 - gridNumber == 131
                    || gridNumber == 67 || 255 - gridNumber == 67 || gridNumber == 195) {
                    triangles.push_back({edges[6], edges[8], edges[9]});
                    triangles.push_back({edges[8], edges[9], edges[10]});
                } 

                // handle case with three adjacent nodes
                if (gridNumber == 224 || 255 - gridNumber == 224 || gridNumber == 225) {
                    triangles.push_back({edges[2], edges[4], edges[6]});
                    triangles.push_back({edges[3], edges[5], edges[6]});
                    triangles.push_back({edges[3], edges[4], edges[6]});
                } else if (gridNumber == 162 || 255 - gridNumber == 162 || gridNumber == 166) {
                    triangles.push_back({edges[0], edges[5], edges[11]});
                    triangles.push_back({edges[0], edges[2], edges[8]});
                    triangles.push_back({edges[0], edges[8], edges[11]});
                } else if (gridNumber == 168 || 255 - gridNumber == 168 || gridNumber == 172) {
                    triangles.push_back({edges[2], edges[6], edges[9]});
                    triangles.push_back({edges[0], edges[2], edges[3]});
                    triangles.push_back({edges[2], edges[3], edges[9]});
                } else if (gridNumber == 208 || 255 - gridNumber == 208 || gridNumber == 209) {
                    triangles.push_back({edges[1], edges[3], edges[8]});
                    triangles.push_back({edges[3], edges[4], edges[7]});
                    triangles.push_back({edges[3], edges[7], edges[8]});
                } else if (gridNumber == 176 || 255 - gridNumber == 176 || gridNumber == 177) {
                    triangles.push_back({edges[0], edges[5], edges[7]});
                    triangles.push_back({edges[5], edges[6], edges[7]});
                    triangles.push_back({edges[6], edges[7], edges[8]});
                } else if (gridNumber == 200 || 255 - gridNumber == 200 || gridNumber == 202) {
                    triangles.push_back({edges[2], edges[4], edges[9]});
                    triangles.push_back({edges[1], edges[2], edges[5]});
                    triangles.push_back({edges[2], edges[5], edges[9]});
                } else if (gridNumber == 196 || 255 - gridNumber == 196 || gridNumber == 198) {
                    triangles.push_back({edges[1], edges[3], edges[10]});
                    triangles.push_back({edges[1], edges[2], edges[7]});
                    triangles.push_back({edges[1], edges[7], edges[10]});
                } else if (gridNumber == 22 || 255 - gridNumber == 22 || gridNumber == 30) {
                    triangles.push_back({edges[2], edges[4], edges[6]});
                    triangles.push_back({edges[4], edges[6], edges[10]});
                    triangles.push_back({edges[6], edges[10], edges[11]});
                } else if (gridNumber == 50 || 255 - gridNumber == 50 || gridNumber == 114) {
                    triangles.push_back({edges[5], edges[7], edges[11]});
                    triangles.push_back({edges[1], edges[2], edges[5]});
                    triangles.push_back({edges[2], edges[5], edges[7]});
                } else if (gridNumber == 35 || 255 - gridNumber == 35 || gridNumber == 99) {
                    triangles.push_back({edges[1], edges[8], edges[10]});
                    triangles.push_back({edges[1], edges[5], edges[9]});
                    triangles.push_back({edges[1], edges[9], edges[10]});
                } else if (gridNumber == 21 || 255 - gridNumber == 21 || gridNumber == 53) {
                    triangles.push_back({edges[2], edges[4], edges[9]});
                    triangles.push_back({edges[2], edges[8], edges[9]});
                    triangles.push_back({edges[8], edges[9], edges[11]});
                } else if (gridNumber == 19 || 255 - gridNumber == 19 || gridNumber == 83) {
                    triangles.push_back({edges[2], edges[6], edges[9]});
                    triangles.push_back({edges[2], edges[7], edges[9]});
                    triangles.push_back({edges[7], edges[9], edges[10]});
                } else if (gridNumber == 13 || 255 - gridNumber == 13 || gridNumber == 141) {
                    triangles.push_back({edges[5], edges[7], edges[11]});
                    triangles.push_back({edges[3], edges[4], edges[5]});
                    triangles.push_back({edges[4], edges[5], edges[7]});
                } else if (gridNumber == 11 || 255 - gridNumber == 11 || gridNumber == 139) {
                    triangles.push_back({edges[3], edges[8], edges[10]});
                    triangles.push_back({edges[3], edges[5], edges[6]});
                    triangles.push_back({edges[3], edges[6], edges[8]});
                } else if (gridNumber == 7 || 255 - gridNumber == 7 || gridNumber == 135) {
                    triangles.push_back({edges[4], edges[6], edges[9]});
                    triangles.push_back({edges[4], edges[6], edges[7]});
                    triangles.push_back({edges[6], edges[7], edges[8]});
                } else if (gridNumber == 42 || 255 - gridNumber == 42 || gridNumber == 46) {
                    triangles.push_back({edges[1], edges[3], edges[8]});
                    triangles.push_back({edges[3], edges[8], edges[9]});
                    triangles.push_back({edges[8], edges[9], edges[11]});
                } else if (gridNumber == 148 || 255 - gridNumber == 148 || gridNumber == 156) {
                    triangles.push_back({edges[1], edges[8], edges[10]});
                    triangles.push_back({edges[0], edges[1], edges[4]});
                    triangles.push_back({edges[1], edges[4], edges[10]});
                } else if (gridNumber == 146 || 255 - gridNumber == 146 || gridNumber == 154) {
                    triangles.push_back({edges[0], edges[7], edges[11]});
                    triangles.push_back({edges[0], edges[1], edges[6]});
                    triangles.push_back({edges[0], edges[6], edges[11]});
                } else if (gridNumber == 104 || 255 - gridNumber == 104 || gridNumber == 120) {
                    triangles.push_back({edges[4], edges[6], edges[9]});
                    triangles.push_back({edges[0], edges[1], edges[4]});
                    triangles.push_back({edges[1], edges[4], edges[6]});
                } else if (gridNumber == 73 || 255 - gridNumber == 73 || gridNumber == 89) {
                    triangles.push_back({edges[0], edges[5], edges[11]});
                    triangles.push_back({edges[0], edges[4], edges[10]});
                    triangles.push_back({edges[0], edges[10], edges[11]});
                } else if (gridNumber == 84 || 255 - gridNumber == 84 || gridNumber == 116) {
                    triangles.push_back({edges[3], edges[8], edges[10]});
                    triangles.push_back({edges[0], edges[2], edges[3]});
                    triangles.push_back({edges[2], edges[3], edges[8]});
                } else if (gridNumber == 69 || 255 - gridNumber == 69 || gridNumber == 101) {
                    triangles.push_back({edges[0], edges[7], edges[11]});
                    triangles.push_back({edges[0], edges[3], edges[9]});
                    triangles.push_back({edges[0], edges[9], edges[11]});
                } else if (gridNumber == 41 || 255 - gridNumber == 41 || gridNumber == 57) {
                    triangles.push_back({edges[1], edges[3], edges[10]});
                    triangles.push_back({edges[1], edges[6], edges[10]});
                    triangles.push_back({edges[6], edges[10], edges[11]});
                } else if (gridNumber == 76 || 255 - gridNumber == 76 || gridNumber == 78) {
                    triangles.push_back({edges[0], edges[5], edges[7]});
                    triangles.push_back({edges[5], edges[7], edges[9]});
                    triangles.push_back({edges[7], edges[9], edges[10]});
                }

                // handle flat sheet 4 nodes
                if (gridNumber == 212 || gridNumber == 43) {
                    triangles.push_back({edges[1], edges[3], edges[8]});
                    triangles.push_back({edges[3], edges[8], edges[10]});
                }
                if (gridNumber == 232 || gridNumber == 23) {
                    triangles.push_back({edges[2], edges[4], edges[6]});
                    triangles.push_back({edges[4], edges[6], edges[9]});
                }
                if (gridNumber == 178 || gridNumber == 77) {
                    triangles.push_back({edges[0], edges[5], edges[7]});
                    triangles.push_back({edges[5], edges[7], edges[11]});
                }

                // super corner 4 nodes
                if (gridNumber == 240 || gridNumber == 15) {
                    triangles.push_back({edges[3], edges[4], edges[7]});
                    triangles.push_back({edges[3], edges[5], edges[7]});
                    triangles.push_back({edges[5], edges[7], edges[8]});
                    triangles.push_back({edges[5], edges[6], edges[8]});
                } else if (gridNumber == 204 || gridNumber == 51) {
                    triangles.push_back({edges[1], edges[2], edges[5]});
                    triangles.push_back({edges[2], edges[5], edges[7]});
                    triangles.push_back({edges[5], edges[7], edges[9]});
                    triangles.push_back({edges[7], edges[9], edges[10]});
                } else if (gridNumber == 170 || gridNumber == 85) {
                    triangles.push_back({edges[0], edges[2], edges[3]});
                    triangles.push_back({edges[2], edges[3], edges[8]});
                    triangles.push_back({edges[3], edges[8], edges[9]});
                    triangles.push_back({edges[8], edges[9], edges[11]});
                } else if (gridNumber == 150 || gridNumber == 105) {
                    triangles.push_back({edges[1], edges[6], edges[11]});
                    triangles.push_back({edges[0], edges[1], edges[10]});
                    triangles.push_back({edges[1], edges[10], edges[11]});
                    triangles.push_back({edges[0], edges[4], edges[10]});
                }

                // line of 4 nodes 
                if (gridNumber == 201 || gridNumber == 54) {
                    triangles.push_back({edges[1], edges[2], edges[4]});
                    triangles.push_back({edges[1], edges[5], edges[11]});
                    triangles.push_back({edges[1], edges[4], edges[11]});
                    triangles.push_back({edges[4], edges[10], edges[11]});
                } else if (gridNumber == 197 || gridNumber == 58) {
                    triangles.push_back({edges[1], edges[2], edges[3]});
                    triangles.push_back({edges[2], edges[7], edges[11]});
                    triangles.push_back({edges[2], edges[3], edges[11]});
                    triangles.push_back({edges[3], edges[9], edges[11]});
                } else if (gridNumber == 149 || gridNumber == 106) {
                    triangles.push_back({edges[0], edges[1], edges[8]});
                    triangles.push_back({edges[8], edges[9], edges[11]});
                    triangles.push_back({edges[0], edges[8], edges[9]});
                    triangles.push_back({edges[0], edges[4], edges[9]});
                } else if (gridNumber == 147 || gridNumber == 108) {
                    triangles.push_back({edges[0], edges[1], edges[7]});
                    triangles.push_back({edges[1], edges[6], edges[9]});
                    triangles.push_back({edges[7], edges[9], edges[10]});
                    triangles.push_back({edges[1], edges[7], edges[9]});
                } else if (gridNumber == 169 || gridNumber == 86) {
                    triangles.push_back({edges[0], edges[2], edges[6]});
                    triangles.push_back({edges[6], edges[10], edges[11]});
                    triangles.push_back({edges[0], edges[3], edges[10]});
                    triangles.push_back({edges[0], edges[6], edges[10]});
                } else if (gridNumber == 163 || gridNumber == 92) {
                    triangles.push_back({edges[0], edges[2], edges[5]});
                    triangles.push_back({edges[2], edges[8], edges[10]});
                    triangles.push_back({edges[5], edges[9], edges[10]});
                    triangles.push_back({edges[2], edges[5], edges[10]});
                } else if (gridNumber == 226 || gridNumber == 29) {
                    triangles.push_back({edges[2], edges[3], edges[4]});
                    triangles.push_back({edges[2], edges[5], edges[11]});
                    triangles.push_back({edges[3], edges[5], edges[11]});
                    triangles.push_back({edges[2], edges[3], edges[11]});
                } else if (gridNumber == 210 || gridNumber == 45) {
                    triangles.push_back({edges[1], edges[6], edges[11]});
                    triangles.push_back({edges[1], edges[3], edges[4]});
                    triangles.push_back({edges[4], edges[7], edges[11]});
                    triangles.push_back({edges[1], edges[4], edges[11]});
                } else if (gridNumber == 75 || gridNumber == 180) {
                    triangles.push_back({edges[4], edges[5], edges[8]});
                    triangles.push_back({edges[5], edges[6], edges[8]});
                    triangles.push_back({edges[4], edges[8], edges[10]});
                    triangles.push_back({edges[0], edges[4], edges[5]});
                } else if (gridNumber == 71 || gridNumber == 184) {
                    triangles.push_back({edges[0], edges[3], edges[7]});
                    triangles.push_back({edges[6], edges[7], edges[8]});
                    triangles.push_back({edges[3], edges[6], edges[9]});
                    triangles.push_back({edges[3], edges[6], edges[7]});
                } else if (gridNumber == 228 || gridNumber == 27) {
                    triangles.push_back({edges[3], edges[5], edges[10]});
                    triangles.push_back({edges[2], edges[5], edges[6]});
                    triangles.push_back({edges[2], edges[7], edges[10]});
                    triangles.push_back({edges[2], edges[5], edges[10]});
                } else if (gridNumber == 39 || gridNumber == 216) {
                    triangles.push_back({edges[4], edges[5], edges[9]});
                    triangles.push_back({edges[1], edges[5], edges[8]});
                    triangles.push_back({edges[4], edges[7], edges[8]});
                    triangles.push_back({edges[4], edges[5], edges[8]});
                }


                
                // neighbors.push_back({i,j,k}); neighbors.push_back({i+1,j,k}); 
                // neighbors.push_back({i,j+1,k}); neighbors.push_back({i,j,k+1});
                // neighbors.push_back({i+1,j+1,k}); neighbors.push_back({i,j+1,k+1}); 
                // neighbors.push_back({i+1,j,k+1}); neighbors.push_back({i+1,j+1,k+1});

                //check each edge for cross boundary
                // std::vector<glm::vec3> vertices;
                // std::vector<glm::vec3> directions;
                // std::vector<PairTriples> boundingTriples;

                // if (triangles.size() > 0) {
                //     std::cout << triangles.size() << " triangles in this neighborhood" << std::endl;
                // }
                for (int n = 0; n < triangles.size(); n++) {
                    std::vector<PairTriples> triangle = triangles[n];
                    std::vector<glm::vec3> vertices;
                    std::vector<glm::vec3> directions;

                    bool e_brake = false;
                    for (int o = 0; o < triangle.size(); o++) {
                        PairTriples p = triangle[o];
                        Triple t1 = p.t1;
                        Triple t2 = p.t2;

                        float phi1 = fluidSystem->grid_ptr->GetPhi(t1);
                        float phi2 = fluidSystem->grid_ptr->GetPhi(t2);
                        glm::vec3 position1, position2;

                        if (phi1 * phi2 > 0.00001f) {
                            //same sign, not at boundary
                            e_brake = true;
                            break;
                            // std::cout << "WE HAVE MAJOR PROBLEMS" << std::endl;
                        } 
                        position1 = fluidSystem->grid_ptr->GetPosition(t1);
                        position2 = fluidSystem->grid_ptr->GetPosition(t2);

                        glm::vec3 interpolated_vertex = position1 + phi1 * cellWidth / (phi1 - phi2) * glm::normalize(position2 - position1);
                        glm::vec3 dir_to_air = phi1 > 0.f ? position1 - position2 : position2 - position1;
                        directions.push_back(glm::normalize(dir_to_air));
                        vertices.push_back(interpolated_vertex);
                    }
                    if (e_brake) {
                        continue;
                    }
                    glm::vec3 p1 = vertices[0], p2 = vertices[1], p3 = vertices[2];
                    glm::vec3 e1 = p1 - p2, e2 = p2 - p3;
                    glm::vec3 cross = glm::cross(e1, e2);

                    glm::vec3 average_dir = directions[0] + directions[1] + directions[2];
                    if (glm::dot(cross, average_dir) > 0.f) {
                        normals->push_back(cross);
                    } else {
                        normals->push_back(-cross);
                    }
                    const PairTriples pair1 = triangle[0];
                    const PairTriples pair2 = triangle[1];
                    const PairTriples pair3 = triangle[2];

                    const TriplePairTriples triplePair = {pair1, pair2, pair3};
                    if (set_surfaces->find(triplePair) != set_surfaces->end()) {
                        std::cout << gridNumber << std::endl;
                        continue;
                    }
                    set_surfaces->insert(triplePair);

                    auto it1 = edge_to_index_ptr->find(pair1);
                    auto it2 = edge_to_index_ptr->find(pair2);
                    auto it3 = edge_to_index_ptr->find(pair3);

                    int index1, index2, index3;
                    if (it1 != edge_to_index_ptr->end()) {
                        index1 = it1->second;
                    } else {
                        index1 = edge_to_index_ptr->size();
                        (*edge_to_index_ptr)[pair1] = index1;
                        positions->push_back(p1); 
                    }

                    if (it2 != edge_to_index_ptr->end()) {
                        index2 = it2->second;
                    } else {
                        index2 = edge_to_index_ptr->size();
                        (*edge_to_index_ptr)[pair2] = index2;
                        positions->push_back(p2); 
                    }

                    if (it3 != edge_to_index_ptr->end()) {
                        index3 = it3->second;
                    } else {
                        index3 = edge_to_index_ptr->size();
                        (*edge_to_index_ptr)[pair3] = index3;
                        positions->push_back(p3);
                    }

                    indices->push_back(index1); 
                    indices->push_back(index2); 
                    indices->push_back(index3);
                }



                // for (int n = 0; n < edgesToCheck.size(); n = n + 2) {
                //     const Triple p1 = neighbors[edgesToCheck[n]];
                //     const Triple p2 = neighbors[edgesToCheck[n+1]];
                //     if (fluidSystem->grid_ptr->GetCellType(p1) == CellType::SOLID || fluidSystem->grid_ptr->GetCellType(p2) == CellType::SOLID) {
                //         continue;
                //     }
                //     float phi1 = fluidSystem->grid_ptr->GetPhi(p1);
                //     float phi2 = fluidSystem->grid_ptr->GetPhi(p2);
                    
                //     glm::vec3 position1, position2;
                //     if (phi1 * phi2 > 0.f) {
                //         //same sign, not at boundary
                //         continue;
                //     } else {
                //         position1 = fluidSystem->grid_ptr->GetPosition(p1);
                //         position2 = fluidSystem->grid_ptr->GetPosition(p2);
                //     }
                //     glm::vec3 interpolated_vertex = position1 + phi1 * cellWidth / (phi1 - phi2) * glm::normalize(position2 - position1);
                //     glm::vec3 dir_to_air = phi1 > 0.f ? position1 - position2 : position2 - position1;
                //     directions.push_back(glm::normalize(dir_to_air));
                //     vertices.push_back(interpolated_vertex);

                //     PairTriples pair = {p1, p2};
                //     boundingTriples.push_back(pair);

                // }

                // if (vertices.size() == 0) {
                //     continue;
                // }
                // // if (m == 16282) {
                // //     std::cout << "made it to triangulation" << std::endl;
                // // }
                // //loop over every set of 3 vertices, and add them to positons and normals and indices
                // int a,b=1,c = 2;
                // for (a = 0; a < b; a++) {
                //     for (b = a+1; b < c; b++) {
                //         for (c = b+1; c < vertices.size(); c++) {
                //             glm::vec3 p1 = vertices[a], p2 = vertices[b], p3 = vertices[c];
                //             glm::vec3 e1 = p1 - p2, e2 = p2 - p3;
                //             glm::vec3 cross = glm::normalize(glm::cross(e1, e2));

                //             glm::vec3 average_dir = directions[a] + directions[b] + directions[c];
                //             if (glm::dot(cross, average_dir) > 0.f) {
                //                 normals->push_back(cross);
                //             } else {
                //                 normals->push_back(-cross);
                //             }
                //             const PairTriples pair1 = boundingTriples[a];
                //             const PairTriples pair2 = boundingTriples[b];
                //             const PairTriples pair3 = boundingTriples[c];

                //             const TriplePairTriples triplePair = {pair1, pair2, pair3};
                //             if (set_surfaces->find(triplePair) != set_surfaces->end()) {
                //                 continue;
                //             }
                //             set_surfaces->insert(triplePair);

                //             auto it1 = edge_to_index_ptr->find(pair1);
                //             auto it2 = edge_to_index_ptr->find(pair2);
                //             auto it3 = edge_to_index_ptr->find(pair3);

                //             int index1, index2, index3;
                //             if (it1 != edge_to_index_ptr->end()) {
                //                 index1 = it1->second;
                //             } else {
                //                 index1 = edge_to_index_ptr->size();
                //                 (*edge_to_index_ptr)[pair1] = index1;
                //                 positions->push_back(p1); 
                //             }

                //             if (it2 != edge_to_index_ptr->end()) {
                //                 index2 = it2->second;
                //             } else {
                //                 index2 = edge_to_index_ptr->size();
                //                 (*edge_to_index_ptr)[pair2] = index2;
                //                 positions->push_back(p2); 
                //             }

                //             if (it3 != edge_to_index_ptr->end()) {
                //                 index3 = it3->second;
                //             } else {
                //                 index3 = edge_to_index_ptr->size();
                //                 (*edge_to_index_ptr)[pair3] = index3;
                //                 positions->push_back(p3);
                //             }

                //             indices->push_back(index1); 
                //             indices->push_back(index2); 
                //             indices->push_back(index3);
                //         }
                //     }
                // }
            }
            // std::cout << "Size of positions: " << positions->size() << std::endl;
            // std::cout << "Size of indices: " << indices->size() << std::endl;

            vtx_obj_ptr->UpdatePositions(std::move(positions));
            vtx_obj_ptr->UpdateNormals(std::move(normals));
            vtx_obj_ptr->UpdateIndices(std::move(indices));
            UpdateNormals();
        }

        std::vector<std::vector<PairTriples>> ComputeTriangles(Triple cellOrigin, std::vector<Triple> corners) {
            int i = cellOrigin.i, j = cellOrigin.j, k = cellOrigin.k;
            //check if on boundary with solid, if so return triangles of face no in solid.
            std::vector<PairTriples> edges = {{corners[0], corners[1]}, {corners[0], corners[2]}, {corners[0], corners[3]},
                                              {corners[1], corners[4]}, {corners[1], corners[5]}, {corners[2], corners[4]},
                                              {corners[2], corners[6]}, {corners[3], corners[5]}, {corners[3], corners[6]},
                                              {corners[4], corners[7]}, {corners[5], corners[7]}, {corners[6], corners[7]}};
            //face above
            if (fluidSystem->grid_ptr->GetCellType({i,j+1,k}) == CellType::SOLID) {
                return std::vector<std::vector<PairTriples>>();
            } /* face behind */ else if (fluidSystem->grid_ptr->GetCellType({i,j,k+1}) == CellType::SOLID) {
                return std::vector<std::vector<PairTriples>>();
            } /* face to right */ else if (fluidSystem->grid_ptr->GetCellType({i+1,j,k}) == CellType::SOLID) {
                return std::vector<std::vector<PairTriples>>();
            }
            int gridNumber = 0;
            int powerOf2 = 128;
            for (int z = 0; z < corners.size(); z++) {
                if (fluidSystem->grid_ptr->GetPhi(corners[z]) < 0.f) {
                    gridNumber += powerOf2;
                }
                powerOf2 /= 2;
            }
            std::vector<std::vector<PairTriples>> triangles;

            // handle no triangles case (all fluid or all air)
            std::set<int> no_triangles = {0};
            if (no_triangles.find(gridNumber) != no_triangles.end() || no_triangles.find(255 - gridNumber) != no_triangles.end()) {
                return triangles;
            }

            // handle case with alone triangle  

            if (gridNumber == 128 || 255 - gridNumber == 128 || gridNumber == 137 || 255 - gridNumber == 137 
             || gridNumber == 133 || 255 - gridNumber == 133 || gridNumber == 131 || 255 - gridNumber == 131
             || gridNumber == 140 || 255 - gridNumber == 140 || gridNumber == 134 || 255 - gridNumber == 134 
             || gridNumber == 138 || 255 - gridNumber == 138 || gridNumber == 136 || 255 - gridNumber == 136
             || gridNumber == 132 || 255 - gridNumber == 132 || gridNumber == 130 || 255 - gridNumber == 130 || gridNumber == 142 
             || gridNumber == 129 || 255 - gridNumber == 129 || gridNumber == 141 || gridNumber == 139 || gridNumber == 135) {
                triangles.push_back({edges[0], edges[1], edges[2]});
            } 
            if (gridNumber == 64 || 255 - gridNumber == 64 || gridNumber == 98 || 255 - gridNumber == 98 
             || gridNumber == 82 || 255 - gridNumber == 82 || gridNumber == 67 || 255 - gridNumber == 67 
             || gridNumber == 112 || 255 - gridNumber == 112 || gridNumber == 81 || 255 - gridNumber == 81 
             || gridNumber == 97 || 255 - gridNumber == 97 || gridNumber == 96 || 255 - gridNumber == 96 
             || gridNumber == 80 || 255 - gridNumber == 80 || gridNumber == 66 || 255 - gridNumber == 66 || gridNumber == 113 
             || gridNumber == 65 || 255 - gridNumber == 65 || gridNumber == 83 || gridNumber == 99 || gridNumber == 114) {
                triangles.push_back({edges[0], edges[3], edges[4]}); 
            } 
            if (gridNumber == 32 || 255 - gridNumber == 32 || gridNumber == 100 || 255 - gridNumber == 100 
             || gridNumber == 52 || 255 - gridNumber == 52 || gridNumber == 37 || 255 - gridNumber == 37 
             || gridNumber == 112 || 255 - gridNumber == 112 || gridNumber == 49 || 255 - gridNumber == 49 
             || gridNumber == 97 || 255 - gridNumber == 97 || gridNumber == 96 || 255 - gridNumber == 96 
             || gridNumber == 48 || 255 - gridNumber == 48 || gridNumber == 36 || 255 - gridNumber == 36 || gridNumber == 113 
             || gridNumber == 33 || 255 - gridNumber == 33 || gridNumber == 116 || gridNumber == 101 || gridNumber == 53) {
                triangles.push_back({edges[1], edges[5], edges[6]}); 
            } 
            if (gridNumber == 16 || 255 - gridNumber == 16 || gridNumber == 88 || 255 - gridNumber == 88 
             || gridNumber == 56 || 255 - gridNumber == 56 || gridNumber == 25 || 255 - gridNumber == 25 
             || gridNumber == 112 || 255 - gridNumber == 112 || gridNumber == 81 || 255 - gridNumber == 81 
             || gridNumber == 49 || 255 - gridNumber == 49 || gridNumber == 80 || 255 - gridNumber == 80 
             || gridNumber == 48 || 255 - gridNumber == 48 || gridNumber == 24 || 255 - gridNumber == 24 || gridNumber == 113 
             || gridNumber == 17 || 255 - gridNumber == 17 || gridNumber == 120 || gridNumber == 89 || gridNumber == 57) {
                triangles.push_back({edges[2], edges[7], edges[8]}); 
            } 
            if (gridNumber == 8 || 255 - gridNumber == 8 || gridNumber == 152 || 255 - gridNumber == 152 
             || gridNumber == 28 || 255 - gridNumber == 28 || gridNumber == 26 || 255 - gridNumber == 26 
             || gridNumber == 140 || 255 - gridNumber == 140 || gridNumber == 138 || 255 - gridNumber == 138 
             || gridNumber == 14 || 255 - gridNumber == 14 || gridNumber == 136 || 255 - gridNumber == 136 
             || gridNumber == 24 || 255 - gridNumber == 24 || gridNumber == 12 || 255 - gridNumber == 12 || gridNumber == 142
             || gridNumber == 10 || 255 - gridNumber == 10 || gridNumber == 30 || gridNumber == 156 || gridNumber == 154) {
                triangles.push_back({edges[3], edges[5], edges[9]}); 
            } 
            if (gridNumber == 4 || 255 - gridNumber == 4 || gridNumber == 164 || 255 - gridNumber == 164 
             || gridNumber == 44 || 255 - gridNumber == 44 || gridNumber == 38 || 255 - gridNumber == 38 
             || gridNumber == 133 || 255 - gridNumber == 133 || gridNumber == 131 || 255 - gridNumber == 131 
             || gridNumber == 133 || 255 - gridNumber == 133 || gridNumber == 132 || 255 - gridNumber == 132
             || gridNumber == 36 || 255 - gridNumber == 36 || gridNumber == 12 || 255 - gridNumber == 12 || gridNumber == 142 
             || gridNumber == 6 || 255 - gridNumber == 6 || gridNumber == 166 || gridNumber == 172 || gridNumber == 46) {
                triangles.push_back({edges[4], edges[7], edges[10]}); 
            } 
            if (gridNumber == 2 || 255 - gridNumber == 2 || gridNumber == 194 || 255 - gridNumber == 194 
             || gridNumber == 74 || 255 - gridNumber == 74 || gridNumber == 70 || 255 - gridNumber == 70 
             || gridNumber == 134 || 255 - gridNumber == 134 || gridNumber == 138 || 255 - gridNumber == 138 
             || gridNumber == 14 || 255 - gridNumber == 14 || gridNumber == 130 || 255 - gridNumber == 130 
             || gridNumber == 66 || 255 - gridNumber == 66 || gridNumber == 10 || 255 - gridNumber == 10 || gridNumber == 142 
             || gridNumber == 6 || 255 - gridNumber == 6 || gridNumber == 202 || gridNumber == 198 || gridNumber == 78) {
                triangles.push_back({edges[6], edges[8], edges[11]}); 
            } 
            if (gridNumber == 1 || 255 - gridNumber == 1 || gridNumber == 193 || 255 - gridNumber == 193 
             || gridNumber == 161 || 255 - gridNumber == 161 || gridNumber == 145 || 255 - gridNumber == 145 
             || gridNumber == 81 || 255 - gridNumber == 81 || gridNumber == 49 || 255 - gridNumber == 49 
             || gridNumber == 97 || 255 - gridNumber == 97 || gridNumber == 129 || 255 - gridNumber == 129 
             || gridNumber == 65 || 255 - gridNumber == 65 || gridNumber == 33 || 255 - gridNumber == 33 || gridNumber == 113 
             || gridNumber == 17 || 255 - gridNumber == 17 || gridNumber == 225 || gridNumber == 209 || gridNumber == 177) {
                triangles.push_back({edges[9], edges[10], edges[11]}); 
            } 

            // handle case of adjacent two nodes
            if (gridNumber == 192 || 255 - gridNumber == 192 || gridNumber == 194 || 255 - gridNumber == 194
                || gridNumber == 193 || 255 - gridNumber == 193 || gridNumber == 195) {
                triangles.push_back({edges[1], edges[2], edges[4]});
                triangles.push_back({edges[1], edges[3], edges[4]});
            } 
            if (gridNumber == 160 || 255 - gridNumber == 160 || gridNumber == 164 || 255 - gridNumber == 164
                || gridNumber == 161 || 255 - gridNumber == 161 || gridNumber == 165) {
                triangles.push_back({edges[0], edges[2], edges[5]});
                triangles.push_back({edges[2], edges[5], edges[6]});
            } 
            if (gridNumber == 144 || 255 - gridNumber == 144 || gridNumber == 152 || 255 - gridNumber == 152
                || gridNumber == 145 || 255 - gridNumber == 145 || gridNumber == 153) {
                triangles.push_back({edges[0], edges[1], edges[8]});
                triangles.push_back({edges[0], edges[7], edges[8]});
            } 
            if (gridNumber == 72 || 255 - gridNumber == 72 || gridNumber == 88 || 255 - gridNumber == 88
                || gridNumber == 74 || 255 - gridNumber == 74 || gridNumber == 90) {
                triangles.push_back({edges[0], edges[4], edges[5]});
                triangles.push_back({edges[4], edges[5], edges[9]});
            } 
            if (gridNumber == 68 || 255 - gridNumber == 68 || gridNumber == 100 || 255 - gridNumber == 100
                || gridNumber == 70 || 255 - gridNumber == 70 || gridNumber == 102) {
                triangles.push_back({edges[0], edges[3], edges[7]});
                triangles.push_back({edges[3], edges[7], edges[10]});
            } 
            if (gridNumber == 40 || 255 - gridNumber == 40 || gridNumber == 56 || 255 - gridNumber == 56
                || gridNumber == 44 || 255 - gridNumber == 44 || gridNumber == 60) {
                triangles.push_back({edges[1], edges[3], edges[6]});
                triangles.push_back({edges[3], edges[6], edges[9]});
            } 
            if (gridNumber == 34 || 255 - gridNumber == 34 || gridNumber == 98 || 255 - gridNumber == 98
                || gridNumber == 38 || 255 - gridNumber == 38 || gridNumber == 102) {
                triangles.push_back({edges[1], edges[5], edges[8]});
                triangles.push_back({edges[5], edges[8], edges[11]});
            }
            if (gridNumber == 20 || 255 - gridNumber == 20 || gridNumber == 52 || 255 - gridNumber == 52
                || gridNumber == 28 || 255 - gridNumber == 28 || gridNumber == 60) {
                triangles.push_back({edges[2], edges[4], edges[8]});
                triangles.push_back({edges[4], edges[8], edges[10]});
            } 
            if (gridNumber == 18 || 255 - gridNumber == 18 || gridNumber == 82 || 255 - gridNumber == 82
                || gridNumber == 26 || 255 - gridNumber == 26 || gridNumber == 90) {
                triangles.push_back({edges[2], edges[6], edges[7]});
                triangles.push_back({edges[6], edges[7], edges[11]});
            } 
            if (gridNumber == 9 || 255 - gridNumber == 9 || gridNumber == 137 || 255 - gridNumber == 137
                || gridNumber == 25 || 255 - gridNumber == 25 || gridNumber == 153) {
                triangles.push_back({edges[3], edges[5], edges[10]});
                triangles.push_back({edges[5], edges[10], edges[11]});
            }
            if (gridNumber == 5 || 255 - gridNumber == 5 || gridNumber == 133 || 255 - gridNumber == 133
                || gridNumber == 37 || 255 - gridNumber == 37 || gridNumber == 165) {
                triangles.push_back({edges[4], edges[7], edges[9]});
                triangles.push_back({edges[7], edges[9], edges[11]});
            } else if (gridNumber == 3 || 255 - gridNumber == 3 || gridNumber == 131 || 255 - gridNumber == 131
                || gridNumber == 67 || 255 - gridNumber == 67 || gridNumber == 195) {
                triangles.push_back({edges[6], edges[8], edges[9]});
                triangles.push_back({edges[8], edges[9], edges[10]});
            } 

            // handle case with three adjacent nodes
            if (gridNumber == 224 || 255 - gridNumber == 224 || gridNumber == 225) {
                triangles.push_back({edges[2], edges[4], edges[6]});
                triangles.push_back({edges[3], edges[5], edges[6]});
                triangles.push_back({edges[3], edges[4], edges[6]});
            } else if (gridNumber == 162 || 255 - gridNumber == 162 || gridNumber == 166) {
                triangles.push_back({edges[0], edges[5], edges[11]});
                triangles.push_back({edges[0], edges[2], edges[8]});
                triangles.push_back({edges[0], edges[8], edges[11]});
            } else if (gridNumber == 168 || 255 - gridNumber == 168 || gridNumber == 172) {
                triangles.push_back({edges[2], edges[6], edges[9]});
                triangles.push_back({edges[0], edges[2], edges[3]});
                triangles.push_back({edges[2], edges[3], edges[9]});
            } else if (gridNumber == 208 || 255 - gridNumber == 208 || gridNumber == 209) {
                triangles.push_back({edges[1], edges[3], edges[8]});
                triangles.push_back({edges[3], edges[4], edges[7]});
                triangles.push_back({edges[3], edges[7], edges[8]});
            } else if (gridNumber == 176 || 255 - gridNumber == 176 || gridNumber == 177) {
                triangles.push_back({edges[0], edges[5], edges[7]});
                triangles.push_back({edges[5], edges[6], edges[7]});
                triangles.push_back({edges[6], edges[7], edges[8]});
            } else if (gridNumber == 200 || 255 - gridNumber == 200 || gridNumber == 202) {
                triangles.push_back({edges[2], edges[4], edges[9]});
                triangles.push_back({edges[1], edges[2], edges[5]});
                triangles.push_back({edges[2], edges[5], edges[9]});
            } else if (gridNumber == 196 || 255 - gridNumber == 196 || gridNumber == 198) {
                triangles.push_back({edges[1], edges[3], edges[10]});
                triangles.push_back({edges[1], edges[2], edges[7]});
                triangles.push_back({edges[1], edges[7], edges[10]});
            } else if (gridNumber == 22 || 255 - gridNumber == 22 || gridNumber == 30) {
                triangles.push_back({edges[2], edges[4], edges[6]});
                triangles.push_back({edges[4], edges[6], edges[10]});
                triangles.push_back({edges[6], edges[10], edges[11]});
            } else if (gridNumber == 50 || 255 - gridNumber == 50 || gridNumber == 114) {
                triangles.push_back({edges[5], edges[7], edges[11]});
                triangles.push_back({edges[1], edges[2], edges[5]});
                triangles.push_back({edges[2], edges[5], edges[7]});
            } else if (gridNumber == 35 || 255 - gridNumber == 35 || gridNumber == 99) {
                triangles.push_back({edges[1], edges[8], edges[10]});
                triangles.push_back({edges[1], edges[5], edges[9]});
                triangles.push_back({edges[1], edges[9], edges[10]});
            } else if (gridNumber == 21 || 255 - gridNumber == 21 || gridNumber == 53) {
                triangles.push_back({edges[2], edges[4], edges[9]});
                triangles.push_back({edges[2], edges[8], edges[9]});
                triangles.push_back({edges[8], edges[9], edges[11]});
            } else if (gridNumber == 19 || 255 - gridNumber == 19 || gridNumber == 83) {
                triangles.push_back({edges[2], edges[6], edges[9]});
                triangles.push_back({edges[2], edges[7], edges[9]});
                triangles.push_back({edges[7], edges[9], edges[10]});
            } else if (gridNumber == 13 || 255 - gridNumber == 13 || gridNumber == 141) {
                triangles.push_back({edges[5], edges[7], edges[11]});
                triangles.push_back({edges[3], edges[4], edges[5]});
                triangles.push_back({edges[4], edges[5], edges[7]});
            } else if (gridNumber == 11 || 255 - gridNumber == 11 || gridNumber == 139) {
                triangles.push_back({edges[3], edges[8], edges[10]});
                triangles.push_back({edges[3], edges[5], edges[6]});
                triangles.push_back({edges[3], edges[6], edges[8]});
            } else if (gridNumber == 7 || 255 - gridNumber == 7 || gridNumber == 135) {
                triangles.push_back({edges[4], edges[6], edges[9]});
                triangles.push_back({edges[4], edges[6], edges[7]});
                triangles.push_back({edges[6], edges[7], edges[8]});
            } else if (gridNumber == 42 || 255 - gridNumber == 42 || gridNumber == 46) {
                triangles.push_back({edges[1], edges[3], edges[8]});
                triangles.push_back({edges[3], edges[8], edges[9]});
                triangles.push_back({edges[8], edges[9], edges[11]});
            } else if (gridNumber == 148 || 255 - gridNumber == 148 || gridNumber == 156) {
                triangles.push_back({edges[1], edges[8], edges[10]});
                triangles.push_back({edges[0], edges[1], edges[4]});
                triangles.push_back({edges[1], edges[4], edges[10]});
            } else if (gridNumber == 146 || 255 - gridNumber == 146 || gridNumber == 154) {
                triangles.push_back({edges[0], edges[7], edges[11]});
                triangles.push_back({edges[0], edges[1], edges[6]});
                triangles.push_back({edges[0], edges[6], edges[11]});
            } else if (gridNumber == 104 || 255 - gridNumber == 104 || gridNumber == 120) {
                triangles.push_back({edges[4], edges[6], edges[9]});
                triangles.push_back({edges[0], edges[1], edges[4]});
                triangles.push_back({edges[1], edges[4], edges[6]});
            } else if (gridNumber == 73 || 255 - gridNumber == 73 || gridNumber == 89) {
                triangles.push_back({edges[0], edges[5], edges[11]});
                triangles.push_back({edges[0], edges[4], edges[10]});
                triangles.push_back({edges[0], edges[10], edges[11]});
            } else if (gridNumber == 84 || 255 - gridNumber == 84 || gridNumber == 116) {
                triangles.push_back({edges[3], edges[8], edges[10]});
                triangles.push_back({edges[0], edges[2], edges[3]});
                triangles.push_back({edges[2], edges[3], edges[8]});
            } else if (gridNumber == 69 || 255 - gridNumber == 69 || gridNumber == 101) {
                triangles.push_back({edges[0], edges[7], edges[11]});
                triangles.push_back({edges[0], edges[3], edges[9]});
                triangles.push_back({edges[0], edges[9], edges[11]});
            } else if (gridNumber == 41 || 255 - gridNumber == 41 || gridNumber == 57) {
                triangles.push_back({edges[1], edges[3], edges[10]});
                triangles.push_back({edges[1], edges[6], edges[10]});
                triangles.push_back({edges[6], edges[10], edges[11]});
            } else if (gridNumber == 76 || 255 - gridNumber == 76 || gridNumber == 78) {
                triangles.push_back({edges[0], edges[5], edges[7]});
                triangles.push_back({edges[5], edges[7], edges[9]});
                triangles.push_back({edges[7], edges[9], edges[10]});
            }

            // handle flat sheet 4 nodes
            if (gridNumber == 212 || gridNumber == 43) {
                triangles.push_back({edges[1], edges[3], edges[8]});
                triangles.push_back({edges[3], edges[8], edges[10]});
            }
            if (gridNumber == 232 || gridNumber == 23) {
                triangles.push_back({edges[2], edges[4], edges[6]});
                triangles.push_back({edges[4], edges[6], edges[9]});
            }
            if (gridNumber == 178 || gridNumber == 77) {
                triangles.push_back({edges[0], edges[5], edges[7]});
                triangles.push_back({edges[5], edges[7], edges[11]});
            }

            // super corner 4 nodes
            if (gridNumber == 240 || gridNumber == 15) {
                triangles.push_back({edges[3], edges[4], edges[7]});
                triangles.push_back({edges[3], edges[5], edges[7]});
                triangles.push_back({edges[5], edges[7], edges[8]});
                triangles.push_back({edges[5], edges[6], edges[8]});
            } else if (gridNumber == 204 || gridNumber == 51) {
                triangles.push_back({edges[1], edges[2], edges[5]});
                triangles.push_back({edges[2], edges[5], edges[7]});
                triangles.push_back({edges[5], edges[7], edges[9]});
                triangles.push_back({edges[7], edges[9], edges[10]});
            } else if (gridNumber == 170 || gridNumber == 85) {
                triangles.push_back({edges[0], edges[2], edges[3]});
                triangles.push_back({edges[2], edges[3], edges[8]});
                triangles.push_back({edges[3], edges[8], edges[9]});
                triangles.push_back({edges[8], edges[9], edges[11]});
            } else if (gridNumber == 150 || gridNumber == 105) {
                triangles.push_back({edges[1], edges[6], edges[11]});
                triangles.push_back({edges[0], edges[1], edges[10]});
                triangles.push_back({edges[1], edges[10], edges[11]});
                triangles.push_back({edges[0], edges[4], edges[10]});
            }

            // line of 4 nodes 
            if (gridNumber == 201 || gridNumber == 54) {
                triangles.push_back({edges[1], edges[2], edges[4]});
                triangles.push_back({edges[1], edges[5], edges[11]});
                triangles.push_back({edges[1], edges[4], edges[11]});
                triangles.push_back({edges[4], edges[10], edges[11]});
            } else if (gridNumber == 197 || gridNumber == 58) {
                triangles.push_back({edges[1], edges[2], edges[3]});
                triangles.push_back({edges[2], edges[7], edges[11]});
                triangles.push_back({edges[2], edges[3], edges[11]});
                triangles.push_back({edges[3], edges[9], edges[11]});
            } else if (gridNumber == 149 || gridNumber == 106) {
                triangles.push_back({edges[0], edges[1], edges[8]});
                triangles.push_back({edges[8], edges[9], edges[11]});
                triangles.push_back({edges[0], edges[8], edges[9]});
                triangles.push_back({edges[0], edges[4], edges[9]});
            } else if (gridNumber == 147 || gridNumber == 108) {
                triangles.push_back({edges[0], edges[1], edges[7]});
                triangles.push_back({edges[1], edges[6], edges[9]});
                triangles.push_back({edges[7], edges[9], edges[10]});
                triangles.push_back({edges[1], edges[7], edges[9]});
            } else if (gridNumber == 169 || gridNumber == 86) {
                triangles.push_back({edges[0], edges[2], edges[6]});
                triangles.push_back({edges[6], edges[10], edges[11]});
                triangles.push_back({edges[0], edges[3], edges[10]});
                triangles.push_back({edges[0], edges[6], edges[10]});
            } else if (gridNumber == 163 || gridNumber == 92) {
                triangles.push_back({edges[0], edges[2], edges[5]});
                triangles.push_back({edges[2], edges[8], edges[10]});
                triangles.push_back({edges[5], edges[9], edges[10]});
                triangles.push_back({edges[2], edges[5], edges[10]});
            } else if (gridNumber == 226 || gridNumber == 29) {
                triangles.push_back({edges[2], edges[3], edges[4]});
                triangles.push_back({edges[2], edges[5], edges[11]});
                triangles.push_back({edges[3], edges[5], edges[11]});
                triangles.push_back({edges[2], edges[3], edges[11]});
            } else if (gridNumber == 210 || gridNumber == 45) {
                triangles.push_back({edges[1], edges[6], edges[11]});
                triangles.push_back({edges[1], edges[3], edges[4]});
                triangles.push_back({edges[4], edges[7], edges[11]});
                triangles.push_back({edges[1], edges[4], edges[11]});
            } else if (gridNumber == 75 || gridNumber == 180) {
                triangles.push_back({edges[4], edges[5], edges[8]});
                triangles.push_back({edges[5], edges[6], edges[8]});
                triangles.push_back({edges[4], edges[8], edges[10]});
                triangles.push_back({edges[0], edges[4], edges[5]});
            } else if (gridNumber == 71 || gridNumber == 184) {
                triangles.push_back({edges[0], edges[3], edges[7]});
                triangles.push_back({edges[6], edges[7], edges[8]});
                triangles.push_back({edges[3], edges[6], edges[9]});
                triangles.push_back({edges[3], edges[6], edges[7]});
            } else if (gridNumber == 228 || gridNumber == 27) {
                triangles.push_back({edges[3], edges[5], edges[10]});
                triangles.push_back({edges[2], edges[5], edges[6]});
                triangles.push_back({edges[2], edges[7], edges[10]});
                triangles.push_back({edges[2], edges[5], edges[10]});
            } else if (gridNumber == 39 || gridNumber == 216) {
                triangles.push_back({edges[4], edges[5], edges[9]});
                triangles.push_back({edges[1], edges[5], edges[8]});
                triangles.push_back({edges[4], edges[7], edges[8]});
                triangles.push_back({edges[4], edges[5], edges[8]});
            }

            if (triangles.size() == 0) {
                std::cout << gridNumber << std::endl;
            }

            return triangles;

        }

        void UpdateNormals() {
            IndexArray indices = vtx_obj_ptr->GetIndices();
            PositionArray vertices = vtx_obj_ptr->GetPositions();
            NormalArray old_normals = vtx_obj_ptr->GetNormals();
            auto normals = make_unique<NormalArray>();

            std::vector<glm::vec3> weighted_normals(vertices.size(), glm::vec3(0));

            for (int i = 0; i < indices.size(); i = i+3) {
                glm::vec3 position1 = vertices[indices[i]];
                glm::vec3 position2 = vertices[indices[i+1]];
                glm::vec3 position3 = vertices[indices[i+2]];

                glm::vec3 edge1 = position2 - position1;
                glm::vec3 edge2 = position3 - position1;
                glm::vec3 normal = old_normals[i/3];
                // glm::vec3 normal = glm::cross(edge2, edge1);

                weighted_normals[indices[i]] += normal;
                weighted_normals[indices[i+1]] += normal;
                weighted_normals[indices[i+2]] += normal;
            }
            for (int k = 0; k < weighted_normals.size(); k++) {
                normals->push_back(glm::normalize(weighted_normals[k]));
            }
            vtx_obj_ptr->UpdateNormals(std::move(normals));
        }  

        bool keyPressedE, keyPressedR;
        float cellWidth;
        Triple gridSize;
        std::unique_ptr<FluidSystem> fluidSystem;
        float current_time = 0.0f;
        float max_time_step = 0.01f;
        std::shared_ptr<VertexObject> vtx_obj_ptr;
        
};

}

#endif