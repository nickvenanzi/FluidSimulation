#include "FluidSystem.hpp"
#include "stdlib.h"

namespace GLOO {

void FluidSystem::AdvectPhi(float time_step) {
    float epsilon = 0.00001f;
    auto newAirCells = std::make_shared<std::vector<Triple>>();
    auto newLiquidCells = std::make_shared<std::vector<Triple>>();

    for (int i = 0; i < grid_ptr->gridWidth; i++) {
        for (int j = 0; j < grid_ptr->gridHeight; j++) {
            for (int k = 0; k < grid_ptr->gridDepth; k++) {
                if (grid_ptr->GetCellType({i,j,k}) == CellType::SOLID) {
                    grid_ptr->SetPhi({i,j,k}, MAXFLOAT);
                    continue;
                } else if (abs(grid_ptr->GetPhi({i,j,k})) > 15.f * grid_ptr->cellWidth) {
                    grid_ptr->SetPhi({i,j,k}, grid_ptr->GetPhi({i,j,k}));
                    continue;
                }

                // now on a cell near the surface

                // if liquid, advect as usual
                // if air, first obtain interpolated U,V,W from nearest surface point.
                float u, v, w;
                
                glm::vec3 curr_position = glm::vec3(i, j, k) * grid_ptr->cellWidth
                                        + grid_ptr->gridPosition;
                float phi_old = grid_ptr->GetPhi({i,j,k});
                int current_i = i, current_j = j, current_k = k;
                // std::cout << "ITERATING: " << i << ", " << j << ", " << k << std::endl;
                int iter = 0;
                while (grid_ptr->GetCellType({current_i, current_j, current_k}) != CellType::LIQUID && iter < 20) {
                    iter++;
                    Triple best_neighbor;
                    float best_phi = phi_old;

                    if (grid_ptr->GetCellType({current_i+1,current_j,current_k}) != CellType::SOLID && grid_ptr->GetPhi({current_i+1,current_j,current_k}) < best_phi) { 
                        best_phi = grid_ptr->GetPhi({current_i+1,current_j,current_k}); 
                        best_neighbor = {current_i+1,current_j,current_k};
                    } 

                    if (grid_ptr->GetCellType({current_i-1,current_j,current_k}) != CellType::SOLID && grid_ptr->GetPhi({current_i-1,current_j,current_k}) < best_phi) { 
                        best_phi = grid_ptr->GetPhi({current_i-1,current_j,current_k}); 
                        best_neighbor = {current_i-1,current_j,current_k};
                    } 

                    if (grid_ptr->GetCellType({current_i,current_j+1,current_k}) != CellType::SOLID && grid_ptr->GetPhi({current_i,current_j+1,current_k}) < best_phi) { 
                        best_phi = grid_ptr->GetPhi({current_i,current_j+1,current_k}); 
                        best_neighbor = {current_i,current_j+1,current_k};
                    } 

                    if (grid_ptr->GetCellType({current_i,current_j-1,current_k}) != CellType::SOLID && grid_ptr->GetPhi({current_i,current_j-1,current_k}) < best_phi) { 
                        best_phi = grid_ptr->GetPhi({current_i,current_j-1,current_k}); 
                        best_neighbor = {current_i,current_j-1,current_k};
                    } 

                    if (grid_ptr->GetCellType({current_i,current_j,current_k+1}) != CellType::SOLID && grid_ptr->GetPhi({current_i,current_j,current_k+1}) < best_phi) { 
                        best_phi = grid_ptr->GetPhi({current_i,current_j,current_k+1}); 
                        best_neighbor = {current_i,current_j,current_k+1};
                    } 

                    if (grid_ptr->GetCellType({current_i,current_j,current_k-1}) != CellType::SOLID && grid_ptr->GetPhi({current_i,current_j,current_k-1}) < best_phi) { 
                        best_phi = grid_ptr->GetPhi({current_i,current_j,current_k-1}); 
                        best_neighbor = {current_i,current_j,current_k-1};
                    } 

                    
                    current_i = best_neighbor.i;
                    current_j = best_neighbor.j;
                    current_k = best_neighbor.k;
                    if (best_phi > phi_old - epsilon) {
                        iter = 20;
                        continue;
                    }
                    // if (best_phi < 0.f && grid_ptr->GetCellType({current_i, current_j, current_k}) != CellType::LIQUID) {
                    //     std::cout << "WAHT THE FUCK" << std::endl;
                    // }
                    // std::cout << "  Moving to: " << current_i << ", " << current_j << ", " << current_k << " -- with phi: " << best_phi << std::endl;

                } 
                if (iter >= 20) {
                    grid_ptr->SetPhi({i,j,k}, grid_ptr->GetPhi({i,j,k}));
                    continue;
                }
                glm::vec3 extended_position = glm::vec3(current_i, current_j, current_k) * grid_ptr->cellWidth + grid_ptr->gridPosition;
                u = grid_ptr->InterpolateU(extended_position);
                v = grid_ptr->InterpolateV(extended_position);
                w = grid_ptr->InterpolateW(extended_position);

                float velocity_mag = sqrt(u*u + v*v + w*w + epsilon);
                int max_cells_traversed = (int)((velocity_mag * time_step)/grid_ptr->cellWidth) + 1;

                glm::vec3 interpolated_position = curr_position - time_step * glm::vec3(u, v, w);
                
                interpolated_position = grid_ptr->BringInboundsP(interpolated_position);
                Triple interpolated_coords = grid_ptr->GetCellCoordinates(interpolated_position);
                for (int h = 0; h < max_cells_traversed ; h++) {
                    if (grid_ptr->GetCellType(interpolated_coords) != CellType::SOLID) {
                        break;
                    }
                    interpolated_position += glm::normalize(glm::vec3(u + epsilon,v + epsilon,w + epsilon)) * grid_ptr->cellWidth;
                    interpolated_position = grid_ptr->BringInboundsP(interpolated_position);
                    interpolated_coords = grid_ptr->GetCellCoordinates(interpolated_position);
                }
                float phi_new = grid_ptr->InterpolatePhi(interpolated_position);
 
                //once we are in a liquid cell, update old U with new U
                grid_ptr->SetPhi({i,j,k}, phi_new);

                if (abs(phi_new) < 15.f && phi_new < 0.f) {
                    newLiquidCells->push_back({i,j,k});
                } else if (abs(phi_new) < 15.f && phi_new > 0.f) {
                    newAirCells->push_back({i,j,k});
                }

                if (phi_new > 0.f) {
                    grid_ptr->SetCellType({i,j,k}, CellType::AIR);
                } else {
                    grid_ptr->SetCellType({i,j,k}, CellType::LIQUID);
                }

                if (phi_old > 0.f && phi_new < 0.f) {
                    grid_ptr->FlipStorage();
                    grid_ptr->SetUMinus({i,j,k}, u);
                    grid_ptr->SetUMinus({i+1,j,k}, u);
                    grid_ptr->SetVMinus({i,j,k}, v);
                    grid_ptr->SetVMinus({i,j+1,k}, v);
                    grid_ptr->SetWMinus({i,j,k}, w);
                    grid_ptr->SetWMinus({i,j,k+1}, w);
                    grid_ptr->SetPressure({i,j,k}, grid_ptr->GetPressure(interpolated_coords));
                    grid_ptr->FlipStorage();
                }

            }
        }
    }

    grid_ptr->airCellsAtSurface = newAirCells;
    grid_ptr->liquidCellsAtSurface = newLiquidCells;
    grid_ptr->FlipPhiStorage();
}

void FluidSystem::AdvectVelocity(float time_step) {
    float epsilon = 0.00001f;
    float delta = grid_ptr->cellWidth;
    //update U's
    for (int i = 0; i <= grid_ptr->gridWidth; i++) {
        for (int j = 0; j < grid_ptr->gridHeight; j++) {
            for (int k = 0; k < grid_ptr->gridDepth; k++) {
                Triple curr_coords = {i, j, k};

                if (grid_ptr->GetCellType(curr_coords) != CellType::LIQUID && grid_ptr->GetCellType({i-1,j,k}) != CellType::LIQUID) {
                    float u_old = grid_ptr->GetUMinus(curr_coords);
                    grid_ptr->SetUMinus(curr_coords, u_old);
                    continue;
                }

                float u_old = grid_ptr->GetUMinus(curr_coords);
                // std::cout << "At coords: " << i << " - 1/2, " << j << ", " << k << std::endl;
                // std::cout << "U: " << u_old << std::endl;
                // std::cout << std::endl;
                
                glm::vec3 curr_position = (glm::vec3(-0.5f, 0.f, 0.f) 
                                        + glm::vec3(i, j, k)) * grid_ptr->cellWidth
                                        + grid_ptr->gridPosition;
                float v_old = grid_ptr->InterpolateV(curr_position);
                float w_old = grid_ptr->InterpolateW(curr_position);

                float velocity_mag = sqrt(u_old*u_old + v_old*v_old + w_old*w_old + epsilon);
                int max_cells_traversed = (int)((velocity_mag * time_step)/grid_ptr->cellWidth) + 1;

                glm::vec3 interpolated_position = curr_position - time_step * glm::vec3(u_old, v_old, w_old);
                
                interpolated_position = grid_ptr->BringInboundsU(interpolated_position);
                Triple interpolated_coords = grid_ptr->GetCellCoordinates(interpolated_position);
                for (int h = 0; h < max_cells_traversed ; h++) {
                    if (grid_ptr->GetCellType(interpolated_coords) == CellType::LIQUID) {
                        break;
                    }
                    interpolated_position += glm::normalize(glm::vec3(u_old + epsilon,v_old + epsilon,w_old + epsilon)) * grid_ptr->cellWidth;
                    interpolated_position = grid_ptr->BringInboundsU(interpolated_position);
                    interpolated_coords = grid_ptr->GetCellCoordinates(interpolated_position);
                }
                float u_new = grid_ptr->InterpolateU(interpolated_position);
                if (isnan(u_new)) {
                    std::cout << "NAN found for x at: " << interpolated_position.x << ", " << interpolated_position.y << ", " << interpolated_position.z << std::endl;
                }
                //once we are in a liquid cell, update old U with new U
                grid_ptr->SetUMinus(curr_coords, u_new);
            }
        }
    }

    //update V's
    for (int i = 0; i < grid_ptr->gridWidth; i++) {
        for (int j = 0; j <= grid_ptr->gridHeight; j++) {
            for (int k = 0; k < grid_ptr->gridDepth; k++) {
                Triple curr_coords = {i, j, k};

                if (grid_ptr->GetCellType(curr_coords) != CellType::LIQUID && grid_ptr->GetCellType({i,j-1,k}) != CellType::LIQUID) {
                    float v_old = grid_ptr->GetVMinus(curr_coords);
                    grid_ptr->SetVMinus(curr_coords, v_old);
                    continue;
                }

                float v_old = grid_ptr->GetVMinus(curr_coords);

                // std::cout << "At coords: " << i << ", " << j << " - 1/2, " << k << std::endl;
                // std::cout << "V: " << v_old << std::endl;
                // std::cout << std::endl;

                glm::vec3 curr_position = (glm::vec3(0.f, -0.5f, 0.f) 
                                        + glm::vec3(i, j, k)) * glm::vec3(grid_ptr->cellWidth)
                                        + grid_ptr->gridPosition;
                float u_old = grid_ptr->InterpolateU(curr_position);
                float w_old = grid_ptr->InterpolateW(curr_position);

                float velocity_mag = sqrt(u_old*u_old + v_old*v_old + w_old*w_old + epsilon);
                int max_cells_traversed = (int)((velocity_mag * time_step)/grid_ptr->cellWidth) + 1;
                // if (i == 12 && k == 12) {
                //     std::cout << "Currently advecting V, at 12, " << j << ", 12" << std::endl;
                //     std::cout << "  initial v: " << v_old << std::endl;

                // }
                // std::cout << "Currently advecting V, at 12, " << j << ", 12" << std::endl;
                glm::vec3 interpolated_position = curr_position - time_step * glm::vec3(u_old, v_old, w_old);
                // if (i == 12 && k == 12) {
                //     std::cout << "  New interpolated position (pre bounding) at " << interpolated_position.x << ", " << interpolated_position.y << ", " << interpolated_position.z<< std::endl;

                //     }
                interpolated_position = grid_ptr->BringInboundsV(interpolated_position);
                Triple interpolated_coords = grid_ptr->GetCellCoordinates(interpolated_position);
                // if (i == 12 && k == 12) {
                //     std::cout << "  New interpolated coords at " << interpolated_coords.i << ", " << interpolated_coords.j << ", " << interpolated_coords.k<< std::endl;
                //     std::cout << "  New interpolated position at " << interpolated_position.x << ", " << interpolated_position.y << ", " << interpolated_position.z<< std::endl;

                //     }
                for (int h = 0; h < max_cells_traversed ; h++) {
                    if (grid_ptr->GetCellType(interpolated_coords) == CellType::LIQUID) {
                        break;
                    }
                    interpolated_position += glm::normalize(glm::vec3(u_old + epsilon,v_old + epsilon,w_old + epsilon)) * grid_ptr->cellWidth;
                    interpolated_position = grid_ptr->BringInboundsV(interpolated_position);
                    interpolated_coords = grid_ptr->GetCellCoordinates(interpolated_position);
                    // if (i == 12 && k == 12) {
                    // std::cout << "  New interpolated coords at " << interpolated_coords.i << ", " << interpolated_coords.j << ", " << interpolated_coords.k<< std::endl;
                    // std::cout << "  New interpolated position at " << interpolated_position.x << ", " << interpolated_position.y << ", " << interpolated_position.z<< std::endl;

                    // }
                } 
                // if (i == 12 && k == 12) {

                //     std::cout << "  interpolating from V, at " << interpolated_coords.i << ", " << interpolated_coords.j << ", " << interpolated_coords.k<< std::endl;
                // }
                float v_new = grid_ptr->InterpolateV(interpolated_position);
                if (isnan(v_new)) {
                    std::cout << "NAN found for y at: " << interpolated_position.x << ", " << interpolated_position.y << ", " << interpolated_position.z << std::endl;
                }
                //once we are in a liquid cell, update old V with new V
                grid_ptr->SetVMinus(curr_coords, v_new);
            }
        }
    }

    //update W's
    for (int i = 0; i < grid_ptr->gridWidth; i++) {
        for (int j = 0; j < grid_ptr->gridHeight; j++) {
            for (int k = 0; k <= grid_ptr->gridDepth; k++) {
                Triple curr_coords = {i, j, k};

                if (grid_ptr->GetCellType(curr_coords) != CellType::LIQUID && grid_ptr->GetCellType({i,j,k-1}) != CellType::LIQUID) {
                    float w_old = grid_ptr->GetWMinus(curr_coords);
                    grid_ptr->SetWMinus(curr_coords, w_old);
                    continue;
                }

                float w_old = grid_ptr->GetWMinus(curr_coords);

                // std::cout << "At coords: " << i << ", " << j << ", " << k << " - 1/2" << std::endl;
                // std::cout << "W: " << w_old << std::endl;
                // std::cout << std::endl;

                glm::vec3 curr_position = (glm::vec3(0.f, 0.f, -0.5f) 
                                        + glm::vec3(i, j, k)) * glm::vec3(grid_ptr->cellWidth)
                                        + grid_ptr->gridPosition;
                float v_old = grid_ptr->InterpolateV(curr_position);
                float u_old = grid_ptr->InterpolateU(curr_position);

                float velocity_mag = sqrt(u_old*u_old + v_old*v_old + w_old*w_old + epsilon);
                int max_cells_traversed = (int)((velocity_mag * time_step)/grid_ptr->cellWidth) + 1;

                glm::vec3 interpolated_position = curr_position - time_step * glm::vec3(u_old, v_old, w_old);
                
                interpolated_position = grid_ptr->BringInboundsW(interpolated_position);
                Triple interpolated_coords = grid_ptr->GetCellCoordinates(interpolated_position);
                for (int h = 0; h < max_cells_traversed ; h++) {
                    if (grid_ptr->GetCellType(interpolated_coords) == CellType::LIQUID) {
                        break;
                    }
                    interpolated_position += glm::normalize(glm::vec3(u_old + epsilon, v_old + epsilon, w_old + epsilon)) * grid_ptr->cellWidth;
                    interpolated_position = grid_ptr->BringInboundsW(interpolated_position);
                    interpolated_coords = grid_ptr->GetCellCoordinates(interpolated_position);
                }
                float w_new = grid_ptr->InterpolateW(interpolated_position);
                if (isnan(w_new)) {
                    std::cout << "NAN found for y at: " << interpolated_position.x << ", " << interpolated_position.y << ", " << interpolated_position.z << std::endl;
                }
                //once we are in a liquid cell, update old W with new W
                grid_ptr->SetWMinus(curr_coords, w_new);
            }
        }
    }
    grid_ptr->FlipStorage();
}

void FluidSystem::UpdateBodyForces(float time_step) {
    //update U's
    for (int i = 0; i <= grid_ptr->gridWidth; i++) {
        for (int j = 0; j < grid_ptr->gridHeight; j++) {
            for (int k = 0; k < grid_ptr->gridDepth; k++) {
                Triple curr_coords = {i, j, k};
                //only advect if currently in liquid
                if (grid_ptr->GetCellType(curr_coords) != CellType::LIQUID && grid_ptr->GetCellType({i-1,j,k}) != CellType::LIQUID) {
                    float u_old = grid_ptr->GetUMinus(curr_coords);
                    grid_ptr->SetUMinus(curr_coords, u_old);
                    continue;
                }
                float u_old = grid_ptr->GetUMinus(curr_coords);         

                grid_ptr->SetUMinus(curr_coords, u_old);
            }
        }
    }
    //update V's
    for (int i = 0; i < grid_ptr->gridWidth; i++) {
        for (int j = 0; j <= grid_ptr->gridHeight; j++) {
            for (int k = 0; k < grid_ptr->gridDepth; k++) {
                Triple curr_coords = {i, j, k};
                //only advect if currently in liquid
                if (grid_ptr->GetCellType(curr_coords) != CellType::LIQUID && grid_ptr->GetCellType({i,j-1,k}) != CellType::LIQUID) {
                    float v_old = grid_ptr->GetVMinus(curr_coords);
                    grid_ptr->SetVMinus(curr_coords, v_old);
                    continue;
                }
                float v_old = grid_ptr->GetVMinus(curr_coords);

                grid_ptr->SetVMinus(curr_coords, v_old + grid_ptr->g * time_step);
            }
        }
    }
    //update W's
    for (int i = 0; i < grid_ptr->gridWidth; i++) {
        for (int j = 0; j < grid_ptr->gridHeight; j++) {
            for (int k = 0; k <= grid_ptr->gridDepth; k++) {
                Triple curr_coords = {i, j, k};
                //only advect if currently in boundary with liquid
                if (grid_ptr->GetCellType(curr_coords) != CellType::LIQUID && grid_ptr->GetCellType({i,j,k-1}) != CellType::LIQUID) {
                    float w_old = grid_ptr->GetWMinus(curr_coords);
                    grid_ptr->SetWMinus(curr_coords, w_old);
                    continue;
                }
                float w_old = grid_ptr->GetWMinus(curr_coords);

                grid_ptr->SetWMinus(curr_coords, w_old);
            }
        }
    }
    grid_ptr->FlipStorage();
}

void FluidSystem::UpdatePressureSOE(float time_step) {
    for (int i = 0; i < grid_ptr->gridWidth; i++) {
        for (int j = 0; j < grid_ptr->gridHeight; j++) {
            for (int k = 0; k < grid_ptr->gridDepth; k++) {
                if (grid_ptr->GetCellType({i,j,k}) != CellType::LIQUID) {
                    continue;
                }
                int non_solid_neighbors = 0;
                float d = 0.f;
                //left neighbor
                if (i > 0) { 
                    CellType neighbor_type = grid_ptr->GetCellType({i - 1, j, k});
                    if (neighbor_type != CellType::SOLID) { 
                        non_solid_neighbors++; 
                        d -= grid_ptr->GetUMinus({i, j, k});
                    }
                    if (neighbor_type == CellType::AIR) {
                        non_solid_neighbors++;
                    }
                    grid_ptr->SetAPlusI({i-1, j, k}, neighbor_type == CellType::LIQUID ? -1 : 0);

                }
                //right neighbor
                if (i < grid_ptr->gridWidth - 1) { 
                    CellType neighbor_type = grid_ptr->GetCellType({i + 1, j, k});
                    if (neighbor_type != CellType::SOLID) { 
                        non_solid_neighbors++; 
                        d += grid_ptr->GetUMinus({i + 1, j, k});
                    }
                    if (neighbor_type == CellType::AIR) {
                        non_solid_neighbors++;
                    }
                    grid_ptr->SetAPlusI({i, j, k}, neighbor_type == CellType::LIQUID ? -1 : 0);
                }

                //bottom neighbor
                if (j > 0) { 
                    CellType neighbor_type = grid_ptr->GetCellType({i, j - 1, k});
                    if (neighbor_type != CellType::SOLID) { 
                        non_solid_neighbors++; 
                        d -= grid_ptr->GetVMinus({i, j, k});
                    }
                    if (neighbor_type == CellType::AIR) {
                        non_solid_neighbors++;
                    }
                    grid_ptr->SetAPlusJ({i, j-1, k}, neighbor_type == CellType::LIQUID ? -1 : 0);
                }
                //top neighbor
                if (j < grid_ptr->gridHeight - 1) { 
                    CellType neighbor_type = grid_ptr->GetCellType({i, j + 1, k});
                    if (neighbor_type != CellType::SOLID) { 
                        non_solid_neighbors++; 
                        d += grid_ptr->GetVMinus({i, j + 1, k});
                    }
                    if (neighbor_type == CellType::AIR) {
                        non_solid_neighbors++;
                    }
                    grid_ptr->SetAPlusJ({i, j, k}, neighbor_type == CellType::LIQUID ? -1 : 0);
                }

                //front neighbor
                if (k > 0) { 
                    CellType neighbor_type = grid_ptr->GetCellType({i, j, k - 1});
                    if (neighbor_type != CellType::SOLID) { 
                        non_solid_neighbors++; 
                        d -= grid_ptr->GetWMinus({i, j, k});
                    }
                    if (neighbor_type == CellType::AIR) {
                        non_solid_neighbors++;
                    }
                    grid_ptr->SetAPlusK({i, j, k-1}, neighbor_type == CellType::LIQUID ? -1 : 0);

                }
                //back neighbor
                if (k < grid_ptr->gridDepth - 1) { 
                    CellType neighbor_type = grid_ptr->GetCellType({i, j, k + 1});
                    if (neighbor_type != CellType::SOLID) { 
                        non_solid_neighbors++;
                        d += grid_ptr->GetWMinus({i, j, k + 1}); 
                    }
                    if (neighbor_type == CellType::AIR) {
                        non_solid_neighbors++;
                    }
                    grid_ptr->SetAPlusK({i, j, k}, neighbor_type == CellType::LIQUID ? -1 : 0);
                }
                // std::cout << "velocity: " << d << std::endl;
                // std::cout << "rho: " << grid_ptr->rho << std::endl;
                // std::cout <<
                grid_ptr->SetADiag({i, j, k}, non_solid_neighbors);
                grid_ptr->SetD({i, j, k}, - (grid_ptr->rho * grid_ptr->cellWidth / time_step) * d);
                // std::cout << "D: " << grid_ptr->GetD({i,j,k}) << std::endl;
            }
        }
    }
}


                /*                                                          |   ...   |
                                                                            |p_i-1,j,k|
                                                                            |   ...   |
                                                                            |p_i,j-1,k|
                                                                            |   ...   |
                                                                            |p_i,j,k-1|
                                                                            |   ...   |
[...A_minusI...A_minusJ...A_minusK...A_diag...A_plusI...AplusJ...AplusK...] | p_i,j,k |
                                                                            |   ...   |
                                                                            |p_i+1,j,k|
                                                                            |   ...   |
                                                                            |p_i,j+1,k|
                                                                            |   ...   |
                                                                            |p_i,j,k+1|
                                                                            |   ...   |
                */
void FluidSystem::UpdatePressure(int max_iter) {
    float epsilon = 0.00001f;
    Triple grid_size = {grid_ptr->gridWidth, grid_ptr->gridHeight, grid_ptr->gridDepth};
    float r_mag_sqrd = 0.f;
    //initialize r to d - Ap and s to d - Ap
    auto residual_ptr = std::make_shared<Grid3D<float>>(Grid3D<float>(grid_size, 0.f));
    auto search_ptr = std::make_shared<Grid3D<float>>(Grid3D<float>(grid_size, 0.f));

    for (int i = 0; i < grid_size.i; i++) {
        for (int j = 0; j < grid_size.j; j++) {
            for (int k = 0; k < grid_size.k; k++) {
                if (grid_ptr->GetCellType({i,j,k}) != CellType::LIQUID) {
                    continue;
                }

                float Ap_ijk = 0.f;
                Ap_ijk += grid_ptr->GetAPlusI({i-1, j, k}) * grid_ptr->GetPressure({i-1,j,k});
                Ap_ijk += grid_ptr->GetAPlusJ({i, j-1, k}) * grid_ptr->GetPressure({i,j-1,k});
                Ap_ijk += grid_ptr->GetAPlusK({i, j, k-1}) * grid_ptr->GetPressure({i,j,k-1});
                Ap_ijk += grid_ptr->GetAPlusI({i, j, k}) * grid_ptr->GetPressure({i+1,j,k});
                Ap_ijk += grid_ptr->GetAPlusJ({i, j, k}) * grid_ptr->GetPressure({i,j+1,k});
                Ap_ijk += grid_ptr->GetAPlusK({i, j, k}) * grid_ptr->GetPressure({i,j,k+1});
                Ap_ijk += grid_ptr->GetADiag({i, j, k}) * grid_ptr->GetPressure({i,j,k});
                
                float d_ijk = grid_ptr->GetD({i,j,k});
                float r_ijk = d_ijk - Ap_ijk;
                residual_ptr->Set({i,j,k}, r_ijk);
                search_ptr->Set({i,j,k}, r_ijk);
                r_mag_sqrd += r_ijk * r_ijk;
            }
        }
    }

    if (r_mag_sqrd < epsilon * epsilon) {
        // std::cout << "########### TOOK 0 ITERATIONS ##################" << std::endl;
        return;
    }

    //replace p_ptr with preconditioned p in future
    float r_next_mag_sqrd = 0.f;
    float alpha, beta = 0.f;
    
    int n = 0;
    while (n < max_iter) {
        n++;
        //compute alpha
        float alpha_denom = 0.f;
        for (int i = 0; i < grid_size.i; i++) {
            for (int j = 0; j < grid_size.j; j++) {
                for (int k = 0; k < grid_size.k; k++) {
                    if (grid_ptr->GetCellType({i,j,k}) != CellType::LIQUID) {
                        continue;
                    }
                    float As_ijk = 0.f;
                    As_ijk += grid_ptr->GetAPlusI({i-1, j, k}) * search_ptr->Get({i-1,j,k});
                    As_ijk += grid_ptr->GetAPlusJ({i, j-1, k}) * search_ptr->Get({i,j-1,k});
                    As_ijk += grid_ptr->GetAPlusK({i, j, k-1}) * search_ptr->Get({i,j,k-1});
                    As_ijk += grid_ptr->GetAPlusI({i, j, k}) * search_ptr->Get({i+1,j,k});
                    As_ijk += grid_ptr->GetAPlusJ({i, j, k}) * search_ptr->Get({i,j+1,k});
                    As_ijk += grid_ptr->GetAPlusK({i, j, k}) * search_ptr->Get({i,j,k+1});
                    As_ijk += grid_ptr->GetADiag({i, j, k}) * search_ptr->Get({i,j,k});
                    
                    alpha_denom += search_ptr->Get({i,j,k}) * As_ijk;
                }
            }
        }
        alpha = r_mag_sqrd/alpha_denom;

        r_next_mag_sqrd = 0.f;
        // compute next p, next r, and next r mag
        for (int i = 0; i < grid_size.i; i++) {
            for (int j = 0; j < grid_size.j; j++) {
                for (int k = 0; k < grid_size.k; k++) {
                    if (grid_ptr->GetCellType({i,j,k}) != CellType::LIQUID) {
                        continue;
                    }
                    float As_ijk = 0.f;
                    As_ijk += grid_ptr->GetAPlusI({i-1, j, k}) * search_ptr->Get({i-1,j,k});
                    As_ijk += grid_ptr->GetAPlusJ({i, j-1, k}) * search_ptr->Get({i,j-1,k});
                    As_ijk += grid_ptr->GetAPlusK({i, j, k-1}) * search_ptr->Get({i,j,k-1});
                    As_ijk += grid_ptr->GetAPlusI({i, j, k}) * search_ptr->Get({i+1,j,k});
                    As_ijk += grid_ptr->GetAPlusJ({i, j, k}) * search_ptr->Get({i,j+1,k});
                    As_ijk += grid_ptr->GetAPlusK({i, j, k}) * search_ptr->Get({i,j,k+1});
                    As_ijk += grid_ptr->GetADiag({i, j, k}) * search_ptr->Get({i,j,k});
                    
                    float old_p_ijk = grid_ptr->GetPressure({i,j,k});
                    float s_ijk = search_ptr->Get({i,j,k});
                    grid_ptr->SetPressure({i,j,k}, old_p_ijk + alpha * s_ijk);

                    float old_r_ijk = residual_ptr->Get({i,j,k});
                    float new_r_ijk = old_r_ijk - alpha * As_ijk;
                    // std::cout << "OLD r_ijk: " << old_r_ijk << ", NEW r_ijk: " << new_r_ijk << ", Alpha: " << alpha << ", As_ijk: " << As_ijk << std::endl;
                    residual_ptr->Set({i,j,k}, new_r_ijk);
                    r_next_mag_sqrd += new_r_ijk * new_r_ijk;
                    // std::cout << "Delta new increased to: " << delta_new << std::endl;
                }
            }
        }

        if (r_next_mag_sqrd < epsilon * epsilon) {
            // std::cout << "########### TOOK " << n << " ITERATIONS ##################" << std::endl;
            return;
        }

        beta = r_next_mag_sqrd / r_mag_sqrd;

        //update new s
        for (int i = 0; i < grid_size.i; i++) {
            for (int j = 0; j < grid_size.j; j++) {
                for (int k = 0; k < grid_size.k; k++) {
                    if (grid_ptr->GetCellType({i,j,k}) != CellType::LIQUID) {
                        continue;
                    }
                    float s_ijk = search_ptr->Get({i,j,k});
                    float r_ijk = residual_ptr->Get({i,j,k});
                    search_ptr->Set({i,j,k}, r_ijk + beta * s_ijk);
                }
            }
        }
        r_mag_sqrd = r_next_mag_sqrd;
        // std::cout << "Magnitude sqrd of Residual: " << r_next_mag_sqrd << std::endl;
    }
    // std::cout << "########### TOOK " << n << " ITERATIONS ##################" << std::endl;
}

void FluidSystem::ProjectVelocity(float time_step) {
    float inv_rho = 1.f / grid_ptr->rho;
    //update Us
    for (int i = 0; i <= grid_ptr->gridWidth; i++) {
        for (int j = 0; j < grid_ptr->gridHeight; j++) {
            for (int k = 0; k < grid_ptr->gridDepth; k++) {
                float u_old = grid_ptr->GetUMinus({i,j,k});

                //if U is on boundary of two non-liquids, no use in calculating velocity
                if (grid_ptr->GetCellType({i,j,k}) != CellType::LIQUID && grid_ptr->GetCellType({i-1,j,k}) != CellType::LIQUID) {
                    grid_ptr->SetUMinus({i,j,k}, u_old);
                    continue;
                }
                float pressure_above = grid_ptr->GetPressure({i,j,k}, {i-1,j,k}, time_step);
                float pressure_below = grid_ptr->GetPressure({i-1,j,k}, {i,j,k}, time_step);

                float u_new = u_old - time_step * grid_ptr->inv_width * inv_rho * (pressure_above - pressure_below);
                grid_ptr->SetUMinus({i,j,k}, u_new);

                float curr_max_time_step = (5.f * grid_ptr->cellWidth)/(abs(u_new) + 0.0001f);
                if (curr_max_time_step < grid_ptr->max_time_step) {
                    grid_ptr->max_time_step = curr_max_time_step;
                }
                // std::cout << "Projecting U from: " << u_old << ", to: " << u_new << std::endl;

            }
        }
    }

    //update Vs
    for (int i = 0; i < grid_ptr->gridWidth; i++) {
        for (int j = 0; j <= grid_ptr->gridHeight; j++) {
            for (int k = 0; k < grid_ptr->gridDepth; k++) {
                float v_old = grid_ptr->GetVMinus({i,j,k});

                //if V is on boundary of two non-liquids, no use in calculating velocity
                if (grid_ptr->GetCellType({i,j,k}) != CellType::LIQUID && grid_ptr->GetCellType({i,j-1,k}) != CellType::LIQUID) {
                    grid_ptr->SetVMinus({i,j,k}, v_old);
                    continue;
                }
                float pressure_above = grid_ptr->GetPressure({i,j,k}, {i,j-1,k}, time_step);
                float pressure_below = grid_ptr->GetPressure({i,j-1,k}, {i,j,k}, time_step);

                float v_new = v_old - time_step * grid_ptr->inv_width * inv_rho * (pressure_above - pressure_below);
                grid_ptr->SetVMinus({i,j,k}, v_new);

                float curr_max_time_step = (5.f * grid_ptr->cellWidth)/(abs(v_new) + 0.0001f);
                if (curr_max_time_step < grid_ptr->max_time_step) {
                    grid_ptr->max_time_step = curr_max_time_step;
                }
                // std::cout << "Projecting V from: " << v_old << ", to: " << v_new << std::endl;
            }
        }
    }

    //update Ws
    for (int i = 0; i < grid_ptr->gridWidth; i++) {
        for (int j = 0; j < grid_ptr->gridHeight; j++) {
            for (int k = 0; k <= grid_ptr->gridDepth; k++) {
                float w_old = grid_ptr->GetWMinus({i,j,k});

                //if W is on boundary of two non-liquids, no use in calculating velocity
                if (grid_ptr->GetCellType({i,j,k}) != CellType::LIQUID && grid_ptr->GetCellType({i,j,k-1}) != CellType::LIQUID) {
                    grid_ptr->SetWMinus({i,j,k}, w_old);
                    continue;
                }
                float pressure_above = grid_ptr->GetPressure({i,j,k}, {i,j,k-1}, time_step);
                float pressure_below = grid_ptr->GetPressure({i,j,k-1}, {i,j,k}, time_step);

                float w_new = w_old - time_step * grid_ptr->inv_width * inv_rho * (pressure_above - pressure_below);
                grid_ptr->SetWMinus({i,j,k}, w_new);

                float curr_max_time_step = (5.f * grid_ptr->cellWidth)/(abs(w_new) + 0.0001f);
                if (curr_max_time_step < grid_ptr->max_time_step) {
                    grid_ptr->max_time_step = curr_max_time_step;
                }
                // std::cout << "Projecting W from: " << w_old << ", to: " << w_new << std::endl;

            }
        }
    }
    grid_ptr->FlipStorage();
}


void FluidSystem::UpdatePhi() {
    //iterate (0,0,1) direction:  
    for (int i = 0; i < grid_ptr->gridWidth; i++) {
        for (int j = 0; j < grid_ptr->gridHeight; j++) {
            for (int k = 0; k < grid_ptr->gridDepth; k++) {
                UpdateIndividualPhi(i, j, k);
            }
        }
    }
    grid_ptr->FlipPhiStorage();
    //iterate (0,1,0) direction:  
    for (int i = 0; i < grid_ptr->gridWidth; i++) {
        for (int k = 0; k < grid_ptr->gridDepth; k++) {
            for (int j = 0; j < grid_ptr->gridHeight; j++) {
                UpdateIndividualPhi(i, j, k);
            }
        }
    }
    grid_ptr->FlipPhiStorage();
    //iterate (1,0,0) direction:  
    for (int k = 0; k < grid_ptr->gridDepth; k++) {
        for (int j = 0; j < grid_ptr->gridHeight; j++) {
            for (int i = 0; i < grid_ptr->gridWidth; i++) {
                UpdateIndividualPhi(i, j, k);
            }
        }
    }
    grid_ptr->FlipPhiStorage();
    //iterate (0,0,-1) direction:  
    for (int i = 0; i < grid_ptr->gridWidth; i++) {
        for (int j = 0; j < grid_ptr->gridHeight; j++) {
            for (int k = grid_ptr->gridDepth - 1; k >= 0; k--) {
                UpdateIndividualPhi(i, j, k);
            }
        }
    }
    grid_ptr->FlipPhiStorage();
    //iterate (0,-1,0) direction:  
    for (int i = 0; i < grid_ptr->gridWidth; i++) {
        for (int k = 0; k < grid_ptr->gridDepth; k++) {
            for (int j = grid_ptr->gridHeight - 1; j >= 0; j--) {
                UpdateIndividualPhi(i, j, k);
            }
        }
    }
    grid_ptr->FlipPhiStorage();
    //iterate (-1,0,0) direction:  
    for (int k = 0; k < grid_ptr->gridDepth; k++) {
        for (int j = 0; j < grid_ptr->gridHeight; j++) {
            for (int i = grid_ptr->gridWidth - 1; i >= 0; i--) {
                UpdateIndividualPhi(i, j, k);
            }
        }
    }
    grid_ptr->FlipPhiStorage();
}

void FluidSystem::UpdateIndividualPhi(int i, int j, int k) {
    if (grid_ptr->GetCellType({i,j,k}) == CellType::SOLID) {
        grid_ptr->SetPhi({i,j,k}, MAXFLOAT);
        return;
    }
    float h = grid_ptr->cellWidth;
    float a1 = MAXFLOAT, a2 = MAXFLOAT, a3 = MAXFLOAT;
    if (abs(grid_ptr->GetPhi({i,j,k})) < h) {
        grid_ptr->SetPhi({i,j,k}, grid_ptr->GetPhi({i,j,k}));
        return;
    }
    for (int l = 0; l < face_directions->size(); l++) {
        Triple delta = face_directions->at(l);
        Triple neighbor_coords = {i + delta.i, j + delta.j, k + delta.k};

        float neighbor_phi = grid_ptr->GetPhi(neighbor_coords);
        // if (abs(neighbor_phi) < abs(a1)) {
        //     a3 = a2;
        //     a2 = a1;
        //     a1 = neighbor_phi;
        // } else if (abs(neighbor_phi) < abs(a2)) {
        //     a3 = a2;
        //     a2 = neighbor_phi;
        // } else if (abs(neighbor_phi) < abs(a3)) {
        //     a3 = neighbor_phi;
        // }
        if (abs(neighbor_phi) < a1) {
            a3 = a2;
            a2 = a1;
            a1 = abs(neighbor_phi);
        } else if (abs(neighbor_phi) < a2) {
            a3 = a2;
            a2 = abs(neighbor_phi);
        } else if (abs(neighbor_phi) < a3) {
            a3 = abs(neighbor_phi);
        }
    }

    // grid_ptr->SetPhi({i,j,k}, a1 + (a1 > 0.f ? h : -h));
    // return;

    float determinant2 = 2 * h * h - (a1 - a2) * (a1 - a2);
    float determinant3 = 3 * h * h - (a1 - a2) * (a1 - a2) - (a1 - a3) * (a1 - a3) - (a2 - a3) * (a2 - a3);

    float resultant2 = determinant2 > 0.f ? sqrt(determinant2) : 0.f;
    float resultant3 = determinant3 > 0.f ? sqrt(determinant3) : 0.f;

    // float s1 = a1 + (a1 > 0.f ? h : -h);
    // float s2 = (a1 + a2 + ((a1 + a2 > 0.f) ? resultant2 : -resultant2)) / 2.f;
    // float s3 = (a1 + a2 + a3 + ((a1 + a2 + a3 > 0.f) ? resultant3 : -resultant3)) / 3.f;
    float s1 = a1 + h;
    float s2 = (a1 + a2 + resultant2) / 2.f;
    float s3 = (a1 + a2 + a3 + resultant3) / 3.f;

    bool positive = grid_ptr->GetCellType({i,j,k}) == CellType::AIR;
    if (abs(s1) > abs(grid_ptr->GetPhi({i,j,k}))) {
        grid_ptr->SetPhi({i,j,k}, grid_ptr->GetPhi({i,j,k}));
    } else if (abs(s1) < a2) {
        grid_ptr->SetPhi({i,j,k}, positive ? s1 : -s1);
    } else if (abs(s2) < a3) {
        grid_ptr->SetPhi({i,j,k}, positive ? s2 : -s2);
    } else {
        grid_ptr->SetPhi({i,j,k}, positive ? s3 : -s3);
    }
}

void FluidSystem::UpdateGridMarkerCells(float time_step) {
    srand(time(NULL)); // (rand() % 5)
    int max_i = grid_ptr->gridWidth, max_j = grid_ptr->gridHeight, max_k = grid_ptr->gridDepth;
    float delta = grid_ptr->cellWidth;
    auto new_cell_types = std::make_shared<Grid3D<CellType>>(Grid3D<CellType>({max_i, max_j, max_k}, CellType::AIR));
    for (int i = 0; i < max_i; i++) {
        for (int j = 0; j < max_j; j++) {
            for (int k = 0; k < max_k; k++) {
                if (grid_ptr->GetCellType({i,j,k}) != CellType::LIQUID) {
                    if (grid_ptr->GetCellType({i,j,k}) == CellType::SOLID) {
                        new_cell_types->Set({i,j,k}, CellType::SOLID);
                    }
                    continue;
                }
                // place marker particle in each of 8 quadrants
                // and project along velocity field
                glm::vec3 center_of_cell = grid_ptr->gridPosition + delta * glm::vec3(i, j, k);
                std::vector<glm::vec3> randomPositions;
                std::vector<float> randNums;
                for (int l = 0; l < 24; l++) {
                    randNums.push_back( (rand() % 100) / 99.f * delta/2.f );
                }
                randomPositions.push_back(glm::vec3(randNums[0], randNums[1], randNums[2]));
                randomPositions.push_back(glm::vec3(-randNums[3], randNums[4], randNums[5]));
                randomPositions.push_back(glm::vec3(randNums[6], -randNums[7], randNums[8]));
                randomPositions.push_back(glm::vec3(randNums[9], randNums[10], -randNums[11]));
                randomPositions.push_back(glm::vec3(-randNums[12], -randNums[13], randNums[14]));
                randomPositions.push_back(glm::vec3(-randNums[15], randNums[16], -randNums[17]));
                randomPositions.push_back(glm::vec3(randNums[18], -randNums[19], -randNums[20]));
                randomPositions.push_back(glm::vec3(-randNums[21], -randNums[22], -randNums[23]));

                for (int m = 0; m < randomPositions.size(); m++) {
                    glm::vec3 starting_position = center_of_cell + randomPositions[m];
                    // std::cout << "Starting at: " << i << ", " << j << ", " << k << std::endl;
                    // std::cout << "  start: " << starting_position.x << ", " << starting_position.y << ", " << starting_position.z << std::endl;
                    // //project forward from this position and mark new cells as liquid
                    float u = grid_ptr->InterpolateU(starting_position);
                    float v = grid_ptr->InterpolateV(starting_position);
                    float w = grid_ptr->InterpolateW(starting_position);
                    float velocity_mag = sqrt(u*u + v*v + w*w + 0.0001f);
                    int max_cells_traversed = (int)(velocity_mag/grid_ptr->cellWidth) + 1;
                    // std::cout << "      u: " << u << ", v: " << v << ", w: " << w << std::endl;
                    glm::vec3 projected_position = starting_position + glm::vec3(u,v,w) * time_step;
                    // std::cout << "  end: " << projected_position.x << ", " << projected_position.y << ", " << projected_position.z << std::endl;
                    projected_position = grid_ptr->BringInboundsP(projected_position);
                    Triple projected_coords = grid_ptr->GetCellCoordinates(projected_position);
                    // std::cout << "  at coords: " << projected_coords.i << ", " << projected_coords.j << ", " << projected_coords.k << std::endl;
                    for (int h = 0; h < max_cells_traversed ; h++) {
                        if (grid_ptr->GetCellType(projected_coords) != CellType::SOLID) {
                            break;
                        }
                        projected_position -= glm::normalize(glm::vec3(u,v,w)) * grid_ptr->cellWidth;
                        projected_position = grid_ptr->BringInboundsP(projected_position);
                        projected_coords = grid_ptr->GetCellCoordinates(projected_position);
                    }

                    new_cell_types->Set(projected_coords, CellType::LIQUID);
                    if (grid_ptr->GetCellType(projected_coords) == CellType::AIR) {
                        //update u's, v's w's
                        int p = projected_coords.i, q = projected_coords.j, r = projected_coords.k;
                        grid_ptr->FlipStorage();
                        grid_ptr->SetUMinus(projected_coords, u);
                        grid_ptr->SetUMinus({p+1,q,r}, u);
                        grid_ptr->SetVMinus(projected_coords, v);
                        grid_ptr->SetVMinus({p,q+1,r}, v);
                        grid_ptr->SetWMinus(projected_coords, w);
                        grid_ptr->SetWMinus({p,q,r+1}, w);
                        grid_ptr->SetPressure(projected_coords, grid_ptr->GetPressure({i,j,k}));
                        grid_ptr->FlipStorage();
                    }
                }


            }
        }
    }
    for (int i = 0; i < max_i; i++) {
        for (int j = 0; j < max_j; j++) {
            for (int k = 0; k < max_k; k++) {
                CellType type = new_cell_types->Get({i,j,k});
                grid_ptr->SetCellType({i,j,k}, type);
            }
        }
    }
}


}

