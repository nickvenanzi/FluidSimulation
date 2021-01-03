#include "FluidGrid.hpp"
#include <iostream>

namespace GLOO {

FluidGrid::FluidGrid(std::vector<Triple> solidCells, 
                     std::vector<Triple> fluidCells, 
                     Triple gridSize, float cellSize, float start_time_step) {

    cellTypes_ptr = std::make_shared<Grid3D<CellType>>(Grid3D<CellType>(gridSize, CellType::AIR));
    max_time_step = MAXFLOAT;
    for (int i = 0; i < solidCells.size(); i++) {
        Triple coords = solidCells[i];
        cellTypes_ptr->Set(coords.i, coords.j, coords.k, CellType::SOLID);
    }
    for (int i = 0; i < fluidCells.size(); i++) {
        Triple coords = fluidCells[i];
        cellTypes_ptr->Set(coords.i, coords.j, coords.k, CellType::LIQUID);
    }

    cellWidth = cellSize;
    gridWidth = gridSize.i;
    gridHeight = gridSize.j;
    gridDepth = gridSize.k;
    storage = StorageSetting::ALPHA;
    phiStorage = StorageSetting::ALPHA;

    u_alpha_ptr = std::make_shared<Grid3D<float>>(Grid3D<float>({gridSize.i + 1, gridSize.j, gridSize.k}));
    u_bravo_ptr = std::make_shared<Grid3D<float>>(Grid3D<float>({gridSize.i + 1, gridSize.j, gridSize.k}));
    
    v_alpha_ptr = std::make_shared<Grid3D<float>>(Grid3D<float>({gridSize.i, gridSize.j + 1, gridSize.k}));
    v_bravo_ptr = std::make_shared<Grid3D<float>>(Grid3D<float>({gridSize.i, gridSize.j + 1, gridSize.k}));

    w_alpha_ptr = std::make_shared<Grid3D<float>>(Grid3D<float>({gridSize.i, gridSize.j, gridSize.k + 1}));
    w_bravo_ptr = std::make_shared<Grid3D<float>>(Grid3D<float>({gridSize.i, gridSize.j, gridSize.k + 1}));

    p_ptr = std::make_shared<Grid3D<float>>(Grid3D<float>(gridSize, 0.f));
    d_ptr = std::make_shared<Grid3D<float>>(Grid3D<float>(gridSize, 0.f));

    A_diag_ptr = std::make_shared<Grid3D<float>>(Grid3D<float>(gridSize));
    A_plusI_ptr = std::make_shared<Grid3D<float>>(Grid3D<float>(gridSize));
    A_plusJ_ptr = std::make_shared<Grid3D<float>>(Grid3D<float>(gridSize));
    A_plusK_ptr = std::make_shared<Grid3D<float>>(Grid3D<float>(gridSize));

    phi_alpha_ptr = std::make_shared<Grid3D<float>>(Grid3D<float>(gridSize, 0.f));
    phi_bravo_ptr = std::make_shared<Grid3D<float>>(Grid3D<float>(gridSize, 0.f));

    airCellsAtSurface = std::make_shared<std::vector<Triple>>();
    liquidCellsAtSurface = std::make_shared<std::vector<Triple>>();

    inv_width = 1.f/(float)cellWidth;
}

glm::vec3 FluidGrid::GetPosition(Triple coords) {
    return gridPosition + cellWidth * glm::vec3(coords.i, coords.j, coords.k);
}

void FluidGrid::SetCellType(Triple coords, CellType type) {
    cellTypes_ptr->Set(coords, type);
}

CellType FluidGrid::GetCellType(Triple coords) {
    return cellTypes_ptr->Get(coords);
}

float FluidGrid::GetUMinus(Triple coords) {
    if (storage == StorageSetting::ALPHA) {
        return u_bravo_ptr->Get(coords);
    } else {
        return u_alpha_ptr->Get(coords);
    }
}

float FluidGrid::GetVMinus(Triple coords) {
    if (storage == StorageSetting::ALPHA) {
        return v_bravo_ptr->Get(coords);
    } else {
        return v_alpha_ptr->Get(coords);
    }
}

float FluidGrid::GetWMinus(Triple coords) {
    if (storage == StorageSetting::ALPHA) {
        return w_bravo_ptr->Get(coords);
    } else {
        return w_alpha_ptr->Get(coords);
    }
}

void FluidGrid::SetUMinus(Triple coords, float u) {
    if (storage == StorageSetting::ALPHA) {
        u_alpha_ptr->Set(coords, u);
    } else {
        u_bravo_ptr->Set(coords, u);
    }
}

void FluidGrid::SetVMinus(Triple coords, float v) {
    if (storage == StorageSetting::ALPHA) {
        v_alpha_ptr->Set(coords, v);
    } else {
        v_bravo_ptr->Set(coords, v);
    }
}

void FluidGrid::SetWMinus(Triple coords, float w) {
    if (storage == StorageSetting::ALPHA) {
        w_alpha_ptr->Set(coords, w);
    } else {
        w_bravo_ptr->Set(coords, w);
    }
}


float FluidGrid::InterpolateU(glm::vec3 position) {
    glm::vec3 inbounds_position = BringInboundsU(position);
    glm::vec3 numCellsFromOrigin = (inbounds_position - gridPosition) * glm::vec3(inv_width) + glm::vec3(0.5f, 0.f, 0.f);
    // std::cout << "Interpolating U from: " << numCellsFromOrigin.x << ", " << numCellsFromOrigin.y << ", " << numCellsFromOrigin.z << std::endl;
    int i = (int) numCellsFromOrigin.x;
    int j = (int) numCellsFromOrigin.y;
    int k = (int) numCellsFromOrigin.z;
    // std::cout << "Interpolating U from: " << i << ", " << j << ", " << k << std::endl;
    float alpha = numCellsFromOrigin.x - (float)i;
    float beta = numCellsFromOrigin.y - (float)j;
    float gamma = numCellsFromOrigin.z - (float)k;
    // std::cout << "numCellsFromOrigin.z: " << numCellsFromOrigin.z << ", K: " << k << ", and GAMMA: " << gamma << std::endl;
    // std::cout << "      alpha: " << alpha << ", beta: " << beta << ", gamma: " << gamma << std::endl;

    float u_1 = interpolate1D(GetUMinus({i, j, k}), GetUMinus({i+1, j, k}), alpha);
    float u_2 = interpolate1D(GetUMinus({i, j, k+1}), GetUMinus({i+1, j, k+1}), alpha);
    float u_3 = interpolate1D(GetUMinus({i, j+1, k}), GetUMinus({i+1, j+1, k}), alpha);
    float u_4 = interpolate1D(GetUMinus({i, j+1, k+1}), GetUMinus({i+1, j+1, k+1}), alpha);

    float v_1 = interpolate1D(u_1, u_3, beta);
    float v_2 = interpolate1D(u_2, u_4, beta);

    return interpolate1D(v_1, v_2, gamma);
}

float FluidGrid::InterpolateV(glm::vec3 position) {
    glm::vec3 inbounds_position = BringInboundsV(position);
    glm::vec3 numCellsFromOrigin = (inbounds_position - gridPosition) * glm::vec3(inv_width) + glm::vec3(0.f, 0.5f, 0.f);
    int i = (int) numCellsFromOrigin.x;
    int j = (int) numCellsFromOrigin.y;
    int k = (int) numCellsFromOrigin.z;

    float alpha = numCellsFromOrigin.x - (float)i;
    float beta = numCellsFromOrigin.y - (float)j;
    float gamma = numCellsFromOrigin.z - (float)k;
   
    float u_1 = interpolate1D(GetVMinus({i, j, k}), GetVMinus({i+1, j, k}), alpha);
    float u_2 = interpolate1D(GetVMinus({i, j, k+1}), GetVMinus({i+1, j, k+1}), alpha);
    float u_3 = interpolate1D(GetVMinus({i, j+1, k}), GetVMinus({i+1, j+1, k}), alpha);
    float u_4 = interpolate1D(GetVMinus({i, j+1, k+1}), GetVMinus({i+1, j+1, k+1}), alpha);

    float v_1 = interpolate1D(u_1, u_3, beta);
    float v_2 = interpolate1D(u_2, u_4, beta);

    return interpolate1D(v_1, v_2, gamma);
}

float FluidGrid::InterpolateW(glm::vec3 position) {
    glm::vec3 inbounds_position = BringInboundsW(position);
    glm::vec3 numCellsFromOrigin = (inbounds_position - gridPosition) * glm::vec3(inv_width) + glm::vec3(0.f, 0.f, 0.5f);
    
    int i = (int) numCellsFromOrigin.x;
    int j = (int) numCellsFromOrigin.y;
    int k = (int) numCellsFromOrigin.z;

    float alpha = numCellsFromOrigin.x - (float)i;
    float beta = numCellsFromOrigin.y - (float)j;
    float gamma = numCellsFromOrigin.z - (float)k;
   
    float u_1 = interpolate1D(GetWMinus({i, j, k}), GetWMinus({i+1, j, k}), alpha);
    float u_2 = interpolate1D(GetWMinus({i, j, k+1}), GetWMinus({i+1, j, k+1}), alpha);
    float u_3 = interpolate1D(GetWMinus({i, j+1, k}), GetWMinus({i+1, j+1, k}), alpha);
    float u_4 = interpolate1D(GetWMinus({i, j+1, k+1}), GetWMinus({i+1, j+1, k+1}), alpha);

    float v_1 = interpolate1D(u_1, u_3, beta);
    float v_2 = interpolate1D(u_2, u_4, beta);

    return interpolate1D(v_1, v_2, gamma);
}

float FluidGrid::InterpolatePhi(glm::vec3 position) {
    glm::vec3 inbounds_position = BringInboundsP(position);
    glm::vec3 numCellsFromOrigin = (inbounds_position - gridPosition) * glm::vec3(inv_width);
    
    int i = (int) numCellsFromOrigin.x;
    int j = (int) numCellsFromOrigin.y;
    int k = (int) numCellsFromOrigin.z;

    float alpha = numCellsFromOrigin.x - (float)i;
    float beta = numCellsFromOrigin.y - (float)j;
    float gamma = numCellsFromOrigin.z - (float)k;
   
    float u_1 = interpolate1D(GetPhi({i, j, k}), GetPhi({i+1, j, k}), alpha);
    float u_2 = interpolate1D(GetPhi({i, j, k+1}), GetPhi({i+1, j, k+1}), alpha);
    float u_3 = interpolate1D(GetPhi({i, j+1, k}), GetPhi({i+1, j+1, k}), alpha);
    float u_4 = interpolate1D(GetPhi({i, j+1, k+1}), GetPhi({i+1, j+1, k+1}), alpha);

    float v_1 = interpolate1D(u_1, u_3, beta);
    float v_2 = interpolate1D(u_2, u_4, beta);

    return interpolate1D(v_1, v_2, gamma);
}

Triple FluidGrid::GetCellCoordinates(glm::vec3 position) {
    glm::vec3 numCellsFromOrigin = (position - gridPosition) * glm::vec3(inv_width) + glm::vec3(0.5f, 0.5f, 0.5f);
    int i = (int) numCellsFromOrigin.x;
    int j = (int) numCellsFromOrigin.y;
    int k = (int) numCellsFromOrigin.z;

    if (i < 0) { i = 0; }
    if (j < 0) { j = 0; }
    if (k < 0) { k = 0; }

    if (i >= gridWidth) { i = gridWidth - 1; }
    if (j >= gridHeight) { j = gridHeight - 1; }
    if (k >= gridDepth) { k = gridDepth - 1; }
 
    return {i, j, k};
}


float FluidGrid::GetPressure(Triple coords) {
    return p_ptr->Get(coords);
}  

float FluidGrid::GetPressure(Triple coords, Triple reference_coords, float time_step) {
    if (GetCellType(coords) == CellType::SOLID) {
        // if boundary on y-z plane
        if (reference_coords.i > coords.i) {
            return GetPressure(reference_coords) + rho * cellWidth * GetUMinus(reference_coords) / time_step;
        } else if (reference_coords.i < coords.i) {
            return GetPressure(reference_coords) + rho * cellWidth * GetUMinus(coords) / time_step;
        }

        // if boundary on x-z plane
        if (reference_coords.j > coords.j) {
            return GetPressure(reference_coords) + rho * cellWidth * GetVMinus(reference_coords) / time_step;
        } else if (reference_coords.j < coords.j) {
            return GetPressure(reference_coords) + rho * cellWidth * GetVMinus(coords) / time_step;
        }

        // if boundary on x-y plane
        if (reference_coords.k > coords.k) {
            return GetPressure(reference_coords) + rho * cellWidth * GetWMinus(reference_coords) / time_step;
        } else if (reference_coords.k < coords.k) {
            return GetPressure(reference_coords) + rho * cellWidth * GetWMinus(coords) / time_step;
        }
    } else if (GetCellType(coords) == CellType::AIR) {
        return -GetPressure(reference_coords);
    } else {
        return GetPressure(coords);
    }
}

float FluidGrid::GetADiag(Triple coords) {
    return A_diag_ptr->Get(coords);
}

float FluidGrid::GetAPlusI(Triple coords) {
    return A_plusI_ptr->Get(coords);
}

float FluidGrid::GetAPlusJ(Triple coords) {
    return A_plusJ_ptr->Get(coords);
}

float FluidGrid::GetAPlusK(Triple coords) {
    return A_plusK_ptr->Get(coords);
}

float FluidGrid::GetD(Triple coords) {
    return d_ptr->Get(coords);
}

float FluidGrid::GetPhi(Triple coords) {
    if (GetCellType(coords) == CellType::SOLID) {
        return MAXFLOAT;
    }
    if (phiStorage == StorageSetting::ALPHA) {
        return phi_bravo_ptr->Get(coords);
    } else {
        return phi_alpha_ptr->Get(coords);
    }
}

void FluidGrid::SetPressure(Triple coords, float val) {
    if (GetCellType(coords) == CellType::SOLID) {
        return;
    } else if (GetCellType(coords) == CellType::AIR) {
        p_ptr->Set(coords, 0.f);
    } else {
        p_ptr->Set(coords, val);
    }
}

void FluidGrid::SetADiag(Triple coords, float val) {
    A_diag_ptr->Set(coords, val);
}

void FluidGrid::SetAPlusI(Triple coords, float val) {
    A_plusI_ptr->Set(coords, val);
}

void FluidGrid::SetAPlusJ(Triple coords, float val) {
    A_plusJ_ptr->Set(coords, val);
}

void FluidGrid::SetAPlusK(Triple coords, float val) {
    A_plusK_ptr->Set(coords, val);
}

void FluidGrid::SetD(Triple coords, float val) {
    d_ptr->Set(coords, val);
}

void FluidGrid::SetPhi(Triple coords, float val) {
    if (phiStorage == StorageSetting::ALPHA) {
        phi_alpha_ptr->Set(coords, val);
    } else {
        phi_bravo_ptr->Set(coords, val);
    }
}

glm::vec3 FluidGrid::BringInboundsU(glm::vec3 position) {
    // std::cout << "BRING U IN BOUND: " << position.x << ", " << position.y << ", " << position.z << std::endl;
    glm::vec3 inbounds_position(position);
    if (position.x < epsilon - cellWidth/2.f) {
        inbounds_position.x = epsilon - cellWidth/2.f;
    } else if (position.x > (gridWidth - 0.5f) * cellWidth - epsilon) {
        inbounds_position.x = (gridWidth - 0.5f) * cellWidth - epsilon;
    }
    if (position.y < epsilon) {
        inbounds_position.y = epsilon;
    } else if (position.y > (gridHeight - 1.f) * cellWidth - epsilon) {
        inbounds_position.y = (gridHeight - 1.f) * cellWidth - epsilon;
    }
    if (position.z < epsilon) {
        inbounds_position.z = epsilon;
    } else if (position.z > (gridDepth - 1.f) * cellWidth - epsilon) {
        inbounds_position.z = (gridDepth - 1.f) * cellWidth - epsilon;
    }
    // std::cout << "OUTPUT: " << inbounds_position.x << ", " << inbounds_position.y << ", " << inbounds_position.z << std::endl;

    return inbounds_position;
}

glm::vec3 FluidGrid::BringInboundsV(glm::vec3 position) {
    glm::vec3 inbounds_position(position);
    if (position.x < epsilon) {
        inbounds_position.x = epsilon;
    } else if (position.x > (gridWidth - 1.f) * cellWidth - epsilon) {
        inbounds_position.x = (gridWidth - 1.f) * cellWidth - epsilon;
    }
    if (position.y < epsilon - cellWidth/2.f) {
        inbounds_position.y = epsilon - cellWidth/2.f;
    } else if (position.y > (gridHeight - 0.5f) * cellWidth - epsilon) {
        inbounds_position.y = (gridHeight - 0.5f) * cellWidth - epsilon;
    }
    if (position.z < epsilon) {
        inbounds_position.z = epsilon;
    } else if (position.z > (gridDepth - 1.f) * cellWidth - epsilon) {
        inbounds_position.z = (gridDepth - 1.f) * cellWidth - epsilon;
    }
    return inbounds_position;
}

glm::vec3 FluidGrid::BringInboundsW(glm::vec3 position) {
    glm::vec3 inbounds_position(position);
    if (position.x < epsilon) {
        inbounds_position.x = epsilon;
    } else if (position.x > (gridWidth - 1.f) * cellWidth - epsilon) {
        inbounds_position.x = (gridWidth - 1.f) * cellWidth - epsilon;
    }
    if (position.y < epsilon) {
        inbounds_position.y = epsilon;
    } else if (position.y > (gridHeight - 1.f) * cellWidth - epsilon) {
        inbounds_position.y = (gridHeight - 1.f) * cellWidth - epsilon;
    }
    if (position.z < epsilon - cellWidth/2.f) {
        inbounds_position.z = epsilon - cellWidth/2.f;
    } else if (position.z > (gridDepth - 0.5f) * cellWidth - epsilon) {
        inbounds_position.z = (gridDepth - 0.5f) * cellWidth - epsilon;
    }
    return inbounds_position;
}

glm::vec3 FluidGrid::BringInboundsP(glm::vec3 position) {
    glm::vec3 inbounds_position(position);
    if (position.x < epsilon) {
        inbounds_position.x = epsilon;
    } else if (position.x > (gridWidth - 1.f) * cellWidth - epsilon) {
        inbounds_position.x = (gridWidth - 1.f) * cellWidth - epsilon;
    }
    if (position.y < epsilon) {
        inbounds_position.y = epsilon;
    } else if (position.y > (gridHeight - 1.f) * cellWidth - epsilon) {
        inbounds_position.y = (gridHeight - 1.f) * cellWidth - epsilon;
    }
    if (position.z < epsilon) {
        inbounds_position.z = epsilon;
    } else if (position.z > (gridDepth - 1.f) * cellWidth - epsilon) {
        inbounds_position.z = (gridDepth - 1.f) * cellWidth - epsilon;
    }
    return inbounds_position;
}

}