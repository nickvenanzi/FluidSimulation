#ifndef FLUID_GRID_H_
#define FLUID_GRID_H_

#include <vector>
#include <glm/glm.hpp>

namespace GLOO {

struct Triple {
    int i;
    int j;
    int k;

    bool operator== (const Triple &other) const {
        return (i == other.i) && (j == other.j) && (k == other.k);
    }
};

template<typename T>
class Grid3D {
    public:
        Grid3D() {
        }

        Grid3D(Triple gridSize) {
            vector_ = std::vector<T>(gridSize.i * gridSize.j * gridSize.k);
            max_i = gridSize.i, max_j = gridSize.j, max_k = gridSize.k;
        }

        Grid3D(Triple gridSize, T defaultValue) {
            vector_ = std::vector<T>(gridSize.i * gridSize.j * gridSize.k, defaultValue);
            max_i = gridSize.i, max_j = gridSize.j, max_k = gridSize.k;
        }

        T& Get(int i, int j, int k) {
            if (i >= max_i) { i = max_i - 1; }
            if (j >= max_j) { j = max_j - 1; }
            if (k >= max_k) { k = max_k - 1; }
            
            if (i < 0) { i = 0; }
            if (j < 0) { j = 0; }
            if (k < 0) { k = 0; }

            return vector_[i * max_j * max_k + j * max_k + k];
        }

        T& Get(Triple coords) {
            return Get(coords.i, coords.j, coords.k);
        }

        void Set(int i, int j, int k, T val) {
            // if (i >= max_i || j >= max_j || k >= max_k) {
            //     throw std::out_of_range("Index out of range");
            // }
            if (i >= max_i) { i = max_i - 1; }
            if (j >= max_j) { j = max_j - 1; }
            if (k >= max_k) { k = max_k - 1; }
            
            if (i < 0) { i = 0; }
            if (j < 0) { j = 0; }
            if (k < 0) { k = 0; }
            vector_[i * max_j * max_k + j * max_k + k] = val;
        }

        void Set(Triple coords, T val) {
            Set(coords.i, coords.j, coords.k, val);
        }

    private:
        int max_i, max_j, max_k;
        std::vector<T> vector_;
};

enum class CellType {
    SOLID, AIR, LIQUID
};

enum class StorageSetting {
    ALPHA, BRAVO
};


class FluidGrid {
    public:
        FluidGrid(std::vector<Triple> solidCells, std::vector<Triple> fluidCells, Triple gridSize, float cellSize, float start_time_step);  
        
        void SetCellType(Triple coords, CellType type); 
        CellType GetCellType(Triple coords); 
        Triple GetCellCoordinates(glm::vec3 position);
        glm::vec3 GetPosition(Triple coords);

        float GetUMinus(Triple coords);
        float GetVMinus(Triple coords);
        float GetWMinus(Triple coords);

        void SetUMinus(Triple coords, float u);
        void SetVMinus(Triple coords, float v);
        void SetWMinus(Triple coords, float w);

        // float InterpolatePressure(glm::vec3 position);
        float InterpolateU(glm::vec3 position);
        float InterpolateV(glm::vec3 position);
        float InterpolateW(glm::vec3 position);
        float InterpolatePhi(glm::vec3 position);

        float GetPressure(Triple coords);
        float GetPressure(Triple coords, Triple reference_coords, float time_step);
        float GetADiag(Triple coords);
        float GetAPlusI(Triple coords);
        float GetAPlusJ(Triple coords);
        float GetAPlusK(Triple coords);
        float GetD(Triple coords);
        float GetPhi(Triple coords);

        void SetPressure(Triple coords, float val);
        void SetADiag(Triple coords, float val);
        void SetAPlusI(Triple coords, float val);
        void SetAPlusJ(Triple coords, float val);
        void SetAPlusK(Triple coords, float val);
        void SetD(Triple coords, float val);
        void SetPhi(Triple coords, float val);

        glm::vec3 BringInboundsU(glm::vec3 position);
        glm::vec3 BringInboundsV(glm::vec3 position);
        glm::vec3 BringInboundsW(glm::vec3 position);
        glm::vec3 BringInboundsP(glm::vec3 position);


        void FlipStorage() {
            if (storage == StorageSetting::ALPHA) {storage = StorageSetting::BRAVO;} 
            else {storage = StorageSetting::ALPHA;}
        }

        void FlipPhiStorage() {
            if (phiStorage == StorageSetting::ALPHA) {phiStorage = StorageSetting::BRAVO;} 
            else {phiStorage = StorageSetting::ALPHA;}
        }

        float cellWidth;
        int gridWidth;
        int gridHeight;
        int gridDepth;
        StorageSetting storage;
        StorageSetting phiStorage;
        glm::vec3 gridPosition = glm::vec3(0.f);
        float g = -9.8f;
        float rho = 1000.f;
        float inv_width;

        float max_time_step;
        std::shared_ptr<Grid3D<CellType>> cellTypes_ptr;
        std::shared_ptr<std::vector<Triple>> airCellsAtSurface;
        std::shared_ptr<std::vector<Triple>> liquidCellsAtSurface;

    
    private:
        static float interpolate1D(float p1, float p2, float alpha) {
            return p1 * (1 - alpha) + p2 * alpha;
        }
        float epsilon = 0.000001f;
        std::shared_ptr<Grid3D<float>> u_alpha_ptr, u_bravo_ptr;
        std::shared_ptr<Grid3D<float>> v_alpha_ptr, v_bravo_ptr;
        std::shared_ptr<Grid3D<float>> w_alpha_ptr, w_bravo_ptr;
        std::shared_ptr<Grid3D<float>> phi_alpha_ptr, phi_bravo_ptr;
        std::shared_ptr<Grid3D<float>> p_ptr;
        std::shared_ptr<Grid3D<float>> A_diag_ptr, A_plusI_ptr, A_plusJ_ptr, A_plusK_ptr, d_ptr;
        
};

}

#endif
