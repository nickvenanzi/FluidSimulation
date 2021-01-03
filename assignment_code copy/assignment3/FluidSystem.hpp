#ifndef FLUID_SYSTEM_H_
#define FLUID_SYSTEM_H_

#include "FluidGrid.hpp"
#include "gloo/SceneNode.hpp"

namespace GLOO {

class FluidSystem {
    public:
        FluidSystem(FluidGrid& grid) {
            grid_ptr = make_unique<FluidGrid>(grid);
            face_directions = make_unique<std::vector<Triple>>();
            face_directions->push_back({1,0,0});
            face_directions->push_back({-1,0,0});
            face_directions->push_back({0,1,0});
            face_directions->push_back({0,-1,0});
            face_directions->push_back({0,0,1});
            face_directions->push_back({0,0,-1});
        }
        void AdvectPhi(float time_step);
        void AdvectVelocity(float time_step);
        void UpdateBodyForces(float time_step);
        void UpdatePressureSOE(float time_step);
        void UpdatePressure(int max_iter);
        void ProjectVelocity(float time_step);
        void UpdatePhi();
        void UpdateGridMarkerCells(float time_step);
        std::unique_ptr<FluidGrid> grid_ptr;

    private:
        void UpdateIndividualPhi(int i, int j, int k);
        std::unique_ptr<std::vector<Triple>> face_directions;

};

}

#endif