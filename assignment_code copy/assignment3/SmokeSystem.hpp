#ifndef SMOKE_SYSTEM_H_
#define SMOKE_SYSTEM_H_

#include "ParticleSystemBase.hpp"
#include "ParticleState.hpp"
#include "gloo/InputManager.hpp"
#include <math.h>


namespace GLOO {

class SmokeSystem: public ParticleSystemBase {
    public:
        SmokeSystem() {
            masses = std::vector<float>();
            // drag_constant = 0;
        }
        SmokeSystem(std::vector<float> masses_/*, glm::vec3 emission_position_*/) {
            masses = masses_;
            // emission_position = emission_position_;
            // drag_constant = drag_constant_;
        }

        void AddMass() {
            masses.push_back(0.044);
        }

        ParticleState ComputeTimeDerivative(const ParticleState& state, float time) const {
            std::vector<glm::vec3> accelerations;
            std::vector<glm::vec3> velocities;
            //loop over and get accelerations from the drag and g force
            for( int i = 0; i < masses.size(); i++) {
                float mass = masses[i];
                velocities.push_back(state.velocities[i]);
                accelerations.push_back(mass*g);
            }
            ParticleState derivative = {velocities, accelerations};
            return derivative;
        }

    private:
        float r_ = 8.314; // J/(mol*K)
        float rho_ = 1.225; // kg/m^3
        // glm::vec3 emission_position;
        // std::vector<int> fixed;
        std::vector<float> masses;
        // std::vector<Spring> springs;
        // float drag_constant;
        glm::vec3 g = glm::vec3(0.0f, 2.8f, 0.0f);

};

}

#endif