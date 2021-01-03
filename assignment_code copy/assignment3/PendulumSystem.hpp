#ifndef PENDULUM_SYSTEM_H_
#define PENDULUM_SYSTEM_H_

#include "ParticleSystemBase.hpp"
#include "ParticleState.hpp"
#include "gloo/InputManager.hpp"
#include <math.h>


namespace GLOO {

struct Spring {
    int start;
    int end;
    float stiffness;
    float restLength;
};

class PendulumSystem: public ParticleSystemBase {
    public:
        PendulumSystem() {
            masses = std::vector<float>();
            drag_constant = 0;
        }
        PendulumSystem(std::vector<float> masses_, float drag_constant_) {
            masses = masses_;
            drag_constant = drag_constant_;
        }

        void AddSpring(Spring spring) {
            springs.push_back(spring);
        }

        void fixParticle(int n) {
            fixed.push_back(n);
        }

        void negateG() {
            g = glm::vec3(0.f, 5.0f, 0.f);
        }

        void positiveG() {
            g = glm::vec3(0.f, -5.0f, 0.f);
        }

        ParticleState ComputeTimeDerivative(const ParticleState& state, float time) const {
            std::vector<glm::vec3> accelerations;
            std::vector<glm::vec3> velocities;
            //loop over and get accelerations from the drag and g force
            for (int i = 0; i < state.positions.size(); i++) {
                glm::vec3 velocity = state.velocities[i];
                float mass = masses[i];
                
                glm::vec3 gravity = mass * g;
                glm::vec3 drag = -drag_constant * velocity;

                accelerations.push_back(drag + gravity);
                velocities.push_back(velocity);
            }
            //loop over springs to get forces attributed from springs
            for (int l = 0; l < springs.size(); l++) {
                Spring spring = springs[l];
                glm::vec3 position1 = state.positions[spring.start];
                glm::vec3 position2 = state.positions[spring.end];
                glm::vec3 d_1 = position1 - position2;
                float magnitude = sqrt(glm::dot(d_1, d_1));
                glm::vec3 force_1 = -spring.stiffness*(magnitude - spring.restLength)*glm::normalize(d_1);
                accelerations[spring.end] -= force_1;

                accelerations[spring.start] += force_1;
            }
            // for all particles fixed, set them to zero v and a
            for (int index: fixed) {
                velocities[index] = glm::vec3(0.f, 0.f, 0.f);
                accelerations[index] = glm::vec3(0.f, 0.f, 0.f);
            }
            ParticleState derivative = {velocities, accelerations};
            return derivative;
        }

    private:
        std::vector<int> fixed;
        std::vector<float> masses;
        std::vector<Spring> springs;
        float drag_constant;
        glm::vec3 g = glm::vec3(0.0f, -5.0f, 0.0f);

};

}

#endif