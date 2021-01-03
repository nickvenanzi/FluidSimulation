#ifndef PENDULUM_NODE_H_
#define PENDULUM_NODE_H_

#include "gloo/SceneNode.hpp"
#include "ParticleSystemBase.hpp"
#include "ParticleState.hpp"
#include "ForwardEulerIntegrator.hpp"
#include "PendulumSystem.hpp"
#include "IntegratorFactory.hpp"
#include "gloo/debug/PrimitiveFactory.hpp"
#include "gloo/components/RenderingComponent.hpp"
#include "gloo/components/ShadingComponent.hpp"
#include "gloo/components/MaterialComponent.hpp"
#include "gloo/shaders/PhongShader.hpp"
#include "gloo/shaders/SimpleShader.hpp"

namespace GLOO {

class PendulumNode: public SceneNode {
    public:
        PendulumNode(IntegratorType integrator_type, double integrator_step) {
            integrator = IntegratorFactory::CreateIntegrator<PendulumSystem, ParticleState>(integrator_type);

            std::vector<float> masses;
            float drag = 0.05f;
            //make node components as normal
            std::shared_ptr<VertexObject> sphere_mesh_ = PrimitiveFactory::CreateSphere(0.1f, 25, 25);
            glm::vec3 color(0.0f, 0.0f, 1.0f);
            std::vector<glm::vec3> positions;
            std::vector<glm::vec3> velocities;
            glm::vec3 initial_position(-3.0f, 4.0f, 0.0f);
            //loop over each of the four nodes we will
            //have in the pendulum
            for (int i = 0; i < 4; i++) {
                masses.push_back(1.0f);
                velocities.push_back(glm::vec3(0));
                auto particleNode = make_unique<SceneNode>();
                auto& rc = particleNode->CreateComponent<RenderingComponent>(sphere_mesh_);
                rc.SetDrawMode(DrawMode::Triangles);
                particleNode->CreateComponent<ShadingComponent>(std::make_shared<PhongShader>());
                particleNode->CreateComponent<MaterialComponent>(std::make_shared<Material>(color, color, color, 0));
                particleNode->GetTransform().SetPosition(initial_position);
                particles.push_back(particleNode.get());
                AddChild(std::move(particleNode));
                positions.push_back(initial_position);
                initial_position += glm::vec3(0.2f, -1.0f, 0.0f);
            }
            system = PendulumSystem(masses, drag);
            Spring s1 = {0, 1, 9.0f, 1.0f};
            Spring s2 = {1, 2, 9.0f, 1.15f};
            Spring s3 = {2, 3, 9.0f, 0.85f};
            system.AddSpring(s1);
            system.AddSpring(s2);
            system.AddSpring(s3);
            system.fixParticle(0);
            state = {positions, velocities};            
        }

        void Update(double delta_time) override {
            ParticleState newState = integrator->Integrate(system, state, current_time, delta_time);
            state = newState;
            current_time += delta_time;
            // for each particle, update their new position
            for (int i = 0; i < state.positions.size(); i++) {
                SceneNode* particleNode = particles[i];
                particleNode->GetTransform().SetPosition(state.positions[i]);
            }
        }
    private:
        PendulumSystem system;
        float current_time = 0.0f;
        std::vector<SceneNode*> particles;
        std::unique_ptr<IntegratorBase<PendulumSystem, ParticleState>> integrator;
        ParticleState state;
};
}

#endif