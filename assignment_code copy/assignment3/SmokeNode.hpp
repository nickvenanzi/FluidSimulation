#ifndef SMOKE_NODE_H_
#define SMOKE_NODE_H_

#include "gloo/SceneNode.hpp"
#include "ParticleSystemBase.hpp"
#include "ParticleState.hpp"
#include "ForwardEulerIntegrator.hpp"
#include "SmokeSystem.hpp"
#include "IntegratorFactory.hpp"
#include "gloo/debug/PrimitiveFactory.hpp"
#include "gloo/components/RenderingComponent.hpp"
#include "gloo/components/ShadingComponent.hpp"
#include "gloo/components/MaterialComponent.hpp"
#include "gloo/shaders/PhongShader.hpp"
#include "gloo/shaders/SimpleShader.hpp"

#include "stdlib.h"
namespace GLOO {

class SmokeNode: public SceneNode {
    public:
        SmokeNode(IntegratorType integrator_type, double integrator_step) {
            integrator = IntegratorFactory::CreateIntegrator<SmokeSystem, ParticleState>(integrator_type);

            std::vector<float> masses;
            //make node components as normal
            std::shared_ptr<VertexObject> sphere_mesh_ = PrimitiveFactory::CreateSphere(0.1f, 25, 25);
            glm::vec3 color(0.3f, 0.3f, 0.3f);
            std::vector<glm::vec3> positions;
            std::vector<glm::vec3> velocities;
            glm::vec3 initial_position(0.f);
  
            system = SmokeSystem(masses);
            state = {positions, velocities};            
        }

        void Update(double delta_time) override {
            srand(time(NULL));

            if ((int)(current_time / delta_time) % 10 == 0) {
                glm::vec3 initial_velocity(((rand() % 5) - 2.5f)/10.f, 0.f, ((rand() % 5) - 2.5f)/10.f);
                system.AddMass();
                state.velocities.push_back(initial_velocity);
                auto particleNode = make_unique<SceneNode>();
                auto& rc = particleNode->CreateComponent<RenderingComponent>(sphere_mesh_);
                rc.SetDrawMode(DrawMode::Triangles);
                particleNode->CreateComponent<ShadingComponent>(std::make_shared<PhongShader>());
                particleNode->CreateComponent<MaterialComponent>(std::make_shared<Material>(smoke_color, smoke_color, smoke_color, 0.5f));
                particleNode->GetTransform().SetPosition(glm::vec3(0.f));
                particles.push_back(particleNode.get());
                AddChild(std::move(particleNode));
                state.positions.push_back(glm::vec3(0.f));
            }
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
        SmokeSystem system;
        float current_time = 0.0f;
        std::vector<SceneNode*> particles;
        std::unique_ptr<IntegratorBase<SmokeSystem, ParticleState>> integrator;
        ParticleState state;
        std::shared_ptr<VertexObject> sphere_mesh_ = PrimitiveFactory::CreateSphere(0.1f, 25, 25);
        glm::vec3 smoke_color = glm::vec3(0.3f, 0.3f, 0.3f);

};
}

#endif