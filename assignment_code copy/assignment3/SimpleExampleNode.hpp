#ifndef SIMPLE_EXAMPLE_NODE_H_
#define SIMPLE_EXAMPLE_NODE_H_

#include "gloo/debug/PrimitiveFactory.hpp"
#include "gloo/components/RenderingComponent.hpp"
#include "gloo/components/ShadingComponent.hpp"
#include "gloo/components/MaterialComponent.hpp"
#include "gloo/shaders/PhongShader.hpp"
#include "gloo/shaders/SimpleShader.hpp"
#include "gloo/SceneNode.hpp"
#include "ParticleSystemBase.hpp"
#include "ParticleState.hpp"
#include "ForwardEulerIntegrator.hpp"
#include "IntegratorFactory.hpp"


namespace GLOO {

class SimpleExampleSystem: public ParticleSystemBase {
    public:
        ParticleState ComputeTimeDerivative(const ParticleState& state, float time) const {
            std::vector<glm::vec3> derivative = std::vector<glm::vec3>();

            glm::vec3 d(-state.positions[0].y, state.positions[0].x, 0.f);
            derivative.push_back(d);
            //doesnt matter what we put in for the velocities
            ParticleState derivativeState = {derivative, derivative};

            return derivativeState;
        }
};

class SimpleExampleNode: public SceneNode {
    public:
        SimpleExampleNode(IntegratorType integrator_type, double integrator_step) {
            integrator = IntegratorFactory::CreateIntegrator<SimpleExampleSystem, ParticleState>(integrator_type);

            //make all the normal rendering stuff for the one node
            auto sphere_mesh_ = PrimitiveFactory::CreateSphere(0.025f, 25, 25);
            auto& rc = CreateComponent<RenderingComponent>(std::move(sphere_mesh_));
            rc.SetDrawMode(DrawMode::Triangles);
            CreateComponent<ShadingComponent>(std::make_shared<SimpleShader>());
            glm::vec3 color(1.0f, 0.0f, 0.0f);
            CreateComponent<MaterialComponent>(std::make_shared<Material>(color, color, color, 0));
            glm::vec3 ip(-1.5f, 1.0f, 1.0f);
            GetTransform().SetPosition(ip);
            std::vector<glm::vec3> position_array = std::vector<glm::vec3>();
            position_array.push_back(ip);
            //doesnt matter what we put for velocities
            state = {position_array, position_array};            
        }

        void Update(double delta_time) override {
            ParticleState newState = integrator->Integrate(system, state, current_time, delta_time);
            state = newState;
            current_time += delta_time;
            // update the position of the one node
            GetTransform().SetPosition(state.positions[0]);
        }
    private:
        ParticleState state;
        SimpleExampleSystem system;
        float current_time = 0.0f;
        std::unique_ptr<IntegratorBase<SimpleExampleSystem, ParticleState>> integrator;
        
};
}

#endif