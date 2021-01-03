#ifndef CLOTH_NODE_H_
#define CLOTH_NODE_H_

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
#include "gloo/InputManager.hpp"

#include <stdlib.h>

namespace GLOO {

class ClothNode: public SceneNode {
    public:
        ClothNode(IntegratorType integrator_type, double integrator_step) {
            integrator = IntegratorFactory::CreateIntegrator<PendulumSystem, ParticleState>(integrator_type);
            // physical size of cloth
            float grid_length = 5.0f;
            // drag constant
            float drag_constant = 0.5f;
            // number of points per side
            int grid_size = 200;
            

            std::shared_ptr<VertexObject> sphere_mesh_ = PrimitiveFactory::CreateSphere(0.1f, 25, 25);
            glm::vec3 color(1.0f, 1.0f, 1.0f);

            std::vector<glm::vec3> positions, velocities;
            std::vector<float> masses;
            auto net = make_unique<SceneNode>();
            vtx_obj = std::make_shared<VertexObject>();
            vtx_obj->UpdatePositions(make_unique<PositionArray>());
            vtx_obj->UpdateIndices(make_unique<IndexArray>());
            auto& rc = net->CreateComponent<RenderingComponent>(vtx_obj);
            rc.SetDrawMode(DrawMode::Triangles);
            net->CreateComponent<ShadingComponent>(std::make_shared<PhongShader>());
            net->CreateComponent<MaterialComponent>(std::make_shared<Material>(color, color, color, 0));
            AddChild(std::move(net));

            //positions of the net
            auto netPositions = make_unique<PositionArray>();
            auto netIndices = make_unique<IndexArray>();
            //loop over the whole grid
            for (int i = 0; i < grid_size; i++) {
                for (int j = 0; j < grid_size; j++) {
                    masses.push_back(1.0f/(float)grid_size);
                    velocities.push_back(glm::vec3(0));
                    glm::vec3 position(3.0f + (float)j * grid_length / (float)grid_size, -(float)i* grid_length / (float)grid_size, 0.f);
                    netPositions->push_back(position);

                    positions.push_back(position);

                }
            }

            system = PendulumSystem(masses, drag_constant);
            //make structural springs and lines
            float structural_rest_length = grid_length/(float)grid_size;
            float shear_rest_length = structural_rest_length * sqrt(2.0f);
            float flex_rest_length = 2.0f * structural_rest_length;
            float stiffness = 99.0f;
            for (int i = 0; i < grid_size; i++) {
                for (int j = 0; j < grid_size; j++) {
                    //add springs in all directions specified in the assignment
                    if (i < grid_size - 1 && j < grid_size - 1) {
                        Spring spring = {j*grid_size + i, (j+1)*grid_size + i + 1, stiffness, shear_rest_length};
                        system.AddSpring(spring);
                        //draw surface on each set of three indices
                        netIndices->push_back(j*grid_size + i);
                        netIndices->push_back(j*grid_size + i + 1);
                        netIndices->push_back((j+1)*grid_size + i + 1);
                        netIndices->push_back(j*grid_size + i);
                        netIndices->push_back((j+1)*grid_size + i);
                        netIndices->push_back((j+1)*grid_size + i + 1);
                    }
                    if (i < grid_size - 1) {
                        Spring spring = {j*grid_size + i, j*grid_size + i + 1, stiffness, structural_rest_length};
                        system.AddSpring(spring);
                    }
                    if (i < grid_size - 2) {
                        Spring spring = {j*grid_size + i, j*grid_size + i + 2, stiffness, flex_rest_length};
                        system.AddSpring(spring);
                    }
                    if (j < grid_size - 2) {
                        Spring spring = {j*grid_size + i, (j+2)*grid_size + i, stiffness, flex_rest_length};
                        system.AddSpring(spring);
                    }
                    if (j < grid_size - 1) {
                        Spring spring = {j*grid_size + i, (j+1)*grid_size + i, stiffness, structural_rest_length};
                        system.AddSpring(spring);
                    }
                    if (i > 0 && j < grid_size - 1) {
                        Spring spring = {j*grid_size + i, (j+1)*grid_size + i - 1, stiffness, shear_rest_length};
                        system.AddSpring(spring);
                    }
                }
            }
            system.fixParticle(0);
            system.fixParticle(grid_size - 1);
            vtx_obj->UpdatePositions(std::move(netPositions));
            vtx_obj->UpdateIndices(std::move(netIndices));
            UpdateNormals();
            state = {positions, velocities};
        }

        void Update(double delta_time) override {
            ParticleState newState = integrator->Integrate(system, state, current_time, delta_time);
            state = newState;
            current_time += delta_time;
            //update positions for all individual particles and their normals
            auto newPositions = make_unique<PositionArray>();
            for (int i = 0; i < state.positions.size(); i++) {
                newPositions->push_back(state.positions[i]);
            }

            vtx_obj->UpdatePositions(std::move(newPositions));
            UpdateNormals();

            // if 'n' is pressed, negate g
            if (InputManager::GetInstance().IsKeyPressed(GLFW_KEY_N)) {
                system.negateG();
            }
            // if 'p' is pressed, positive g
            if (InputManager::GetInstance().IsKeyPressed(GLFW_KEY_P)) {
                system.positiveG();
            }

        }

        // FROM ASSIGNMENT 2
        void UpdateNormals() {
            IndexArray indices = vtx_obj->GetIndices();
            PositionArray vertices = vtx_obj->GetPositions();
            auto normals = make_unique<NormalArray>();
            std::vector<glm::vec3> weighted_normals(vertices.size(), glm::vec3(0));
            for (int i = 0; i < indices.size(); i = i+3) {
                glm::vec3 position1 = vertices[indices[i]];
                glm::vec3 position2 = vertices[indices[i+1]];
                glm::vec3 position3 = vertices[indices[i+2]];
                glm::vec3 edge1 = position2 - position1;
                glm::vec3 edge2 = position3 - position1;
                glm::vec3 normal = glm::cross(edge2, edge1);
                weighted_normals[indices[i]] += normal;
                weighted_normals[indices[i+1]] += normal;
                weighted_normals[indices[i+2]] += normal;
            }
            for (int k = 0; k < weighted_normals.size(); k++) {
                normals->push_back(glm::normalize(weighted_normals[k]));
            }
            vtx_obj->UpdateNormals(std::move(normals));
        }
    
    private:
        float current_time = 0.0f;
        std::vector<SceneNode*> particles;
        std::shared_ptr<VertexObject> vtx_obj;
        std::unique_ptr<IntegratorBase<PendulumSystem, ParticleState>> integrator;
        ParticleState state;
        PendulumSystem system;
};
}

#endif