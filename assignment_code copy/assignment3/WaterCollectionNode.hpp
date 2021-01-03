#ifndef WATER_COLLECTION_NODE_H_
#define WATER_COLLECTION_NODE_H_

#include "gloo/SceneNode.hpp"
#include "ParticleState.hpp"
#include "gloo/debug/PrimitiveFactory.hpp"
#include "gloo/components/RenderingComponent.hpp"
#include "gloo/components/ShadingComponent.hpp"
#include "gloo/components/MaterialComponent.hpp"
#include "gloo/shaders/PhongShader.hpp"
#include "gloo/shaders/SimpleShader.hpp"
#include "gloo/InputManager.hpp"
#include "FluidGrid.hpp"
#include "FluidSystem.hpp"

namespace GLOO {
const int p = 2;

class WaterCollectionNode: public SceneNode {
    public:
        WaterCollectionNode(std::shared_ptr<Grid3D<CellType>> cellGrid, float cellWidth, Triple gridSize) {
            glm::vec3 color(0.f, 0.f, 0.5f);
            // std::shared_ptr<VertexObject> sphere_mesh_ = PrimitiveFactory::CreateSphere(cellWidth/2.f, 4, 4);
            std::shared_ptr<VertexObject> cube_mesh_ = PrimitiveFactory::CreateCube(cellWidth);
            for (int i = 0; i < gridSize.i; i++) {
                for (int j = 0; j < gridSize.j; j++) {
                    for (int k = 0; k < gridSize.k; k++) {
                        if (cellGrid->Get({i,j,k}) == CellType::LIQUID) {
                            auto waterNode = make_unique<SceneNode>();
                            auto& rc = waterNode->CreateComponent<RenderingComponent>(cube_mesh_);
                            rc.SetDrawMode(DrawMode::Triangles);
                            waterNode->CreateComponent<ShadingComponent>(std::make_shared<PhongShader>());
                            waterNode->CreateComponent<MaterialComponent>(std::make_shared<Material>(color, color, color, 0));
                            waterNode->GetTransform().SetPosition(cellWidth * glm::vec3(i,j,k));
                            AddChild(std::move(waterNode));
                        }
                        if (cellGrid->Get({i,j,k}) == CellType::SOLID && !(i == 0 || j == 0 || k == 0 || i == gridSize.i - 1 || j == gridSize.j - 1 || k == gridSize.k - 1) ) {
                            auto solidNode = make_unique<SceneNode>();
                            glm::vec3 solidColor(1.f);
                            auto& rc = solidNode->CreateComponent<RenderingComponent>(cube_mesh_);
                            rc.SetDrawMode(DrawMode::Triangles);
                            solidNode->CreateComponent<ShadingComponent>(std::make_shared<PhongShader>());
                            solidNode->CreateComponent<MaterialComponent>(std::make_shared<Material>(solidColor, solidColor, solidColor, 0));
                            solidNode->GetTransform().SetPosition(cellWidth * glm::vec3(i,j,k));
                            AddChild(std::move(solidNode));
                        }
                    }
                }

            }
        }
        
};

}

#endif