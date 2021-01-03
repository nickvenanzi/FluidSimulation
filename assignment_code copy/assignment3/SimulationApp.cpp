#include "SimulationApp.hpp"

#include "glm/gtx/string_cast.hpp"

#include "gloo/shaders/PhongShader.hpp"
#include "gloo/components/RenderingComponent.hpp"
#include "gloo/components/ShadingComponent.hpp"
#include "gloo/components/CameraComponent.hpp"
#include "gloo/components/LightComponent.hpp"
#include "gloo/components/MaterialComponent.hpp"
#include "gloo/MeshLoader.hpp"
#include "gloo/lights/DirectionalLight.hpp"
#include "gloo/lights/AmbientLight.hpp"
#include "gloo/cameras/ArcBallCameraNode.hpp"
#include "gloo/debug/AxisNode.hpp"

#include "IntegratorBase.hpp"
#include "SimpleExampleNode.hpp"
#include "PendulumNode.hpp"
#include "ClothNode.hpp"
// #include "SmokeNode.hpp"
#include "FluidGrid.hpp"
#include "FluidSystem.hpp"
#include "FluidNode.hpp"

namespace GLOO {
SimulationApp::SimulationApp(const std::string& app_name,
                             glm::ivec2 window_size,
                             IntegratorType integrator_type,
                             float integration_step)
    : Application(app_name, window_size),
      integrator_type_(integrator_type),
      integration_step_(integration_step) {

}

void SimulationApp::SetupScene() {
  SceneNode& root = scene_->GetRootNode();

  auto camera_node = make_unique<ArcBallCameraNode>(40.f, -3.0f, 10.0f);
  scene_->ActivateCamera(camera_node->GetComponentPtr<CameraComponent>());
  root.AddChild(std::move(camera_node));

  root.AddChild(make_unique<AxisNode>('A'));

  auto ambient_light = std::make_shared<AmbientLight>();
  ambient_light->SetAmbientColor(glm::vec3(0.2f));
  root.CreateComponent<LightComponent>(ambient_light);

  auto directional_light = std::make_shared<DirectionalLight>();
  directional_light->SetDiffuseColor(glm::vec3(0.8f, 0.8f, 0.8f));
  directional_light->SetSpecularColor(glm::vec3(1.0f, 1.0f, 1.0f));
  directional_light->SetDirection(glm::normalize(glm::vec3(1.f, -1.f, 0.25f)));
  // directional_light->SetAttenuation(glm::vec3(1.0f, 0.09f, 0.032f));
  auto directional_light_node = make_unique<SceneNode>();
  directional_light_node->CreateComponent<LightComponent>(directional_light);
  // directional_light_node->GetTransform().SetPosition(glm::vec3(0.0f, 2.0f, 4.f));
  root.AddChild(std::move(directional_light_node));
  // root.AddChild(make_unique<SimpleExampleNode>(integrator_type_, integration_step_));
  // root.AddChild(make_unique<PendulumNode>(integrator_type_, integration_step_));
  // root.AddChild(make_unique<ClothNode>(integrator_type_, integration_step_));
  // root.AddChild(make_unique<SmokeNode>(integrator_type_, integration_step_));
  auto fluidNode = make_unique<FluidNode>(integration_step_);
  root.AddChild(std::move(fluidNode));
}
}  // namespace GLOO
