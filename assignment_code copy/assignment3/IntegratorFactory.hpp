#ifndef INTEGRATOR_FACTORY_H_
#define INTEGRATOR_FACTORY_H_

#include "IntegratorBase.hpp"
#include "ForwardEulerIntegrator.hpp"
#include "ForwardTrapezoidalIntegrator.hpp"
#include "ForwardRKIntegrator.hpp"
#include <stdexcept>
#include "SimpleExampleNode.hpp"
#include "gloo/utils.hpp"

#include "IntegratorType.hpp"

namespace GLOO {
class IntegratorFactory {
 public:
  template <class TSystem, class TState>
  static std::unique_ptr<IntegratorBase<TSystem, TState>> CreateIntegrator(
      IntegratorType type) {
    if (type == IntegratorType::Trapezoidal) {
      return make_unique<ForwardTrapezoidalIntegrator<TSystem, TState>>();
    } else if (type == IntegratorType::Euler) {
      return make_unique<ForwardEulerIntegrator<TSystem, TState>>();
    } else if (type == IntegratorType::RK4) {
      return make_unique<ForwardRKIntegrator<TSystem, TState>>();
    }

  }
};
}  // namespace GLOO

#endif
