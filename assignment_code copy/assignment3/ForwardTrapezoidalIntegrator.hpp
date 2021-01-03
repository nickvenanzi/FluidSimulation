#ifndef FORWARD_TRAPEZOIDAL_INTEGRATOR_H_
#define FORWARD_TRAPEZOIDAL_INTEGRATOR_H_

#include "IntegratorBase.hpp"

namespace GLOO {
template <class TSystem, class TState>
class ForwardTrapezoidalIntegrator : public IntegratorBase<TSystem, TState> {
  TState Integrate(const TSystem& system,
                   const TState& state,
                   float start_time,
                   float dt) const override {

    TState f0 = system.ComputeTimeDerivative(state, start_time);
    TState mid = state + system.ComputeTimeDerivative(state, start_time)*dt;

    TState f1 = system.ComputeTimeDerivative(mid, start_time + dt);
    return state+dt/2.f*(f1 + f0);
  }
};
}  // namespace GLOO

#endif