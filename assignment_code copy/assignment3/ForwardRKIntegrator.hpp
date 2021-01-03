#ifndef FORWARD_RK_INTEGRATOR_H_
#define FORWARD_RK_INTEGRATOR_H_

#include "IntegratorBase.hpp"

namespace GLOO {
template <class TSystem, class TState>
class ForwardRKIntegrator : public IntegratorBase<TSystem, TState> {
    TState Integrate(const TSystem& system,
        const TState& state,
        float start_time,
        float dt) const override {
    TState k_1 = system.ComputeTimeDerivative(state, start_time);
    TState m1 = state + dt/2.0 * k_1;

    TState k_2 = system.ComputeTimeDerivative(m1, start_time + dt/2);
    TState m2 = state + dt/2.0 * k_2;

    TState k_3 = system.ComputeTimeDerivative(m2, start_time + dt/2);
    TState m3 = state + dt*k_3;

    TState k_4 = system.ComputeTimeDerivative(m3, start_time + dt);
    return state + (1.0/6.0 * dt * (k_1 + 2.0*k_2 + 2.0*k_3 + k_4));
}
};
}// namespace GLOO

#endif