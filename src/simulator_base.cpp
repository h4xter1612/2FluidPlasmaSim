#include "simulator_base.hh"

SimulatorBase::SimulatorBase(const PlasmaParams& params)
    : params_(params), dispersion_(params) {}
