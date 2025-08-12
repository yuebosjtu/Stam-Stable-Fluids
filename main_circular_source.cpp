#include "fluid_simulator.h"
#include <iostream>
#include <vector>
#include <cmath>

int main()
{
    // Create simulation parameters with circular source
    FluidSimulator::SimulationParams params;
    // params.width = 64;
    // params.height = 64;
    params.width = 128;
    params.height = 128;
    params.dt = 0.016f;
    params.total_time = 6.0f;
    params.enable_gravity = true;  // Enable gravity to see pure source effect
    params.pressure_iterations = 500;

    // Configure circular source in lower part of domain
    // params.source_center_x = 32.0f;    // Center horizontally
    // params.source_center_y = 16.0f;    // Lower third of domain
    params.source_center_x = 64.0f;
    params.source_center_y = 32.0f;
    params.source_radius = 8.0f;       // 8 grid units radius
    params.source_velocity = 5.0f;     // Upward velocity magnitude
    
    // Create fluid simulator
    FluidSimulator simulator(params);
    
    // Create initial velocity field with circular source
    std::vector<float> initial_u_field, initial_v_field;
    simulator.CreateCircularSourceField(initial_u_field, initial_v_field);
    
    // Initialize simulation with custom fields
    simulator.Initialize("fluid_output", &initial_u_field, &initial_v_field, nullptr);

    // Run simulation
    simulator.RunSimulation();
    
    std::cout << "\n=== Simulation completed! ===" << std::endl;
    std::cout << "The simulation shows:" << std::endl;
    std::cout << "- A circular source at (" << params.source_center_x << ", " << params.source_center_y << ")" << std::endl;
    std::cout << "- Source radius: " << params.source_radius << " grid units" << std::endl;
    std::cout << "- Source velocity: " << params.source_velocity << " units/s upward" << std::endl;
    std::cout << "- No-slip boundary conditions on all walls" << std::endl;
    std::cout << "- Constant velocity maintained inside the source region" << std::endl;
    
    return 0;
}
