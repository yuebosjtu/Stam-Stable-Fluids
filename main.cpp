#include "fluid_simulator.h"
#include "utils.h"
#include <iostream>
#include <vector>
#include <cmath>

int main()
{
    // Create simulation parameters with circular source
    FluidSimulator::SimulationParams params;
    params.width = 256;
    params.height = 512;
    params.dt = 0.1f;
    params.total_time = 60.0f;
    params.enable_gravity = false;  // Enable gravity to see pure source effect
    params.pressure_iterations = 800;

    // Configure circular source in lower part of domain
    params.source_center_x = 128.0f;
    params.source_center_y = 64.0f;
    params.source_radius = 16.0f;       // 16 grid units radius
    params.source_velocity = 5.0f;    // Upward velocity magnitude

    params.advect_method_ = FluidSimulator::AdvectMethod::MacCormack;

    // Create fluid simulator
    FluidSimulator simulator(params);
    
    // Create initial velocity field with circular source
    std::vector<float> initial_u_field, initial_v_field;
    CreateCircularSourceField(params.width, params.height,
        params.source_center_x, params.source_center_y,
        params.source_radius, params.source_velocity,
        initial_u_field, initial_v_field);

    std::string output_dir = "fluid_output_mac";
    std::string img_dir = "velocity_field_mac";

    // Initialize simulation with custom fields
    simulator.Initialize(output_dir, img_dir, &initial_u_field, &initial_v_field, nullptr);

    // Run simulation
    simulator.RunSimulation(img_dir);

    std::cout << "\n=== Simulation completed! ===" << std::endl;
    std::cout << "The simulation shows:" << std::endl;
    std::cout << "- A circular source at (" << params.source_center_x << ", " << params.source_center_y << ")" << std::endl;
    std::cout << "- Source radius: " << params.source_radius << " grid units" << std::endl;
    std::cout << "- Source velocity: " << params.source_velocity << " units/s upward" << std::endl;
    std::cout << "- No-slip boundary conditions on all walls" << std::endl;
    std::cout << "- Constant velocity maintained inside the source region" << std::endl;
    
    return 0;
}
