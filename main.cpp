#include "fluid_simulator.h"
#include <iostream>
#include <chrono>

int main()
{
    std::cout << "=== Stam Stable Fluid Simulation ===" << std::endl;
    
    // Configure simulation parameters
    FluidSimulator::SimulationParams params;
    params.width = 128;
    params.height = 128;
    params.dt = 0.016f;
    params.total_time = 5.0f;             // 5 seconds simulation
    params.viscosity = 0.0001f;
    params.diffusion = 0.0005f;           // Diffusion Constant
    params.dissipation = 0.001f;          // Dye dissipation rate
    params.vorticity_strength = 0.2f;
    params.pressure_iterations = 30;
    params.enable_vorticity = false;
    params.enable_gravity = true;        // Disable gravity for this demo
    
    // Create fluid simulator
    FluidSimulator simulator(params);
    
    // Initialize the simulation
    simulator.Initialize("fluid_output");
    
    // Add some example source events
    Eigen::Vector2f right_dir(1.0f, 0.0f);
    Eigen::Vector2f up_dir(0.0f, 1.0f);
    Eigen::Vector2f diagonal_dir(0.707f, 0.707f);
    
    // Add sources at different times and locations
    simulator.AddSourceEvent(0.5f, 32, 64, 10, right_dir, 60.0f, 2.0f);        // Initial horizontal jet from left
    simulator.AddSourceEvent(1.5f, 96, 64, 8, -right_dir, 50.0f, 1.8f);        // Counter-flow from right
    simulator.AddSourceEvent(2.5f, 64, 32, 8, up_dir, 45.0f, 2.2f);            // Upward flow from bottom
    simulator.AddSourceEvent(3.5f, 64, 64, 12, Eigen::Vector2f::Zero(), 0.0f, 3.0f); // Pure dye injection at center
    simulator.AddSourceEvent(4.0f, 48, 80, 6, diagonal_dir, 35.0f, 1.5f);      // Final diagonal push
    
    std::cout << "Added " << 5 << " source events to the simulation." << std::endl;
    
    // Record start time for performance measurement
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Run the simulation
    try 
    {
        simulator.RunSimulation();
    }
    catch (const std::exception& e)
    {
        std::cerr << "Simulation error: " << e.what() << std::endl;
        return -1;
    }
    
    // Calculate and display performance statistics
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    std::cout << "\n=== Simulation Statistics ===" << std::endl;
    std::cout << "Grid size: " << params.width << " x " << params.height << std::endl;
    std::cout << "Total frames: " << simulator.GetCurrentFrame() << std::endl;
    std::cout << "Simulation time: " << simulator.GetCurrentTime() << " seconds" << std::endl;
    std::cout << "Wall clock time: " << duration.count() << " ms" << std::endl;
    std::cout << "Average FPS: " << (simulator.GetCurrentFrame() * 1000.0f) / duration.count() << std::endl;
    
    // Display information about output files
    std::cout << "\n=== Output Files ===" << std::endl;
    std::cout << "Velocity data: fluid_output/velocity_data.txt" << std::endl;
    std::cout << "Pressure data: fluid_output/pressure_data.txt" << std::endl;
    std::cout << "Dye data: fluid_output/dye_data.txt" << std::endl;
    
    std::cout << "\nData format:" << std::endl;
    std::cout << "velocity_data.txt: frame time x y vel_x vel_y" << std::endl;
    std::cout << "pressure_data.txt: frame time x y pressure" << std::endl;
    std::cout << "dye_data.txt: frame time x y dye_concentration" << std::endl;
    
    std::cout << "Simulation completed successfully!" << std::endl;
    
    return 0;
}
