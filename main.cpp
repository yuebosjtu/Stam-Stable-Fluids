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
    params.pressure_iterations = 30;
    params.enable_gravity = true;
    
    // Create fluid simulator
    FluidSimulator simulator(params);
    
    // Initialize the simulation
    simulator.Initialize("fluid_output");
    
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
