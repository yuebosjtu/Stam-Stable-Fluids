#include "fluid_simulator.h"
#include <iostream>
#include <vector>
#include <cmath>

// Helper function to create a vortex field
void CreateVortexField(int width, int height, float center_x, float center_y, float strength,
                      std::vector<float>& u_field, std::vector<float>& v_field)
{
    // Resize fields to proper MAC grid sizes
    int u_total_cells = (width + 1) * height;
    int v_total_cells = width * (height + 1);
    
    u_field.resize(u_total_cells, 0.0f);
    v_field.resize(v_total_cells, 0.0f);
    
    // Create vortex in u field (stored at (i, j+0.5))
    for (int i = 0; i < width + 1; ++i)
    {
        for (int j = 0; j < height; ++j)
        {
            float x = static_cast<float>(i);
            float y = static_cast<float>(j) + 0.5f;
            
            float dx = x - center_x;
            float dy = y - center_y;
            float r_sq = dx*dx + dy*dy;
            
            if (r_sq > 0.0f)
            {
                // Vortex velocity: tangential component
                u_field[IXY(i, j, width + 1)] = -strength * dy * std::exp(-r_sq / 10.0f);
            }
        }
    }
    
    // Create vortex in v field (stored at (i+0.5, j))
    for (int i = 0; i < width; ++i)
    {
        for (int j = 0; j < height + 1; ++j)
        {
            float x = static_cast<float>(i) + 0.5f;
            float y = static_cast<float>(j);
            
            float dx = x - center_x;
            float dy = y - center_y;
            float r_sq = dx*dx + dy*dy;
            
            if (r_sq > 0.0f)
            {
                // Vortex velocity: tangential component
                v_field[IXY(i, j, width)] = strength * dx * std::exp(-r_sq / 10.0f);
            }
        }
    }
}

// Helper function to create a dye blob
void CreateDyeBlob(int width, int height, float center_x, float center_y, 
                   float radius, float concentration, std::vector<float>& dye_field)
{
    int p_total_cells = width * height;
    dye_field.resize(p_total_cells, 0.0f);
    
    // Create circular blob
    for (int i = 0; i < width; ++i)
    {
        for (int j = 0; j < height; ++j)
        {
            float x = static_cast<float>(i) + 0.5f;  // Cell center
            float y = static_cast<float>(j) + 0.5f;  // Cell center
            
            float dx = x - center_x;
            float dy = y - center_y;
            float r = std::sqrt(dx*dx + dy*dy);
            
            if (r <= radius)
            {
                // Smooth falloff
                float factor = 1.0f - (r / radius);
                dye_field[IXY(i, j, width)] = concentration * factor * factor;
            }
        }
    }
}

int main()
{
    // Create simulation parameters
    FluidSimulator::SimulationParams params;
    params.width = 64;
    params.height = 64;
    params.dt = 0.016f;
    params.total_time = 5.0f;
    params.enable_gravity = false;  // Disable gravity for vortex example
    
    // Create fluid simulator
    FluidSimulator simulator(params);
    
    // Example 1: Default initialization (zero fields)
    std::cout << "=== Running simulation with zero initial fields ===" << std::endl;
    simulator.Initialize("output_zero");
    simulator.RunSimulation();
    
    std::cout << "\n=== Running simulation with custom vortex field ===" << std::endl;
    
    // Example 2: Create custom initial fields
    std::vector<float> initial_u_field, initial_v_field;
    CreateVortexField(params.width, params.height, 32.0f, 32.0f, 2.0f, initial_u_field, initial_v_field);
    
    std::vector<float> initial_dye_field;
    CreateDyeBlob(params.width, params.height, 20.0f, 20.0f, 5.0f, 1.0f, initial_dye_field);
    
    // Create new simulator for second example
    FluidSimulator simulator2(params);
    
    // Initialize with custom fields
    simulator2.Initialize("output_vortex", &initial_u_field, &initial_v_field, &initial_dye_field);
    simulator2.RunSimulation();
    
    return 0;
}
