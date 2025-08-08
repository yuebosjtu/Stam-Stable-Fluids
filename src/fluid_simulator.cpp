#include "fluid_simulator.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <direct.h>  // For _mkdir on Windows

FluidSimulator::FluidSimulator(const SimulationParams& params)
    : params_(params), current_time_(0.0f), current_frame_(0), next_source_index_(0)
{
    InitializeFields();
}

FluidSimulator::~FluidSimulator()
{
    if (velocity_file_.is_open()) velocity_file_.close();
    if (pressure_file_.is_open()) pressure_file_.close();
    if (dye_file_.is_open()) dye_file_.close();
    
    // Clean up allocated memory
    delete velocity_pair_;
    delete pressure_pair_;
    delete dye_pair_;
}

void FluidSimulator::InitializeFields()
{
    int total_cells = params_.width * params_.height;
    
    // Initialize velocity fields
    velocity_field_.resize(total_cells, Vector2f::Zero());
    velocity_temp_.resize(total_cells, Vector2f::Zero());
    velocity_pair_ = new TexPair<std::vector<Vector2f> >(velocity_field_, velocity_temp_);
    
    // Initialize pressure fields
    pressure_field_.resize(total_cells, 0.0f);
    pressure_temp_.resize(total_cells, 0.0f);
    pressure_pair_ = new TexPair<std::vector<float> >(pressure_field_, pressure_temp_);
    
    // Initialize dye fields
    dye_field_.resize(total_cells, 0.0f);
    dye_temp_.resize(total_cells, 0.0f);
    dye_pair_ = new TexPair<std::vector<float> >(dye_field_, dye_temp_);
    
    // Initialize helper fields
    divergence_field_.resize(total_cells, 0.0f);
    curl_field_.resize(total_cells, 0.0f);
    curl_force_.resize(total_cells, Vector2f::Zero());
}

void FluidSimulator::Initialize(const std::string& output_dir)
{
    CreateOutputDirectory(output_dir);
    
    // Open output files
    velocity_file_.open((output_dir + "/velocity_data.txt").c_str());
    pressure_file_.open((output_dir + "/pressure_data.txt").c_str());
    dye_file_.open((output_dir + "/dye_data.txt").c_str());
    
    if (!velocity_file_.is_open() || !pressure_file_.is_open() || !dye_file_.is_open())
    {
        std::cerr << "Error: Failed to open output files!" << std::endl;
        return;
    }
    
    // Write headers
    velocity_file_ << "# Velocity field data: frame time x y vel_x vel_y" << std::endl;
    pressure_file_ << "# Pressure field data: frame time x y pressure" << std::endl;
    dye_file_ << "# Dye field data: frame time x y dye_concentration" << std::endl;
    
    std::cout << "Fluid simulation initialized." << std::endl;
    std::cout << "Grid size: " << params_.width << "x" << params_.height << std::endl;
    std::cout << "Time step: " << params_.dt << "s" << std::endl;
    std::cout << "Total time: " << params_.total_time << "s" << std::endl;
    std::cout << "Output directory: " << output_dir << std::endl;
}

void FluidSimulator::Step()
{
    // 1. Apply external forces (gravity, user input, etc.)
    ApplyForces();
    
    // Process any source events at current time
    ProcessSources();
    
    // 2. Advect velocity field
    AdvectVelocity();
    
    // 3. Diffuse velocity field (viscosity)
    DiffuseVelocity();
    
    // 4. Project - Solve pressure to ensure incompressibility
    SolvePressure();
    
    // Apply vorticity confinement if enabled (after pressure projection)
    if (params_.enable_vorticity)
    {
        ApplyVorticityConfinement();
    }
    
    // Advect and diffuse dye field
    AdvectDye();
    DiffuseDye();
    
    // Apply dissipation to dye field
    DissipateDye();
    
    // Save frame data
    SaveFrameData();
    
    // Update time and frame
    current_time_ += params_.dt;
    current_frame_++;
    
    // Progress output
    if (current_frame_ % 60 == 0) // Print every second (assuming 60 FPS)
    {
        std::cout << "Frame " << current_frame_ 
                  << ", Time: " << std::fixed << std::setprecision(2) << current_time_ 
                  << "s" << std::endl;
    }
}

void FluidSimulator::RunSimulation()
{
    std::cout << "Starting fluid simulation..." << std::endl;
    
    while (!IsComplete())
    {
        Step();
    }
    
    std::cout << "Simulation completed!" << std::endl;
    std::cout << "Total frames: " << current_frame_ << std::endl;
    std::cout << "Final time: " << current_time_ << "s" << std::endl;
}

void FluidSimulator::AddSourceEvent(float time, int x, int y, int radius, 
                                   const Vector2f& direction, float strength, float dye_amount)
{
    source_events_.push_back(SourceEvent(time, x, y, radius, direction, strength, dye_amount));
    
    // Keep source events sorted by time (simple bubble sort for compatibility)
    for (size_t i = 0; i < source_events_.size(); ++i) {
        for (size_t j = 0; j < source_events_.size() - 1 - i; ++j) {
            if (source_events_[j].time > source_events_[j + 1].time) {
                SourceEvent temp = source_events_[j];
                source_events_[j] = source_events_[j + 1];
                source_events_[j + 1] = temp;
            }
        }
    }
}

void FluidSimulator::ApplyForces()
{
    if (params_.enable_gravity)
    {
        apply_force<Vector2f>(params_.width, params_.height, velocity_pair_->cur, -9.8f, params_.dt);
    }
}

void FluidSimulator::ApplyVorticityConfinement()
{
    // Calculate curl of velocity field
    get_curl<Vector2f>(params_.width, params_.height, velocity_pair_->cur, curl_field_);
    
    // Apply vorticity confinement force
    vorticity_confinement<Vector2f>(params_.width, params_.height, velocity_pair_->cur, 
                                   curl_field_, curl_force_, params_.vorticity_strength, params_.dt);
}

void FluidSimulator::SolvePressure()
{
    // Calculate divergence of velocity field
    get_divergence<Vector2f>(params_.width, params_.height, velocity_pair_->cur, divergence_field_);
    
    // Clear pressure field
    std::fill(pressure_pair_->cur.begin(), pressure_pair_->cur.end(), 0.0f);
    
    // Solve pressure using Gauss-Seidel iterations
    for (int iter = 0; iter < params_.pressure_iterations; ++iter)
    {
        pressure_gauss_sidel<Vector2f>(params_.width, params_.height, divergence_field_, 
                                      pressure_pair_->cur, pressure_pair_->nxt);
        pressure_pair_->Swap();
    }
    
    // Subtract pressure gradient from velocity
    subtract_gradient<Vector2f>(params_.width, params_.height, velocity_pair_->cur, pressure_pair_->cur);
}

void FluidSimulator::AdvectVelocity()
{
    advection_velocity<Vector2f>(params_.width, params_.height, velocity_pair_->cur, 
                                velocity_pair_->nxt, params_.dt, 0.999f);
    velocity_pair_->Swap();
}

void FluidSimulator::AdvectDye()
{
    advection<float, Vector2f>(params_.width, params_.height, velocity_pair_->cur, 
                              dye_pair_->cur, dye_pair_->nxt, params_.dt, 0.999f);
    dye_pair_->Swap();
}

void FluidSimulator::DiffuseVelocity()
{
    if (params_.viscosity <= 0.0f) return;
    
    // Copy current velocity to temp for diffusion source
    velocity_temp_ = velocity_pair_->cur;
    
    // Perform Gauss-Seidel iterations for viscosity diffusion
    for (int iter = 0; iter < params_.pressure_iterations; ++iter)
    {
        diffuse_gauss_seidel<Vector2f, float>(params_.width, params_.height, 
                                             params_.viscosity, params_.dt,
                                             velocity_temp_, velocity_pair_->cur);
    }
}

void FluidSimulator::DiffuseDye()
{
    if (params_.diffusion <= 0.0f) return;
    
    // Copy current dye to temp for diffusion source
    dye_temp_ = dye_pair_->cur;
    
    // Perform Gauss-Seidel iterations for dye diffusion
    for (int iter = 0; iter < params_.pressure_iterations; ++iter)
    {
        diffuse_gauss_seidel<float, float>(params_.width, params_.height, 
                                          params_.diffusion, params_.dt,
                                          dye_temp_, dye_pair_->cur);
    }
}

void FluidSimulator::DissipateDye()
{
    // Apply dissipation (decay) to dye field
    dissipate<float>(params_.width, params_.height, dye_pair_->cur, params_.dissipation, params_.dt);
}

void FluidSimulator::ProcessSources()
{
    while (next_source_index_ < source_events_.size() && 
           source_events_[next_source_index_].time <= current_time_)
    {
        const SourceEvent& event = source_events_[next_source_index_];
        
        // Check bounds
        if (event.x >= 0 && event.x < params_.width && 
            event.y >= 0 && event.y < params_.height)
        {
            Vector2f normalized_dir = event.direction.normalized();
            Vector2f scaled_velocity = normalized_dir * event.strength;
            add_source<Vector2f>(params_.width, event.x, event.y, event.radius, 
                               event.dye_amount, scaled_velocity,
                               dye_pair_->cur, velocity_pair_->cur);
            
            std::cout << "Applied source at frame " << current_frame_ 
                      << ", position (" << event.x << ", " << event.y << ")" << std::endl;
        }
        
        next_source_index_++;
    }
}

void FluidSimulator::SaveFrameData()
{
    if (!velocity_file_.is_open() || !pressure_file_.is_open() || !dye_file_.is_open())
        return;
    
    for (int j = 0; j < params_.height; ++j)
    {
        for (int i = 0; i < params_.width; ++i)
        {
            int idx = IXY(i, j, params_.width);
            
            // Save velocity data
            velocity_file_ << current_frame_ << " " << current_time_ << " " 
                          << i << " " << j << " " 
                          << velocity_pair_->cur[idx](0) << " " 
                          << velocity_pair_->cur[idx](1) << std::endl;
            
            // Save pressure data
            pressure_file_ << current_frame_ << " " << current_time_ << " " 
                          << i << " " << j << " " 
                          << pressure_pair_->cur[idx] << std::endl;
            
            // Save dye data
            dye_file_ << current_frame_ << " " << current_time_ << " " 
                       << i << " " << j << " " 
                       << dye_pair_->cur[idx] << std::endl;
        }
    }
}

void FluidSimulator::CreateOutputDirectory(const std::string& dir)
{
    _mkdir(dir.c_str());
    std::cout << "Output directory: " << dir << std::endl;
}
