#include "fluid_simulator.h"
#include "../visualize/colorramp.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <direct.h>  // For _mkdir on Windows

FluidSimulator::FluidSimulator(const SimulationParams& params)
    : params_(params), current_time_(0.0f), save_frame_(0), current_frame_(0), 
      poisson_matrix_(nullptr), amg_initialized_(false)
{
    InitializeFields();
}

FluidSimulator::~FluidSimulator()
{
    if (velocity_file_.is_open()) velocity_file_.close();
    if (pressure_file_.is_open()) pressure_file_.close();
    if (dye_file_.is_open()) dye_file_.close();
    
    delete u_pair_;
    delete v_pair_;
    delete pressure_pair_;
    delete dye_pair_;
    
    delete poisson_matrix_;
    for (auto* matrix : A_L_) delete matrix;
    for (auto* matrix : R_L_) delete matrix;
    for (auto* matrix : P_L_) delete matrix;
}

void FluidSimulator::InitializeFields(const std::vector<float>* initial_u_field,
                                      const std::vector<float>* initial_v_field,
                                      const std::vector<float>* initial_dye_field)
{
    int u_total_cells = (params_.width + 1) * params_.height;
    int v_total_cells = params_.width * (params_.height + 1);
    int p_total_cells = params_.width * params_.height;
    
    // Initialize u velocity component field
    u_field_.resize(u_total_cells, 0.0f);
    u_temp_.resize(u_total_cells, 0.0f);
    
    // Set initial u velocity field if provided
    if (initial_u_field != nullptr)
    {
        if (initial_u_field->size() == u_total_cells)
        {
            u_field_ = *initial_u_field;
            std::cout << "Using custom initial u velocity field." << std::endl;
        }
        else
        {
            std::cerr << "Warning: Initial u field size mismatch. Using zero initialization." << std::endl;
        }
    }
    
    u_pair_ = new TexPair<std::vector<float> >(u_field_, u_temp_);
    
    // Initialize v velocity component field
    v_field_.resize(v_total_cells, 0.0f);
    v_temp_.resize(v_total_cells, 0.0f);
    
    // Set initial v velocity field if provided
    if (initial_v_field != nullptr)
    {
        if (initial_v_field->size() == v_total_cells)
        {
            v_field_ = *initial_v_field;
            std::cout << "Using custom initial v velocity field." << std::endl;
        }
        else
        {
            std::cerr << "Warning: Initial v field size mismatch. Using zero initialization." << std::endl;
        }
    }
    
    v_pair_ = new TexPair<std::vector<float> >(v_field_, v_temp_);
    
    // Initialize pressure fields
    pressure_field_.resize(p_total_cells, 0.0f);
    pressure_temp_.resize(p_total_cells, 0.0f);
    pressure_pair_ = new TexPair<std::vector<float> >(pressure_field_, pressure_temp_);
    
    // Initialize dye fields
    dye_field_.resize(p_total_cells, 0.0f);
    dye_temp_.resize(p_total_cells, 0.0f);
    
    // Set initial dye field if provided
    if (initial_dye_field != nullptr)
    {
        if (initial_dye_field->size() == p_total_cells)
        {
            dye_field_ = *initial_dye_field;
            std::cout << "Using custom initial dye field." << std::endl;
        }
        else
        {
            std::cerr << "Warning: Initial dye field size mismatch. Using zero initialization." << std::endl;
        }
    }
    
    dye_pair_ = new TexPair<std::vector<float> >(dye_field_, dye_temp_);
    
    divergence_field_.resize(p_total_cells, 0.0f);
}

void FluidSimulator::Initialize(const std::string& output_dir,
                                const std::string& image_dir,
                               const std::vector<float>* initial_u_field,
                               const std::vector<float>* initial_v_field,
                               const std::vector<float>* initial_dye_field)
{
    InitializeFields(initial_u_field, initial_v_field, initial_dye_field);
    
    // Create output directory
    _mkdir(output_dir.c_str());
    _mkdir(image_dir.c_str());
    std::cout << "Output directory: " << output_dir << std::endl;
    std::cout << "image directory: " << image_dir << std::endl;
    
    // Open output files
    velocity_file_.open((output_dir + "/velocity_data.txt").c_str());
    pressure_file_.open((output_dir + "/pressure_data.txt").c_str());
    dye_file_.open((output_dir + "/dye_data.txt").c_str());
    
    // Write headers
    velocity_file_ << "# Velocity field data: frame time x y vel_x vel_y" << std::endl;
    pressure_file_ << "# Pressure field data: frame time x y pressure" << std::endl;
    dye_file_ << "# Dye field data: frame time x y dye_concentration" << std::endl;
    
    std::cout << "Grid size: " << params_.width << "x" << params_.height << std::endl;
    std::cout << "Time step: " << params_.dt << "s" << std::endl;
    std::cout << "Total time: " << params_.total_time << "s" << std::endl;
    std::cout << "Output directory: " << output_dir << std::endl;
}

void FluidSimulator::Step(const std::string& image_dir)
{
    // 1. Apply external forces (gravity, user input, etc.)
    ApplyForces();
    
    // 2. Advect velocity field
    AdvectVelocity(params_.advect_method_);
    
    // 3. TODO: Diffuse velocity field (viscosity)
    
    // 4. Project - Solve pressure to ensure incompressibility
    SolvePressure();
    
    // Advect and diffuse dye field
    AdvectDye();
    // DiffuseDye();
    
    // Apply dissipation to dye field
    DissipateDye();
    
    // 5. Apply boundary conditions and inside source constraints
    ApplyBoundaryConditions();
    ApplyInsideSource();
    
    // Update time and frame
    current_time_ += params_.dt;
    current_frame_++;
    
    // Progress output
    if (current_frame_ % 10 == 0)
    {
        SaveFrameData(image_dir);
        std::cout << "Frame " << save_frame_ 
                  << ", Time: " << std::fixed << std::setprecision(2) << current_time_ 
                  << "s" << std::endl;
        save_frame_++;
    }
}

void FluidSimulator::RunSimulation(const std::string& image_dir)
{
    std::cout << "Starting fluid simulation..." << std::endl;
    
    while (!IsComplete())
    {
        Step(image_dir);
    }
    
    std::cout << "Simulation completed!" << std::endl;
    std::cout << "Total frames: " << current_frame_ << std::endl;
    std::cout << "Final time: " << current_time_ << "s" << std::endl;
}

void FluidSimulator::ApplyForces()
{
    if (params_.enable_gravity)
    {
        apply_force<float>(params_.width, params_.height, u_pair_->cur, v_pair_->cur, -9.8f, params_.dt);
    }
}

void FluidSimulator::InitializeAMGSolver()
{
    if (amg_initialized_) return;
    
    // Create Poisson matrix for pressure solve
    int n = params_.width * params_.height;
    poisson_matrix_ = new FixedSparseMatrix<float>(n, n);
    setupPoissonMatrix2D(*poisson_matrix_, params_.width, params_.height);
    
    // Initialize AMG hierarchy (simplified - just allocate placeholders)
    A_L_.clear();
    R_L_.clear();
    P_L_.clear();
    S_L_.clear();
    
    // For now, just use the original matrix as the only level
    A_L_.push_back(poisson_matrix_);
    S_L_.push_back(Vec2i(params_.width, params_.height));
    
    amg_initialized_ = true;
}

void FluidSimulator::SolvePressure()
{
    // Initialize AMG solver if not already done
    if (!amg_initialized_) {
        InitializeAMGSolver();
    }
    
    // Calculate divergence of velocity field
    get_divergence<float>(params_.width, params_.height, u_pair_->cur, v_pair_->cur, divergence_field_);
    
    // Clear pressure field
    std::fill(pressure_pair_->cur.begin(), pressure_pair_->cur.end(), 0.0f);
    
    // Solve pressure using AMG PCG solver
    float residual_out;
    int iterations_out;
    float tolerance = 1e-6f;
    int max_iterations = params_.pressure_iterations;
    
    // Copy divergence to RHS
    std::vector<float> rhs = divergence_field_;
    for (auto& val : rhs) {
        val = -val;
    }
    
    bool converged = AMGPCGSolvePrebuilt2D<float>(
        *poisson_matrix_,                   // System matrix
        rhs,                                // Right-hand side
        pressure_pair_->cur,                // Solution vector (pressure field)
        A_L_,                               // AMG level matrices
        R_L_,                               // Restriction matrices
        P_L_,                               // Prolongation matrices
        S_L_,                               // Grid sizes for each level
        1,                                  // Total AMG levels
        tolerance,                          // Tolerance factor
        max_iterations,                     // Maximum iterations
        residual_out,                       // Output residual
        iterations_out,                     // Output iteration count
        params_.width,                      // Grid width
        params_.height,                     // Grid height
        false                               // Not pure Neumann
    );
    
    if (!converged) {
        std::cout << "Warning: Pressure solver did not converge. Residual: " << residual_out 
                  << ", Iterations: " << iterations_out << std::endl;
    }
    
    // Subtract pressure gradient from velocity
    subtract_gradient<float>(params_.width, params_.height, u_pair_->cur, v_pair_->cur, pressure_pair_->cur);
}

void FluidSimulator::AdvectVelocity(AdvectMethod method)
{
    switch (method) {
        case AdvectMethod::SemiLagrange:
            advection_velocity<float>(params_.width, params_.height,
                                      u_pair_->cur, v_pair_->cur,
                                      u_pair_->nxt, v_pair_->nxt, params_.dt);
            break;
        case AdvectMethod::MacCormack:
            macCormackVelocity<float>(params_.width, params_.height,
                                      u_pair_->cur, v_pair_->cur,
                                      u_pair_->nxt, v_pair_->nxt, params_.dt);
            break;
        default:
            std::cerr << "Unknown advection method!" << std::endl;
            break;
    }
    u_pair_->Swap();
    v_pair_->Swap();
}

void FluidSimulator::AdvectDye()
{
    advection_dye<float>(params_.width, params_.height, u_pair_->cur, v_pair_->cur,
                              dye_pair_->cur, dye_pair_->nxt, params_.dt);
    dye_pair_->Swap();
}

// void FluidSimulator::DiffuseDye()
// {
//     if (params_.diffusion <= 0.0f) return;
    
//     // Copy current dye to temp for diffusion source
//     dye_temp_ = dye_pair_->cur;
    
//     // Perform Gauss-Seidel iterations for dye diffusion
//     for (int iter = 0; iter < params_.pressure_iterations; ++iter)
//     {
//         diffuse_gauss_seidel<float, float>(params_.width, params_.height, 
//                                           params_.diffusion, params_.dt,
//                                           dye_temp_, dye_pair_->cur);
//     }
// }

void FluidSimulator::DissipateDye()
{
    // Apply dissipation (decay) to dye field
    dissipate<float>(params_.width, params_.height, dye_pair_->cur, params_.dissipation, params_.dt);
}

void FluidSimulator::SaveFrameData(const std::string& image_dir)
{
    if (!velocity_file_.is_open() || !pressure_file_.is_open() || !dye_file_.is_open())
        return;

    const std::vector<Vector2f>& velocity_field = GetVelocityField();
    const std::vector<float>& pressure_field = pressure_pair_->cur;
    const std::vector<float>& dye_field = dye_pair_->cur;
    for (int j = 0; j < params_.height; ++j)
    {
        for (int i = 0; i < params_.width; ++i)
        {
            int id = IXY(i, j, params_.width);
            
            // Save velocity data
            velocity_file_ << current_frame_ << " " << current_time_ << " " 
                          << i << " " << j << " " 
                          << velocity_field[id](0) << " " 
                          << velocity_field[id](1) << std::endl;

            // Save pressure data
            pressure_file_ << current_frame_ << " " << current_time_ << " " 
                          << i << " " << j << " " 
                          << pressure_field[id] << std::endl;

            // Save dye data
            dye_file_ << current_frame_ << " " << current_time_ << " " 
                       << i << " " << j << " " 
                       << dye_field[id] << std::endl;
        }
    }
    
    // Save velocity field as color-mapped image
    SaveVelocityFieldImage(image_dir, velocity_field, current_frame_);
}

void FluidSimulator::ApplyBoundaryConditions()
{
    // Apply no-slip boundary conditions
    
    // Left and right boundaries for u component
    for (int j = 0; j < params_.height; ++j)
    {
        u_pair_->cur[IXY(0, j, params_.width + 1)] = 0.0f;  // Left boundary
        u_pair_->cur[IXY(params_.width, j, params_.width + 1)] = 0.0f;  // Right boundary
    }
    
    // Top and bottom boundaries for v component
    for (int i = 0; i < params_.width; ++i)
    {
        v_pair_->cur[IXY(i, 0, params_.width)] = 0.0f;  // Bottom boundary
        v_pair_->cur[IXY(i, params_.height, params_.width)] = 0.0f;  // Top boundary
    }
}

void FluidSimulator::ApplyInsideSource()
{
    // Apply constant velocity inside the circular source region
    
    for (int i = 0; i < params_.width + 1; ++i)
    {
        for (int j = 0; j < params_.height; ++j)
        {
            // u component is stored at (i, j+0.5)
            float x = static_cast<float>(i);
            float y = static_cast<float>(j) + 0.5f;
            
            float dx = x - params_.source_center_x;
            float dy = y - params_.source_center_y;
            float distance = std::sqrt(dx*dx + dy*dy);
            
            if (distance <= params_.source_radius)
            {
                u_pair_->cur[IXY(i, j, params_.width + 1)] = 0.0f;
            }
        }
    }
    
    for (int i = 0; i < params_.width; ++i)
    {
        for (int j = 0; j < params_.height + 1; ++j)
        {
            // v component is stored at (i+0.5, j)
            float x = static_cast<float>(i) + 0.5f;
            float y = static_cast<float>(j);
            
            float dx = x - params_.source_center_x;
            float dy = y - params_.source_center_y;
            float distance = std::sqrt(dx*dx + dy*dy);
            
            if (distance <= params_.source_radius)
            {
                v_pair_->cur[IXY(i, j, params_.width)] = params_.source_velocity;
            }
        }
    }
}

void FluidSimulator::SaveVelocityFieldImage(const std::string& image_dir, const std::vector<Vector2f>& velocity_field, int current_frame)
{
    // Create RGB data array for velocity magnitude visualization
    std::vector<float> rgb_data(params_.width * params_.height * 3);
    
    // Initialize color ramp for velocity magnitude mapping
    Mfree::ColorRamp color_ramp;
    
    // Find maximum velocity magnitude for normalization
    float max_velocity = 0.0f;
    for (int j = 0; j < params_.height; ++j)
    {
        for (int i = 0; i < params_.width; ++i)
        {
            int id = IXY(i, j, params_.width);
            float magnitude = std::sqrt(velocity_field[id](0) * velocity_field[id](0) + 
                                       velocity_field[id](1) * velocity_field[id](1));
            max_velocity = std::max(max_velocity, magnitude);
        }
    }
    
    // Avoid division by zero
    if (max_velocity < 1e-6f) max_velocity = 1.0f;
    
    // Generate color-mapped image data
    for (int j = 0; j < params_.height; ++j)
    {
        for (int i = 0; i < params_.width; ++i)
        {
            int id = IXY(i, j, params_.width);
            int rgb_id = ((params_.height - 1 - j) * params_.width + i) * 3; // Flip Y for image output
            
            // Calculate velocity magnitude
            float magnitude = std::sqrt(velocity_field[id](0) * velocity_field[id](0) + 
                                       velocity_field[id](1) * velocity_field[id](1));
            
            // Normalize velocity magnitude to [0, 1]
            float normalized_velocity = magnitude / max_velocity;
            
            // Get color from color ramp using JET colormap
            Mfree::vec3 color(0.0f, 0.0f, 0.0f);
            color_ramp.set_GLcolor(normalized_velocity, Mfree::COLOR_JET, color, false);
            
            // Store RGB values
            rgb_data[rgb_id] = color.x;
            rgb_data[rgb_id + 1] = color.y;
            rgb_data[rgb_id + 2] = color.z;
        }
    }
    
    // Save as PPM image with frame number in filename
    char filename[512];
    sprintf(filename, "%s/velocity_magnitude_frame_%05d.ppm", image_dir.c_str(), current_frame);
    
    std::ofstream ppm_file(filename, std::ios::out | std::ios::binary);
    if (ppm_file.is_open())
    {
        ppm_file << "P6\n" << params_.width << " " << params_.height << "\n255\n";
        for (int i = 0; i < params_.width * params_.height * 3; ++i) {
            ppm_file << (unsigned char)(rgb_data[i] * 255);
        }
        ppm_file.close();
    }
    else
    {
        std::cerr << "Error: Could not create velocity image file: " << filename << std::endl;
    }
}
