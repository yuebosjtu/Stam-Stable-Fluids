#ifndef FLUID_SIMULATOR_H
#define FLUID_SIMULATOR_H

#include "stable_fluid.h"
#include <Eigen/Dense>
#include <string>
#include <fstream>
#include <memory>
#include <vector>

/**
 * @brief Fluid simulation class implementing Jos Stam's stable fluid method
 */
class FluidSimulator
{
public:
    typedef Eigen::Vector2f Vector2f;
    typedef Eigen::VectorXf VectorXf;
    
    struct SimulationParams
    {
        int width;                    // Grid width
        int height;                   // Grid height
        float dt;                     // Time step
        float viscosity;              // Fluid viscosity
        float diffusion;              // Dye diffusion constant
        float dissipation;            // Dye dissipation rate
        int pressure_iterations;      // Gauss-Seidel iterations for pressure
        float total_time;             // Total simulation time
        bool enable_gravity;          // Enable gravity
        
        // Constructor with default values
        SimulationParams() : 
            width(64), height(64), dt(0.016f), viscosity(0.0001f), diffusion(0.0001f),
            dissipation(0.001f), pressure_iterations(80), total_time(10.0f), enable_gravity(true) {}
    };

private:
    SimulationParams params_;
    float current_time_;
    int current_frame_;
    
    // Velocity components fields
    std::vector<float> u_field_;
    std::vector<float> u_temp_;
    std::vector<float> v_field_;
    std::vector<float> v_temp_;

    // Scalar fields stored at cell centers 
    std::vector<float> pressure_field_;
    std::vector<float> pressure_temp_;
    std::vector<float> dye_field_;
    std::vector<float> dye_temp_;
    std::vector<float> divergence_field_;
    
    // TexPair wrappers for easy swapping
    TexPair<std::vector<float> >* u_pair_;
    TexPair<std::vector<float> >* v_pair_;
    TexPair<std::vector<float> >* pressure_pair_;
    TexPair<std::vector<float> >* dye_pair_;
    
    // Output files for saving data
    std::ofstream velocity_file_;
    std::ofstream pressure_file_;
    std::ofstream dye_file_;

public:
    explicit FluidSimulator(const SimulationParams& params = SimulationParams());
    
    ~FluidSimulator();
    
    /**
     * @brief Initialize the simulation
     * @param output_dir Directory to save output files
     * @param initial_u_field Optional initial u velocity field (size: (width+1) x height)
     * @param initial_v_field Optional initial v velocity field (size: width x (height+1))
     * @param initial_dye_field Optional initial dye field (size: width x height)
     */
    void Initialize(const std::string& output_dir = "output",
                   const std::vector<float>* initial_u_field = nullptr,
                   const std::vector<float>* initial_v_field = nullptr,
                   const std::vector<float>* initial_dye_field = nullptr);
    
    /**
     * @brief Run one simulation step
     */
    void Step();
    
    /**
     * @brief Run the complete simulation
     */
    void RunSimulation();
    
    float GetCurrentTime() const { return current_time_; }
    
    int GetCurrentFrame() const { return current_frame_; }
    
    bool IsComplete() const { return current_time_ >= params_.total_time; }
    
    const SimulationParams& GetParams() const { return params_; }
    
    const std::vector<float>& GetVelocityField_u() const { return u_field_; }

    const std::vector<float>& GetVelocityField_v() const { return v_field_; }

    const std::vector<Vector2f>& GetVelocityField() const
    {
        static std::vector<Vector2f> velocity_field;
        velocity_field.resize(params_.width * params_.height);
        
        for (int j = 0; j < params_.height; ++j)
        {
            for (int i = 0; i < params_.width; ++i)
            {
                int idx = IXY(i, j, params_.width);
                
                // Average u velocities at left and right faces
                float u_left = u_field_[IXY(i, j, params_.width + 1)];
                float u_right = u_field_[IXY(i + 1, j, params_.width + 1)];
                float u_avg = 0.5f * (u_left + u_right);
                
                // Average v velocities at bottom and top faces
                float v_bottom = v_field_[IXY(i, j, params_.width)];
                float v_top = v_field_[IXY(i, j + 1, params_.width)];
                float v_avg = 0.5f * (v_bottom + v_top);
                
                velocity_field[idx] = Vector2f(u_avg, v_avg);
            }
        }
        
        return velocity_field;
    }

    const std::vector<float>& GetPressureField() const { return pressure_field_; }

    const std::vector<float>& GetDyeField() const { return dye_field_; }

private:
    /**
     * @brief Initialize all fields
     * @param initial_u_field Optional initial u velocity field (size: (width+1) x height). If nullptr, zero initialized.
     * @param initial_v_field Optional initial v velocity field (size: width x (height+1)). If nullptr, zero initialized.
     * @param initial_dye_field Optional initial dye field (size: width x height). If nullptr, zero initialized.
     */
    void InitializeFields(const std::vector<float>* initial_u_field = nullptr,
                         const std::vector<float>* initial_v_field = nullptr,
                         const std::vector<float>* initial_dye_field = nullptr);
    
    /**
     * @brief Apply external forces (gravity, user input, etc.)
     */
    void ApplyForces();
    
    /**
     * @brief Solve for pressure to make velocity field divergence-free
     */
    void SolvePressure();
    
    /**
     * @brief Advect velocity field
     */
    void AdvectVelocity();
    
    /**
     * @brief Advect dye field
     */
    void AdvectDye();
    
    /**
     * @brief Diffuse dye field
     */
    // void DiffuseDye();
    
    /**
     * @brief Dissipate (decay) dye field
     */
    void DissipateDye();
    
    /**
     * @brief Save current frame data to files
     */
    void SaveFrameData();
    
    /**
     * @brief Create output directory if it doesn't exist
     */
    void CreateOutputDirectory(const std::string& dir);
};

#endif // FLUID_SIMULATOR_H
