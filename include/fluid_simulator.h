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
    std::vector<Vector2f> u_field_;
    std::vector<Vector2f> u_temp_;
    std::vector<Vector2f> v_field_;
    std::vector<Vector2f> v_temp_;

    // Scalar fields stored at cell centers 
    std::vector<float> pressure_field_;
    std::vector<float> pressure_temp_;
    std::vector<float> dye_field_;
    std::vector<float> dye_temp_;
    std::vector<float> divergence_field_;
    
    // TexPair wrappers for easy swapping
    TexPair<std::vector<Vector2f> >* u_pair_;
    TexPair<std::vector<Vector2f> >* v_pair_;
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
     */
    void Initialize(const std::string& output_dir = "output");
    
    /**
     * @brief Run one simulation step
     */
    void Step();
    
    /**
     * @brief Run the complete simulation
     */
    void RunSimulation();
    
    /**
     * @brief Add a source event to the simulation
     * @param time Time when the source should be activated
     * @param x X position
     * @param y Y position
     * @param radius Source radius
     * @param direction Velocity direction
     * @param strength Source strength
     * @param dye_amount Amount of dye to add
     */
    void AddSourceEvent(float time, int x, int y, int radius, 
                       const Vector2f& direction, float strength, float dye_amount = 1.0f);
    
    float GetCurrentTime() const { return current_time_; }
    
    int GetCurrentFrame() const { return current_frame_; }
    
    bool IsComplete() const { return current_time_ >= params_.total_time; }
    
    const SimulationParams& GetParams() const { return params_; }
    
    const std::vector<Vector2f>& GetVelocityField() const { return velocity_pair_->cur; }
    
    const std::vector<float>& GetPressureField() const { return pressure_pair_->cur; }
    
    const std::vector<float>& GetDyeField() const { return dye_pair_->cur; }

private:
    /**
     * @brief Initialize all fields
     */
    void InitializeFields();
    
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
    void DiffuseDye();
    
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
