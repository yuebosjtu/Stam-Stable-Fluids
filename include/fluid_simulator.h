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
        float vorticity_strength;     // Vorticity confinement strength
        int pressure_iterations;      // Gauss-Seidel iterations for pressure
        float total_time;             // Total simulation time
        bool enable_vorticity;        // Enable vorticity confinement
        bool enable_gravity;          // Enable gravity
        
        // Constructor with default values
        SimulationParams() : 
            width(64), height(64), dt(0.016f), viscosity(0.0001f), diffusion(0.0001f),
            dissipation(0.001f), vorticity_strength(0.1f), pressure_iterations(20), total_time(10.0f),
            enable_vorticity(true), enable_gravity(true) {}
    };

    struct SourceEvent
    {
        float time;              // Time to add source
        int x, y;               // Position
        int radius;             // Source radius
        Vector2f direction;     // Velocity direction
        float strength;         // Source strength
        float dye_amount;       // Amount of dye to add
        
        SourceEvent(float t, int px, int py, int r, Vector2f dir, float s, float dye = 1.0f)
            : time(t), x(px), y(py), radius(r), direction(dir), strength(s), dye_amount(dye) {}
    };

private:
    SimulationParams params_;
    float current_time_;
    int current_frame_;
    
    // Fluid fields
    std::vector<Vector2f> velocity_field_;
    std::vector<Vector2f> velocity_temp_;
    std::vector<float> pressure_field_;
    std::vector<float> pressure_temp_;
    std::vector<float> dye_field_;
    std::vector<float> dye_temp_;
    std::vector<float> divergence_field_;
    std::vector<float> curl_field_;
    std::vector<Vector2f> curl_force_;
    
    // TexPair wrappers for easy swapping
    TexPair<std::vector<Vector2f> >* velocity_pair_;
    TexPair<std::vector<float> >* pressure_pair_;
    TexPair<std::vector<float> >* dye_pair_;
    
    // Source events
    std::vector<SourceEvent> source_events_;
    size_t next_source_index_;
    
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
     * @brief Apply vorticity confinement
     */
    void ApplyVorticityConfinement();
    
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
     * @brief Diffuse velocity field
     */
    void DiffuseVelocity();
    
    /**
     * @brief Diffuse dye field
     */
    void DiffuseDye();
    
    /**
     * @brief Dissipate (decay) dye field
     */
    void DissipateDye();
    
    /**
     * @brief Process source events at current time
     */
    void ProcessSources();
    
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
