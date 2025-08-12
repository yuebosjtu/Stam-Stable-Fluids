# ifndef UTILS_H
# define UTILS_H

# include <iostream>
# include <vector>
# include <cmath>

void CreateCircularSourceField(float width, float height, int source_center_x, int source_center_y,
                               float source_radius, float source_velocity,
                               std::vector<float>& u_field, std::vector<float>& v_field)
{
    int u_total_cells = (width + 1) * height;
    int v_total_cells = width * (height + 1);

    u_field.resize(u_total_cells, 0.0f);
    v_field.resize(v_total_cells, 0.0f);
    
    // Set v velocity in circular region (upward flow)
    for (int i = 0; i < width; ++i)
    {
        for (int j = 0; j < height + 1; ++j)
        {
            // v component is stored at (i+0.5, j)
            float x = static_cast<float>(i) + 0.5f;
            float y = static_cast<float>(j);
            
            float dx = x - source_center_x;
            float dy = y - source_center_y;
            float distance = std::sqrt(dx*dx + dy*dy);
            
            if (distance <= source_radius)
            {
                v_field[width * j + i] = source_velocity;
            }
        }
    }

    std::cout << "Created circular source field at (" << source_center_x << ", " 
              << source_center_y << ") with radius " << source_radius 
              << " and velocity " << source_velocity << std::endl;
}

# endif // UTILS_H