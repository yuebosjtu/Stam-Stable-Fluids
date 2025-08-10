#ifndef STABLE_FLUID
#define STABLE_FLUID

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

/**
 * @brief Data aggregation of current and next frame
 * @tparam T the type of data stored in aggregation
 */
template <typename T> struct TexPair
{
  public:
    T &cur;
    T &nxt;
    TexPair(T &a, T &b) : cur(a), nxt(b) {}
    void Swap()
    {
        std::swap(cur, nxt);
    }
};

template <typename T> inline T min(T a, T b)
{
    return a < b ? a : b;
}
template <typename T> inline T max(T a, T b)
{
    return a > b ? a : b;
}
/**
 * @brief Get index of 1d array out of 2d index (i,j)
 *
 * @param i index of x-axis
 * @param j index of y-axis
 * @param N width of 2d array
 */
inline int IXY(int i, int j, int N)
{
    return N * j + i;
}

template <typename T, typename SCALAR> inline T lerp(T l, T r, SCALAR t)
{
    return l + t * (r - l);
}

/**
 * @brief Get qf[u, v]
 */
template <typename T, typename SCALAR>
T sample(const int N, const int M, const std::vector<T> &qf, const SCALAR u, const SCALAR v)
{
    int x = static_cast<int>(u);
    x = max<int>(0, min<int>(N - 1, x));
    int y = static_cast<int>(v);
    y = max<int>(0, min<int>(M - 1, y));
    return qf[IXY(x, y, N)];
}

/**
 * @brief Bilinear interpolation
 *
 * @param N width
 * @param M height
 * @param p position to interpolate at
 */
template <typename T, typename VEC, typename SCALAR = typename VEC::Scalar>
T bilerp(const int N, const int M, const std::vector<T> &qf, const VEC &p)
{
    SCALAR s = p(0) - 0.5f;
    SCALAR t = p(1) - 0.5f;
    SCALAR iu = floor(s);
    SCALAR iv = floor(t);
    SCALAR fu = s - iu;
    SCALAR fv = t - iv;
    T a = sample<T, SCALAR>(N, M, qf, iu, iv);
    T b = sample<T, SCALAR>(N, M, qf, iu + 1, iv);
    T c = sample<T, SCALAR>(N, M, qf, iu, iv + 1);
    T d = sample<T, SCALAR>(N, M, qf, iu + 1, iv + 1);
    return lerp<T, SCALAR>(lerp<T, SCALAR>(a, b, fu), lerp<T, SCALAR>(c, d, fu), fv);
}

template <typename VEC, typename SCALAR = typename VEC::Scalar>
SCALAR bilerp_velocity(const int N, const int M, const std::vector<VEC> &qf, const VEC &p, int index)
{
    SCALAR s, t;
    if (index == 0) // Sample horizontal velocity
    {
        s = p(0);
        t = p(1) - 0.5f;
    }
    else if (index == 1) // Sample vertival velocity
    {
        s = p(0) - 0.5f;
        t = p(1);
    }
    else
    {
        printf("Error! No accessible index\n");
        assert(false);
        return 0.0f;
    }
    SCALAR iu = floor(s);
    SCALAR iv = floor(t);
    SCALAR fu = s - iu;
    SCALAR fv = t - iv;
    VEC a = sample<VEC, SCALAR>(N, M, qf, iu, iv);
    VEC b = sample<VEC, SCALAR>(N, M, qf, iu + 1, iv);
    VEC c = sample<VEC, SCALAR>(N, M, qf, iu, iv + 1);
    VEC d = sample<VEC, SCALAR>(N, M, qf, iu + 1, iv + 1);
    return (lerp<VEC, SCALAR>(lerp<VEC, SCALAR>(a, b, fu), lerp<VEC, SCALAR>(c, d, fu), fv))(index);
}

/**
 * @brief Locate which point will move to position p in next dt time
 *
 * @param N width
 * @param M height
 * @param vel velocity field
 */
template <typename VEC, typename SCALAR = typename VEC::Scalar>
VEC backtrace(const int N, const int M, VEC p, SCALAR dt, const std::vector<VEC> &vel)
{
    VEC v1(bilerp<VEC, VEC>(N, M, vel, p));
    VEC p1(p(0) - 0.5 * dt * v1(0), p(1) - 0.5 * dt * v1(1));
    VEC v2(bilerp<VEC, VEC>(N, M, vel, p1));
    VEC p2(p(0) - 0.75 * dt * v2(0), p(1) - 0.75 * dt * v2(1));
    VEC v3(bilerp<VEC, VEC>(N, M, vel, p2));
    p = p + (-1.f) * dt * ((2.f / 9.f) * v1 + (1.f / 3.f) * v2 + (4.f / 9.f) * v3);
    return p * 0.9999f;
}

template <typename SCALAR>
void backtrace_u(const int N, const int M, SCALAR x, SCALAR y, SCALAR dt,
                     const std::vector<SCALAR> &u_vel, const std::vector<SCALAR> &v_vel,
                     SCALAR &back_x, SCALAR &back_y)
{
    // RK3 integration for u component at position (x, y+0.5)
    // Stage 1
    SCALAR u1 = (x > 0.5f && x < N - 0.5f) ? 0.5f * (u_vel[IXY(static_cast<int>(x-0.5f), static_cast<int>(y), N+1)] + 
                                                     u_vel[IXY(static_cast<int>(x+0.5f), static_cast<int>(y), N+1)]) : 0.0f;
    SCALAR v1 = 0.0f;
    if (y > 0.5f && y < M - 0.5f) {
        int ix = static_cast<int>(x);
        int iy = static_cast<int>(y);
        v1 = 0.25f * (v_vel[IXY(max(0, ix-1), iy, N)] + v_vel[IXY(min(N-1, ix), iy, N)] + 
                      v_vel[IXY(max(0, ix-1), iy+1, N)] + v_vel[IXY(min(N-1, ix), iy+1, N)]);
    }
    
    SCALAR x1 = x - 0.5f * dt * u1;
    SCALAR y1 = y - 0.5f * dt * v1;
    
    // Stage 2
    SCALAR u2 = (x1 > 0.5f && x1 < N - 0.5f) ? 0.5f * (u_vel[IXY(static_cast<int>(x1-0.5f), static_cast<int>(y1), N+1)] + 
                                                        u_vel[IXY(static_cast<int>(x1+0.5f), static_cast<int>(y1), N+1)]) : 0.0f;
    SCALAR v2 = 0.0f;
    if (y1 > 0.5f && y1 < M - 0.5f) {
        int ix = static_cast<int>(x1);
        int iy = static_cast<int>(y1);
        v2 = 0.25f * (v_vel[IXY(max(0, ix-1), iy, N)] + v_vel[IXY(min(N-1, ix), iy, N)] + 
                      v_vel[IXY(max(0, ix-1), iy+1, N)] + v_vel[IXY(min(N-1, ix), iy+1, N)]);
    }
    
    SCALAR x2 = x - 0.75f * dt * u2;
    SCALAR y2 = y - 0.75f * dt * v2;
    
    // Stage 3
    SCALAR u3 = (x2 > 0.5f && x2 < N - 0.5f) ? 0.5f * (u_vel[IXY(static_cast<int>(x2-0.5f), static_cast<int>(y2), N+1)] + 
                                                        u_vel[IXY(static_cast<int>(x2+0.5f), static_cast<int>(y2), N+1)]) : 0.0f;
    SCALAR v3 = 0.0f;
    if (y2 > 0.5f && y2 < M - 0.5f) {
        int ix = static_cast<int>(x2);
        int iy = static_cast<int>(y2);
        v3 = 0.25f * (v_vel[IXY(max(0, ix-1), iy, N)] + v_vel[IXY(min(N-1, ix), iy, N)] + 
                      v_vel[IXY(max(0, ix-1), iy+1, N)] + v_vel[IXY(min(N-1, ix), iy+1, N)]);
    }
    
    // Final position
    back_x = x - dt * ((2.0f/9.0f) * u1 + (1.0f/3.0f) * u2 + (4.0f/9.0f) * u3);
    back_y = y - dt * ((2.0f/9.0f) * v1 + (1.0f/3.0f) * v2 + (4.0f/9.0f) * v3);
}

template <typename SCALAR>
void backtrace_v(const int N, const int M, SCALAR x, SCALAR y, SCALAR dt,
                     const std::vector<SCALAR> &u_vel, const std::vector<SCALAR> &v_vel,
                     SCALAR &back_x, SCALAR &back_y)
{
    // RK3 integration for v component at position (x+0.5, y)
    // Stage 1
    SCALAR u1 = 0.0f;
    if (x > 0.5f && x < N - 0.5f) {
        int ix = static_cast<int>(x);
        int iy = static_cast<int>(y);
        u1 = 0.25f * (u_vel[IXY(ix, max(0, iy-1), N+1)] + u_vel[IXY(ix+1, max(0, iy-1), N+1)] + 
                      u_vel[IXY(ix, min(M-1, iy), N+1)] + u_vel[IXY(ix+1, min(M-1, iy), N+1)]);
    }
    SCALAR v1 = (y > 0.5f && y < M - 0.5f) ? 0.5f * (v_vel[IXY(static_cast<int>(x), static_cast<int>(y-0.5f), N)] + 
                                                     v_vel[IXY(static_cast<int>(x), static_cast<int>(y+0.5f), N)]) : 0.0f;
    
    SCALAR x1 = x - 0.5f * dt * u1;
    SCALAR y1 = y - 0.5f * dt * v1;
    
    // Stage 2
    SCALAR u2 = 0.0f;
    if (x1 > 0.5f && x1 < N - 0.5f) {
        int ix = static_cast<int>(x1);
        int iy = static_cast<int>(y1);
        u2 = 0.25f * (u_vel[IXY(ix, max(0, iy-1), N+1)] + u_vel[IXY(ix+1, max(0, iy-1), N+1)] + 
                      u_vel[IXY(ix, min(M-1, iy), N+1)] + u_vel[IXY(ix+1, min(M-1, iy), N+1)]);
    }
    SCALAR v2 = (y1 > 0.5f && y1 < M - 0.5f) ? 0.5f * (v_vel[IXY(static_cast<int>(x1), static_cast<int>(y1-0.5f), N)] + 
                                                        v_vel[IXY(static_cast<int>(x1), static_cast<int>(y1+0.5f), N)]) : 0.0f;
    
    SCALAR x2 = x - 0.75f * dt * u2;
    SCALAR y2 = y - 0.75f * dt * v2;
    
    // Stage 3
    SCALAR u3 = 0.0f;
    if (x2 > 0.5f && x2 < N - 0.5f) {
        int ix = static_cast<int>(x2);
        int iy = static_cast<int>(y2);
        u3 = 0.25f * (u_vel[IXY(ix, max(0, iy-1), N+1)] + u_vel[IXY(ix+1, max(0, iy-1), N+1)] + 
                      u_vel[IXY(ix, min(M-1, iy), N+1)] + u_vel[IXY(ix+1, min(M-1, iy), N+1)]);
    }
    SCALAR v3 = (y2 > 0.5f && y2 < M - 0.5f) ? 0.5f * (v_vel[IXY(static_cast<int>(x2), static_cast<int>(y2-0.5f), N)] + 
                                                        v_vel[IXY(static_cast<int>(x2), static_cast<int>(y2+0.5f), N)]) : 0.0f;
    
    // Final position
    back_x = x - dt * ((2.0f/9.0f) * u1 + (1.0f/3.0f) * u2 + (4.0f/9.0f) * u3);
    back_y = y - dt * ((2.0f/9.0f) * v1 + (1.0f/3.0f) * v2 + (4.0f/9.0f) * v3);
}

/**
 * @brief Bilinear interpolation for u component in MAC grid
 */
template <typename SCALAR>
SCALAR interpolate_u(const int N, const int M, const std::vector<SCALAR> &u_vel, SCALAR x, SCALAR y)
{
    // u component is stored at (i, j+0.5), so y needs adjustment
    SCALAR adjusted_y = y - 0.5f;
    
    int ix = static_cast<int>(x);
    int iy = static_cast<int>(adjusted_y);
    SCALAR fx = x - ix;
    SCALAR fy = adjusted_y - iy;
    
    ix = max<int>(0, min<int>(N, ix));
    iy = max<int>(0, min<int>(M-1, iy));
    
    SCALAR u00 = (ix <= N && iy < M) ? u_vel[IXY(ix, iy, N+1)] : 0.0f;
    SCALAR u10 = (ix+1 <= N && iy < M) ? u_vel[IXY(ix+1, iy, N+1)] : 0.0f;
    SCALAR u01 = (ix <= N && iy+1 < M) ? u_vel[IXY(ix, iy+1, N+1)] : 0.0f;
    SCALAR u11 = (ix+1 <= N && iy+1 < M) ? u_vel[IXY(ix+1, iy+1, N+1)] : 0.0f;
    
    SCALAR u_bottom = u00 * (1.0f - fx) + u10 * fx;
    SCALAR u_top = u01 * (1.0f - fx) + u11 * fx;
    
    return u_bottom * (1.0f - fy) + u_top * fy;
}

/**
 * @brief Bilinear interpolation for v component in MAC grid
 */
template <typename SCALAR>
SCALAR interpolate_v(const int N, const int M, const std::vector<SCALAR> &v_vel, SCALAR x, SCALAR y)
{
    // v component is stored at (i+0.5, j), so x needs adjustment
    SCALAR adjusted_x = x - 0.5f;
    
    int ix = static_cast<int>(adjusted_x);
    int iy = static_cast<int>(y);
    SCALAR fx = adjusted_x - ix;
    SCALAR fy = y - iy;
    
    ix = max<int>(0, min<int>(N-1, ix));
    iy = max<int>(0, min<int>(M, iy));
    
    SCALAR v00 = (ix < N && iy <= M) ? v_vel[IXY(ix, iy, N)] : 0.0f;
    SCALAR v10 = (ix+1 < N && iy <= M) ? v_vel[IXY(ix+1, iy, N)] : 0.0f;
    SCALAR v01 = (ix < N && iy+1 <= M) ? v_vel[IXY(ix, iy+1, N)] : 0.0f;
    SCALAR v11 = (ix+1 < N && iy+1 <= M) ? v_vel[IXY(ix+1, iy+1, N)] : 0.0f;
    
    SCALAR v_bottom = v00 * (1.0f - fx) + v10 * fx;
    SCALAR v_top = v01 * (1.0f - fx) + v11 * fx;
    
    return v_bottom * (1.0f - fy) + v_top * fy;
}

/**
 * @brief update advection of fluid property qf
 * @param N width
 * @param M height
 * @param vel velocity
 * @param dt time interval
 */
template <typename T, typename VEC, typename SCALAR = typename VEC::Scalar>
void advection(const int N, const int M, const std::vector<VEC> &vel, const std::vector<T> &qf,
               std::vector<T> &new_qf, SCALAR dt, SCALAR damping = 0.9999f)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            VEC p(i + .5f, j + .5f);
            p = backtrace<VEC>(N, M, p, dt, vel) * damping;

            new_qf[IXY(i, j, N)] = bilerp<T, VEC>(N, M, qf, p);
        }
    }
}

/**
 * @brief Velocity advection for MAC grid
 * @param N width
 * @param M height
 * @param u_vel current u component field
 * @param v_vel current v component field
 * @param new_u_vel new u component field after advection
 * @param new_v_vel new v component field after advection
 * @param dt time step
 * @param damping damping factor
 */
template <typename SCALAR>
void advection_velocity(const int N, const int M, 
                       const std::vector<SCALAR> &u_vel, const std::vector<SCALAR> &v_vel,
                       std::vector<SCALAR> &new_u_vel, std::vector<SCALAR> &new_v_vel,
                       SCALAR dt, SCALAR damping = 0.9999f)
{
    // Advect u component (stored at cell faces i+0.5, j)
    for (int i = 0; i < N + 1; i++)
        for (int j = 0; j < M; j++)
        {
            // Boundary conditions for u component
            if (i == 0 || i == N)
            {
                new_u_vel[IXY(i, j, N + 1)] = 0.0f;  // No-slip at left/right boundaries
            }
            else
            {
                SCALAR curr_x = static_cast<SCALAR>(i);
                SCALAR curr_y = static_cast<SCALAR>(j) + 0.5f;
                SCALAR back_x, back_y;
                backtrace_u<SCALAR>(N, M, curr_x, curr_y, dt, u_vel, v_vel, back_x, back_y);
                
                back_x = max<SCALAR>(0.0f, min<SCALAR>(static_cast<SCALAR>(N), back_x));
                back_y = max<SCALAR>(0.5f, min<SCALAR>(static_cast<SCALAR>(M) - 0.5f, back_y));
                new_u_vel[IXY(i, j, N + 1)] = interpolate_u<SCALAR>(N, M, u_vel, back_x, back_y) * damping;
            }
        }
    
    // Advect v component (stored at cell faces i, j+0.5)
    for (int i = 0; i < N; i++)
        for (int j = 0; j < M + 1; j++)
        {
            // Boundary conditions for v component
            if (j == 0 || j == M)
            {
                new_v_vel[IXY(i, j, N)] = 0.0f;  // No-slip at top/bottom boundaries
            }
            else
            {
                SCALAR curr_x = static_cast<SCALAR>(i) + 0.5f;
                SCALAR curr_y = static_cast<SCALAR>(j);
                SCALAR back_x, back_y;
                backtrace_v<SCALAR>(N, M, curr_x, curr_y, dt, u_vel, v_vel, back_x, back_y);
                
                back_x = max<SCALAR>(0.5f, min<SCALAR>(static_cast<SCALAR>(N) - 0.5f, back_x));
                back_y = max<SCALAR>(0.0f, min<SCALAR>(static_cast<SCALAR>(M), back_y));
                new_v_vel[IXY(i, j, N)] = interpolate_v<SCALAR>(N, M, v_vel, back_x, back_y) * damping;
            }
        }
}

/**
 * @brief apply global forces and damping
 * @param vel velocity field
 * @param g scale of gravity(-9.8f as default)
 * @param dt time interval
 */
template <typename SCALAR>
void apply_force(const int N, const int M, std::vector<SCALAR> &u_vel, std::vector<SCALAR> &v_vel,
                 const SCALAR g = -9.8f, const SCALAR dt = 0.03f)
{
    int v_total = N * (M + 1);
    for (int i = 0; i < v_total; i++)
    {
        v_vel[i] += g * 1.0 * dt;
    }
}

/**
 * @brief Get the divergence of velocity field for MAC grid
 * @param N width
 * @param M height  
 * @param u_vel u component field (size: (N+1) x M)
 * @param v_vel v component field (size: N x (M+1))
 * @param divergence divergence of velocity field (size: N x M)
 */
template <typename SCALAR>
void get_divergence(const int N, const int M, const std::vector<SCALAR> &u_vel, 
                   const std::vector<SCALAR> &v_vel, std::vector<SCALAR> &divergence)
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
        {            
            SCALAR u_right = u_vel[IXY(i + 1, j, N + 1)];
            SCALAR u_left = u_vel[IXY(i, j, N + 1)];
            SCALAR v_top = v_vel[IXY(i, j + 1, N)];
            SCALAR v_bottom = v_vel[IXY(i, j, N)];
            
            divergence[IXY(i, j, N)] = (u_right - u_left) + (v_top - v_bottom);
        }
}

/**
 * @brief Single iteration of gauss-sidel method for diffusion
 * @param N width
 * @param M height
 * @param diffusion_rate diffusion coefficient
 * @param dt time step
 * @param old_field previous field values
 * @param new_field new field values after diffusion step
 */
template <typename T, typename SCALAR>
void diffuse_gauss_seidel(const int N, const int M, const SCALAR diffusion_rate, const SCALAR dt,
                         const std::vector<T> &old_field, std::vector<T> &new_field)
{
    SCALAR a = dt * diffusion_rate * N * M;
    SCALAR c = 1.0f + 4.0f * a;
    
    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
        {
            T left = (i > 0) ? new_field[IXY(i-1, j, N)] : T();
            T right = (i < N-1) ? new_field[IXY(i+1, j, N)] : T();
            T bottom = (j > 0) ? new_field[IXY(i, j-1, N)] : T();
            T top = (j < M-1) ? new_field[IXY(i, j+1, N)] : T();
            
            new_field[IXY(i, j, N)] = (old_field[IXY(i, j, N)] + a * (left + right + bottom + top)) / c;
        }
}

/**
 * @brief Single iteration of gauss-sidel method
 * @param N width
 * @param M height
 * @param pressure pressure field
 * @param divergence divergence of velocity field
 */
template <typename VEC, typename SCALAR = typename VEC::Scalar>
void pressure_gauss_sidel(const int N, const int M, const std::vector<SCALAR> &divergence,
                          const std::vector<SCALAR> &pressure, std::vector<SCALAR> &new_pressure)
{
    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
        {
            SCALAR pl = sample<SCALAR, SCALAR>(N, M, pressure, i - 1, j);
            SCALAR pr = sample<SCALAR, SCALAR>(N, M, pressure, i + 1, j);
            SCALAR pb = sample<SCALAR, SCALAR>(N, M, pressure, i, j - 1);
            SCALAR pt = sample<SCALAR, SCALAR>(N, M, pressure, i, j + 1);
            SCALAR diver = divergence[IXY(i, j, N)];
            new_pressure.at(IXY(i, j, N)) = (pl + pr + pb + pt + (-1.f) * diver) * 0.25f;
        }
}

/**
 * @brief Apply pressure gradient to velocity field for MAC grid
 * @param N width
 * @param M height
 * @param u_vel u component field
 * @param v_vel v component field
 * @param pressure pressure field
 */
template <typename SCALAR>
void subtract_gradient(const int N, const int M, std::vector<SCALAR> &u_vel, 
                      std::vector<SCALAR> &v_vel, const std::vector<SCALAR> &pressure)
{
    // Update u component: u[i,j] -= (pressure[i,j] - pressure[i-1,j])
    for (int i = 1; i < N; i++)  // Start from 1, skip boundary at i=0
        for (int j = 0; j < M; j++)
        {
            SCALAR p_right = pressure[IXY(i, j, N)];
            SCALAR p_left = pressure[IXY(i-1, j, N)];
            u_vel[IXY(i, j, N+1)] -= (p_right - p_left);
        }
    
    // Update v component: v[i,j] -= (pressure[i,j] - pressure[i,j-1])
    for (int i = 0; i < N; i++)
        for (int j = 1; j < M; j++)  // Start from 1, skip boundary at j=0
        {
            SCALAR p_top = pressure[IXY(i, j, N)];
            SCALAR p_bottom = pressure[IXY(i, j-1, N)];
            v_vel[IXY(i, j, N)] -= (p_top - p_bottom);
        }
}

template <typename SCALAR> SCALAR smooth_step(SCALAR a, SCALAR x)
{
    SCALAR y = (a - x) / a;
    if (y < 0.0)
        y = 0.0;
    if (y > 1.0)
        y = 1.0;
    SCALAR rst = y * y;
    return rst;
}

/**
 * @brief Apply dissipation to scalar field
 * @param N width
 * @param M height
 * @param field scalar field to dissipate
 * @param dissipation_rate dissipation coefficient
 * @param dt time step
 */
template <typename SCALAR>
void dissipate(const int N, const int M, std::vector<SCALAR> &field, const SCALAR dissipation_rate, const SCALAR dt)
{
    SCALAR decay_factor = 1.0f - dissipation_rate * dt;
    decay_factor = max<SCALAR>(0.0f, decay_factor); // Ensure non-negative
    
    for (int i = 0; i < N * M; ++i)
    {
        field[i] *= decay_factor;
    }
}

/**
 * @brief add round-shape fluid momentum at point [x,y]
 *
 * @param N width
 * @param x x-position
 * @param y y-position
 * @param r radius
 * @param dir direction of added velocity
 * @param dye color
 * @param vel velocity field
 * @param value scale of added color and velocity
 */
template <typename VEC, typename SCALAR = typename VEC::Scalar>
void add_source(const int N, int x, int y, int r, SCALAR value, VEC dir,
                std::vector<SCALAR> &dye, std::vector<VEC> &vel)
{
    for (int i = -r; i <= r; i++)
        for (int j = -r; j <= r; j++)
        {
            int index = IXY(x + i, y + j, N);
            SCALAR smooth = smooth_step<SCALAR>(r * r, i * i + j * j);
            smooth *= value;
            if (index < 0 || index >= dye.size())
                printf("Error info: index out of range {%d, %d, %d, %d, %d}\n", x, y, i, j, r);
            dye[index] = min(smooth + dye[index], 3.0f);
            vel[index] += dir * smooth * 100.0f;
        }
}

#endif