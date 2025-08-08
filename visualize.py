import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import Normalize
import os
import argparse

class FluidVisualizer:
    """
    Visualizer for Jos Stam's Stable Fluid simulation output data
    
    Reads velocity, pressure, and dye data from txt files and creates
    various visualizations including animations and static plots.
    """
    
    def __init__(self, data_dir="fluid_output"):
        self.data_dir = data_dir
        self.velocity_data = None
        self.pressure_data = None
        self.dye_data = None
        self.frames = []
        self.grid_width = 0
        self.grid_height = 0
        self.times = []
        
    def load_data(self):
        """Load all data files"""
        print("Loading simulation data...")
        
        # Load velocity data
        velocity_file = os.path.join(self.data_dir, "velocity_data.txt")
        if os.path.exists(velocity_file):
            print(f"Loading velocity data from {velocity_file}")
            self.velocity_data = self._load_velocity_data(velocity_file)
        else:
            print(f"Warning: {velocity_file} not found")
            
        # Load pressure data  
        pressure_file = os.path.join(self.data_dir, "pressure_data.txt")
        if os.path.exists(pressure_file):
            print(f"Loading pressure data from {pressure_file}")
            self.pressure_data = self._load_scalar_data(pressure_file)
        else:
            print(f"Warning: {pressure_file} not found")
            
        # Load dye data
        dye_file = os.path.join(self.data_dir, "dye_data.txt")
        if os.path.exists(dye_file):
            print(f"Loading dye data from {dye_file}")
            self.dye_data = self._load_scalar_data(dye_file)
        else:
            print(f"Warning: {dye_file} not found")
            
        print(f"Data loading complete. Grid size: {self.grid_width}x{self.grid_height}, Frames: {len(self.frames)}")
    
    def _load_velocity_data(self, filename):
        """Load velocity data from file
        Format: frame time x y vel_x vel_y
        """
        data = []
        frames_set = set()
        
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('#') or line.strip() == '':
                    continue
                    
                parts = line.strip().split()
                if len(parts) >= 6:
                    frame = int(parts[0])
                    time = float(parts[1])
                    x = int(parts[2])
                    y = int(parts[3])
                    vel_x = float(parts[4])
                    vel_y = float(parts[5])
                    
                    data.append((frame, time, x, y, vel_x, vel_y))
                    frames_set.add(frame)
                    
                    # Update grid dimensions
                    self.grid_width = max(self.grid_width, x + 1)
                    self.grid_height = max(self.grid_height, y + 1)
        
        self.frames = sorted(list(frames_set))
        return data
    
    def _load_scalar_data(self, filename):
        """Load scalar data (pressure or dye) from file
        Format: frame time x y value
        """
        data = []
        
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('#') or line.strip() == '':
                    continue
                    
                parts = line.strip().split()
                if len(parts) >= 5:
                    frame = int(parts[0])
                    time = float(parts[1])
                    x = int(parts[2])
                    y = int(parts[3])
                    value = float(parts[4])
                    
                    data.append((frame, time, x, y, value))
        
        return data
    
    def get_frame_data(self, frame_num, data_type='dye'):
        """Extract data for a specific frame and reshape to grid"""
        if data_type == 'velocity' and self.velocity_data:
            vel_x = np.zeros((self.grid_height, self.grid_width))
            vel_y = np.zeros((self.grid_height, self.grid_width))
            
            for frame, time, x, y, vx, vy in self.velocity_data:
                if frame == frame_num:
                    vel_x[y, x] = vx
                    vel_y[y, x] = vy
                    
            return vel_x, vel_y
            
        elif data_type == 'pressure' and self.pressure_data:
            field = np.zeros((self.grid_height, self.grid_width))
            
            for frame, time, x, y, value in self.pressure_data:
                if frame == frame_num:
                    field[y, x] = value
                    
            return field
            
        elif data_type == 'dye' and self.dye_data:
            field = np.zeros((self.grid_height, self.grid_width))
            
            for frame, time, x, y, value in self.dye_data:
                if frame == frame_num:
                    field[y, x] = value
                    
            return field
        
        return None
    
    def plot_frame(self, frame_num, save_path=None):
        """Create a multi-panel plot for a specific frame"""
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle(f'Fluid Simulation - Frame {frame_num}', fontsize=16)
        
        # Dye concentration
        if self.dye_data:
            dye_field = self.get_frame_data(frame_num, 'dye')
            if dye_field is not None:
                im1 = axes[0, 0].imshow(dye_field, cmap='hot', origin='lower')
                axes[0, 0].set_title('Dye Concentration')
                axes[0, 0].set_xlabel('X')
                axes[0, 0].set_ylabel('Y')
                plt.colorbar(im1, ax=axes[0, 0])
        
        # Pressure field
        if self.pressure_data:
            pressure_field = self.get_frame_data(frame_num, 'pressure')
            if pressure_field is not None:
                im2 = axes[0, 1].imshow(pressure_field, cmap='coolwarm', origin='lower')
                axes[0, 1].set_title('Pressure Field')
                axes[0, 1].set_xlabel('X')
                axes[0, 1].set_ylabel('Y')
                plt.colorbar(im2, ax=axes[0, 1])
        
        # Velocity magnitude
        if self.velocity_data:
            vel_x, vel_y = self.get_frame_data(frame_num, 'velocity')
            if vel_x is not None and vel_y is not None:
                vel_mag = np.sqrt(vel_x**2 + vel_y**2)
                im3 = axes[1, 0].imshow(vel_mag, cmap='viridis', origin='lower')
                axes[1, 0].set_title('Velocity Magnitude')
                axes[1, 0].set_xlabel('X')
                axes[1, 0].set_ylabel('Y')
                plt.colorbar(im3, ax=axes[1, 0])
                
                # Velocity field (streamlines)
                step = max(1, self.grid_width // 20)  # Subsample for clarity
                x_grid = np.arange(0, self.grid_width, step)
                y_grid = np.arange(0, self.grid_height, step)
                X, Y = np.meshgrid(x_grid, y_grid)
                
                U = vel_x[::step, ::step]
                V = vel_y[::step, ::step]
                
                axes[1, 1].streamplot(X, Y, U.T, V.T, density=1.5, color='blue', linewidth=1)
                axes[1, 1].set_title('Velocity Field (Streamlines)')
                axes[1, 1].set_xlabel('X')
                axes[1, 1].set_ylabel('Y')
                axes[1, 1].set_xlim(0, self.grid_width)
                axes[1, 1].set_ylim(0, self.grid_height)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Frame {frame_num} saved to {save_path}")
        else:
            plt.show()
        
        plt.close()
    
    def create_animation(self, field_type='dye', save_path=None, fps=15, max_frames=None):
        """Create an animation of the specified field"""
        if not self.frames:
            print("No frame data available for animation")
            return
            
        frames_to_animate = self.frames[:max_frames] if max_frames else self.frames
        
        fig, ax = plt.subplots(figsize=(8, 8))
        
        # Get data for first frame to set up plot
        if field_type == 'dye':
            first_data = self.get_frame_data(frames_to_animate[0], 'dye')
            title = 'Dye Concentration Animation'
            cmap = 'hot'
        elif field_type == 'pressure':
            first_data = self.get_frame_data(frames_to_animate[0], 'pressure')
            title = 'Pressure Field Animation'
            cmap = 'coolwarm'
        elif field_type == 'velocity':
            vel_x, vel_y = self.get_frame_data(frames_to_animate[0], 'velocity')
            first_data = np.sqrt(vel_x**2 + vel_y**2) if vel_x is not None else None
            title = 'Velocity Magnitude Animation'
            cmap = 'viridis'
        else:
            print(f"Unknown field type: {field_type}")
            return
            
        if first_data is None:
            print(f"No data available for field type: {field_type}")
            return
        
        # Set up the plot
        im = ax.imshow(first_data, cmap=cmap, origin='lower', animated=True)
        ax.set_title(title)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        
        # Frame text
        frame_text = ax.text(0.02, 0.98, '', transform=ax.transAxes, 
                            fontsize=12, verticalalignment='top',
                            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        def animate(frame_idx):
            frame_num = frames_to_animate[frame_idx]
            
            if field_type == 'dye':
                data = self.get_frame_data(frame_num, 'dye')
            elif field_type == 'pressure':
                data = self.get_frame_data(frame_num, 'pressure')
            elif field_type == 'velocity':
                vel_x, vel_y = self.get_frame_data(frame_num, 'velocity')
                data = np.sqrt(vel_x**2 + vel_y**2) if vel_x is not None else None
            
            if data is not None:
                im.set_data(data)
                im.set_clim(vmin=data.min(), vmax=data.max())
                frame_text.set_text(f'Frame: {frame_num}')
            
            return [im, frame_text]
        
        # Create animation
        anim = animation.FuncAnimation(fig, animate, frames=len(frames_to_animate),
                                      interval=1000/fps, blit=True, repeat=True)
        
        if save_path:
            print(f"Saving animation to {save_path}...")
            anim.save(save_path, writer='pillow', fps=fps)
            print(f"Animation saved to {save_path}")
        else:
            plt.show()
        
        plt.close()
        return anim
    
    def plot_statistics(self):
        """Plot various statistics of the simulation"""
        if not self.frames:
            print("No data available for statistics")
            return
            
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle('Simulation Statistics', fontsize=16)
        
        # Calculate statistics for each frame
        total_dye = []
        max_velocity = []
        avg_pressure = []
        frame_numbers = []
        
        for frame_num in self.frames[::5]:  # Sample every 5th frame for performance
            frame_numbers.append(frame_num)
            
            # Dye statistics
            if self.dye_data:
                dye_field = self.get_frame_data(frame_num, 'dye')
                if dye_field is not None:
                    total_dye.append(np.sum(dye_field))
                else:
                    total_dye.append(0)
            
            # Velocity statistics
            if self.velocity_data:
                vel_x, vel_y = self.get_frame_data(frame_num, 'velocity')
                if vel_x is not None and vel_y is not None:
                    vel_mag = np.sqrt(vel_x**2 + vel_y**2)
                    max_velocity.append(np.max(vel_mag))
                else:
                    max_velocity.append(0)
            
            # Pressure statistics
            if self.pressure_data:
                pressure_field = self.get_frame_data(frame_num, 'pressure')
                if pressure_field is not None:
                    avg_pressure.append(np.mean(pressure_field))
                else:
                    avg_pressure.append(0)
        
        # Plot statistics
        if total_dye:
            axes[0, 0].plot(frame_numbers, total_dye, 'r-', linewidth=2)
            axes[0, 0].set_title('Total Dye Amount')
            axes[0, 0].set_xlabel('Frame')
            axes[0, 0].set_ylabel('Total Dye')
            axes[0, 0].grid(True)
        
        if max_velocity:
            axes[0, 1].plot(frame_numbers, max_velocity, 'b-', linewidth=2)
            axes[0, 1].set_title('Maximum Velocity')
            axes[0, 1].set_xlabel('Frame')
            axes[0, 1].set_ylabel('Max Velocity')
            axes[0, 1].grid(True)
        
        if avg_pressure:
            axes[1, 0].plot(frame_numbers, avg_pressure, 'g-', linewidth=2)
            axes[1, 0].set_title('Average Pressure')
            axes[1, 0].set_xlabel('Frame')
            axes[1, 0].set_ylabel('Avg Pressure')
            axes[1, 0].grid(True)
        
        # Energy plot (kinetic energy proxy)
        if max_velocity and total_dye:
            energy_proxy = [v * d for v, d in zip(max_velocity, total_dye)]
            axes[1, 1].plot(frame_numbers, energy_proxy, 'm-', linewidth=2)
            axes[1, 1].set_title('Energy Proxy (Vel × Dye)')
            axes[1, 1].set_xlabel('Frame')
            axes[1, 1].set_ylabel('Energy Proxy')
            axes[1, 1].grid(True)
        
        plt.tight_layout()
        plt.show()


def main():
    parser = argparse.ArgumentParser(description='Visualize Stam Stable Fluid simulation data')
    parser.add_argument('--data-dir', default='fluid_output', 
                       help='Directory containing simulation data files')
    parser.add_argument('--mode', choices=['frame', 'animation', 'stats'], default='frame',
                       help='Visualization mode')
    parser.add_argument('--frame', type=int, default=0,
                       help='Frame number to visualize (for frame mode)')
    parser.add_argument('--field', choices=['dye', 'pressure', 'velocity'], default='dye',
                       help='Field type for animation')
    parser.add_argument('--save', help='Save output to file')
    parser.add_argument('--fps', type=int, default=15,
                       help='Frames per second for animation')
    parser.add_argument('--max-frames', type=int,
                       help='Maximum number of frames to animate')
    
    args = parser.parse_args()
    
    # Create visualizer and load data
    visualizer = FluidVisualizer(args.data_dir)
    visualizer.load_data()
    
    if not visualizer.frames:
        print("No simulation data found. Make sure to run the simulation first.")
        return
    
    # Execute based on mode
    if args.mode == 'frame':
        frame_num = min(args.frame, max(visualizer.frames))
        print(f"Plotting frame {frame_num}")
        visualizer.plot_frame(frame_num, args.save)
        
    elif args.mode == 'animation':
        print(f"Creating {args.field} animation...")
        visualizer.create_animation(args.field, args.save, args.fps, args.max_frames)
        
    elif args.mode == 'stats':
        print("Plotting simulation statistics...")
        visualizer.plot_statistics()


if __name__ == "__main__":
    main()