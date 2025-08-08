# Stam-Stable-Fluids

An implementation of stable fluids simulation (refer to Stam's paper in 1999)

## Command Line Builds
```bash
mkdir build
cd build
cmake ..
cmake --build . --config Release
cd ..
.\main.exe
```

## Visualize by python scripts
Ensure that the following Python packages are installed:

```bash
pip install numpy matplotlib
```

### Basic Usage

```bash
# Visualize a single frame (default frame 0)
python visualize.py --mode frame

# Visualize a specific frame
python visualize.py --mode frame --frame 100

# Create dye concentration animation
python visualize.py --mode animation --field dye

# Create velocity field animation
python visualize.py --mode animation --field velocity

# Create pressure field animation
python visualize.py --mode animation --field pressure

# Display simulation statistics
python visualize.py --mode stats
```

### Advanced Options

```bash
# Specify data directory
python visualize.py --data-dir "custom_output" --mode frame

# Save image to file
python visualize.py --mode frame --frame 50 --save "frame_50.png"

# Save animation to GIF file
python visualize.py --mode animation --field dye --save "dye_animation.gif" --fps 20

# Limit animation frames (improve performance)
python visualize.py --mode animation --field velocity --max-frames 100
```