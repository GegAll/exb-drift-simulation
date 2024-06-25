# ExB Drift Simulation

![ExBDrift_E_straight_pot_2](https://github.com/GegAll/exb-drift-simulation/assets/170819708/748a55d4-1d52-49fc-8bd2-63769d0a12ec)

This repository contains a Poisson solver, and more importantly, an electric field and ExB drift simulator for a grid emulating a plasma reactor chamber with a mask on the bottom side where a bias voltage is set. A magnetic field perpendicular to the plane is applied, showing (for a positive value) into the grid’s plane as illustrated below.

![Maske_straight_ExBdrift](https://github.com/GegAll/exb-drift-simulation/assets/170819708/5998f53a-2379-43ea-a33f-2881fa7e27fc)

The interaction between the electric field created by the potential difference and the magnetic field results in an ExB drift of the particles on the surface of the mask. The mask can take three different configurations: inwards, outwards, and straight. A pictorial representation of each mask is shown below:

![Maske_straight](https://github.com/GegAll/exb-drift-simulation/assets/170819708/79795bc4-0f4e-48c4-820d-845dadb1455c)
Straight

![Maske_inwards_better](https://github.com/GegAll/exb-drift-simulation/assets/170819708/0fcf6f5f-27ad-44e8-be59-4614a996692c)
Inwards

![Maske_outwards_better](https://github.com/GegAll/exb-drift-simulation/assets/170819708/55a8cee2-ef33-42c5-b0dc-37c6ebd655f6)
Outwards

## Features

Each configuration can be selected directly after compiling the program by choosing it in the terminal. Other adjustable parameters include the bias voltage at the bottom, the length and number of grid points for each coordinate, and the magnitude and direction of the magnetic field (negative value to reverse direction). Note that increasing the number of grid points to enhance resolution also increases compilation time, which can significantly impact running time.

## Output

After compiling and running the program, it displays the number of iterations the Poisson solver needed before convergence. It saves the potential, electric field, and ExB drift data for each point into `.txt` files. The paths for these files can be changed in the program's last lines to suit your directory structure. The following files are saved after the simulation:
- `potential_data.txt`: Contains the potential of the grid.
- `Ex_data.txt` and `Ey_data.txt`: Contain the x and y components of the electric field.
- `Vd_x_data.txt` and `Vd_y_data.txt`: Contain the x and y components of the ExB drift.

## Visualization

The repository includes a visualizer for all three magnitudes. It plots the potential data as a filled contour plot, and both the electric field and ExB drift as 2D arrow fields, with a filled contour in the background showing the magnitude and direction of these vectors at each grid point. The `plt.quiver()` function is used for the arrow plots, and its features such as `headwidth=`, `headlength=`, and `scale=` can be adjusted for better visualization. The density of the arrows is controlled by the variables `step` and `step_drift`, where higher values result in fewer arrows, improving clarity.

### Examples

Some examples of the potential, electric field, and ExB drift for different mask structures are shown below:

#### Inwards
![potential_inwards_pot_2](https://github.com/GegAll/exb-drift-simulation/assets/170819708/f7d0ea1d-3de7-4795-a298-f16cc4403979)
![ExBDrift_E_inwards_pot_2](https://github.com/GegAll/exb-drift-simulation/assets/170819708/82da2989-b5a9-4e62-b453-1a771f124d03)

#### Outwards
![potential_outwards_pot_2](https://github.com/GegAll/exb-drift-simulation/assets/170819708/ab05e20c-a72a-46cd-8469-9836481e3ca0)
![ExBDrift_E_outwards_pot_2](https://github.com/GegAll/exb-drift-simulation/assets/170819708/3e540ee7-7ed9-47d6-a0f6-c8a75e4a5620)

## How to Use (Simulation)

1. **Compile the Simulator (C++)**:
   ```sh
   g++ -o exb-drift-simulation potential_calculator.cpp
   ./exb_simulator
   ```
2. **Adjust Parameters**:
   ```sh
   Choose the mask configuration and other parameters in the terminal after running the simulator.
   ```

## How to Use (Visualiizer)

1. **Clone the repository**:
   ```sh
   git clone https://github.com/GegAll/exb-drift-simulation.git
   ```

2. **Navigate to the project directory**:
   ```sh
   cd exb_simulator
   ```
   
3. **Install the required dependencies**:
   ```sh
   pip install -r requirements.txt
   ```
   
4. **Run the Jupyter notebook to see the visualizer**:
   ```sh
   jupyter notebook academic_success.ipynb
   ```

5. **Change the filepaths**:
   ```sh
   If necessary change the filepahts to match where you saved the simulator and the visualizer in your PC
   ```
   
## License
This project is licensed under the MIT License. See the LICENSE file for details.
