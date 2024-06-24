#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
using namespace std;

// Constants
const int Ly = 6; // y length
const int Ny = 120; // y grid points
const double Lx = 4;  // x length
const int Nx = 80;  // x grid points

const double dx = Lx / (Nx - 1);  // Grid spacing
const double dy = dx;  // Grid spacing (assuming square grid cells)

const double omega = 1.9;  // SOR relaxation parameter
const int max_iter = 1000000;  // Maximum number of iterations
const double tolerance = 1e-6;  // Convergence criterion
const double B = 1; // Constant magnetic field magnitude
const double gr_pot = -2.0; // ground potential

// Initialize potential
std::vector<std::vector<double>> phi(Nx, std::vector<double>(Ny, 0.0));

// Set boundary conditions
void setBoundaryConditions(std::vector<std::vector<double>>& phi) {
    for (int i = 0; i < Nx; ++i) {
        phi[i][0] = gr_pot;  // Grounded electrode at the bottom
        phi[i][Ny - 1] = 0;  // Higher potential electrode at the top
    }
}
// Save the potential to a file
void savePotential(const std::vector<std::vector<double>>& phi, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        for (const auto& row : phi) {
            for (const auto& val : row) {
                file << val << " ";
            }
            file << "\n";
        }
        file.close();
    } else {
        std::cerr << "Unable to open file for writing: " << filename << std::endl;
    }
}

std::vector<std::vector<double>> createPotential (int mask) {
    // Main loop for SOR method with periodic boundary conditions
    for (int iter = 0; iter < max_iter; ++iter) {

        double width = Lx; // grid's width [mm]
        double x_r = 3*width/8; // mask length [mm]
        double right_boundary = x_r/dx; // mask length [grid points]
        double Gap = 0.5; // gap between the mask and the substrate [mm]
        double Mask_heigth = 1; // height of the mask [mm]
        double upper_boundary = Mask_heigth / dy; // height of the mask [grid points]
        double lower_boundary = Gap / dx; // gap between the mask and the substrate [mm]
        double left_boundary = x_r/dx - (upper_boundary - lower_boundary); // fase length [grid points]
        double m_pot = gr_pot; // Mask potential
        std::vector<std::vector<double>> phi_old = phi;

        switch (mask) // calculate the potential for each geometry
        {

        case 0: // Outwards
            for (int i = 1; i < Nx - 1; ++i) {
                for (int j = 1; j < Ny - 1; ++j) {
                    if (j > lower_boundary && j < lower_boundary + upper_boundary) {
                        if (i > left_boundary && i < right_boundary) {
                            if(-2*i >= j - lower_boundary - 2*right_boundary) {
                                phi[i][j] = m_pot;
                            } else { // solving poisson equation
                                phi[i][j] = (1 - omega) * phi[i][j] + omega * 0.25 * (phi[i + 1][j] + phi[i - 1][j] + phi[i][j + 1] + phi[i][j - 1]);
                            }
                        } else if (i < left_boundary) {
                            phi[i][j] = m_pot;
                        } else if (i > Nx - right_boundary && i < Nx - left_boundary) {
                            if (2*i >= j - lower_boundary + 2*Nx - 2*right_boundary){
                                phi[i][j] = m_pot;
                            } else {
                                phi[i][j] = (1 - omega) * phi[i][j] + omega * 0.25 * (phi[i + 1][j] + phi[i - 1][j] + phi[i][j + 1] + phi[i][j - 1]);
                            }
                        } else if (i > Nx - left_boundary) {
                            phi[i][j] = m_pot;
                        } else {
                            phi[i][j] = (1 - omega) * phi[i][j] + omega * 0.25 * (phi[i + 1][j] + phi[i - 1][j] + phi[i][j + 1] + phi[i][j - 1]);
                        }
                    } else {
                        phi[i][j] = (1 - omega) * phi[i][j] + omega * 0.25 * (phi[i + 1][j] + phi[i - 1][j] + phi[i][j + 1] + phi[i][j - 1]);
                    }     
                }
            }
            break;

        case 1: // Inwards 
            for (int i = 1; i < Nx - 1; ++i) {
                for (int j = 1; j < Ny - 1; ++j) {
                    if (j > lower_boundary && j < lower_boundary + upper_boundary) {
                        if (i > left_boundary && i < right_boundary) {
                            if(2*i <= j + 2*left_boundary - lower_boundary) {
                                phi[i][j] = m_pot;
                            } else { 
                                phi[i][j] = (1 - omega) * phi[i][j] + omega * 0.25 * (phi[i + 1][j] + phi[i - 1][j] + phi[i][j + 1] + phi[i][j - 1]);
                            }
                        } else if (i < left_boundary) {
                            phi[i][j] = m_pot;
                        } else if (i > Nx - right_boundary && i < Nx - left_boundary) {
                            if (2*i >= -j + lower_boundary + 2*Nx - 2*left_boundary){
                                phi[i][j] = m_pot;
                            } else {
                                phi[i][j] = (1 - omega) * phi[i][j] + omega * 0.25 * (phi[i + 1][j] + phi[i - 1][j] + phi[i][j + 1] + phi[i][j - 1]);
                            }
                        } else if (i > Nx - left_boundary) {
                            phi[i][j] = m_pot;
                        } else {
                            phi[i][j] = (1 - omega) * phi[i][j] + omega * 0.25 * (phi[i + 1][j] + phi[i - 1][j] + phi[i][j + 1] + phi[i][j - 1]);
                        }
                    } else {
                        phi[i][j] = (1 - omega) * phi[i][j] + omega * 0.25 * (phi[i + 1][j] + phi[i - 1][j] + phi[i][j + 1] + phi[i][j - 1]);
                    }     
                }
            }
            break;
            
        case 2: // Straight
            for (int i = 1; i < Nx - 1; ++i) {
                for (int j = 1; j < Ny - 1; ++j) {
                    if (j < upper_boundary + lower_boundary && j > lower_boundary && i < right_boundary) {
                        phi[i][j] = m_pot; // left mask potential
                    } else if (j < upper_boundary + lower_boundary && j > lower_boundary && i > Nx - right_boundary) {
                        phi[i][j] = m_pot; // right mask potential
                    } else { // calcualte the potential by solving the poisson equation
                        phi[i][j] = (1 - omega) * phi[i][j] + omega * 0.25 * (
                            phi[i + 1][j] + phi[i - 1][j] + phi[i][j + 1] + phi[i][j - 1]);
                    }
                }
            }
            break;
        
        default:
            cerr << "wrong mask" << endl;
            break;
        }

        // Treat boundary conditions
        for (int j = 0; j < Ny; ++j) {
            phi[0][j] = phi[Nx - 1][j];  // Left boundary
            phi[Nx - 1][j] = phi[1][j];  // Right boundary
        }

        // Check for convergence
        bool converged = true;
        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {
                if (abs(phi[i][j] - phi_old[i][j]) > tolerance) {
                    converged = false;
                    break;
                }
            }
            if (!converged) break;
        }

        if (converged) {
            cout << "Convergence achieved after " << iter + 1 << " iterations." << endl; // print the number of iterations needed
            break;
        }
    }

    return phi;
}

// Calculate and save the electric field components
void calculateAndSaveElectricField(const std::vector<std::vector<double>>& phi, double dx, double dy, const std::string& Ex_filename, const std::string& Ey_filename) {
    std::ofstream Ex_file(Ex_filename);
    std::ofstream Ey_file(Ey_filename);
    if (Ex_file.is_open() && Ey_file.is_open()) {
        // two grid points less for each direction because the electric field on the boundaries cannot be calculated
        for (int i = 1; i < Nx - 1; ++i) { 
            for (int j = 1; j < Ny - 1; ++j) {
                // calculate the electric field from the potential by solving the differential equation
                double Ex = -(phi[i+1][j] - phi[i-1][j]) / (2 * dx);
                double Ey = -(phi[i][j+1] - phi[i][j-1]) / (2 * dy);
                Ex_file << Ex << " ";
                Ey_file << Ey << " ";
            }
            Ex_file << "\n";
            Ey_file << "\n";
        }
        // save the files
        Ex_file.close();
        Ey_file.close();
    } else {
        std::cerr << "Unable to open files for writing electric field components." << std::endl;
    }
}

// Calculate and save the E x B drift components
void calculateAndSaveExBDrift(const std::string& Ex_filename, const std::string& Ey_filename, const std::string& Vd_x_filename, const std::string& Vd_y_filename) {
    // Open/Create the electric field and drift files
    std::ifstream Ex_file(Ex_filename);
    std::ifstream Ey_file(Ey_filename);
    std::ofstream Vd_x_file(Vd_x_filename);
    std::ofstream Vd_y_file(Vd_y_filename);

    if (Ex_file.is_open() && Ey_file.is_open() && Vd_x_file.is_open() && Vd_y_file.is_open()) {
        double Ex, Ey;
        while (Ex_file >> Ex && Ey_file >> Ey) {
            // Calculate the ExB drift for this specific structure (B pointing into the z direction)
            double Vd_x = Ey / B; 
            double Vd_y = -Ex / B;
            Vd_x_file << Vd_x << " ";
            Vd_y_file << Vd_y << " ";
        }
        // save the files
        Ex_file.close();
        Ey_file.close();
        Vd_x_file.close();
        Vd_y_file.close();
    } else {
        std::cerr << "Unable to open files for writing ExB drift components." << std::endl;
    }
}


int main() {
    setBoundaryConditions(phi);
    
    // Choose the desired mask configuration
    cout << "Outwards : 0" << endl;
    cout << "Inwards : 1" << endl;
    cout << "Straight : 2" << endl;
    cout << "Enter the type of mask from which the potential will be calculated: ";
    int mask;
    cin >> mask;

    // Create the Potential
    createPotential(mask);

    // Save the final potential to a file
    savePotential(phi, "potential_data\\potential_data.txt");

    // Calculate and save the electric field components
    calculateAndSaveElectricField(phi, dx, dy, "potential_data\\Ex_data.txt", "C:\\Users\\Asus\\Desktop\\python_profile\\potential_data\\Ey_data.txt");
    
    // Calculate and save the E x B drift components
    calculateAndSaveExBDrift("potential_data\\Ex_data.txt", "C:\\Users\\Asus\\Desktop\\python_profile\\potential_data\\Ey_data.txt", "C:\\Users\\Asus\\Desktop\\python_profile\\potential_data\\Vd_x_data.txt", "C:\\Users\\Asus\\Desktop\\python_profile\\potential_data\\Vd_y_data.txt");

    return 0;
}
