#include <fstream>
#include <iostream>
#include <numeric> // For std::iota

#include "mpi.h"

#include <Kokkos_Core.hpp>

#include <CabanaPD.hpp>
#include <CabanaPD_Output.hpp> 

// Simulate cold spray particle impact on a substrate.
void coldSprayImpactExample( const std::string filename )
{
    // ====================================================
    //               Choose Kokkos spaces
    // ====================================================
    using exec_space = Kokkos::DefaultExecutionSpace;
    using memory_space = typename exec_space::memory_space;

    // ====================================================
    //                   Read inputs
    // ====================================================
    CabanaPD::Inputs inputs( filename );

    // ====================================================
    //                Material parameters (Copper)
    // ====================================================
    // 假设铜的参数，这些值需要根据实际情况调整和查阅文献
    double rho0 = 8960.0; // Copper density (kg/m^3)
    double E = 110.0e9;   // Copper Young's Modulus (Pa)
    double nu = 0.34;     // Copper Poisson's ratio
    double K = E / ( 3 * ( 1 - 2 * nu ) ); // Bulk modulus
    double G0 = 1000.0;   // Fracture energy (J/m^2) - Peridynamics parameter, needs calibration
    double sigma_y = 70.0e6; // Copper Yield Strength (Pa) - For elastic-perfectly-plastic model
    double delta = inputs["horizon"]; // Horizon, must be larger than dx
    delta += 1e-10;

    // ====================================================
    //                  Discretization
    // ====================================================
    std::array<double, 3> low_corner = inputs["low_corner"];   // Domain min coordinates
    std::array<double, 3> high_corner = inputs["high_corner"]; // Domain max coordinates
    std::array<int, 3> num_cells = inputs["num_cells"];       // Number of cells in each dimension

    // Ensure horizon is correctly related to grid spacing for halo width calculation
    double dx = (high_corner[0] - low_corner[0]) / num_cells[0];
    int m = std::ceil( delta / dx ); // number of cells to cover horizon
    int halo_width = m + 1;          // Halo width for ghost particles

    // ====================================================
    //                  Force model type
    // ====================================================
    // PMB (Peridynamic Micropolar Bond) is a bond-based model with shear.
    // For more advanced plasticity or State-based PD, you'd typically need
    // to implement a custom force model and mechanics type.
    using model_type = CabanaPD::PMB;
    // We will use ElasticPerfectlyPlastic provided by CabanaPD.
    // For true crystal plasticity, this would need to be replaced with
    // a custom implementation within CabanaPD's framework.
    using mechanics_type = CabanaPD::ElasticPerfectlyPlastic;

    // ====================================================
    //    Custom particle generation and initialization
    // ====================================================
    // Define geometry for the particle and the substrate.
    // Assuming a 2D or 3D simulation domain for simplicity
    // Example: Substrate at the bottom, particle impacting from top

    // Substrate dimensions
    // Let's assume substrate occupies bottom half of the domain.
    // For example, from low_corner[1] to mid-point in Y.
    double substrate_height = (high_corner[1] - low_corner[1]) / 2.0;
    std::array<double, 3> substrate_low_corner = {low_corner[0], low_corner[1], low_corner[2]};
    std::array<double, 3> substrate_high_corner = {high_corner[0], low_corner[1] + substrate_height, high_corner[2]};

    // Particle dimensions
    // Let's assume a single spherical particle for simplicity,
    // positioned above the substrate.
    double particle_radius = inputs["particle_radius"]; // Read from input file
    double particle_initial_y_pos = substrate_high_corner[1] + particle_radius + dx; // Just above the substrate
    std::array<double, 3> particle_center = {
        0.5 * (low_corner[0] + high_corner[0]), // Center in X
        particle_initial_y_pos,                 // Above substrate in Y
        0.5 * (low_corner[2] + high_corner[2])  // Center in Z
    };

    // Initial impact velocity
    double impact_velocity = 500.0; // 500 m/s in negative Y direction

    // Lambda function to initialize particles based on geometry
    auto init_op = KOKKOS_LAMBDA( const int, const double x[3] )
    {
        // Check if particle belongs to substrate
        if ( x[0] >= substrate_low_corner[0] && x[0] <= substrate_high_corner[0] &&
             x[1] >= substrate_low_corner[1] && x[1] <= substrate_high_corner[1] &&
             x[2] >= substrate_low_corner[2] && x[2] <= substrate_high_corner[2] )
        {
            return true; // Is part of substrate
        }
        // Check if particle belongs to the impacting particle
        else
        {
            double dist_sq = ( x[0] - particle_center[0] ) * ( x[0] - particle_center[0] ) +
                             ( x[1] - particle_center[1] ) * ( x[1] - particle_center[1] ) +
                             ( x[2] - particle_center[2] ) * ( x[2] - particle_center[2] );
            if ( dist_sq <= particle_radius * particle_radius )
            {
                return true; // Is part of impacting particle
            }
        }
        return false; // Neither substrate nor particle
    };

    // Create particles based on the defined geometry
    CabanaPD::Particles particles(
        memory_space{}, model_type{}, low_corner, high_corner, num_cells,
        halo_width, Cabana::InitRandom{}, init_op, exec_space{} );

    // Get slices for density and velocity
    auto rho = particles.sliceDensity();
    auto v = particles.sliceVelocity();

    // Initialize density and initial velocity for impact particle
    auto init_functor = KOKKOS_LAMBDA( const int pid )
    {
        // Set density for all created particles (both substrate and impactor)
        rho( pid ) = rho0;

        // Check if particle belongs to the impacting particle (based on initial position)
        double x_pid[3];
        particles.sliceReferencePosition()(pid, 0) = x_pid[0];
        particles.sliceReferencePosition()(pid, 1) = x_pid[1];
        particles.sliceReferencePosition()(pid, 2) = x_pid[2];

        double dist_sq = ( x_pid[0] - particle_center[0] ) * ( x_pid[0] - particle_center[0] ) +
                         ( x_pid[1] - particle_center[1] ) * ( x_pid[1] - particle_center[1] ) +
                         ( x_pid[2] - particle_center[2] ) * ( x_pid[2] - particle_center[2] );

        if ( dist_sq <= particle_radius * particle_radius )
        {
            // Set initial downward velocity for the impacting particle
            v( pid, 0 ) = 0.0;
            v( pid, 1 ) = -impact_velocity; // Negative Y direction for downward impact
            v( pid, 2 ) = 0.0;
        }
        else
        {
            // Substrate particles start with zero velocity
            v( pid, 0 ) = 0.0;
            v( pid, 1 ) = 0.0;
            v( pid, 2 ) = 0.0;
        }
    };
    particles.updateParticles( exec_space{}, init_functor );

    // ====================================================
    //                    Force model
    // ====================================================
    CabanaPD::ForceModel force_model( model_type{}, mechanics_type{},
                                      memory_space{}, delta, K, G0, sigma_y );

    // ====================================================
    //                   Create solver
    // ====================================================
    CabanaPD::Solver solver( inputs, particles, force_model );

    // ====================================================
    //                  Boundary conditions
    // ====================================================
    // Define a region for the fixed bottom of the substrate
    // Assuming a small fixed layer at the very bottom of the substrate
    double fixed_layer_thickness = 2.0 * dx; // Fix a few layers of particles
    CabanaPD::Region<CabanaPD::RectangularPrism> fixed_bottom_region(
        low_corner[0], high_corner[0],
        low_corner[1], low_corner[1] + fixed_layer_thickness,
        low_corner[2], high_corner[2]
    );

    // Create BC last to ensure ghost particles are included.
    auto x_ref = solver.particles.sliceReferencePosition();
    auto u = solver.particles.sliceDisplacement();
    auto disp_func = KOKKOS_LAMBDA( const int pid, const double t )
    {
        // Fix the displacement of particles in the bottom fixed region
        if ( fixed_bottom_region.inside( x_ref, pid ) )
        {
            u( pid, 0 ) = 0.0;
            u( pid, 1 ) = 0.0;
            u( pid, 2 ) = 0.0;
        }
    };
    // Note: The 'false' argument indicates that this BC is not for a moving boundary.
    // We are setting initial velocity for the particle, and fixing the substrate boundary.
    auto bc = CabanaPD::createBoundaryCondition( disp_func, exec_space{},
                                                 solver.particles, false,
                                                 fixed_bottom_region ); // Only one region for fixing

    // ====================================================
    //                      Outputs
    // ====================================================
    // You'll want to output more detailed particle data for deformation analysis
    // For example, individual particle positions, velocities, damage, etc.
    // The following are examples based on the dogbone test.
    // You might need to add specific output for particle properties (e.g., damage, plastic strain if implemented)

    // Example: Output average position of the impactor particle over time (for center of mass tracking)
    // First, identify the IDs of the impactor particles
    Kokkos::View<int*, memory_space> impactor_pids("impactor_pids", particles.size());
    Kokkos::View<int, memory_space> num_impactor_pids("num_impactor_pids");
    auto x_ref_local = particles.sliceReferencePosition();

    Kokkos::parallel_reduce( "count_impactor_pids", particles.size(), 
        KOKKOS_LAMBDA( const int pid, int& count ) {
            double x_pid[3];
            x_ref_local(pid, 0) = x_pid[0];
            x_ref_local(pid, 1) = x_pid[1];
            x_ref_local(pid, 2) = x_pid[2];

            double dist_sq = ( x_pid[0] - particle_center[0] ) * ( x_pid[0] - particle_center[0] ) +
                             ( x_pid[1] - particle_center[1] ) * ( x_pid[1] - particle_center[1] ) +
                             ( x_pid[2] - particle_center[2] ) * ( x_pid[2] - particle_center[2] );

            if ( dist_sq <= particle_radius * particle_radius ) {
                impactor_pids(count++) = pid;
            }
        }, num_impactor_pids);

    // This is a more complex output as it requires averaging over a dynamic set of particles.
    // For simplicity, let's just output the current position of the first particle in the simulation.
    // For real analysis, you'd extract data from VTK files for visualization.
    auto y_current = solver.particles.sliceCurrentPosition();
    auto first_particle_pos_y = KOKKOS_LAMBDA( const int p ) {
        // Output y-position of the first particle in the global list (not necessarily the impactor)
        // For specific particle tracking, you'd need to identify particle IDs.
        if (p == 0) return y_current(p, 1);
        return 0.0; // Return dummy for other particles
    };
    auto output_first_particle_y = CabanaPD::createOutputTimeSeries(
        "output_first_particle_pos_y.txt", inputs, exec_space{}, solver.particles,
        first_particle_pos_y, CabanaPD::Region<CabanaPD::RectangularPrism>(low_corner[0], high_corner[0], low_corner[1], high_corner[1], low_corner[2], high_corner[2]) ); // Whole domain to check all particles

    // Output all particle data to VTK files for visualization
    // This is crucial for seeing the deformation.
//    !CabanaPD::VTKOutput vtk_output( "cold_spray_impact", inputs, solver.particles );

    // ====================================================
    //                   Simulation run
    // ====================================================
    solver.init( bc ); // Apply initial boundary conditions (fixed substrate)
    // solver.run( bc, vtk_output, output_first_particle_y ); // Pass VTK output to solver.run
    solver.run( bc, output_first_particle_y ); // 尝试这个版本
}

// Initialize MPI+Kokkos.
int main( int argc, char* argv[] )
{
    MPI_Init( &argc, &argv );
    Kokkos::initialize( argc, argv );

    if ( argc < 2 ) {
        std::cerr << "Usage: " << argv[0] << " <input_file.json>" << std::endl;
        return 1;
    }

    coldSprayImpactExample( argv[1] );

    Kokkos::finalize();
    MPI_Finalize();
}