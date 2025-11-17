/****************************************************************************
 * Copyright (c) 2022-2023 by Oak Ridge National Laboratory                 *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of CabanaPD. CabanaPD is distributed under a           *
 * BSD 3-clause license. For the licensing terms see the LICENSE file in    *
 * the top-level directory.                                                 *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#include <fstream>
#include <iostream>

#include "mpi.h"

#include <Kokkos_Core.hpp>

#include <CabanaPD.hpp>

// Simulate cold spray.
void coldspray( const std::string filename )
{
      std::cout << "Running bi-material cold spray example with input file: "
                << filename << std::endl;

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
    //                Material parameters
    // ====================================================
// --- Al particle ---
    double rho_Al     = inputs["density"][0];
    double E_Al       = inputs["elastic_modulus"][0];
    double nu_Al      = inputs["Poisson's_ratio"][0];
    double K_Al       = E_Al / (3.0 * (1.0 - 2.0 * nu_Al));
    double G_Al       = E_Al / (2.0 * (1.0 + nu_Al));
    double G0_Al      = inputs["fracture_energy"][0];
    double sigma_y_Al = inputs["yield_stress"][0];

    // --- Cu substrate ---
    double rho_Cu     = inputs["density"][1];
    double E_Cu       = inputs["elastic_modulus"][1];
    double nu_Cu      = inputs["Poisson's_ratio"][1];
    double K_Cu       = E_Cu / (3.0 * (1.0 - 2.0 * nu_Cu));
    double G_Cu       = E_Cu / (2.0 * (1.0 + nu_Cu));
    double K= (K_Cu+K_Al)/2;
    double G0_Cu      = inputs["fracture_energy"][1];
    double sigma_y_Cu = inputs["yield_stress"][1];

    double e = inputs["coefficient_of_restitution"];
    double gamma = inputs["surface_energy"];
    double delta = inputs["horizon"];

    delta += 1e-10;
   
    // ====================================================
    //                  Discretization
    // ====================================================
    std::array<double, 3> low_corner = inputs["low_corner"];
    std::array<double, 3> high_corner = inputs["high_corner"];
    std::array<int, 3> num_cells = inputs["num_cells"];
    int m = std::floor( delta /
                        ( ( high_corner[0] - low_corner[0] ) / num_cells[0] ) );
    int halo_width = m + 1; // Just to be safe.

    // ====================================================
    //                    Force model
    // ====================================================
    using model_type = CabanaPD::PMB;
    using mechanics_type = CabanaPD::ElasticPerfectlyPlastic;

    CabanaPD::ForceModel force_model_Al( model_type{}, mechanics_type{},
                                      memory_space{}, delta,  K_Al, G0_Al, sigma_y_Al);

    CabanaPD::ForceModel force_model_Cu(model_type{}, mechanics_type{},
                                     memory_space{}, delta, K_Cu, G0_Cu, sigma_y_Cu);
    
    // using model_type = CabanaPD::LPS;
    // CabanaPD::ForceModel force_model( model_type{}, delta, K, G, G0 );   

    // using model_type = CabanaPD::PMB;
    // CabanaPD::ForceModel force_model( model_type{}, delta, K, G0 );
    // ====================================================
    //    Custom particle generation and initialization
    // ====================================================
    double x_center = 0.5 * ( low_corner[0] + high_corner[0] );
    double y_center = 0.5 * ( low_corner[1] + high_corner[1] );
    double z_center = 0.5 * ( low_corner[2] + high_corner[2] );



    std::array<double, 3> ball_center = inputs["ball_center"];
    double ball_radius = inputs["ball_radius"];
    double plate_z_min = inputs["plate_z_min"];
    double plate_z_max = inputs["plate_z_max"];
    double vz_ball = inputs["ball_initial_velocity"];

    // Do not create particles outside given cylindrical region
    auto init_op = KOKKOS_LAMBDA(const int, const double x[3])
    {
        // check if the particle is within the ball or plate region
        double rsq = (x[0] - ball_center[0]) * (x[0] - ball_center[0]) +
                     (x[1] - ball_center[1]) * (x[1] - ball_center[1]) +
                     (x[2] - ball_center[2]) * (x[2] - ball_center[2]);

        bool is_ball = rsq <= ball_radius * ball_radius;
        bool is_plate = (x[2] >= plate_z_min) && (x[2] <= plate_z_max);

        return is_ball || is_plate;
    };

    // ====================================================
    //  Simulation run with contact physics
    // ====================================================
    if ( inputs["use_contact"] )
    {
        std::cout << "contact activated"<< std::endl;
        using contact_type = CabanaPD::NormalRepulsionModel;

        CabanaPD::Particles particles(
            memory_space{}, contact_type{}, low_corner, high_corner, num_cells,
            halo_width, Cabana::InitRandom{}, init_op, exec_space{} );

        auto rho = particles.sliceDensity();
        auto x = particles.sliceReferencePosition();
        auto v = particles.sliceVelocity();
        auto f = particles.sliceForce();
        auto dx = particles.dx ;
        auto itype = particles.sliceType(); // particle type: ball or plate    // material type:1=Al_ball, 0=Cu_substrate
        auto nofail= particles.sliceNoFail();   //  make  ball particles not fail

        auto init_functor = KOKKOS_LAMBDA( const int pid )
        {
             
            if ((x(pid, 2) - ball_center[2]) * (x(pid, 2) - ball_center[2]) +
                    (x(pid, 0) - ball_center[0]) * (x(pid, 0) - ball_center[0]) +
                    (x(pid, 1) - ball_center[1]) * (x(pid, 1) - ball_center[1]) <=
                ball_radius * ball_radius)
            {
                v(pid, 2) = -vz_ball; // impact velocity downwards
                itype(pid) = 1; // ball
                rho(pid) = rho_Al;
            }
            else
            {
                v(pid, 2) = 0.0;
                itype(pid) = 0; // plate
                rho(pid) = rho_Cu;
            }
            nofail(pid) = (itype(pid) != 0); // make ball particles not fail
        };
        particles.updateParticles( exec_space{}, init_functor );

        // Use contact radius and extension relative to particle spacing.
        double r_c = inputs["contact_horizon_factor"];
        double r_extend = inputs["contact_horizon_extend_factor"];
        double LJ_r0 = inputs["LJr0"];
        double r0 = LJ_r0 * dx[0]; // lj potential width sigma 1.05dx,
        double beta = inputs["LJbeta"]; 
        double alpha = inputs["LJalpha"]; 
      std::cout << "beta: "
                << beta << std::endl;
      std::cout << "alpha: "
                << alpha << std::endl;
         // NOTE: dx/2 is when particles first touch.
        r_c *= r0;
        r_extend *= dx[0];
 
        double c_czm=inputs["CZM_cohesive_scaling"];  // CZM contact parameters,
        double sy = inputs["CZM_yield_stretch"]; 
        double m_czm = inputs["CZM_degradation_rate"]; 

        // std::cout << "c_czm: "
        //         << c_czm << std::endl;
        // std::cout << "sy: "
        //         << sy << std::endl;
        // std::cout << "m_czm: "
        //         << m_czm << std::endl;
        //NonRepulsiveLJModel
        contact_type contact_model(delta, r_c, r_extend, K, r0, beta, alpha,c_czm,sy,m_czm);

        //HertzianModel contact_model
        //contact_type contact_model( r_c, r_extend, nu, E, e );
        //JKRHertzianModel contact_model
        //contact_type contact_model( r_c, r_extend, nu, E, e, gamma );
    //  Multi-material force model
    // ====================================================
        auto models = CabanaPD::createMultiForceModel(
        particles, CabanaPD::AverageTag{}, force_model_Al, force_model_Cu);

        
        CabanaPD::Solver solver( inputs, particles, models,
                                 contact_model );

        double boundary_layer_thickness = 4*dx[0]; 
        double bottom_z_thickness = 4*dx[0]; 

        double z_bc = low_corner[2]; 
        CabanaPD::Region<CabanaPD::RectangularPrism> bottom_region(
            low_corner[0], high_corner[0],
            low_corner[1], high_corner[1],
            z_bc - bottom_z_thickness, z_bc + bottom_z_thickness );

        CabanaPD::Region<CabanaPD::RectangularPrism> x_min_side_region(
            low_corner[0] - boundary_layer_thickness, low_corner[0] + boundary_layer_thickness,
            low_corner[1], high_corner[1],                           
            low_corner[2], high_corner[2]                            
        );
        CabanaPD::Region<CabanaPD::RectangularPrism> x_max_side_region(
            high_corner[0] - boundary_layer_thickness, high_corner[0] + boundary_layer_thickness,
            low_corner[1], high_corner[1],
            low_corner[2], high_corner[2]
        );
        CabanaPD::Region<CabanaPD::RectangularPrism> y_min_side_region(
            low_corner[0], high_corner[0],
            low_corner[1] - boundary_layer_thickness, low_corner[1] + boundary_layer_thickness,
            low_corner[2], high_corner[2]
        );
        CabanaPD::Region<CabanaPD::RectangularPrism> y_max_side_region(
            low_corner[0], high_corner[0],
            high_corner[1] - boundary_layer_thickness, high_corner[1] + boundary_layer_thickness,
            low_corner[2], high_corner[2]
        );

     
        auto x0 = particles.sliceReferencePosition();
        auto v_slice = particles.sliceVelocity();
        auto f_slice = particles.sliceForce();

        auto fix_velocity_to_zero_op = KOKKOS_LAMBDA(const int pid, const double /*time*/) {
            v_slice(pid,0)=v_slice(pid,1)=v_slice(pid,2)=0.0;
            f_slice(pid,0)=f_slice(pid,1)=f_slice(pid,2)=0.0;
            x(pid,0)=x0(pid,0);
            x(pid,1)=x0(pid,1);
            x(pid,2)=x0(pid,2);
        };

      
        auto bc_combined_fixed = createBoundaryCondition
        (
            fix_velocity_to_zero_op, // 
            exec_space{},            // 
            solver.particles,               // 
            true,                   
            bottom_region
        );        

        // --- 
        solver.init(bc_combined_fixed); 
        solver.run(bc_combined_fixed);

    }
    // ====================================================
    //  Simulation run without contact
    // ====================================================
    else
    {
        std::cout << "no contact"<< std::endl;
        CabanaPD::Particles particles(
            memory_space{}, model_type{}, low_corner, high_corner, num_cells,
            halo_width, Cabana::InitRandom{}, init_op, exec_space{} );

        auto rho = particles.sliceDensity();
        auto x = particles.sliceReferencePosition();
        auto v = particles.sliceVelocity();
        auto f = particles.sliceForce();


        auto init_functor = KOKKOS_LAMBDA( const int pid )
        {
             
            if ((x(pid, 2) - ball_center[2]) * (x(pid, 2) - ball_center[2]) +
                    (x(pid, 0) - ball_center[0]) * (x(pid, 0) - ball_center[0]) +
                    (x(pid, 1) - ball_center[1]) * (x(pid, 1) - ball_center[1]) <=
                ball_radius * ball_radius)
            {
                v(pid, 2) = -vz_ball; // impact velocity downwards
                    
                rho(pid) = rho_Al;
            }
            else
            {
                v(pid, 2) = 0.0; 
                 rho(pid) = rho_Cu;
            }
        };
        particles.updateParticles( exec_space{}, init_functor );

        CabanaPD::Solver solver( inputs, particles, force_model_Cu );
        solver.init();
        solver.run();
    }
}

// Initialize MPI+Kokkos.
int main( int argc, char* argv[] )
{
    MPI_Init( &argc, &argv );
    Kokkos::initialize( argc, argv );

    coldspray( argv[1] );

    Kokkos::finalize();
    MPI_Finalize();
}