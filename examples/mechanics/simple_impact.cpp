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
      std::cout << "Running cold spray example with input file: "
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
    double rho0 = inputs["density"];
    double E = inputs["elastic_modulus"];
    double nu = 0.25; // Use bond-based model

    double K = E / ( 3 * ( 1 - 2 * nu ) );
    double G0 = inputs["fracture_energy"];
    double sigma_y = inputs["yield_stress"];


    double e = inputs["coefficient_of_restitution"];
    double gamma = inputs["surface_energy"];

    double sc = inputs["critical_stretch"];
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
    CabanaPD::ForceModel force_model( model_type{}, mechanics_type{},
                                      memory_space{}, delta, K, G0, sigma_y );

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
        // 球体区域
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
        using contact_type = CabanaPD::HertzianJKRModel;
        CabanaPD::Particles particles(
            memory_space{}, contact_type{}, low_corner, high_corner, num_cells,
            halo_width, Cabana::InitRandom{}, init_op, exec_space{} );

        auto rho = particles.sliceDensity();
        auto x = particles.sliceReferencePosition();
        auto v = particles.sliceVelocity();
        auto f = particles.sliceForce();
        auto dx = particles.dx ;



        auto init_functor = KOKKOS_LAMBDA( const int pid )
        {
            rho(pid) = rho0;
            if ((x(pid, 2) - ball_center[2]) * (x(pid, 2) - ball_center[2]) +
                    (x(pid, 0) - ball_center[0]) * (x(pid, 0) - ball_center[0]) +
                    (x(pid, 1) - ball_center[1]) * (x(pid, 1) - ball_center[1]) <=
                ball_radius * ball_radius)
            {
                v(pid, 2) = -vz_ball; // 冲击速度向下
            }
            else
            {
                v(pid, 2) = 0.0; // 板静止
            }
        };
        particles.updateParticles( exec_space{}, init_functor );

        // Use contact radius and extension relative to particle spacing.
        double r_c = inputs["contact_horizon_factor"];
        double r_extend = inputs["contact_horizon_extend_factor"];
        // NOTE: dx/2 is when particles first touch.
        r_c *= dx[0] / 2.0;
        r_extend *= dx[0];

        contact_type contact_model( r_c, r_extend, nu, E, e, gamma );;

        CabanaPD::Solver solver( inputs, particles, force_model,
                                 contact_model );

             double boundary_layer_thickness = 0.003; 
        double bottom_z_thickness = 0.003; 

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

        // --- 2. 获取速度切片，并定义一个 UserFunctor 来将速度设置为 0.0 (实现固定位移) ---
        // 关键：获取速度切片
        auto v_slice = particles.sliceVelocity(); 

        auto fix_velocity_to_zero_op = KOKKOS_LAMBDA( const int pid, const double /*time*/ )
        {
            // 将粒子的 X, Y, Z 三个方向的速度都设置为 0.0
            v_slice( pid, 0 ) = 0.0; // 固定 X 方向的位移 (通过固定速度)
            v_slice( pid, 1 ) = 0.0; // 固定 Y 方向的位移 (通过固定速度)
            v_slice( pid, 2 ) = 0.0; // 固定 Z 方向的位移 (通过固定速度)
        };

        // --- 3. 创建一个单一的边界条件对象，将 UserFunctor 和所有区域作为参数传入 ---
        auto bc_combined_fixed = createBoundaryCondition(
            fix_velocity_to_zero_op, // 第一个参数：UserFunctor
            exec_space{},            // 第二个参数：执行空间
            particles,               // 第三个参数：粒子数据
            true,                    // 第四个参数：布尔值 (表示每次都强制更新)
            bottom_region
        //     ,           // 区域 1
        //    x_min_side_region,       // 区域 2
        //    x_max_side_region,       // 区域 3
        //    y_min_side_region,       // 区域 4
        //    y_max_side_region        // 区域 5
        );        

        // --- 4. 将这一个边界条件对象传入求解器 ---
        solver.init(bc_combined_fixed); 
        solver.run(bc_combined_fixed);

    }
    // ====================================================
    //  Simulation run without contact
    // ====================================================
    else
    {
        CabanaPD::Particles particles(
            memory_space{}, model_type{}, low_corner, high_corner, num_cells,
            halo_width, Cabana::InitRandom{}, init_op, exec_space{} );

        auto rho = particles.sliceDensity();
        auto x = particles.sliceReferencePosition();
        auto v = particles.sliceVelocity();
        auto f = particles.sliceForce();


        auto init_functor = KOKKOS_LAMBDA( const int pid )
        {
            rho(pid) = rho0;
            if ((x(pid, 2) - ball_center[2]) * (x(pid, 2) - ball_center[2]) +
                    (x(pid, 0) - ball_center[0]) * (x(pid, 0) - ball_center[0]) +
                    (x(pid, 1) - ball_center[1]) * (x(pid, 1) - ball_center[1]) <=
                ball_radius * ball_radius)
            {
                v(pid, 2) = -vz_ball; // 冲击速度向下
            }
            else
            {
                v(pid, 2) = 0.0; // 板静止
            }
        };
        particles.updateParticles( exec_space{}, init_functor );

        CabanaPD::Solver solver( inputs, particles, force_model );
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
