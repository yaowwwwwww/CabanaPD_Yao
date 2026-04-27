/****************************************************************************
 * Copyright (c) 2022 by Oak Ridge National Laboratory                      *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of CabanaPD. CabanaPD is distributed under a           *
 * BSD 3-clause license. For the licensing terms see the LICENSE file in    *
 * the top-level directory.                                                 *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#ifndef CONTACTMODELS_H
#define CONTACTMODELS_H

#include <cmath>

#include <CabanaPD_Force.hpp>
#include <CabanaPD_Input.hpp>
#include <CabanaPD_Output.hpp>

namespace CabanaPD
{
/******************************************************************************
  Contact model
******************************************************************************/
struct ContactModel
{
    using base_model = Contact;
    using material_type = SingleMaterial;

    // Contact neighbor search radius.
    double radius;
    // Extend neighbor search radius to reuse lists.
    double radius_extend;

    ContactModel() {}

    // PD horizon
    // Contact radius
    ContactModel( const double _radius, const double _radius_extend )
        : radius( _radius )
        , radius_extend( _radius_extend ){};
};

/* Normal repulsion */
struct NormalRepulsionModel : public ContactModel
{
    using base_type = ContactModel;
    using base_model = base_type::base_model;
    using model_type = NormalRepulsionModel;
    using fracture_type = NoFracture;
    using thermal_type = TemperatureIndependent;

    double delta;
    using ContactModel::radius;
    using ContactModel::radius_extend;

    double c;
    double K;
    
     
    double r0;      // lj potential width sigma 1.05dx
    double beta;    //   β   
    double alpha1;    //     α  
    double particle_volume;
    double lj_linearize_force_density_cap;
    double lj_linearize_r;
    double lj_linearize_value;
    double lj_linearize_slope;
    double lj_linearize_force_density_max;
    double lj_linearize_force_density_at_zero;
    
        // parameters for CZM cohesive law
    double c_czm;   // cohesive scaling
    double sy;      // yield stretch
    double m_czm;   // exponential decay rate

    NormalRepulsionModel() {}
    NormalRepulsionModel( const double _delta, 
                          const double _radius,
                          const double radius_extend, 
                          const double _K,
                          const double _r0,
                          const double _beta,
                          const double _alpha,
                          const double _particle_volume,
                          const double _lj_linearize_force_density_cap,
                          const double _lj_linearize_force_density_max,
                          const double _c_czm,
                          const double _sy,
                          const double _m_czm )
        : ContactModel(_radius, radius_extend )
        , delta( _delta )
        , K( _K )
        , r0( _r0 )
        , beta( _beta )
        , alpha1( _alpha )
        , particle_volume( _particle_volume )
        , lj_linearize_force_density_cap( _lj_linearize_force_density_cap )
        , lj_linearize_r( 0.0 )
        , lj_linearize_value( 0.0 )
        , lj_linearize_slope( 0.0 )
        , lj_linearize_force_density_max( _lj_linearize_force_density_max )
        , lj_linearize_force_density_at_zero( 0.0 )
        , c_czm( _c_czm )
        , sy( _sy )
        , m_czm( _m_czm )
    {
        K = _K;
        // This could inherit from PMB (same c)
        c = 18.0 * K / ( 3.1415926  * delta * delta * delta * delta * alpha1);

        if ( lj_linearize_force_density_cap > 0.0 )
        {
            const double r_zero = r0 / std::pow( beta, 1.0 / 6.0 );
            double lo = 1.0e-14;
            double hi = r_zero;
            const double raw_lo = rawLJForceDensity( lo, particle_volume );
            const bool repulsion_is_negative = raw_lo < 0.0;
            const double cap =
                repulsion_is_negative ? -lj_linearize_force_density_cap
                                      : lj_linearize_force_density_cap;

            if ( repulsion_is_negative )
                lj_linearize_force_density_at_zero =
                    -lj_linearize_force_density_max;
            else
                lj_linearize_force_density_at_zero =
                    lj_linearize_force_density_max;

            const bool cap_is_reached =
                repulsion_is_negative ? raw_lo < cap : raw_lo > cap;
            if ( cap_is_reached )
            {
                for ( int iter = 0; iter < 100; ++iter )
                {
                    const double mid = 0.5 * ( lo + hi );
                    const double raw_mid =
                        rawLJForceDensity( mid, particle_volume );
                    const bool above_cap =
                        repulsion_is_negative ? raw_mid < cap : raw_mid > cap;
                    if ( above_cap )
                        lo = mid;
                    else
                        hi = mid;
                }
                lj_linearize_r = hi;
                lj_linearize_value =
                    rawLJForceDensity( lj_linearize_r, particle_volume );

                if ( lj_linearize_r > 0.0 )
                    lj_linearize_slope =
                        ( lj_linearize_value -
                          lj_linearize_force_density_at_zero ) /
                        lj_linearize_r;
            }
        }

    }

    KOKKOS_INLINE_FUNCTION
    auto forceCoeff( const double r, const double vol ) const
    {
        if ( r <= 1e-14 ) return 0.0;
         if ( r > radius ) return 0.0;
        // Contact "stretch"
        //const double sc = ( r - radius ) / delta;
        const double lj_force_density = linearizedLJForceDensity( r, vol );
        // const double Fc = lj_force_density * vol;
        
        //  CZM attraction (tensile)
        // double s = (r - 2.0e-6) / 2.0e-6;
        // double F_czm = 0.0; 

        // if ( s < 0.0 && s >= -sy )
        //     F_czm = c_czm * (-s); // linear elastic
        // else if ( s < -sy )
        //     F_czm = c_czm * sy * exp( -m_czm * ( -s - sy ) ); // exponential softening

        // combine: repulsion (compressive) + CZM attraction (tensile)
        // double Fc_total = Fc- F_czm;  // note minus: CZM acts in opposite (tensile) direction

        // Normal repulsion uses a 15 factor compared to the PMB force
        return lj_force_density;
    }

  private:
    KOKKOS_INLINE_FUNCTION
    double rawLJForceDensity( const double r, const double vol ) const
    {
        const double alpha =
            c * r0 * r0 * vol * vol / 72 / pow( beta, 7.0 / 3.0 );
        const double term13 = pow( r0 / r, 13 );
        const double term7 = pow( r0 / r, 7 );
        return ( ( 12.0 * alpha ) / r0 ) * ( beta * term7 - term13 ) / vol;
    }

    KOKKOS_INLINE_FUNCTION
    double rawLJForceDensityDerivative( const double r, const double vol ) const
    {
        const double alpha =
            c * r0 * r0 * vol * vol / 72 / pow( beta, 7.0 / 3.0 );
        const double pref = ( 12.0 * alpha ) / r0;
        return pref / vol *
               ( 13.0 * pow( r0, 13 ) / pow( r, 14 ) -
                 7.0 * beta * pow( r0, 7 ) / pow( r, 8 ) );
    }

    KOKKOS_INLINE_FUNCTION
    double linearizedLJForceDensity( const double r, const double vol ) const
    {
        const double raw = rawLJForceDensity( r, vol );
        if ( lj_linearize_force_density_cap <= 0.0 || lj_linearize_r <= 0.0 ||
             r >= lj_linearize_r )
            return raw;

        if ( lj_linearize_force_density_at_zero < 0.0 )
        {
            if ( raw >= -lj_linearize_force_density_cap ) return raw;
        }
        else if ( raw <= lj_linearize_force_density_cap )
            return raw;

        return lj_linearize_force_density_at_zero + lj_linearize_slope * r;
    }
};

} // namespace CabanaPD

#endif
