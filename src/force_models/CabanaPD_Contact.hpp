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
                          const double _c_czm,
                          const double _sy,
                          const double _m_czm )
        : ContactModel(_radius, radius_extend )
        , delta( _delta )
        , K( _K )
        , r0( _r0 )
        , beta( _beta )
        , alpha1( _alpha )
        , c_czm( _c_czm )
        , sy( _sy )
        , m_czm( _m_czm )
    {
        K = _K;
        // This could inherit from PMB (same c)
        c = 18.0 * K / ( 3.1415926  * delta * delta * delta * delta * alpha1);

    }

    KOKKOS_INLINE_FUNCTION
    auto forceCoeff( const double r, const double vol ) const
    {
        if ( r <= 1e-14 ) return 0.0;
         if ( r > radius ) return 0.0;
        // Contact "stretch"
        //const double sc = ( r - radius ) / delta;
        double alpha = c * r0 * r0 * vol * vol / 72 / pow(beta, 7.0/3.0);
        double term13 = pow( r0 / r, 13 );
        double term7  = pow( r0 / r, 7 );

        double Fc = ( (12.0 * alpha)/r0   ) * ( term13 - beta * term7 );
        
        //  CZM attraction (tensile)
        double s = (r - 2.0e-6) / 2.0e-6;
        double F_czm = 0.0; 

        if ( s < 0.0 && s >= -sy )
            F_czm = c_czm * (-s); // linear elastic
        else if ( s < -sy )
            F_czm = c_czm * sy * exp( -m_czm * ( -s - sy ) ); // exponential softening

        // combine: repulsion (compressive) + CZM attraction (tensile)
        double Fc_total = Fc- F_czm;  // note minus: CZM acts in opposite (tensile) direction

        // Normal repulsion uses a 15 factor compared to the PMB force
        return Fc_total/vol;
    }
};

} // namespace CabanaPD

#endif
