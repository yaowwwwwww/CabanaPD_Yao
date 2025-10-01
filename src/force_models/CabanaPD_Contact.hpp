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
    NormalRepulsionModel() {}
    NormalRepulsionModel( const double _delta, 
                          const double _radius,
                          const double radius_extend, 
                          const double _K,
                          const double _r0,
                          const double _beta,
                          const double _alpha  )
        : ContactModel(_radius, radius_extend )
        , delta( _delta )
        , K( _K )
        , r0( _r0 )
        , beta( _beta )
        , alpha1( _alpha )
    {
        K = _K;
        // This could inherit from PMB (same c)
        c = 18.0 * K / ( 3.1415926  * delta * delta * delta * alpha1 );

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
        
        // Normal repulsion uses a 15 factor compared to the PMB force
        return Fc/vol;
    }
};

} // namespace CabanaPD

#endif
