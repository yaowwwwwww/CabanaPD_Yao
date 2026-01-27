/****************************************************************************
 * Copyright (c) 2024 by Oak Ridge National Laboratory                      *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of CabanaPD. CabanaPD is distributed under a           *
 * BSD 3-clause license. For the licensing terms see the LICENSE file in    *
 * the top-level directory.                                                 *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#ifndef FORCE_MODELS_PMB_JOHNSONCOOK_H
#define FORCE_MODELS_PMB_JOHNSONCOOK_H

#include <cstdio>

#include <Kokkos_Core.hpp>

#include <CabanaPD_Constants.hpp>
#include <CabanaPD_ForceModels.hpp>
#include <CabanaPD_Types.hpp>
#include <force_models/CabanaPD_PMB.hpp>

namespace CabanaPD
{

// Johnson-Cook mechanics tag (kept local to this header).
struct JohnsonCook
{
};

// Johnson-Cook plasticity for PMB (bond-based), hardening only.
template <class MemorySpace>
struct BaseForceModelPMB<JohnsonCook, MemorySpace>
    : public BaseForceModelPMB<Elastic>, public BasePlasticity<MemorySpace>
{
    using base_type = BaseForceModelPMB<Elastic>;
    using base_plasticity_type = BasePlasticity<MemorySpace>;
    using mechanics_type = JohnsonCook;

    using base_type::K;
    using base_type::c;
    using base_plasticity_type::_s_p;
    using base_plasticity_type::updateBonds;

    double A;
    double B;
    double n;
    int sample_pid;
    Kokkos::View<int*, MemorySpace> _print_counter;

    BaseForceModelPMB( PMB model, mechanics_type, MemorySpace,
                       const double delta, const double _K, const double _A,
                       const double _B, const double _n,
                       const int _sample_pid = -1 )
        : base_type( model, NoFracture{}, delta, _K )
        , base_plasticity_type()
        , A( _A )
        , B( _B )
        , n( _n )
        , sample_pid( _sample_pid )
    {
    }

    // Constructor to average from existing models.
    template <typename ModelType1, typename ModelType2>
    BaseForceModelPMB( const ModelType1& model1, const ModelType2& model2 )
        : base_type( model1, model2 )
        , base_plasticity_type()
    {
        A = ( model1.A + model2.A ) / 2.0;
        B = ( model1.B + model2.B ) / 2.0;
        n = ( model1.n + model2.n ) / 2.0;
        sample_pid = -1;
    }

    void updateBonds( const int num_local, const int max_neighbors )
    {
        base_plasticity_type::updateBonds( num_local, max_neighbors );
        Kokkos::realloc( _print_counter, 1 );
        Kokkos::deep_copy( _print_counter, 0 );
    }

    KOKKOS_INLINE_FUNCTION
    double hardeningYieldStress( const double eps_p ) const
    {
        return A + B * Kokkos::pow( eps_p, n );
    }

    KOKKOS_INLINE_FUNCTION
    double yieldStretch( const double eps_p ) const
    {
        return hardeningYieldStress( eps_p ) / ( 3.0 * K );
    }

    KOKKOS_INLINE_FUNCTION
    auto operator()( ForceCoeffTag, const int i, const int, const double s,
                     const double vol, const int n_id ) const
    {
        // Update bond plastic stretch.
        const double s_p = _s_p( i, n_id );
        const double eps_p = Kokkos::abs( s_p );
        const double sigma_y = hardeningYieldStress( eps_p );
        const double s_Y = sigma_y / ( 3.0 * K );

        // Yield in tension.
        if ( s >= s_p + s_Y )
            _s_p( i, n_id ) = s - s_Y;
        // Yield in compression.
        else if ( s <= s_p - s_Y )
            _s_p( i, n_id ) = s + s_Y;
        // else: Elastic (in between), do not modify.

        if ( sample_pid >= 0 && i == sample_pid && n_id == 0 )
        {
            const int step =
                Kokkos::atomic_fetch_add( &_print_counter( 0 ), 1 );
            if ( step % 50 == 0 )
            {
                const double s_p_new = _s_p( i, n_id );
                printf(
                    "JC hardening step=%d eps_p=%e sigma_y=%e s_Y=%e s=%e s_p=%e s_p_new=%e\n",
                    step, eps_p, sigma_y, s_Y, s, s_p, s_p_new );
            }
        }

        const double s_eff = s - _s_p( i, n_id );
        return c * s_eff * vol;
    }

    // This energy calculation is only valid for pure tension or pure
    // compression.
    KOKKOS_INLINE_FUNCTION
    auto operator()( EnergyTag, const int i, const int, const double s,
                     const double xi, const double vol, const int n_id ) const
    {
        const double s_p = _s_p( i, n_id );
        const double eps_p = Kokkos::abs( s_p );
        const double s_Y = yieldStretch( eps_p );
        double stretch_term;

        // Yield in tension.
        if ( s >= s_p + s_Y )
            stretch_term = s_Y * ( 2.0 * s - s_Y );
        // Yield in compression.
        else if ( s <= s_p - s_Y )
            stretch_term = s_Y * ( -2.0 * s - s_Y );
        else
            // Elastic (in between).
            stretch_term = s * s;

        // 0.25 factor is due to 1/2 from outside the integral and 1/2 from
        // the integrand (pairwise potential).
        return 0.25 * c * stretch_term * xi * vol;
    }
};

// Base PMB + Johnson-Cook, no fracture/temperature.
template <typename MemorySpace>
struct ForceModel<PMB, JohnsonCook, NoFracture, TemperatureIndependent,
                  MemorySpace>
    : public BaseForceModelPMB<JohnsonCook, MemorySpace>,
      BaseNoFractureModel,
      BaseTemperatureModel<TemperatureIndependent>
{
    using base_type = BaseForceModelPMB<JohnsonCook, MemorySpace>;

    using base_type::operator();
    using base_type::updateBonds;
    using base_type::base_type;

    // Constructor to average from existing models.
    template <typename ModelType1, typename ModelType2>
    ForceModel( const ModelType1& model1, const ModelType2& model2 )
        : base_type( model1, model2 )
    {
    }
};

template <typename MemorySpace>
struct ForceModel<PMB, JohnsonCook, Fracture, TemperatureIndependent,
                  MemorySpace>
    : public BaseForceModelPMB<JohnsonCook, MemorySpace>,
      public BaseFractureModel,
      public BaseTemperatureModel<TemperatureIndependent>
{
    using base_type = BaseForceModelPMB<JohnsonCook, MemorySpace>;
    using base_fracture_type = BaseFractureModel;
    using base_temperature_type = BaseTemperatureModel<TemperatureIndependent>;

    using base_type::operator();
    using base_fracture_type::operator();
    using base_temperature_type::operator();
    using base_type::updateBonds;

    ForceModel( PMB model, JohnsonCook mechanics, MemorySpace space,
                const double delta, const double K, const double G0,
                const double A, const double B, const double n,
                const int sample_id = -1 )
        : base_type( model, mechanics, space, delta, K, A, B, n, sample_id )
        , base_fracture_type( delta, K, G0 )
        , base_temperature_type()
    {
    }

    // Constructor to average from existing models.
    template <typename ModelType1, typename ModelType2>
    ForceModel( const ModelType1& model1, const ModelType2& model2 )
        : base_type( model1, model2 )
        , base_fracture_type( model1, model2 )
    {
    }
};

template <typename TemperatureType>
struct ForceModel<PMB, JohnsonCook, NoFracture, DynamicTemperature,
                  TemperatureType>
    : public BaseForceModelPMB<JohnsonCook, typename TemperatureType::memory_space>,
      BaseNoFractureModel,
      BaseTemperatureModel<TemperatureDependent, TemperatureType>,
      BaseDynamicTemperatureModel
{
    using base_type =
        BaseForceModelPMB<JohnsonCook, typename TemperatureType::memory_space>;
    using base_temperature_type =
        BaseTemperatureModel<TemperatureDependent, TemperatureType>;
    using base_heat_transfer_type = BaseDynamicTemperatureModel;

    using thermal_type = DynamicTemperature;

    using base_type::operator();
    using base_temperature_type::operator();
    using base_type::updateBonds;

    ForceModel( PMB model, JohnsonCook mechanics, const double delta,
                const double K, const double A, const double B,
                const double n, const TemperatureType& temp,
                const double kappa, const double cp, const double alpha,
                const double temp0 = 0.0,
                const bool constant_microconductivity = true )
        : base_type( model, mechanics, typename TemperatureType::memory_space{},
                     delta, K, A, B, n )
        , base_temperature_type( temp, alpha, temp0 )
        , base_heat_transfer_type( delta, kappa, cp,
                                   constant_microconductivity )
    {
    }
};

template <typename TemperatureType>
struct ForceModel<PMB, JohnsonCook, Fracture, DynamicTemperature,
                  TemperatureType>
    : public BaseForceModelPMB<JohnsonCook, typename TemperatureType::memory_space>,
      public ThermalFractureModel<TemperatureType>,
      BaseDynamicTemperatureModel
{
    using base_type =
        BaseForceModelPMB<JohnsonCook, typename TemperatureType::memory_space>;
    using base_temperature_type = ThermalFractureModel<TemperatureType>;
    using base_heat_transfer_type = BaseDynamicTemperatureModel;

    using thermal_type = DynamicTemperature;

    using base_type::operator();
    using base_temperature_type::operator();
    using base_type::updateBonds;

    ForceModel( PMB model, JohnsonCook mechanics, const double delta,
                const double K, const double G0, const double A,
                const double B, const double n, const TemperatureType& temp,
                const double kappa, const double cp, const double alpha,
                const double temp0 = 0.0,
                const bool constant_microconductivity = true )
        : base_type( model, mechanics, typename TemperatureType::memory_space{},
                     delta, K, A, B, n )
        , base_temperature_type( delta, K, G0, temp, alpha, temp0 )
        , base_heat_transfer_type( delta, kappa, cp,
                                   constant_microconductivity )
    {
    }
};

template <typename ModelType, typename MemorySpace>
ForceModel( ModelType, JohnsonCook, MemorySpace, const double delta,
            const double K, const double A, const double B, const double n )
    -> ForceModel<ModelType, JohnsonCook, NoFracture, TemperatureIndependent,
                  MemorySpace>;

template <typename ModelType, typename MemorySpace>
ForceModel( ModelType, JohnsonCook, MemorySpace, const double delta,
            const double K, const double G0, const double A, const double B,
            const double n, const int sample_id = -1 )
    -> ForceModel<ModelType, JohnsonCook, Fracture, TemperatureIndependent,
                  MemorySpace>;

template <typename ModelType, typename TemperatureType>
ForceModel( ModelType, JohnsonCook, const double delta, const double K,
            const double A, const double B, const double n,
            const TemperatureType& temp, const double kappa, const double cp,
            const double alpha, const double temp0 = 0.0,
            const bool constant_microconductivity = true )
    -> ForceModel<ModelType, JohnsonCook, NoFracture, DynamicTemperature,
                  TemperatureType>;

template <typename ModelType, typename TemperatureType>
ForceModel( ModelType, JohnsonCook, const double delta, const double K,
            const double G0, const double A, const double B, const double n,
            const TemperatureType& temp, const double kappa, const double cp,
            const double alpha, const double temp0 = 0.0,
            const bool constant_microconductivity = true )
    -> ForceModel<ModelType, JohnsonCook, Fracture, DynamicTemperature,
                  TemperatureType>;

} // namespace CabanaPD

#endif
