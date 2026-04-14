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

// Johnson-Cook plasticity for PMB:
// sigma = (A + B*eps_p^n) * rateFactor(epsdot) * (1 - (T*)^m)
// rateFactor (piecewise):
//   epsdot <= epsdot_u:
//     1 + C1*ln(epsdot/epsdot0)
//   epsdot > epsdot_u:
//     1 + C1*ln(epsdot_u/epsdot0) + C2*ln(epsdot/epsdot_u)
// T* = (T - T_ref) / (T_melt - T_ref), clamped to [0,1].
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
    // Rate term coefficients.
    double C;   // C1
    double C2;  // C2 (defaults to C1 when not provided)
    double eps_dot0;
    double eps_dot_u; // critical strain-rate for piecewise split
    // Dislocation drag parameters for sigma_d = (B_d / (rho_m b^2)) * epsdot_p.
    double drag_Bd;
    double drag_rho_mobile;
    double drag_burgers;
    // Adiabatic heating parameters.
    double cp_adiabatic;
    double taylor_quinney;
    double T_ref;
    double T_melt;
    double m_thermal;
    int sample_pid;
    double dt;
    Kokkos::View<int*, MemorySpace> _print_counter;
    // Point-level equivalent plastic strain (lagged) for JC hardening.
    Kokkos::View<double*, MemorySpace> _eps_p;
    Kokkos::View<double*, MemorySpace> _eps_p_prev;
    Kokkos::View<double*, MemorySpace> _temp_adiabatic;

    BaseForceModelPMB( PMB model, mechanics_type, MemorySpace,
                       const double delta, const double _K, const double _A,
                       const double _B, const double _n,
                       const double _C = 0.0, const double _eps_dot0 = 1.0,
                       const int _sample_pid = -1, const double _dt = 0.0,
                       const double _cp_adiabatic = 0.0,
                       const double _taylor_quinney = 0.9,
                       const double _T_ref = 298.0,
                       const double _T_melt = 1356.0,
                       const double _m_thermal = 1.09,
                       const double _C2 = -1.0,
                       const double _eps_dot_u = -1.0,
                       const double _drag_Bd = 0.0,
                       const double _drag_rho_mobile = 0.0,
                       const double _drag_burgers = 0.0 )
        : base_type( model, NoFracture{}, delta, _K )
        , base_plasticity_type()
        , A( _A )
        , B( _B )
        , n( _n )
        , C( _C )
        , C2( _C2 >= 0.0 ? _C2 : _C )
        , eps_dot0( _eps_dot0 )
        , eps_dot_u( _eps_dot_u )
        , drag_Bd( _drag_Bd )
        , drag_rho_mobile( _drag_rho_mobile )
        , drag_burgers( _drag_burgers )
        , cp_adiabatic( _cp_adiabatic )
        , taylor_quinney( _taylor_quinney )
        , T_ref( _T_ref )
        , T_melt( _T_melt )
        , m_thermal( _m_thermal )
        , sample_pid( _sample_pid )
        , dt( _dt )
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
        C = ( model1.C + model2.C ) / 2.0;
        C2 = ( model1.C2 + model2.C2 ) / 2.0;
        eps_dot0 = ( model1.eps_dot0 + model2.eps_dot0 ) / 2.0;
        eps_dot_u = ( model1.eps_dot_u + model2.eps_dot_u ) / 2.0;
        drag_Bd = ( model1.drag_Bd + model2.drag_Bd ) / 2.0;
        drag_rho_mobile =
            ( model1.drag_rho_mobile + model2.drag_rho_mobile ) / 2.0;
        drag_burgers = ( model1.drag_burgers + model2.drag_burgers ) / 2.0;
        cp_adiabatic = ( model1.cp_adiabatic + model2.cp_adiabatic ) / 2.0;
        taylor_quinney =
            ( model1.taylor_quinney + model2.taylor_quinney ) / 2.0;
        T_ref = ( model1.T_ref + model2.T_ref ) / 2.0;
        T_melt = ( model1.T_melt + model2.T_melt ) / 2.0;
        m_thermal = ( model1.m_thermal + model2.m_thermal ) / 2.0;
        sample_pid = -1;
        dt = 0.5 * ( model1.dt + model2.dt );
    }

    void updateBonds( const int num_local, const int max_neighbors )
    {
        base_plasticity_type::updateBonds( num_local, max_neighbors );
        Kokkos::realloc( _print_counter, 1 );
        Kokkos::deep_copy( _print_counter, 0 );
        Kokkos::realloc( _eps_p, num_local );
        Kokkos::realloc( _eps_p_prev, num_local );
        Kokkos::realloc( _temp_adiabatic, num_local );
        Kokkos::deep_copy( _eps_p, 0.0 );
        Kokkos::deep_copy( _eps_p_prev, 0.0 );
        Kokkos::deep_copy( _temp_adiabatic, 0.0 );
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
    double rateFactorLow( const double eps_p_dot ) const
    {
        if ( C == 0.0 || eps_dot0 <= 0.0 )
            return 1.0;

        const double eps_eff = Kokkos::fmax( eps_p_dot, eps_dot0 );
        return 1.0 + C * Kokkos::log( eps_eff / eps_dot0 );
    }

    KOKKOS_INLINE_FUNCTION
    double rateFactor( const double eps_p_dot ) const
    {
        if ( C == 0.0 || eps_dot0 <= 0.0 )
            return 1.0;

        const double eps_eff = Kokkos::fmax( eps_p_dot, eps_dot0 );

        // Legacy branch (or disabled piecewise split).
        if ( eps_dot_u <= eps_dot0 )
        {
            const double ratio = eps_eff / eps_dot0;
            return 1.0 + C * Kokkos::log( ratio );
        }

        // Piecewise JC rate law.
        if ( eps_eff <= eps_dot_u )
        {
            return 1.0 + C * Kokkos::log( eps_eff / eps_dot0 );
        }

        return 1.0 + C * Kokkos::log( eps_dot_u / eps_dot0 ) +
               C2 * Kokkos::log( eps_eff / eps_dot_u );
    }

    KOKKOS_INLINE_FUNCTION
    bool useDragStrength() const
    {
        return drag_Bd > 0.0 && drag_rho_mobile > 0.0 &&
               drag_burgers > 0.0;
    }

    // Accessors for point-level plasticity aggregation.
    auto plasticStretch() const { return _s_p; }
    auto pointPlasticStrain() { return _eps_p; }
    auto pointPlasticStrainPrev() { return _eps_p_prev; }
    auto pointPlasticStrain() const { return _eps_p; }
    auto pointPlasticStrainPrev() const { return _eps_p_prev; }
    auto pointAdiabaticTemperature() { return _temp_adiabatic; }
    auto pointAdiabaticTemperature() const { return _temp_adiabatic; }

    KOKKOS_INLINE_FUNCTION
    double pointPlasticStrainRate( const int i ) const
    {
        if ( dt <= 0.0 )
            return 0.0;
        const double rate = ( _eps_p( i ) - _eps_p_prev( i ) ) / dt;
        return ( rate > 0.0 ) ? rate : 0.0;
    }

    KOKKOS_INLINE_FUNCTION
    double thermalFactor( const int i ) const
    {
        if ( m_thermal <= 0.0 || T_melt <= T_ref )
            return 1.0;

        const double T = T_ref + _temp_adiabatic( i );
        double T_star = ( T - T_ref ) / ( T_melt - T_ref );
        T_star = Kokkos::fmax( 0.0, Kokkos::fmin( 1.0, T_star ) );
        const double soft = 1.0 - Kokkos::pow( T_star, m_thermal );
        return Kokkos::fmax( 0.0, soft );
    }

    KOKKOS_INLINE_FUNCTION
    double dragStress( const double eps_p_dot ) const
    {
        if ( !useDragStrength() )
            return 0.0;

        const double eps_eff = Kokkos::fmax( eps_p_dot, 0.0 );
        if ( eps_eff <= 0.0 )
            return 0.0;

        return drag_Bd * eps_eff /
               ( drag_rho_mobile * drag_burgers * drag_burgers );
    }

    KOKKOS_INLINE_FUNCTION
    double pointYieldStress( const int i ) const
    {
        const double eps_p = _eps_p( i );
        const double eps_p_dot = pointPlasticStrainRate( i );
        const double sigma_low_jc =
            hardeningYieldStress( eps_p ) * rateFactorLow( eps_p_dot ) *
            thermalFactor( i );

        return sigma_low_jc + dragStress( eps_p_dot );
    }

    KOKKOS_INLINE_FUNCTION
    auto operator()( ForceCoeffTag, const int i, const int j, const double s,
                     const double vol, const int n_id ) const
    {
        // Update bond plastic stretch.
        const double s_p = _s_p( i, n_id );
        const double eps_p_i = _eps_p( i );
        const int num_local = static_cast<int>( _eps_p.extent( 0 ) );
        const bool j_local = ( j < num_local );
        const double eps_p_j = j_local ? _eps_p( j ) : eps_p_i;
        const double sigma_y_i =
            pointYieldStress( i );
        const double sigma_y_j =
            j_local ? pointYieldStress( j ) : sigma_y_i;
        const double s_Y = 0.5 * ( sigma_y_i + sigma_y_j ) / ( 3.0 * K );

        // Yield in tension.
        if ( s >= s_p + s_Y )
            _s_p( i, n_id ) = s - s_Y;
        // Yield in compression.
        else if ( s <= s_p - s_Y )
            _s_p( i, n_id ) = s + s_Y;
        // else: Elastic (in between), do not modify.

        const double s_eff = s - _s_p( i, n_id );
        return c * s_eff * vol;
    }

    // This energy calculation is only valid for pure tension or pure
    // compression.
    KOKKOS_INLINE_FUNCTION
    auto operator()( EnergyTag, const int i, const int j, const double s,
                     const double xi, const double vol, const int n_id ) const
    {
        const double s_p = _s_p( i, n_id );
        const double eps_p_i = _eps_p( i );
        const int num_local = static_cast<int>( _eps_p.extent( 0 ) );
        const bool j_local = ( j < num_local );
        const double eps_p_j = j_local ? _eps_p( j ) : eps_p_i;
        const double sigma_y_i =
            pointYieldStress( i );
        const double sigma_y_j =
            j_local ? pointYieldStress( j ) : sigma_y_i;
        const double s_Y = 0.5 * ( sigma_y_i + sigma_y_j ) / ( 3.0 * K );
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
                const double C = 0.0, const double eps_dot0 = 1.0,
                const int sample_id = -1, const double dt = 0.0,
                const double cp_adiabatic = 0.0,
                const double taylor_quinney = 0.9,
                const double T_ref = 298.0,
                const double T_melt = 1356.0,
                const double m_thermal = 1.09,
                const double C2 = -1.0,
                const double eps_dot_u = -1.0,
                const double drag_Bd = 0.0,
                const double drag_rho_mobile = 0.0,
                const double drag_burgers = 0.0 )
        : base_type( model, mechanics, space, delta, K, A, B, n, C, eps_dot0,
                     sample_id, dt, cp_adiabatic, taylor_quinney, T_ref,
                     T_melt, m_thermal, C2, eps_dot_u, drag_Bd,
                     drag_rho_mobile, drag_burgers )
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
                const double n, const TemperatureType& temp, const double kappa,
                const double cp, const double alpha,
                const double temp0 = 0.0,
                const bool constant_microconductivity = true )
        : base_type( model, mechanics, typename TemperatureType::memory_space{},
                     delta, K, A, B, n )
        , base_temperature_type( temp, alpha, temp0 )
        , base_heat_transfer_type( delta, kappa, cp,
                                   constant_microconductivity )
    {
    }

    ForceModel( PMB model, JohnsonCook mechanics, const double delta,
                const double K, const double A, const double B, const double n,
                const TemperatureType& temp, const double kappa,
                const double cp, const double alpha, const double C,
                const double eps_dot0, const double dt,
                const double temp0 = 0.0,
                const bool constant_microconductivity = true )
        : base_type( model, mechanics, typename TemperatureType::memory_space{},
                     delta, K, A, B, n, C, eps_dot0, -1, dt, cp )
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

    ForceModel( PMB model, JohnsonCook mechanics, const double delta,
                const double K, const double G0, const double A,
                const double B, const double n, const TemperatureType& temp,
                const double kappa, const double cp, const double alpha,
                const double C, const double eps_dot0, const double dt,
                const double temp0 = 0.0,
                const bool constant_microconductivity = true )
        : base_type( model, mechanics, typename TemperatureType::memory_space{},
                     delta, K, A, B, n, C, eps_dot0, -1, dt, cp )
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

template <typename ModelType, typename MemorySpace>
ForceModel( ModelType, JohnsonCook, MemorySpace, const double delta,
            const double K, const double G0, const double A, const double B,
            const double n, const double C, const double eps_dot0,
            const int sample_id, const double dt,
            const double cp_adiabatic = 0.0,
            const double taylor_quinney = 0.9,
            const double T_ref = 298.0, const double T_melt = 1356.0,
            const double m_thermal = 1.09, const double C2 = -1.0,
            const double eps_dot_u = -1.0,
            const double drag_Bd = 0.0,
            const double drag_rho_mobile = 0.0,
            const double drag_burgers = 0.0 )
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
            const double A, const double B, const double n,
            const TemperatureType& temp, const double kappa, const double cp,
            const double alpha, const double C, const double eps_dot0,
            const double dt, const double temp0 = 0.0,
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

template <typename ModelType, typename TemperatureType>
ForceModel( ModelType, JohnsonCook, const double delta, const double K,
            const double G0, const double A, const double B, const double n,
            const TemperatureType& temp, const double kappa, const double cp,
            const double alpha, const double C, const double eps_dot0,
            const double dt, const double temp0 = 0.0,
            const bool constant_microconductivity = true )
    -> ForceModel<ModelType, JohnsonCook, Fracture, DynamicTemperature,
                  TemperatureType>;

} // namespace CabanaPD

#endif
