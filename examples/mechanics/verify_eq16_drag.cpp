#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <Kokkos_Core.hpp>
#include <nlohmann/json.hpp>

#include <CabanaPD_Force.hpp>
#include <force_models/CabanaPD_PMB_JohnsonCook.hpp>

namespace
{
using json = nlohmann::json;

double scalarValue( const json& j, const char* key )
{
    return j.at( key ).at( "value" ).get<double>();
}

double arrayValue( const json& j, const char* key, const int idx )
{
    return j.at( key ).at( "value" ).at( idx ).get<double>();
}

struct MaterialInputs
{
    std::string name;
    double A;
    double B;
    double n;
    double C;
    double C2;
    double epsdot0;
    double epsdot_u;
    double T_ref;
    double T_melt;
    double m_thermal;
    double drag_Bd;
    double drag_rho_mobile;
    double drag_burgers;
};

template <class ModelType>
void verifyMaterial( const MaterialInputs& mat, const ModelType& model,
                     std::ostream& out )
{
    const std::vector<double> rates = { 1.01 * mat.epsdot_u, 1.50 * mat.epsdot_u,
                                        2.00 * mat.epsdot_u, 1.0e7, 4.0e7 };

    out << "# material=" << mat.name << " epsdot_u=" << mat.epsdot_u << "\n";
    out << "material\tepsdot\tsigma_drag_code\tsigma_drag_manual\tabs_diff\trel_diff\n";
    for ( const double rate : rates )
    {
        const double sigma_code = model.dragStress( rate );
        const double sigma_manual =
            mat.drag_Bd * rate /
            ( mat.drag_rho_mobile * mat.drag_burgers * mat.drag_burgers );
        const double abs_diff = std::abs( sigma_code - sigma_manual );
        const double rel_diff =
            ( std::abs( sigma_manual ) > 0.0 ) ? abs_diff / std::abs( sigma_manual )
                                               : abs_diff;

        out << mat.name << "\t" << std::setprecision( 10 ) << rate << "\t"
            << sigma_code << "\t" << sigma_manual << "\t" << abs_diff << "\t"
            << rel_diff << "\n";
    }
}

} // namespace

int main( int argc, char* argv[] )
{
    if ( argc < 2 || argc > 3 )
    {
        std::cerr << "usage: VerifyEq16Drag INPUT_JSON [OUT_TSV]\n";
        return 2;
    }

    Kokkos::initialize( argc, argv );

    int rc = 0;
    {
        std::ifstream input( argv[1] );
        if ( !input )
        {
            std::cerr << "failed to open input json: " << argv[1] << "\n";
            Kokkos::finalize();
            return 1;
        }

        json j;
        input >> j;

        std::ostream* out_ptr = &std::cout;
        std::ofstream output_file;
        if ( argc == 3 )
        {
            output_file.open( argv[2] );
            if ( !output_file )
            {
                std::cerr << "failed to open output tsv: " << argv[2] << "\n";
                Kokkos::finalize();
                return 1;
            }
            out_ptr = &output_file;
        }
        auto& out = *out_ptr;

        const double E = arrayValue( j, "elastic_modulus", 0 );
        const double nu = arrayValue( j, "Poisson's_ratio", 0 );
        const double K = E / ( 3.0 * ( 1.0 - 2.0 * nu ) );
        const double G0 = arrayValue( j, "fracture_energy", 0 );
        const double delta = scalarValue( j, "horizon" );
        const double dt = scalarValue( j, "timestep" );
        const double cp = arrayValue( j, "specific_heat_capacity", 0 );
        const double tq =
            j.contains( "taylor_quinney" ) ? scalarValue( j, "taylor_quinney" ) : 0.9;

        using memory_space = Kokkos::HostSpace;

        const std::vector<MaterialInputs> mats = {
            { "Al", arrayValue( j, "yield_stress", 0 ), arrayValue( j, "jc_B", 0 ),
              arrayValue( j, "jc_n", 0 ), arrayValue( j, "jc_C", 0 ),
              j.contains( "jc_C2" ) ? arrayValue( j, "jc_C2", 0 )
                                     : arrayValue( j, "jc_C", 0 ),
              arrayValue( j, "jc_epsdot0", 0 ), j.contains( "jc_epsdot_u" )
                                      ? arrayValue( j, "jc_epsdot_u", 0 )
                                      : -1.0,
              arrayValue( j, "reference_temperature", 0 ),
              j.contains( "jc_Tmelt" ) ? arrayValue( j, "jc_Tmelt", 0 )
                                       : arrayValue( j, "reference_temperature", 0 ) +
                                             1.0,
              j.contains( "jc_m" ) ? arrayValue( j, "jc_m", 0 ) : 0.0,
              j.contains( "drag_Bd" ) ? arrayValue( j, "drag_Bd", 0 ) : 0.0,
              j.contains( "drag_mobile_dislocation_density" )
                  ? arrayValue( j, "drag_mobile_dislocation_density", 0 )
                  : 0.0,
              j.contains( "drag_burgers_vector" )
                  ? arrayValue( j, "drag_burgers_vector", 0 )
                  : 0.0 },
            { "Cu", arrayValue( j, "yield_stress", 1 ), arrayValue( j, "jc_B", 1 ),
              arrayValue( j, "jc_n", 1 ), arrayValue( j, "jc_C", 1 ),
              j.contains( "jc_C2" ) ? arrayValue( j, "jc_C2", 1 )
                                     : arrayValue( j, "jc_C", 1 ),
              arrayValue( j, "jc_epsdot0", 1 ), j.contains( "jc_epsdot_u" )
                                      ? arrayValue( j, "jc_epsdot_u", 1 )
                                      : -1.0,
              arrayValue( j, "reference_temperature", 1 ),
              j.contains( "jc_Tmelt" ) ? arrayValue( j, "jc_Tmelt", 1 )
                                       : arrayValue( j, "reference_temperature", 1 ) +
                                             1.0,
              j.contains( "jc_m" ) ? arrayValue( j, "jc_m", 1 ) : 0.0,
              j.contains( "drag_Bd" ) ? arrayValue( j, "drag_Bd", 1 ) : 0.0,
              j.contains( "drag_mobile_dislocation_density" )
                  ? arrayValue( j, "drag_mobile_dislocation_density", 1 )
                  : 0.0,
              j.contains( "drag_burgers_vector" )
                  ? arrayValue( j, "drag_burgers_vector", 1 )
                  : 0.0 } };

        out << "# VerifyEq16Drag\n";
        out << "# note: current Eq.16 implementation does not switch on at epsdot_u;\n";
        out << "# this check evaluates dragStress(rate) for rates chosen above epsdot_u.\n";

        for ( const auto& mat : mats )
        {
            CabanaPD::ForceModel model(
                CabanaPD::PMB{}, CabanaPD::JohnsonCook{}, memory_space{}, delta, K,
                G0, mat.A, mat.B, mat.n, mat.C, mat.epsdot0, -1, dt, cp, tq,
                mat.T_ref, mat.T_melt, mat.m_thermal, mat.C2, mat.epsdot_u,
                mat.drag_Bd, mat.drag_rho_mobile, mat.drag_burgers );

            verifyMaterial( mat, model, out );
        }
    }

    Kokkos::finalize();
    return rc;
}
