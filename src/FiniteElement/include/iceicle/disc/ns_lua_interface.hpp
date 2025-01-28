#pragma once
#include "iceicle/anomaly_log.hpp"
#include "iceicle/string_utils.hpp"
#include <iceicle/disc/navier_stokes.hpp>
#include <iceicle/tmp_utils.hpp>
#include <optional>
#include <sol/sol.hpp>
#include <type_traits>

namespace iceicle {
namespace navier_stokes {
namespace lua {

/// @brief parse free stream variables
/// from a lua table
///
/// The free stream table is a table that contains named values for
/// quantities that can be used to determine the free stream state
///
/// Case 1: Specify primitive
/// "rho" - density
/// "u" - u
/// "T" - temperature
/// "mu" - viscosity : Optional (determined from T)
/// "l" - reference length : Optional
///
/// Case 2: Reynolds number and mach number
/// "rho" - density
/// "T" - temperature
/// "Re" - Reynolds number
/// "mach" - mach number
/// "l" - reference length : Optional
///
template <class T, int ndim>
[[nodiscard]]
auto parse_free_stream(sol::table free_stream_table, T gamma, T Rgas)
-> std::optional<FreeStream<T, ndim>> 
{
    auto all_has_value = [free_stream_table](
            std::vector<const char *> required_params) -> bool 
    {
        for (std::string_view key : required_params) {
            sol::optional<T> val = free_stream_table[key];
            if (!val)
                return false;
        }

        return true;
    };

    if (all_has_value(std::vector{"rho", "u", "T"})) {
        T temp = free_stream_table["T"];
        T mu = free_stream_table.get_or("mu", Sutherlands<T>{}(temp));
        return std::optional{
            FreeStream<T, ndim>{.rho_inf = free_stream_table["rho"],
                                .u_inf = free_stream_table["u"],
                                .temp_inf = temp,
                                .mu_inf = mu,
                                .l_ref = free_stream_table["l"].get_or((T)1.0)}};
    }

    if (all_has_value(std::vector{"rho", "T", "Re", "mach"})) {
        T temp = free_stream_table["T"];
        T csound = std::sqrt(gamma * Rgas * temp);
        T mach = free_stream_table["mach"];
        T u = mach / csound;
        T l_ref = free_stream_table["l"].get_or((T)1.0);
        T rho = free_stream_table["rho"];
        T Re = free_stream_table["Re"];
        T mu = rho * u * l_ref / Re;

        return std::optional{navier_stokes::FreeStream<T, ndim>{
            .rho_inf = rho,
            .u_inf = u,
            .temp_inf = temp,
            .mu_inf = mu,
            .l_ref = l_ref}
        };
  }

  return std::nullopt;
}

/// @brief get the reference quantities for nondimensionalization 
/// from the conservation_law table
template <class real, int ndim>
[[nodiscard]] inline
auto ref_parameters(sol::table cons_law_tbl)
-> std::optional<ReferenceParameters<real>>
{
    using namespace util;

    // Option 1: keyword for reference parameters
    sol::optional<std::string> ref_param_name_opt =
      cons_law_tbl["reference_quantities"];

    if(eq_icase("free_stream", ref_param_name_opt.value_or(std::string{""}))) {

        // TODO: maybe use EoS instead of direct gamma and Rgas
        sol::table eos_tbl = cons_law_tbl["eos"];
        real gamma = eos_tbl.get_or("gamma", 1.4);
        real Rgas = eos_tbl.get_or("Rgas", 287.052874);

        // get the free stream quantities with error handling
        sol::optional<sol::table> fs_tbl_opt = cons_law_tbl["free_stream"];
        if(!fs_tbl_opt){
          AnomalyLog::log_anomaly("Free stream reference quantities"
                  "specified but no free stream quantities provided");
          return std::nullopt;
        }
        std::optional<FreeStream<real, ndim>> fs_opt 
          = parse_free_stream<real, ndim>(fs_tbl_opt.value(), gamma, Rgas);

        if(!fs_opt){
          AnomalyLog::log_anomaly("free_stream quantities could not be parsed."
                  "check that the specified values describe a complete free stream state"
                  " as outlined in the documentation");
          return std::nullopt;
        }

        // construct the reference quantities based on the free stream
        FreeStream<real, ndim> fs = fs_opt.value();
        return FreeStreamReference(fs.rho_inf, fs.u_inf, fs.temp_inf, fs.mu_inf, fs.l_ref);
    }

    // Option 2: manually specify the reference parameters
    sol::optional<sol::table> ref_param_tbl_opt =
        cons_law_tbl["reference_quantities"];
    if(ref_param_tbl_opt.has_value()){
        sol::table ref_param_tbl = ref_param_tbl_opt.value();

        auto all_has_value = [ref_param_tbl](
                std::vector<const char *> required_params) -> bool 
        {
            for (std::string_view key : required_params) {
                sol::optional<real> val = ref_param_tbl[key];
                if (!val)
                    return false;
            }
            return true;
        };

        if(all_has_value(std::vector{"t", "l", "rho", "u", "e", "p", "T", "mu"})) {
            // user specifies each reference quantity individually
            real t = ref_param_tbl["t"];
            real l = ref_param_tbl["l"];
            real rho = ref_param_tbl["rho"];
            real u = ref_param_tbl["u"];
            real e = ref_param_tbl["e"];
            real p = ref_param_tbl["p"];
            real T = ref_param_tbl["T"];
            real mu = ref_param_tbl["mu"];
            return ReferenceParameters<real>{t, l, rho, u, e, p, T, mu};
        } else {
            AnomalyLog::log_anomaly("reference_quantities table does not contain all "
                "the necessary parameters");
            return std::nullopt;
        }
    }

    // Default: unit reference quantities
    return ReferenceParameters<real>{};
}

/// @brief given the lua table that describes the conservation law
/// select the appropriate viscosity function and set it up
///
/// @param cons_law_tbl the lua table that describes the conservation law
/// @param ref reference quantities
template <class real>
[[nodiscard]]
auto select_viscosity(sol::table cons_law_tbl, ReferenceParameters<real> ref)
-> std::optional<std::function<real(real)>>
{
    sol::optional<sol::table> visc_tbl_opt = cons_law_tbl["viscosity"];
    if(!visc_tbl_opt.has_value()){
        util::AnomalyLog::log_anomaly("specify the viscosity in the 'viscosity' table "
                "of conservation_law");
        return std::nullopt;
    }

    sol::table visc_tbl = visc_tbl_opt.value();
    std::string visc_name = visc_tbl.get_or("name", std::string{""});

    // First case: constant viscosity
    if(util::eq_icase_any(visc_name, "constant", "const")){
        sol::optional<real> mu = visc_tbl["mu"];
        if(!mu.has_value()){
            util::AnomalyLog::log_anomaly(
                    "Must specify a value for 'mu' in the 'viscosity' table");
            return std::nullopt;
        }
        return constant_viscosity<real>{.mu = mu.value()};
    }

    // Second case: Dimensionless Sutherlands (preferred)
    if(util::eq_icase_any(visc_name, "sutherlands", "dimensionless sutherlands")) {
        DimensionlessSutherlands<real> visc{ref};
        return visc;
    }

    // Third case: Dimensional form of Sutherlands Law 
    // This allows manual specification of mu_0, T_0, and T_S
    // NOTE: the reference viscosity must be at the reference temperature 
    if(util::eq_icase_any(visc_name, "sutherlands dimensional")){
        Sutherlands<real> visc{};
        visc.mu_0 = visc_tbl.get_or("mu_0", visc.mu_0);
        visc.T_0 = visc_tbl.get_or("T_0", visc.T_0);
        visc.T_s = visc_tbl.get_or("T_s", visc.T_s);
        return visc;
    }

    util::AnomalyLog::log_anomaly(
        "name: '" + visc_name +"' in 'viscosity' unrecognized, could not parse");
    return std::nullopt;
}

/// @brief variant of the equations of state 
template<class real, int ndim>
using eos_options = std::variant<CaloricallyPerfectEoS<real, ndim>>;

/// @brief the the equation of state to use
template<class real, int ndim>
[[nodiscard]] inline 
auto select_eos(sol::table cons_law_tbl)
-> std::optional<eos_options<real, ndim>>
{
    sol::optional<sol::table> eos_tbl_opt = cons_law_tbl["eos"];
    if(eos_tbl_opt){
        sol::table eos_tbl = eos_tbl_opt.value();
        std::string name = eos_tbl.get_or("name", std::string{""});
        if(util::eq_icase_any(name, "ideal", "ideal gas", "ideal_gas", 
                    "calorically perfect", "calorically_perfect")){
            CaloricallyPerfectEoS<real, ndim> eos{};
            eos.gamma = eos_tbl.get_or("gamma", eos.gamma);
            eos.Rgas = eos_tbl.get_or("Rgas", eos.Rgas);
            return eos;
        } else {
            util::AnomalyLog::log_anomaly(
                "EoS name: " + name + "not recognized");
            return std::nullopt;
        }
    }

    // Default: CaloricallyPerfectEoS
    return CaloricallyPerfectEoS<real, ndim>{};
}

/// @brief variant of all the Physics options 
template< class real, int ndim >
using physics_options = std::variant<
    Physics<real, ndim, CaloricallyPerfectEoS<real, ndim>, VARSET::CONSERVATIVE>,
    Physics<real, ndim, CaloricallyPerfectEoS<real, ndim>, VARSET::RHO_U_P>,
    Physics<real, ndim, CaloricallyPerfectEoS<real, ndim>, VARSET::RHO_U_T>
>;

template< class real, int ndim >
[[nodiscard]] inline 
auto get_physics(sol::table cons_law_tbl) 
-> std::optional<physics_options<real, ndim>>
{
    // Get the reference quantities
    auto ref_opt = ref_parameters<real, ndim>(cons_law_tbl);
    if(!ref_opt){
        util::AnomalyLog::log_anomaly("Error initializing reference quantities");
        return std::nullopt;
    }
    ReferenceParameters<real> ref = ref_opt.value();

    // get the equation of state 
    auto eos_opt = select_eos<real, ndim>(cons_law_tbl);
    if(!eos_opt){
        util::AnomalyLog::log_anomaly("Error initializing Equation of State");
        return std::nullopt;
    }

    auto visc_opt = select_viscosity<real>(cons_law_tbl, ref);
    if(!visc_opt){
        util::AnomalyLog::log_anomaly("Error initializing viscosity function");
        return std::nullopt;
    }

    eos_opt.value() >> tmp::select_fcn{
        [&](const auto& eos) 
        -> std::optional<physics_options<real, ndim>> 
        {
            using eos_t = std::remove_cvref_t<decltype(eos)>;
            sol::optional<std::string> varset_name_opt = cons_law_tbl["varset"];
            if(varset_name_opt){
                std::string varset_name = varset_name_opt.value();
                if(util::eq_icase(varset_name, "conservative"))
                    return Physics<real, ndim, eos_t, VARSET::CONSERVATIVE>{ref, eos, visc_opt.value()};
                if(util::eq_icase_any(varset_name, "rho_u_p", "primitive pressure"))
                    return Physics<real, ndim, eos_t, VARSET::RHO_U_P>{ref, eos, visc_opt.value()};
                if(util::eq_icase_any(varset_name, "rho_u_t", "primitive temperature"))
                    return Physics<real, ndim, eos_t, VARSET::RHO_U_T>{ref, eos, visc_opt.value()};

                // should have found varset by now 
                util::AnomalyLog::log_anomaly("unrecognized varset: " + varset_name);
                return std::nullopt;
            }

            // Default: Conservative variables
            return Physics<real, ndim, eos_t, VARSET::CONSERVATIVE>{ref, eos, visc_opt.value()};
        }
    };

    return std::nullopt;
}

} // namespace lua
} // namespace navier_stokes

//
//         /// @brief parse the physical quantities from
//         /// the lua table describing the conservation law
//         template<class T>
//         auto parse_physical_quantities(sol::table cons_law_tbl)
//         -> std::tuple<T, T, T>
//         {
//             T gamma = cons_law_tbl.get_or("gamma", 1.4);
//             T Rgas = cons_law_tbl.get_or("Rgas", 287.052874);
//             T Pr = cons_law_tbl.get_or("Pr", 0.72);
//             return std::tuple{gamma, Rgas, Pr};
//         }
//
//         /// @brief parse the physics for NS from
//         /// the lua table describing the conservation law
//         template< class T, int ndim >
//         auto parse_physics(sol::table cons_law_tbl){
//             using namespace util;
//             // get physical quantities
//             auto [gamma, Rgas, Pr] =
//             parse_physical_quantities<T>(cons_law_tbl);
//
//             // get the free stream quantities
//             sol::optional<sol::table> free_stream_table_opt
//                 = cons_law_tbl["free_stream"];
//             std::optional<FreeStream<T, ndim>> free_stream_opt =
//             std::nullopt; if(free_stream_table_opt){
//                 free_stream_opt =
//                 parse_free_stream(free_stream_table_opt.value(), gamma,
//                 Rgas);
//             }
//
//             // set up the reference parameters
//             ReferenceParameters<T> ref;
//             sol::optional<std::string> ref_param_name_opt =
//             cons_law_tbl["reference_quantities"]; if(ref_param_name_opt){
//                 std::string ref_param_name = ref_param_name_opt.value();
//
//                 // Case 1: use free stream
//                 if(eq_icase(ref_param_name, "free_stream")){
//                     if(free_stream_opt){
//                         FreeStream<T, ndim> fs = free_stream_opt.value();
//                         ref = ReferenceParameters{
//                             fs.rho_inf, fs.u_inf, fs.temp_inf, fs.mu_inf,
//                             fs.l_ref};
//                     } else {
//                         util::AnomalyLog::log_anomaly("Free stream reference
//                         parameterizaton specified "
//                                 "but free stream state was not successfully
//                                 parsed.");
//                     }
//                 }
//             } else {
//                 ref = ReferenceParameters<T>{};
//             }
//
//             Physics<T, ndim> physics{ref, gamma, Rgas, Pr};
//             if(free_stream_opt)
//                 physics.free_stream = free_stream_opt.value();
//             return physics;
//         }
//     }
} // namespace iceicle
