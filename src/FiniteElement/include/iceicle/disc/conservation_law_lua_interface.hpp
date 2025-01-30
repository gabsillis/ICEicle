#pragma once

#include <iceicle/disc/burgers.hpp>
#include <iceicle/disc/navier_stokes.hpp>
#include "iceicle/anomaly_log.hpp"
#include "iceicle/disc/conservation_law.hpp"
#include "iceicle/disc/ns_lua_interface.hpp"
#include "iceicle/string_utils.hpp"
#include <sol/sol.hpp>


/// @brief namespace for the Lua API
namespace iceicle::lua {

    /// @brief given a conservation law table for the burgers equation 
    /// return a set up discretization for burgers equation
    template<class T, int ndim>
    auto set_up_burgers(sol::table cons_law_tbl)
    {
          // get the coefficients for burgers equation
          BurgersCoefficients<T, ndim> burgers_coeffs{};
          sol::optional<T> mu_input = cons_law_tbl["mu"];
          if (mu_input)
            burgers_coeffs.mu = mu_input.value();

          sol::optional<sol::table> b_adv_input =
              cons_law_tbl["b_adv"];
          sol::optional<sol::table> a_adv_input =
              cons_law_tbl["a_adv"];
          if (a_adv_input.has_value()) {
              for (int idim = 0; idim < ndim; ++idim)
                  burgers_coeffs.a[idim] =
                      cons_law_tbl["a_adv"][idim + 1];
          }
          if (b_adv_input.has_value()) {
              for (int idim = 0; idim < ndim; ++idim)
                  burgers_coeffs.b[idim] =
                      cons_law_tbl["b_adv"][idim + 1];
          }

//          std::cout << burgers_coeffs.mu 
//                    << " " << burgers_coeffs.a[0] 
//                    << " " << burgers_coeffs.b[0] 
//                << std::endl;
//
          // create the discretization
          BurgersFlux physical_flux{burgers_coeffs};
          BurgersUpwind convective_flux{burgers_coeffs};
          BurgersDiffusionFlux diffusive_flux{burgers_coeffs};
          ConservationLawDDG disc{std::move(physical_flux),
                                  std::move(convective_flux),
                                  std::move(diffusive_flux)};
          disc.field_names = std::vector<std::string>{"u"};
          disc.residual_names = std::vector<std::string>{"residual"};

          return disc;
    }

    /// @brief given a conservation law table for the spacetime burgers equation 
    /// return a set up discretization for burgers equation
    template<class T, int ndim>
    auto set_up_st_burgers(sol::table cons_law_tbl){

          static constexpr int ndim_space = ndim - 1;
          // get the coefficients for burgers equation
          BurgersCoefficients<T, ndim_space> burgers_coeffs{};
          sol::optional<T> mu_input = cons_law_tbl["mu"];
          if (mu_input)
              burgers_coeffs.mu = mu_input.value();

          sol::optional<sol::table> b_adv_input = cons_law_tbl["b_adv"];
          sol::optional<sol::table> a_adv_input = cons_law_tbl["a_adv"];
          if (a_adv_input.has_value()) {
              for (int idim = 0; idim < ndim_space; ++idim)
                  burgers_coeffs.a[idim] = cons_law_tbl["a_adv"][idim + 1];
          }
          if (b_adv_input.has_value()) {
              for (int idim = 0; idim < ndim_space; ++idim)
                  burgers_coeffs.b[idim] = cons_law_tbl["b_adv"][idim + 1];
          }
//          std::cout << burgers_coeffs.b[0] << std::endl;

          // create the discretization
          SpacetimeBurgersFlux physical_flux{burgers_coeffs};
          SpacetimeBurgersUpwind convective_flux{burgers_coeffs};
          SpacetimeBurgersDiffusion diffusive_flux{burgers_coeffs};
          ConservationLawDDG disc{std::move(physical_flux),
                                  std::move(convective_flux),
                                  std::move(diffusive_flux)};
          disc.field_names = std::vector<std::string>{"u"};
          disc.residual_names = std::vector<std::string>{"residual"};

          return disc;
    }
} // namespace iceicle::lua
