#pragma once 
#include <functional>
#include <vector>
#include <span>
#include <string>
#include <iceicle/anomaly_log.hpp>

namespace iceicle::io {

    // @brief a representation of a function that takes in a span over the solution variables 
    // and fills a span over the variables to output
    template<class T, int neq>
    struct output_field_function {

        /// @brief the names of all the output fields
        std::vector<std::string> field_names;

        /// @brief the function that takes the span over solution variables 
        /// and fills a span over th
        std::function<void(std::span<T, neq>, std::span<T>)> field_func;

        /// @brief get the number of fields outputted
        [[nodiscard]] inline constexpr 
        auto size() const 
        -> std::size_t 
        { return field_names.size(); }

        /// @brief calling this as a function invokes the field_func
        void operator()(std::span<T, neq> u, std::span<T> fields) const
        { 
            if(fields.size() != size()){
                util::AnomalyLog::log_anomaly(
                        "fields span size does not match the size of the fields function... skipping evaluation.");
                return;
            }
            field_func(u, fields); 
        }
    };

}
