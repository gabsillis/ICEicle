#pragma once
#include <cmath>

namespace iceicle {
    
    /// @brief safely scale the finite difference epsilon 
    /// bounded below by epsilon
    /// @param epsilon the epsilon for finite difference 
    /// @param scale a zero or positive scale value 
    template<class T>
    auto scale_fd_epsilon(T epsilon, T scale) -> T {
//        return epsilon; // for debugging
        return std::max(epsilon, scale * epsilon);
    }

    namespace finite_difference {
        
        template< class func_t, class... args_t>
        auto forward_first_order(func_t, args_t...);

    }
}
