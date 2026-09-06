#ifndef HAMPR_CORE_EXCEPTION_HPP
#define HAMPR_CORE_EXCEPTION_HPP

#include <stdexcept>
#include <string>

namespace hampr {

class HamprException : public std::runtime_error {
public:
    explicit HamprException(const std::string& msg)
        : std::runtime_error(msg) {}
};

} // namespace hampr

#endif // HAMPR_CORE_EXCEPTION_HPP
