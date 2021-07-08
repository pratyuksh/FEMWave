#include "../../include/waveO2/test_cases_factory.hpp"
#include <fmt/format.h>


//! Makes different test cases
std::shared_ptr<WaveO2TestCases>
make_waveO2_test_case(const nlohmann::json& config)
{
    const std::string problem_type
            = config["problem_type"];

    if (problem_type == "unit_square_test1") {
        return std::make_shared
                <WaveO2TestCase<UnitSquare_Test1>> (config);
    }

    throw std::runtime_error
            (fmt::format("Unknown problem type for "
                         "second-order Wave equation. "
                         "[{}]", problem_type));
}
