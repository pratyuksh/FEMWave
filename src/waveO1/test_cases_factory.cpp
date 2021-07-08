#include "../../include/waveO1/test_cases_factory.hpp"
#include <fmt/format.h>


//! Makes different test cases
std::shared_ptr<WaveO1TestCases>
make_waveO1_test_case(const nlohmann::json& config)
{
    const std::string problem_type
            = config["problem_type"];

    if (problem_type == "unit_square_test1") {
        return std::make_shared
                <WaveO1TestCase<UnitSquare_Test1>> (config);
    }
    else if (problem_type == "unit_square_test2") {
        return std::make_shared
                <WaveO1TestCase<UnitSquare_Test2>> (config);
    }
    else if (problem_type == "unit_square_test3") {
        return std::make_shared
                <WaveO1TestCase<UnitSquare_Test3>> (config);
    }
    else if (problem_type == "unit_square_test4") {
        return std::make_shared
                <WaveO1TestCase<UnitSquare_Test4>> (config);
    }
    else if (problem_type == "gammaShaped_test1") {
        return std::make_shared
                <WaveO1TestCase<GammaShaped_Test1>> (config);
    }
    else if (problem_type == "lShaped_test1") {
        return std::make_shared
                <WaveO1TestCase<LShaped_Test1>> (config);
    }
    else if (problem_type == "lShaped_test2") {
        return std::make_shared
                <WaveO1TestCase<LShaped_Test2>> (config);
    }
    /*else if (problem_type == "lShaped_test3") {
        return std::make_shared
                <WaveO1TestCase<LShaped_Test3>> (config);
    }
    else if (problem_type == "lShaped_test4") {
        return std::make_shared
                <WaveO1TestCase<LShaped_Test4>> (config);
    }
    else if (problem_type == "lShaped_test5") {
        return std::make_shared
                <WaveO1TestCase<LShaped_Test5>> (config);
    }*/
    else if (problem_type == "transmission_test1") {
        return std::make_shared
                <WaveO1TestCase<SquareTwoPiece_Test1>> (config);
    }
    else if (problem_type == "transmission_test2") {
        return std::make_shared
                <WaveO1TestCase<SquareTwoPiece_Test2>> (config);
    }

    throw std::runtime_error(fmt::format(
        "Unknown problem type for first-order Wave equation. "
        "[{}]", problem_type));
}
