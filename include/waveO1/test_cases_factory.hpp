#ifndef WAVEO1_TEST_CASES_FACTORY_HPP
#define WAVEO1_TEST_CASES_FACTORY_HPP

#include "test_cases.hpp"


//! Makes different test cases
std::shared_ptr<WaveO1TestCases>
make_waveO1_test_case(const nlohmann::json& config);


#endif /// WAVEO1_TEST_CASES_FACTORY_HPP
