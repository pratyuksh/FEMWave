#ifndef WAVEO2_TEST_CASES_FACTORY_HPP
#define WAVEO2_TEST_CASES_FACTORY_HPP

#include "test_cases.hpp"


//! Makes different test cases
std::shared_ptr<WaveO2TestCases>
make_waveO2_test_case(const nlohmann::json& config);


#endif /// WAVEO2_TEST_CASES_FACTORY_HPP
