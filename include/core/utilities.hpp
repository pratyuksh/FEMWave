#ifndef CORE_UTILITIES_HPP
#define CORE_UTILITIES_HPP

#include <nlohmann/json.hpp>
#include <fmt/format.h>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <Eigen/Dense>

#include "../core/config.hpp"


void write_json_file (std::string system_type,
                      const nlohmann::json& config,
                      Eigen::VectorXd &data_x1,
                      Eigen::VectorXi &data_x2,
                      Eigen::MatrixXd &data_y);


#endif // CORE_UTILITIES_HPP
