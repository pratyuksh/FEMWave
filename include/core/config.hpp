#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <nlohmann/json.hpp>

/// Returns the config stored at `config.json`.
// Note: this is easily misused as essentially a global variable.
nlohmann::json get_global_config(int argc, char* const argv[]);

nlohmann::json get_global_config(const std::string& fileName);

#endif // CONFIG_HPP
