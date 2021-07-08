#include "../../include/core/utilities.hpp"
#include <assert.h>

namespace fs = std::filesystem;


void write_json_file (std::string system_type,
                      const nlohmann::json& config,
                      Eigen::VectorXd &data_x1,
                      Eigen::VectorXi &data_x2,
                      Eigen::MatrixXd &data_y)
{
    assert(data_x1.size() == data_x2.size());

    std::string problem_type = config["problem_type"];

    std::string base_out_dir = "../output";
    if (config.contains("base_out_dir")) {
        base_out_dir = config["base_out_dir"];
    }

    std::string sub_out_dir = "wave_"+problem_type;
    if (config.contains("sub_out_dir")) {
        sub_out_dir = config["sub_out_dir"];
    }

    std::string out_dir = base_out_dir+"/"+sub_out_dir+"/";

    fs::create_directories(out_dir);

    int deg_t = config["deg_t"];
    int deg1_x = config["deg1_x"];
    int deg2_x = config["deg2_x"];

    std::string error_type = "relL2";
    if (config.contains("error_type")) {
        error_type = config["error_type"];
    }

    int stabParamsType = 1;
    if (config.contains("stab_parameters_type")) {
        stabParamsType = config["stab_parameters_type"];
    }

    std::string outfile;
    if (system_type == "waveFG")
    {
        outfile = out_dir+"/convergenceFG"
                +"_"+problem_type
                +"_"+error_type
                +"_degt"+std::to_string(deg_t)
                +"_deg1x"+std::to_string(deg1_x)
                +"_deg2x"+std::to_string(deg2_x)
                +"_stab"+std::to_string(stabParamsType)
                +".json";
    }
    else if (system_type == "waveSG")
    {
        outfile = out_dir+"/convergenceSG"
                +"_"+problem_type
                +"_"+error_type
                +"_degt"+std::to_string(deg_t)
                +"_deg1x"+std::to_string(deg1_x)
                +"_deg2x"+std::to_string(deg2_x)
                +"_stab"+std::to_string(stabParamsType)
                +".json";
    }

    // set data names
    std::string data_x1_name;
    std::string data_x2_name;
    std::vector<std::string> data_y_names;
    if (system_type == "waveFG"
            || system_type == "waveSG")
    {
        data_x1_name = "h_max";
        data_x2_name = "ndofs";
        data_y_names.push_back("pressure");
        data_y_names.push_back("velocity");
    }

    auto json = nlohmann::json{};
    for (unsigned int i=0; i<data_x1.size(); i++) {
        json[data_x1_name][i] = data_x1(i);
        json[data_x2_name][i] = data_x2(i);
        for (unsigned int j=0; j<data_y.rows(); j++) {
            json[data_y_names[j]][i] = data_y(j,i);
        }
    }

    auto file = std::ofstream(outfile);
    assert(file.good());
    file << json.dump(2);
}


// End of file
