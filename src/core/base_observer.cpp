#include "../../include/core/base_observer.hpp"

#include <fstream>
#include <iostream>
#include <filesystem>

namespace fs = std::filesystem;
using namespace std;


//! Constructor for BaseObserver class
BaseObserver
:: BaseObserver (const nlohmann::json& config, int lx)
{
    m_bool_visualize = false;
    if (config.contains("visualization")) {
        m_bool_visualize = config["visualization"];
    }

    m_bool_dumpOut = false;
    if (config.contains("dump_output")) {
        m_bool_dumpOut = config["dump_output"];
    }

    m_meshName_suffix = "_lx"+std::to_string(lx);
}

//! Visualizes grid function in GLVis
void BaseObserver :: visualize (GridFunction& u) const
{
    socketstream sout;
    if (m_bool_visualize)
    {
        char vishost[] = "localhost";
        int  visport   = 19916;
        sout.open(vishost, visport);
        if (!sout)
        {
            cout << "Unable to connect to GLVis server at "
                 << vishost << ':' << visport << endl;
            m_bool_visualize = false;
            cout << "GLVis visualization disabled.\n";
        }
        else
        {
            sout.precision(m_precision);
            sout << "solution\n"
                 << *(u.FESpace()->GetMesh()) << u;
            sout << "pause\n";
            sout << flush;
        }
    }
}

//! Dumps mesh
void BaseObserver
:: dump_mesh (const std::shared_ptr<Mesh>& mesh) const
{
    if (m_bool_dumpOut)
    {
        std::string mesh_name
                = m_output_dir+"mesh"+m_meshName_suffix;
        std::cout << mesh_name << std::endl;

        std::ofstream mesh_ofs(mesh_name.c_str());
        mesh_ofs.precision(m_precision);
        mesh->Print(mesh_ofs);
        mesh_ofs.close();
    }
}


// End of file
