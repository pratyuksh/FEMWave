#ifndef BASE_OBSERVER_HPP
#define BASE_OBSERVER_HPP

#include "mfem.hpp"
using namespace mfem;

#include "../core/config.hpp"


// Base Observer
class BaseObserver
{
public:
    BaseObserver () {}

    BaseObserver (const nlohmann::json&, int);
    
    void operator () (GridFunction& u) const {
        visualize(u);
    }

    void operator() (std::shared_ptr<GridFunction>& u) const {
        visualize(u);
    }

    void visualize (std::shared_ptr<GridFunction>& u) const{
        visualize(*u);
    }

    void visualize (GridFunction&) const;
    void dump_mesh (const std::shared_ptr<Mesh>&) const;
    
protected:
    int m_precision = 15;

    mutable bool m_bool_visualize;
    
    bool m_bool_dumpOut;
    mutable std::string m_output_dir;
    std::string m_meshName_prefix;
    std::string m_meshName_suffix;
};


#endif /// BASE_OBSERVER_HPP
