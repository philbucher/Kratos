//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Lorenzo Gracia
//
//

// External includes

// Project includes
#include "includes/table.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "includes/kratos_parameters.h"

#include "custom_utilities/streamlines_output_3D_utilities.hpp"
#include "custom_utilities/global_joint_stress_utility.hpp"
#include "custom_utilities/transfer_selfweight_stress_utility.hpp"
#include "custom_utilities/construction_utility.hpp"
#include "custom_utilities/mapping_variables_2D_utilities.hpp"
#include "custom_utilities/mapping_variables_3D_utilities.hpp"


namespace Kratos
{

namespace Python
{

inline
void AssignTimeActivation(
        ConstructionUtility& rThisUtil,
        std::string ThermalSubModelPartName,
        int phase, double time_activation)
{
    rThisUtil.AssignTimeActivation(ThermalSubModelPartName, phase, time_activation);
}

inline
void InitializeSolutionStep(
        ConstructionUtility& rThisUtil,
        std::string ThermalSubModelPartName,
        std::string MechanicalSubModelPartName,
        int phase)
{
    rThisUtil.InitializeSolutionStep(ThermalSubModelPartName, MechanicalSubModelPartName, phase);
}


inline
void ActiveHeatFluxNoorzai(
        ConstructionUtility& rThisUtil,
        Parameters& NoorzaiParameters)
{
    rThisUtil.ActiveHeatFluxNoorzai(NoorzaiParameters);
}

inline
void ActiveHeatFluxAzenha(
        ConstructionUtility& rThisUtil,
        Parameters& AzenhaParameters)
{
    rThisUtil.ActiveHeatFluxAzenha(AzenhaParameters);
}


void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    using namespace pybind11;

    typedef Table<double,double> TableType;

    class_< StreamlinesOutput3DUtilities >
    (m, "StreamlinesOutput3DUtilities")
    .def(init<>())
    .def("ComputeOutputStep",&StreamlinesOutput3DUtilities::ComputeOutputStep);

    class_< GlobalJointStressUtility >
    (m, "GlobalJointStressUtility")
    .def(init<ModelPart&, Parameters>())
    .def("ComputingGlobalStress",&GlobalJointStressUtility::ComputingGlobalStress);

    class_< TransferSelfweightStressUtility >
    (m ,"TransferSelfweightStressUtility")
    .def(init<>())
    .def("Transfer",&TransferSelfweightStressUtility::Transfer);

    class_< ConstructionUtility >
    (m, "ConstructionUtility")
    .def(init<ModelPart&, ModelPart&, TableType&, Parameters&>())
    .def("Initialize",&ConstructionUtility::Initialize)
    .def("AssignTimeActivation", AssignTimeActivation)
    .def("InitializeSolutionStep",InitializeSolutionStep)
    .def("SearchingFluxes",&ConstructionUtility::SearchingFluxes)
    .def("ActiveHeatFluxNoorzai",ActiveHeatFluxNoorzai)
    .def("ActiveHeatFluxAzenha",ActiveHeatFluxAzenha)
    .def("AfterOutputStep",&ConstructionUtility::AfterOutputStep);

    class_< MappingVariables2DUtilities >
    (m, "MappingVariables2DUtilities")
    .def(init<>())
    .def("MappingModelParts",&MappingVariables2DUtilities::MappingModelParts);

    class_< MappingVariables3DUtilities >
    (m, "MappingVariables3DUtilities")
    .def(init<>())
    .def("MappingModelParts",&MappingVariables3DUtilities::MappingModelParts);
}

}  // namespace Python.
} // Namespace Kratos
