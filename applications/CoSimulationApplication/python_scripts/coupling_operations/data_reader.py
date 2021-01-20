# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

import json
from pathlib import Path

def Create(*args):
    return DataReader(*args)

class JsonInput:
    def __init__(self, model_part, settings):
        self.model_part = model_part

        nodal_solution_step_data_variable_names = settings["nodal_solution_step_data_variables"].GetStringArray()
        self.nodal_solution_step_data_variables = []

        for variable_name in nodal_solution_step_data_variable_names:
            variable = KM.KratosGlobals.GetVariable(variable_name)
            variable_type = KM.KratosGlobals.GetVariableType(variable_name)
            if variable_type == "Double":
                self.nodal_solution_step_data_variables.append(variable)
            elif variable_type == "Array":
                self.nodal_solution_step_data_variables.append(KM.KratosGlobals.GetVariable(variable_name+"_X"))
                self.nodal_solution_step_data_variables.append(KM.KratosGlobals.GetVariable(variable_name+"_Y"))
                self.nodal_solution_step_data_variables.append(KM.KratosGlobals.GetVariable(variable_name+"_Z"))
            else:
                raise Exception("Wrong variable type!")

    def ReadAndAssignData(self, input_file_name):
        with open(input_file_name, 'r') as data_file:
            data = json.load(data_file)

        variable_indices = {var.Name() : data["variable_names"].index(var.Name()) for var in self.nodal_solution_step_data_variables}

        for node in self.model_part.Nodes:
            nodal_data = data[str(node.Id)]
            for var in self.nodal_solution_step_data_variables:
                value = nodal_data[variable_indices[var.Name()]]
                node.SetSolutionStepValue(var, 0, value)



class DataReader(CoSimulationCouplingOperation):
    """This operation is used to read and assing values from files
    Useful for e.g. one-way coupling when the data was previously saved
    TODO:
    - add support for json, hdf5, gid
    - add tests
    - more cleanup
    """
    def __init__(self, settings, solver_wrappers, process_info):
        super().__init__(settings, process_info)
        self.model = solver_wrappers[settings["solver"].GetString()].model
        self.execution_point = settings["execution_point"].GetString()

        self.input_path = Path(settings["input_path"].GetString())
        self.base_file_name = settings["file_name"].GetString()

        available_execution_points = [
            "initialize_solution_step",
            "finalize_solution_step",
            "initialize_coupling_iteration",
            "finalize_coupling_iteration"
        ]

        if self.execution_point not in available_execution_points:
            err_msg  = 'Execution point "{}" is not available, only the following options are available:\n    '.format(self.execution_point)
            err_msg += "\n    ".join(available_execution_points)
            raise Exception(err_msg)

        self.step = 0

    def Initialize(self):
        data_format = self.settings["data_format"].GetString()
        if data_format == "json":
            self.input = JsonInput(self.model[self.settings["model_part_name"].GetString()], self.settings["input_parameters"])
        else:
            raise Exception('Currently only and "json" is supported!')

    def InitializeSolutionStep(self):
        self.step += 1
        self.coupling_iteration = 0

        if self.execution_point == "initialize_solution_step":
            file_name = self.input_path / self.base_file_name.format(self.step)
            self.input.ReadAndAssignData(file_name)

    def FinalizeSolutionStep(self):
        if self.execution_point == "finalize_solution_step":
            file_name = self.input_path / self.base_file_name.format(self.step)
            self.input.ReadAndAssignData(input_file_name)

    def InitializeCouplingIteration(self):
        self.coupling_iteration += 1

        if self.execution_point == "initialize_coupling_iteration":
            raise NotImplementedError
            output_file_name = self.base_output_file_name + "{}_{}".format(self.step, self.coupling_iteration)
            self.input.PrintOutput(output_file_name)

    def FinalizeCouplingIteration(self):
        if self.execution_point == "finalize_coupling_iteration":
            raise NotImplementedError
            output_file_name = self.base_output_file_name + "{}_{}".format(self.step, self.coupling_iteration)
            self.input.PrintOutput(output_file_name)

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "solver"            : "UNSPECIFIED",
            "execution_point"   : "UNSPECIFIED",
            "model_part_name"   : "UNSPECIFIED",
            "data_format"       : "json",
            "input_path"        : "",
            "file_name"         : "UNSPECIFIED",
            "input_parameters"  : {}
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults
