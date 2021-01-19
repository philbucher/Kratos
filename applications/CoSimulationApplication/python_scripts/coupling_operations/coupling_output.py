# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

import json
from pathlib import Path

def Create(*args):
    return CouplingOutput(*args)

class JsonOutput:
    def __init__(self, model, settings):
        self.model = model
        self.model_part_name = settings["model_part_name"].GetString()
        self.output_path = Path(settings["output_path"].GetString())
        self.output_path.mkdir(exist_ok=True)

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

    def Initialize(self):
        self.model_part = self.model[self.model_part_name]

    def PrintOutput(self, output_file_name):
        output = {"variable_names" : [var.Name() for var in self.nodal_solution_step_data_variables]}
        for node in self.model_part.Nodes:
            output[node.Id] = [node.GetSolutionStepValue(var) for var in self.nodal_solution_step_data_variables]

        file_path = self.output_path / (output_file_name + ".json")

        with open(file_path, 'w') as out_file:
            json.dump(output, out_file)


class CouplingOutput(CoSimulationCouplingOperation):
    """This operation is used to output at different points in the coupling.
    TODO:
    - add support for json, hdf5, gid
    - add tests
    - more cleanup
    """
    def __init__(self, settings, solver_wrappers, process_info):
        super().__init__(settings, process_info)
        self.model = solver_wrappers[self.settings["solver"].GetString()].model
        self.execution_point = self.settings["execution_point"].GetString()
        model_part_name = self.settings["output_parameters"]["model_part_name"].GetString()
        self.base_output_file_name = "{}_{}_{}_".format(self.settings["solver"].GetString(), model_part_name, self.execution_point)

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

        self.step = 0 # this should come from self.process_info
        # TODO check if restarted. If not delete the folder => check self.process_info

        output_format = self.settings["output_format"].GetString()
        if output_format == "vtk":
            self.output = KM.VtkOutput(self.model[model_part_name], self.settings["output_parameters"])
        elif output_format == "json":
            self.output = JsonOutput(self.model, self.settings["output_parameters"])
        else:
            raise Exception('Currently only "vtk" and "json" are supported!')

    def Initialize(self):
        if hasattr(self.output, "Initialize"):
            self.output.Initialize()

    def InitializeSolutionStep(self):
        self.step += 1
        self.coupling_iteration = 0

        if self.execution_point == "initialize_solution_step":
            output_file_name = self.base_output_file_name + str(self.step)
            self.output.PrintOutput(output_file_name)

    def FinalizeSolutionStep(self):
        if self.execution_point == "finalize_solution_step":
            output_file_name = self.base_output_file_name + str(self.step)
            self.output.PrintOutput(output_file_name)

    def InitializeCouplingIteration(self):
        self.coupling_iteration += 1

        if self.execution_point == "initialize_coupling_iteration":
            output_file_name = self.base_output_file_name + "{}_{}".format(self.step, self.coupling_iteration)
            self.output.PrintOutput(output_file_name)

    def FinalizeCouplingIteration(self):
        if self.execution_point == "finalize_coupling_iteration":
            output_file_name = self.base_output_file_name + "{}_{}".format(self.step, self.coupling_iteration)
            self.output.PrintOutput(output_file_name)


    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "solver"            : "UNSPECIFIED",
            "execution_point"   : "UNSPECIFIED",
            "output_format"     : "vtk",
            "output_parameters" : { }
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults
