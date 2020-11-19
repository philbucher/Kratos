# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics import kratos_utilities
default_data_comm = KM.DataCommunicator.GetDefault()
if default_data_comm.IsDistributed():
    import KratosMultiphysics.mpi as KratosMPI

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_io import CoSimulationIO

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication as KratosCoSim
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
from KratosMultiphysics.CoSimulationApplication.utilities import model_part_utilities

# Other imports
import os

def Create(settings, model, solver_name):
    return EmpireIO(settings, model, solver_name)

communication_folder = ".EmpireIO" # hardcoded in C++

class EmpireIO(CoSimulationIO):
    """IO for the legacy EMPIRE_API
    """
    def __init__(self, settings, model, solver_name):
        super().__init__(settings, model, solver_name)
        self.rank = default_data_comm.Rank()
        self.is_distributed = default_data_comm.IsDistributed()
        if self.rank == 0:
            # Note: calling "EMPIRE_API_Connect" is NOT necessary, it is replaced by the next two lines
            KratosCoSim.EMPIRE_API.EMPIRE_API_SetEchoLevel(self.echo_level)
            KratosCoSim.EMPIRE_API.EMPIRE_API_PrintTiming(self.settings["api_print_timing"].GetBool())

            # delete and recreate communication folder to avoid leftover files
            kratos_utilities.DeleteDirectoryIfExisting(communication_folder)
            os.mkdir(communication_folder)

        self.aux_model = KM.Model()
        if self.is_distributed:
            # creating auxiliar Model for gathering ModelParts on rank 0
            self.gather_utilities = {}

        default_data_comm.Barrier()

    def Finalize(self):
        kratos_utilities.DeleteDirectoryIfExisting(communication_folder)

    def __del__(self):
        # make sure no communication files are left even if simulation is terminated prematurely
        if os.path.isdir(communication_folder) and self.rank == 0:
            kratos_utilities.DeleteDirectoryIfExisting(communication_folder)
            if self.echo_level > 0:
                cs_tools.cs_print_info(self._ClassName(), "Deleting Communication folder in destructor")

    def ImportCouplingInterface(self, interface_config):
        model_part_name = interface_config["model_part_name"]
        comm_name = interface_config["comm_name"]

        if not self.model.HasModelPart(model_part_name):
            main_model_part_name, *sub_model_part_names = model_part_name.split(".")
            model_part_utilities.RecursiveCreateModelParts(self.model[main_model_part_name], ".".join(sub_model_part_names))

        model_part = self.model[model_part_name]
        KratosCoSim.EMPIRE_API.EMPIRE_API_recvMesh(model_part, comm_name)

    def ExportCouplingInterface(self, interface_config):
        model_part_name = interface_config["model_part_name"]
        comm_name = interface_config["comm_name"]
        model_part_to_export = self.model[model_part_name]

        if model_part_to_export.IsDistributed():
            gathered_model_part = self.aux_model.CreateModelPart(model_part_to_export.FullName()+"_gather_mesh_export")
            KratosMPI.GatherModelPartUtility(0, model_part_to_export, 0, gathered_model_part)
            model_part_for_api = gathered_model_part
        else:
            model_part_for_api = model_part_to_export

        if self.rank == 0:
            KratosCoSim.EMPIRE_API.EMPIRE_API_sendMesh(model_part_for_api, comm_name)

        # cleanup, the gathered ModelPart is no longer needed
        if model_part_to_export.IsDistributed():
            self.aux_model.DeleteModelPart(model_part_to_export.FullName()+"_gather_mesh_export")

    def ImportData(self, data_config):
        data_type = data_config["type"]
        if data_type == "coupling_interface_data":
            interface_data = data_config["interface_data"]
            KratosCoSim.EMPIRE_API.EMPIRE_API_recvDataField(interface_data.GetModelPart(), self.solver_name+"_"+interface_data.name, interface_data.variable)
        else:
            raise NotImplementedError('Importing interface data of type "{}" is not implemented for this IO: "{}"'.format(data_type, self._ClassName()))

    def ExportData(self, data_config):
        data_type = data_config["type"]
        if data_type == "coupling_interface_data":
            interface_data = data_config["interface_data"]
            KratosCoSim.EMPIRE_API.EMPIRE_API_sendDataField(interface_data.GetModelPart(), self.solver_name+"_"+interface_data.name, interface_data.variable)
        elif data_type == "convergence_signal":
            KratosCoSim.EMPIRE_API.EMPIRE_API_sendConvergenceSignal(data_config["is_converged"], self.solver_name)
        else:
            raise NotImplementedError('Exporting interface data of type "{}" is not implemented for this IO: "{}"'.format(data_type, self._ClassName()))

    def PrintInfo(self):
        print("This is the EMPIRE-IO")

    def Check(self):
        pass

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "api_print_timing" : false
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())

        return this_defaults
