# Importing the Kratos Library
import KratosMultiphysics as KM
default_data_comm = KM.DataCommunicator.GetDefault()
from KratosMultiphysics.mpi import GatherModelPartUtility

# Mapping imports
from KratosMultiphysics.MappingApplication.python_mapper import PythonMapper
from KratosMultiphysics.MappingApplication.python_mapper_factory import CreateMapper

def Create(model_part_origin, model_part_destination, mapper_settings):
    return DistributedMapperWrapper(model_part_origin, model_part_destination, mapper_settings)


class DistributedMapperWrapper(PythonMapper):

    def __init__(self, model_part_origin, model_part_destination, mapper_settings):
        super().__init__(model_part_origin, model_part_destination, mapper_settings)

        self.aux_model = KM.Model()
        self.gathered_model_part_origin = self.aux_model.CreateModelPart(model_part_origin.Name+"_origin_gathered")
        self.gathered_model_part_dest = self.aux_model.CreateModelPart(model_part_destination.Name+"_destination_gathered")

        self.gather_util_origin = GatherModelPartUtility(0, model_part_origin, 0, self.gathered_model_part_origin)
        self.gather_util_dest = GatherModelPartUtility(0, model_part_destination, 0, self.gathered_model_part_dest)


        if default_data_comm.Rank() == 0:
            self._mapper = CreateMapper(self.gathered_model_part_origin, self.gathered_model_part_dest, mapper_settings["inner_mapper_settings"])
            stop

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "inner_mapper_settings" : { }
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults

    def UpdateInterface(self):
        raise NotImplementedError


    def _MapInternal(self, variable_origin, variable_destination, mapper_flags):
        self.gather_util_origin.GatherOnMaster(variable_origin)
        if default_data_comm.Rank() == 0:
            self._mapper.Map(variable_origin, variable_destination, mapper_flags)
        self.gather_util_dest.ScatterFromMaster(variable_destination)

    def _InverseMapInternal(self, variable_origin, variable_destination, mapper_flags):
        self.gather_util_dest.GatherOnMaster(variable_destination)
        if default_data_comm.Rank() == 0:
            self._mapper.InverseMap(variable_origin, variable_destination, mapper_flags)
        self.gather_util_origin.ScatterFromMaster(variable_origin)