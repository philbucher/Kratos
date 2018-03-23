from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import contact_structural_mechanics_analysis

# Other imports
import os

# This utility will control the execution scope in case we need to access files or we depend
# on specific relative locations of the files.

# TODO: Should we move this to KratosUnittest?
class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)


class ContactStructuralMechanicsTestFactory(KratosUnittest.TestCase):
    def setUp(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):

            # Reading the ProjectParameters
            with open(self.file_name + "_parameters.json",'r') as parameter_file:
                ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

        # Checking if frictionless_by_components is defined
        try:
            self.frictionless_by_components
        except AttributeError:
            self.frictionless_by_components = False

        # Setting for frictionless by components
        if (self.frictionless_by_components is True):
            if (ProjectParameters.Has("output_configuration") is True):
                list_nodal_var = ProjectParameters["output_configuration"]["result_file_configuration"]["nodal_results"]
                for i in range(0, list_nodal_var.size()):
                    if (list_nodal_var[i].GetString() == "NORMAL_CONTACT_STRESS"):
                        list_nodal_var[i].SetString("VECTOR_LAGRANGE_MULTIPLIER")
                new_list = list_nodal_var.Clone()
                ProjectParameters["output_configuration"]["result_file_configuration"].RemoveValue("nodal_results")
                ProjectParameters["output_configuration"]["result_file_configuration"].AddValue("nodal_results", new_list)
            ProjectParameters["solver_settings"]["contact_settings"]["mortar_type"].SetString("ALMContactFrictionlessComponents")
            for i in range(ProjectParameters["contact_process_list"].size()):
                ProjectParameters["contact_process_list"][i]["Parameters"]["contact_type"].SetString("FrictionlessComponents")

        # To avoid many prints
        echo_level = ProjectParameters["problem_data"]["echo_level"].GetInt()
        if (echo_level == 0):
            KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

        # Creating the test
        self.test = contact_structural_mechanics_analysis.ContactStructuralMechanicsAnalysis(ProjectParameters)
        self.test.Initialize()

    def test_execution(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.RunMainTemporalLoop()

    def tearDown(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.Finalize()