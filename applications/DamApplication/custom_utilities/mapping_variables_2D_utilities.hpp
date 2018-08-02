
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Joaquin Irazabal
//


#if !defined(KRATOS_MAPPING_VARIABLES_2D_UTILITIES )
#define  KRATOS_MAPPING_VARIABLES_2D_UTILITIES

// System includes
#include <fstream>
#include <iostream>
#include <cmath>

// Project includes
#include "geometries/geometry.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "dam_application_variables.h"

namespace Kratos
{

class MappingVariables2DUtilities
{

protected:

/// Basic Structs for the utility ---------------------------------------------------------------------------------------------------------------------------------------------

    struct UtilityVariables
    {
        double X_max, X_min, Y_max, Y_min;
        int NRows, NColumns;
        double RowSize, ColumnSize;
    };


public:

    KRATOS_CLASS_POINTER_DEFINITION( MappingVariables2DUtilities );

    /// Constructor
    MappingVariables2DUtilities() {}

    ///------------------------------------------------------------------------------------

    /// Destructor
    virtual ~MappingVariables2DUtilities() {}

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void MappingModelParts (Parameters& rParameters, ModelPart& rModelPartOld, ModelPart& rModelPartNew, const bool& move_mesh_flag)
    {
        // Define necessary variables
        UtilityVariables AuxVariables;

        this->InitializeMapping(AuxVariables,rModelPartNew, move_mesh_flag);

        this->NodalVariablesMapping(AuxVariables,rParameters,rModelPartOld,rModelPartNew);

        if(move_mesh_flag==true)
            this->UpdateMeshPosition(rModelPartNew);
    }

protected:

    /// Member Variables

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/// Mapping -------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void InitializeMapping(
        UtilityVariables& rAuxVariables,
        ModelPart& rModelPartNew,
        const bool& move_mesh_flag)
    {
        // Move mesh to the original position to work in the reference state
        if(move_mesh_flag==true)
            this->ResetMeshPosition(rModelPartNew);

        this->ComputeCellMatrixDimensions(rAuxVariables,rModelPartNew);
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void NodalVariablesMapping(
        const UtilityVariables& AuxVariables,
        Parameters& rParameters,
        ModelPart& rModelPartOld,
        ModelPart& rModelPartNew)
    {
        // Define ElementOld Cell matrix
        std::vector< std::vector< std::vector<Element::Pointer> > > ElementOldCellMatrix;
        ElementOldCellMatrix.resize(AuxVariables.NRows);
        for(int i = 0; i < AuxVariables.NRows; i++) ElementOldCellMatrix[i].resize(AuxVariables.NColumns);

        // Locate Old Elments in cells
        double X_me;
        double Y_me;
        int PointsNumber;

        unsigned int NumBodySubModelParts = rParameters["fracture_data"]["body_domain_sub_model_part_list"].size();

        // Loop through all BodySubModelParts
        for(unsigned int m = 0; m < NumBodySubModelParts; m++)
        {
            ModelPart& SubModelPart = rModelPartOld.GetSubModelPart(rParameters["fracture_data"]["body_domain_sub_model_part_list"][m].GetString());

            int NElems = static_cast<int>(SubModelPart.Elements().size());
            ModelPart::ElementsContainerType::iterator el_begin = SubModelPart.ElementsBegin();

            #pragma omp parallel for private(X_me,Y_me,PointsNumber)
            for(int i = 0; i < NElems; i++)
            {
                ModelPart::ElementsContainerType::iterator itElemOld = el_begin + i;

                double X_left = itElemOld->GetGeometry().GetPoint(0).X0();
                double X_right = X_left;
                double Y_top = itElemOld->GetGeometry().GetPoint(0).Y0();
                double Y_bot = Y_top;
                PointsNumber = itElemOld->GetGeometry().PointsNumber();

                for(int j = 1; j < PointsNumber; j++)
                {
                    X_me = itElemOld->GetGeometry().GetPoint(j).X0();
                    Y_me = itElemOld->GetGeometry().GetPoint(j).Y0();

                    if(X_me > X_right) X_right = X_me;
                    else if(X_me < X_left) X_left = X_me;
                    if(Y_me > Y_top) Y_top = Y_me;
                    else if(Y_me < Y_bot) Y_bot = Y_me;
                }

                int Column_left = int((X_left-AuxVariables.X_min)/AuxVariables.ColumnSize);
                int Column_right = int((X_right-AuxVariables.X_min)/AuxVariables.ColumnSize);
                int Row_top = int((AuxVariables.Y_max-Y_top)/AuxVariables.RowSize);
                int Row_bot = int((AuxVariables.Y_max-Y_bot)/AuxVariables.RowSize);

                if(Column_left < 0) Column_left = 0;
                else if(Column_left >= AuxVariables.NColumns) Column_left = AuxVariables.NColumns-1;
                if(Column_right >= AuxVariables.NColumns) Column_right = AuxVariables.NColumns-1;
                else if(Column_right < 0) Column_right = 0;

                if(Row_top < 0) Row_top = 0;
                else if(Row_top >= AuxVariables.NRows) Row_top = AuxVariables.NRows-1;
                if(Row_bot >= AuxVariables.NRows) Row_bot = AuxVariables.NRows-1;
                else if(Row_bot < 0) Row_bot = 0;

                for(int k = Row_top; k <= Row_bot; k++)
                {
                    for(int l = Column_left; l <= Column_right; l++)
                    {
                        #pragma omp critical
                        {
                            ElementOldCellMatrix[k][l].push_back((*(itElemOld.base())));
                        }
                    }
                }
            }
        }

        unsigned int NumInterfaceSubModelPartsOld = rParameters["fracture_data"]["interface_domain_sub_model_part_old_list"].size();

        // Loop through all InterfaceSubModelParts
        for(unsigned int m = 0; m < NumInterfaceSubModelPartsOld; m++)
        {
            ModelPart& SubModelPart = rModelPartOld.GetSubModelPart(rParameters["fracture_data"]["interface_domain_sub_model_part_old_list"][m].GetString());

            int NElems = static_cast<int>(SubModelPart.Elements().size());
            ModelPart::ElementsContainerType::iterator el_begin = SubModelPart.ElementsBegin();

            #pragma omp parallel for private(X_me,Y_me,PointsNumber)
            for(int i = 0; i < NElems; i++)
            {
                ModelPart::ElementsContainerType::iterator itElemOld = el_begin + i;

                double X_left = itElemOld->GetGeometry().GetPoint(0).X0();
                double X_right = X_left;
                double Y_top = itElemOld->GetGeometry().GetPoint(0).Y0();
                double Y_bot = Y_top;
                PointsNumber = itElemOld->GetGeometry().PointsNumber();

                for(int j = 1; j < PointsNumber; j++)
                {
                    X_me = itElemOld->GetGeometry().GetPoint(j).X0();
                    Y_me = itElemOld->GetGeometry().GetPoint(j).Y0();

                    if(X_me > X_right) X_right = X_me;
                    else if(X_me < X_left) X_left = X_me;
                    if(Y_me > Y_top) Y_top = Y_me;
                    else if(Y_me < Y_bot) Y_bot = Y_me;
                }

                int Column_left = int((X_left-AuxVariables.X_min)/AuxVariables.ColumnSize);
                int Column_right = int((X_right-AuxVariables.X_min)/AuxVariables.ColumnSize);
                int Row_top = int((AuxVariables.Y_max-Y_top)/AuxVariables.RowSize);
                int Row_bot = int((AuxVariables.Y_max-Y_bot)/AuxVariables.RowSize);

                if(Column_left < 0) Column_left = 0;
                else if(Column_left >= AuxVariables.NColumns) Column_left = AuxVariables.NColumns-1;
                if(Column_right >= AuxVariables.NColumns) Column_right = AuxVariables.NColumns-1;
                else if(Column_right < 0) Column_right = 0;

                if(Row_top < 0) Row_top = 0;
                else if(Row_top >= AuxVariables.NRows) Row_top = AuxVariables.NRows-1;
                if(Row_bot >= AuxVariables.NRows) Row_bot = AuxVariables.NRows-1;
                else if(Row_bot < 0) Row_bot = 0;

                for(int k = Row_top; k <= Row_bot; k++)
                {
                    for(int l = Column_left; l <= Column_right; l++)
                    {
                        #pragma omp critical
                        {
                            ElementOldCellMatrix[k][l].push_back((*(itElemOld.base())));
                        }
                    }
                }
            }
        }

        // Locate new nodes inside old elements and interpolate nodal variables
        const int NNodes = static_cast<int>(rModelPartNew.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = rModelPartNew.NodesBegin();

        array_1d<double,3> GlobalCoordinates;
        array_1d<double,3> LocalCoordinates;

        #pragma omp parallel for private(X_me,Y_me,PointsNumber,GlobalCoordinates,LocalCoordinates)
        for(int i = 0; i < NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator itNodeNew = node_begin + i;

            X_me = itNodeNew->X0();
            Y_me = itNodeNew->Y0();

            int Column = int((X_me-AuxVariables.X_min)/AuxVariables.ColumnSize);
            int Row = int((AuxVariables.Y_max-Y_me)/AuxVariables.RowSize);

            if(Column >= AuxVariables.NColumns) Column = AuxVariables.NColumns-1;
            else if(Column < 0) Column = 0;
            if(Row >= AuxVariables.NRows) Row = AuxVariables.NRows-1;
            else if(Row < 0) Row = 0;

            noalias(GlobalCoordinates) = itNodeNew->Coordinates(); //Coordinates of new nodes are still in the original position
            bool IsInside = false;
            Element::Pointer pElementOld;

            for(unsigned int m = 0; m < (ElementOldCellMatrix[Row][Column]).size(); m++)
            {
                pElementOld = ElementOldCellMatrix[Row][Column][m];
                IsInside = pElementOld->GetGeometry().IsInside(GlobalCoordinates,LocalCoordinates); //Checks whether the global coordinates fall inside the original old element
                if(IsInside) break;
            }
            if(IsInside==false)
            {
                for(unsigned int m = 0; m < (ElementOldCellMatrix[Row][Column]).size(); m++)
                {
                    pElementOld = ElementOldCellMatrix[Row][Column][m];
                    IsInside = pElementOld->GetGeometry().IsInside(GlobalCoordinates,LocalCoordinates,1.0e-5);
                    if(IsInside) break;
                }
            }
            if(IsInside == false)
                std::cout << "ERROR!!, NONE OF THE OLD ELEMENTS CONTAINS NODE: " << itNodeNew->Id() << std::endl;

            PointsNumber = pElementOld->GetGeometry().PointsNumber();
            Vector ShapeFunctionsValuesVector(PointsNumber);
            Vector NodalVariableVector(PointsNumber);

            pElementOld->GetGeometry().ShapeFunctionsValues(ShapeFunctionsValuesVector,LocalCoordinates);

            // Interpolation of nodal variables
            if( itNodeNew->IsFixed(DISPLACEMENT_X)==false )
            {
                for(int j = 0; j < PointsNumber; j++)
                {
                    NodalVariableVector[j] = pElementOld->GetGeometry().GetPoint(j).FastGetSolutionStepValue(DISPLACEMENT_X);
                }
                itNodeNew->FastGetSolutionStepValue(DISPLACEMENT_X) = inner_prod(ShapeFunctionsValuesVector,NodalVariableVector);
            }
            if( itNodeNew->IsFixed(VELOCITY_X)==false )
            {
                for(int j = 0; j < PointsNumber; j++)
                {
                    NodalVariableVector[j] = pElementOld->GetGeometry().GetPoint(j).FastGetSolutionStepValue(VELOCITY_X);
                }
                itNodeNew->FastGetSolutionStepValue(VELOCITY_X) = inner_prod(ShapeFunctionsValuesVector,NodalVariableVector);
            }
            if( itNodeNew->IsFixed(ACCELERATION_X)==false )
            {
                for(int j = 0; j < PointsNumber; j++)
                {
                    NodalVariableVector[j] = pElementOld->GetGeometry().GetPoint(j).FastGetSolutionStepValue(ACCELERATION_X);
                }
                itNodeNew->FastGetSolutionStepValue(ACCELERATION_X) = inner_prod(ShapeFunctionsValuesVector,NodalVariableVector);
            }
            if( itNodeNew->IsFixed(DISPLACEMENT_Y)==false )
            {
                for(int j = 0; j < PointsNumber; j++)
                {
                    NodalVariableVector[j] = pElementOld->GetGeometry().GetPoint(j).FastGetSolutionStepValue(DISPLACEMENT_Y);
                }
                itNodeNew->FastGetSolutionStepValue(DISPLACEMENT_Y) = inner_prod(ShapeFunctionsValuesVector,NodalVariableVector);
            }
            if( itNodeNew->IsFixed(VELOCITY_Y)==false )
            {
                for(int j = 0; j < PointsNumber; j++)
                {
                    NodalVariableVector[j] = pElementOld->GetGeometry().GetPoint(j).FastGetSolutionStepValue(VELOCITY_Y);
                }
                itNodeNew->FastGetSolutionStepValue(VELOCITY_Y) = inner_prod(ShapeFunctionsValuesVector,NodalVariableVector);
            }
            if( itNodeNew->IsFixed(ACCELERATION_Y)==false )
            {
                for(int j = 0; j < PointsNumber; j++)
                {
                    NodalVariableVector[j] = pElementOld->GetGeometry().GetPoint(j).FastGetSolutionStepValue(ACCELERATION_Y);
                }
                itNodeNew->FastGetSolutionStepValue(ACCELERATION_Y) = inner_prod(ShapeFunctionsValuesVector,NodalVariableVector);
            }
            if( itNodeNew->IsFixed(WATER_PRESSURE)==false )
            {
                for(int j = 0; j < PointsNumber; j++)
                {
                    NodalVariableVector[j] = pElementOld->GetGeometry().GetPoint(j).FastGetSolutionStepValue(WATER_PRESSURE);
                }
                itNodeNew->FastGetSolutionStepValue(WATER_PRESSURE) = inner_prod(ShapeFunctionsValuesVector,NodalVariableVector);
                for(int j = 0; j < PointsNumber; j++)
                {
                    NodalVariableVector[j] = pElementOld->GetGeometry().GetPoint(j).FastGetSolutionStepValue(DT_WATER_PRESSURE);
                }
                itNodeNew->FastGetSolutionStepValue(DT_WATER_PRESSURE) = inner_prod(ShapeFunctionsValuesVector,NodalVariableVector);
            }
        }
    }

/// Common --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ResetMeshPosition(ModelPart& rModelPart)
    {
        // Move mesh to the original position
        const int NNodes = static_cast<int>(rModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = rModelPart.NodesBegin();

        #pragma omp parallel for
        for(int i = 0; i < NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator itNode = node_begin + i;

            itNode->X() = itNode->X0();
            itNode->Y() = itNode->Y0();
            itNode->Z() = itNode->Z0();
        }
    }

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void UpdateMeshPosition(ModelPart& rModelPart)
    {
        // Move mesh to the current position
        const int NNodes = static_cast<int>(rModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = rModelPart.NodesBegin();

        #pragma omp parallel for
        for(int i = 0; i < NNodes; i++)
        {
            ModelPart::NodesContainerType::iterator itNode = node_begin + i;

            itNode->X() = itNode->X0() + itNode->FastGetSolutionStepValue(DISPLACEMENT_X);
            itNode->Y() = itNode->Y0() + itNode->FastGetSolutionStepValue(DISPLACEMENT_Y);
            itNode->Z() = itNode->Z0() + itNode->FastGetSolutionStepValue(DISPLACEMENT_Z);
        }
    }

private:

/// Common --------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    void ComputeCellMatrixDimensions(
        UtilityVariables& rAuxVariables,
        ModelPart& rModelPart)
    {
        // Compute X and Y limits of the current geometry
        unsigned int NumThreads = OpenMPUtils::GetNumThreads();
        std::vector<double> X_max_partition(NumThreads);
        std::vector<double> X_min_partition(NumThreads);
        std::vector<double> Y_max_partition(NumThreads);
        std::vector<double> Y_min_partition(NumThreads);

        const int NNodes = static_cast<int>(rModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = rModelPart.NodesBegin();

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            X_max_partition[k] = node_begin->X();
            X_min_partition[k] = X_max_partition[k];
            Y_max_partition[k] = node_begin->Y();
            Y_min_partition[k] = Y_max_partition[k];

            double X_me, Y_me;

            #pragma omp for
            for(int i = 0; i < NNodes; i++)
            {
                ModelPart::NodesContainerType::iterator itNode = node_begin + i;

                X_me = itNode->X();
                Y_me = itNode->Y();

                if( X_me > X_max_partition[k] ) X_max_partition[k] = X_me;
                else if( X_me < X_min_partition[k] ) X_min_partition[k] = X_me;

                if( Y_me > Y_max_partition[k] ) Y_max_partition[k] = Y_me;
                else if( Y_me < Y_min_partition[k] ) Y_min_partition[k] = Y_me;
            }
        }

        rAuxVariables.X_max = X_max_partition[0];
        rAuxVariables.X_min = X_min_partition[0];
        rAuxVariables.Y_max = Y_max_partition[0];
        rAuxVariables.Y_min = Y_min_partition[0];

        for(unsigned int i=1; i < NumThreads; i++)
        {
            if(X_max_partition[i] > rAuxVariables.X_max) rAuxVariables.X_max = X_max_partition[i];
            if(X_min_partition[i] < rAuxVariables.X_min) rAuxVariables.X_min = X_min_partition[i];
            if(Y_max_partition[i] > rAuxVariables.Y_max) rAuxVariables.Y_max = Y_max_partition[i];
            if(Y_min_partition[i] < rAuxVariables.Y_min) rAuxVariables.Y_min = Y_min_partition[i];
        }

        // Calculate Average Element Length
        double AverageElementLength = 0.0;

        int NElems = static_cast<int>(rModelPart.Elements().size());
        ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();

        #pragma omp parallel for reduction(+:AverageElementLength)
        for(int i = 0; i < NElems; i++)
        {
            ModelPart::ElementsContainerType::iterator itElem = el_begin + i;

            AverageElementLength += itElem->GetGeometry().Length();
        }

        AverageElementLength = AverageElementLength/NElems;

        // Compute FracturePoints CellMatrix dimensions
        rAuxVariables.NRows = int((rAuxVariables.Y_max-rAuxVariables.Y_min)/(AverageElementLength));
        rAuxVariables.NColumns = int((rAuxVariables.X_max-rAuxVariables.X_min)/(AverageElementLength));
        if(rAuxVariables.NRows <= 0) rAuxVariables.NRows = 1;
        if(rAuxVariables.NColumns <= 0) rAuxVariables.NColumns = 1;
        rAuxVariables.RowSize = (rAuxVariables.Y_max-rAuxVariables.Y_min)/rAuxVariables.NRows;
        rAuxVariables.ColumnSize = (rAuxVariables.X_max-rAuxVariables.X_min)/rAuxVariables.NColumns;
    }

}; // Class MappingVariables2DUtilities

} // namespace Kratos.

#endif /* KRATOS_MAPPING_VARIABLES_2D_UTILITIES defined */
