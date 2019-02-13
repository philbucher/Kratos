//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Veronika Singer
//


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/particle_based_conditions/mpm_particle_point_load_condition.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
    //******************************* CONSTRUCTOR ****************************************
    //************************************************************************************

    MPMParticlePointLoadCondition::MPMParticlePointLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry )
        : MPMParticleBaseLoadCondition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    MPMParticlePointLoadCondition::MPMParticlePointLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : MPMParticleBaseLoadCondition( NewId, pGeometry, pProperties )
    {
    }

    //********************************* CREATE *******************************************
    //************************************************************************************

    Condition::Pointer MPMParticlePointLoadCondition::Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const
    {
        return Kratos::make_shared<MPMParticlePointLoadCondition>(NewId, pGeom, pProperties);
    }

    //************************************************************************************
    //************************************************************************************

    Condition::Pointer MPMParticlePointLoadCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
    {
        return Kratos::make_shared<MPMParticlePointLoadCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
    }

    //******************************* DESTRUCTOR *****************************************
    //************************************************************************************

    MPMParticlePointLoadCondition::~MPMParticlePointLoadCondition()
    {
    }

    void MPMParticlePointLoadCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag
        )
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int block_size = this->GetBlockSize();

        // Resizing as needed the LHS
        const unsigned int matrix_size = number_of_nodes * block_size;

        if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
        {
            if ( rLeftHandSideMatrix.size1() != matrix_size )
            {
                rLeftHandSideMatrix.resize( matrix_size, matrix_size, false );
            }

            noalias( rLeftHandSideMatrix ) = ZeroMatrix(matrix_size,matrix_size); //resetting LHS
        }

        // Resizing as needed the RHS
        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
        {
            if ( rRightHandSideVector.size( ) != matrix_size )
            {
                rRightHandSideVector.resize( matrix_size, false );
            }

            noalias( rRightHandSideVector ) = ZeroVector( matrix_size ); //resetting RHS
        }

        // Get imposed displacement and normal vector
        const array_1d<double, 3 > & xg_c = this->GetValue(MPC_COORD);
        const array_1d<double, 3 > & Point_Load = this->GetValue (POINT_LOAD);
        Matrix Nodal_Force = ZeroMatrix(3,number_of_nodes);


        // Prepare variables
        GeneralVariables Variables;

        // Calculating shape function
        Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg_c);

        // Here MP contribution in terms of force are added
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            for (unsigned int j = 0; j < dimension; j++)
            {
                Nodal_Force(j,i) = Variables.N[i] * Point_Load[j];
            }
        }

        for (unsigned int ii = 0; ii < number_of_nodes; ++ii)
        {
            const unsigned int base = ii*dimension;

            for(unsigned int k = 0; k < dimension; ++k)
            {
                rRightHandSideVector[base + k] += GetPointLoadIntegrationWeight() * Nodal_Force(k,ii);
            }
        }
        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    double MPMParticlePointLoadCondition::GetPointLoadIntegrationWeight()
    {
        return 1.0;
    }

    void MPMParticlePointLoadCondition::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
    {
        GeneralVariables Variables;

        Variables.CurrentDisp = CalculateCurrentDisp(Variables.CurrentDisp, rCurrentProcessInfo);

        const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        const array_1d<double,3> & xg = this->GetValue(MPC_COORD);

        array_1d<double,3> delta_xg = ZeroVector(3);

        Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg);

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            if (Variables.N[i] > 1e-16)
            {
                for ( unsigned int j = 0; j < dimension; j++ )
                {
                    delta_xg[j] += Variables.N[i] * Variables.CurrentDisp(i,j);
                }
            }
        }

        // Update the Material Point Condition Position
        const array_1d<double,3>& new_xg = xg + delta_xg ;
        this -> SetValue(MPC_COORD,new_xg);

        //mFinalizedStep = true;

    }

} // Namespace Kratos