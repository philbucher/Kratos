// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors: Klaus B. Sautter
//
//
//
// System includes

// External includes

// Project includes
#include "custom_elements/beam_element_linear_3D2N.hpp"
#include "includes/define.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos
{
BeamElementLinear3D2N::BeamElementLinear3D2N(IndexType NewId,
                                     GeometryType::Pointer pGeometry)
    : BeamElement3D2N(NewId, pGeometry) {}

BeamElementLinear3D2N::BeamElementLinear3D2N(IndexType NewId,
                                     GeometryType::Pointer pGeometry,
                                     PropertiesType::Pointer pProperties)
    : BeamElement3D2N(NewId, pGeometry, pProperties) {}

Element::Pointer
BeamElementLinear3D2N::Create(IndexType NewId, NodesArrayType const& rThisNodes,
                          PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = GetGeometry();
    return Kratos::make_intrusive<BeamElementLinear3D2N>(NewId, rGeom.Create(rThisNodes),
            pProperties);
}

Element::Pointer
BeamElementLinear3D2N::Create(IndexType NewId, GeometryType::Pointer pGeom,
                          PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<BeamElementLinear3D2N>(NewId, pGeom,
            pProperties);
}

BeamElementLinear3D2N::~BeamElementLinear3D2N() {}


void BeamElementLinear3D2N::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    KRATOS_CATCH("")
}


Matrix BeamElementLinear3D2N::CreateElementStiffnessMatrixIntermediate() const
{
    const double numerical_limit = std::numeric_limits<double>::epsilon();
    Matrix b_a = BMatrixIntermediateA();
    Matrix local_stiffness_matrix = prod(CreateElementStiffnessMatrix_Material(),b_a);
    local_stiffness_matrix = prod(trans(b_a),local_stiffness_matrix);
    return local_stiffness_matrix;
}

Matrix BeamElementLinear3D2N::GlobalTangentStiffnessMatrix() const
{
    Matrix b_g = BMatrixGlobal();
    Matrix tangent_stiffness_matrix = prod(CreateElementStiffnessMatrixIntermediate(),b_g);
    tangent_stiffness_matrix = prod(trans(b_g),tangent_stiffness_matrix);
    return tangent_stiffness_matrix;
}


Matrix BeamElementLinear3D2N::CoRotatingCS() const
{
    return CalculateInitialLocalCS();
}


Vector BeamElementLinear3D2N::CalculateGlobalNodalForces() const
{
    Vector global_deformation = ZeroVector(msElementSize);

    for (int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msDimension * 2;
        const auto& disp =
            GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, 0);
        const auto& rot =
            GetGeometry()[i].FastGetSolutionStepValue(ROTATION, 0);

        global_deformation[index] = disp[0];
        global_deformation[index + 1] = disp[1];
        global_deformation[index + 2] = disp[2];

        global_deformation[index + 3] = rot[0];
        global_deformation[index + 4] = rot[1];
        global_deformation[index + 5] = rot[2];
    }

    return prod(GlobalTangentStiffnessMatrix(),global_deformation);
}


BoundedVector<double, BeamElementLinear3D2N::msLocalSize>
BeamElementLinear3D2N::GetCurrentNodalPosition() const
{
    BoundedVector<double, msLocalSize> current_nodal_position =
        ZeroVector(msLocalSize);
    for (unsigned int i = 0; i < msNumberOfNodes; ++i) {
        int index = i * msDimension;
        current_nodal_position[index] = GetGeometry()[i].X0();
        current_nodal_position[index + 1] = GetGeometry()[i].Y0();
        current_nodal_position[index + 2] = GetGeometry()[i].Z0();
    }

    return current_nodal_position;
}


void BeamElementLinear3D2N::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY
    // element with two nodes can only represent results at one node
    const unsigned int& write_points_number =
        GetGeometry().IntegrationPointsNumber(Kratos::GeometryData::GI_GAUSS_3);
    if (rOutput.size() != write_points_number) {
        rOutput.resize(write_points_number);
    }


    // rOutput[GP 1,2,3][x,y,z]

    if (rVariable == MOMENT) {
        Vector nodal_forces_global_qe = CalculateGlobalNodalForces();
        Vector nodal_forces_local_qe = prod(trans(EMatrix()),nodal_forces_global_qe);

        rOutput[0][0] = -1.0 * nodal_forces_local_qe[3] * 0.75 + nodal_forces_local_qe[9] * 0.25;
        rOutput[1][0] = -1.0 * nodal_forces_local_qe[3] * 0.50 + nodal_forces_local_qe[9] * 0.50;
        rOutput[2][0] = -1.0 * nodal_forces_local_qe[3] * 0.25 + nodal_forces_local_qe[9] * 0.75;

        rOutput[0][1] = -1.0 * (nodal_forces_local_qe[4] - (nodal_forces_local_qe[4]+nodal_forces_local_qe[10]) * 0.25);
        rOutput[1][1] = -1.0 * (nodal_forces_local_qe[4] - (nodal_forces_local_qe[4]+nodal_forces_local_qe[10]) * 0.50);
        rOutput[2][1] = -1.0 * (nodal_forces_local_qe[4] - (nodal_forces_local_qe[4]+nodal_forces_local_qe[10]) * 0.75);

        rOutput[0][2] = (nodal_forces_local_qe[5] - (nodal_forces_local_qe[5]+nodal_forces_local_qe[11]) * 0.25);
        rOutput[1][2] = (nodal_forces_local_qe[5] - (nodal_forces_local_qe[5]+nodal_forces_local_qe[11]) * 0.50);
        rOutput[2][2] = (nodal_forces_local_qe[5] - (nodal_forces_local_qe[5]+nodal_forces_local_qe[11]) * 0.75);
    } else if (rVariable == FORCE) {
        Vector nodal_forces_global_qe = CalculateGlobalNodalForces();
        Vector nodal_forces_local_qe = prod(trans(EMatrix()),nodal_forces_global_qe);

        rOutput[0][0] = -1.0 * nodal_forces_local_qe[0] * 0.75 + nodal_forces_local_qe[6] * 0.25;
        rOutput[1][0] = -1.0 * nodal_forces_local_qe[0] * 0.50 + nodal_forces_local_qe[6] * 0.50;
        rOutput[2][0] = -1.0 * nodal_forces_local_qe[0] * 0.25 + nodal_forces_local_qe[6] * 0.75;

        rOutput[0][1] = -1.0 * nodal_forces_local_qe[1] * 0.75 + nodal_forces_local_qe[7] * 0.25;
        rOutput[1][1] = -1.0 * nodal_forces_local_qe[1] * 0.50 + nodal_forces_local_qe[7] * 0.50;
        rOutput[2][1] = -1.0 * nodal_forces_local_qe[1] * 0.25 + nodal_forces_local_qe[7] * 0.75;

        rOutput[0][2] = -1.0 * nodal_forces_local_qe[2] * 0.75 + nodal_forces_local_qe[8] * 0.25;
        rOutput[1][2] = -1.0 * nodal_forces_local_qe[2] * 0.50 + nodal_forces_local_qe[8] * 0.50;
        rOutput[2][2] = -1.0 * nodal_forces_local_qe[2] * 0.25 + nodal_forces_local_qe[8] * 0.75;
    }

    else if (rVariable == LOCAL_AXIS_1) {
        BoundedMatrix<double, msElementSize, msElementSize> rotation_matrix = CoRotatingCS();
        for (SizeType i =0; i<msDimension; ++i) {
            rOutput[1][i] = column(rotation_matrix, 0)[i];
        }
    } else if (rVariable == LOCAL_AXIS_2) {
        BoundedMatrix<double, msElementSize, msElementSize> rotation_matrix = CoRotatingCS();
        for (SizeType i =0; i<msDimension; ++i) {
            rOutput[1][i] = column(rotation_matrix, 1)[i];
        }
    } else if (rVariable == LOCAL_AXIS_3) {
        BoundedMatrix<double, msElementSize, msElementSize> rotation_matrix = CoRotatingCS();
        for (SizeType i =0; i<msDimension; ++i) {
            rOutput[1][i] = column(rotation_matrix, 2)[i];
        }
    }


    KRATOS_CATCH("")
}


} // namespace Kratos.
