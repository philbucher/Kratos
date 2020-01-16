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

#if !defined(KRATOS_BEAM_ELEMENT_LINEAR_3D2N_H_INCLUDED )
#define  KRATOS_BEAM_ELEMENT_LINEAR_3D2N_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/beam_element_3D2N.hpp"
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/serializer.h"

namespace Kratos
{
/**
 * @class BeamElementLinear3D2N
 *
 * @brief This is a 3D-2node beam element with 3 translational dofs and 3 rotational dof per node
 *
 * @author Klaus B Sautter
 */

class BeamElementLinear3D2N : public BeamElement3D2N
{

protected:
    //const values
    static constexpr int msNumberOfNodes = 2;
    static constexpr int msDimension = 3;
    static constexpr unsigned int msLocalSize = msNumberOfNodes * msDimension;
    static constexpr unsigned int msElementSize = msLocalSize * 2;

public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(BeamElementLinear3D2N);


    typedef Element BaseType;
    typedef BaseType::GeometryType GeometryType;
    typedef BaseType::NodesArrayType NodesArrayType;
    typedef BaseType::PropertiesType PropertiesType;
    typedef BaseType::IndexType IndexType;
    typedef BaseType::SizeType SizeType;
    typedef BaseType::MatrixType MatrixType;
    typedef BaseType::VectorType VectorType;
    typedef BaseType::EquationIdVectorType EquationIdVectorType;
    typedef BaseType::DofsVectorType DofsVectorType;

    BeamElementLinear3D2N() {};
    BeamElementLinear3D2N(IndexType NewId, GeometryType::Pointer pGeometry);
    BeamElementLinear3D2N(IndexType NewId, GeometryType::Pointer pGeometry,
                      PropertiesType::Pointer pProperties);


    ~BeamElementLinear3D2N() override;

    /**
    * @brief Creates a new element
    * @param NewId The Id of the new created element
    * @param pGeom The pointer to the geometry of the element
    * @param pProperties The pointer to property
    * @return The pointer to the created element
    */
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const override;

    /**
    * @brief Creates a new element
    * @param NewId The Id of the new created element
    * @param ThisNodes The array containing nodes
    * @param pProperties The pointer to property
    * @return The pointer to the created element
    */
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override;

    /**
     * @brief This function calculates the elastic part of the total stiffness matrix
     */
    Matrix CreateElementStiffnessMatrixIntermediate() const override;
    Matrix GlobalTangentStiffnessMatrix() const override;

    /**
     * @brief This function calculates the current nodal position
     */
    BoundedVector<double,msLocalSize> GetCurrentNodalPosition() const override;


    Matrix CoRotatingCS() const override;


    void CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo) override;


    Vector CalculateGlobalNodalForces() const override;

    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;


private:


    Vector mDeformationCurrentIteration = ZeroVector(msElementSize);
    Vector mDeformationPreviousIteration = ZeroVector(msElementSize);
    Matrix mGlobalRotationNode1 = IdentityMatrix(msDimension);
    Matrix mGlobalRotationNode2 = IdentityMatrix(msDimension);
};


}

#endif
