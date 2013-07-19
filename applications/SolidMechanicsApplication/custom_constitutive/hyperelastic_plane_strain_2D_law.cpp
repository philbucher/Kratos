//   
//   Project Name:        KratosSolidMechanicsApplication $      
//   Last modified by:    $Author:            JMCarbonell $ 
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes
#include <iostream>

// External includes
#include<cmath>

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/hyperelastic_plane_strain_2D_law.hpp"

#include "solid_mechanics_application.h"

namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  HyperElasticPlaneStrain2DLaw::HyperElasticPlaneStrain2DLaw()
  : HyperElastic3DLaw()
  {
  }

  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  HyperElasticPlaneStrain2DLaw::HyperElasticPlaneStrain2DLaw(const HyperElasticPlaneStrain2DLaw& rOther)
  : HyperElastic3DLaw()
  {
  }
  
  //********************************CLONE***********************************************
  //************************************************************************************

  ConstitutiveLaw::Pointer HyperElasticPlaneStrain2DLaw::Clone() const
  {
    HyperElasticPlaneStrain2DLaw::Pointer p_clone(new HyperElasticPlaneStrain2DLaw(*this));
    return p_clone;
  }
  
  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************

  HyperElasticPlaneStrain2DLaw::~HyperElasticPlaneStrain2DLaw()
  {
  }

  

  //***********************COMPUTE TOTAL STRAIN*****************************************
  //************************************************************************************

  void HyperElasticPlaneStrain2DLaw::CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
							Vector& rStrainVector )
  {

    //E= 0.5*(FT*F-1)
    rStrainVector[0] = 0.5 * ( rRightCauchyGreen( 0, 0 ) - 1.00 );
    rStrainVector[1] = 0.5 * ( rRightCauchyGreen( 1, 1 ) - 1.00 );
    rStrainVector[2] = rRightCauchyGreen( 0, 1 );

  }



  //***********************COMPUTE TOTAL STRAIN*****************************************
  //************************************************************************************

  void HyperElasticPlaneStrain2DLaw::CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
						  Vector& rStrainVector )
  {

    // e= 0.5*(1-invbT*invb)   
    Matrix InverseLeftCauchyGreen ( 2 , 2 );
    double det_b=0;
    MathUtils<double>::InvertMatrix( rLeftCauchyGreen, InverseLeftCauchyGreen, det_b);

    rStrainVector.clear();
    rStrainVector[0] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 0, 0 ) );
    rStrainVector[1] = 0.5 * ( 1.0 - InverseLeftCauchyGreen( 1, 1 ) );
    rStrainVector[2] = InverseLeftCauchyGreen( 0, 1 );


  }



  //***********************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX**********************
  //************************************************************************************

  void HyperElasticPlaneStrain2DLaw::CalculateConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
						       Matrix& rConstitutiveMatrix)
  {
    
    rConstitutiveMatrix.clear();
		
    static const unsigned int msIndexVoigt2D [6][2] = { {0, 0}, {1, 1}, {0, 1} };		

    for(unsigned int i=0; i<3; i++)
      {
    	for(unsigned int j=0; j<3; j++)
    	  {
    	    rConstitutiveMatrix( i, j ) = ConstitutiveComponent(rConstitutiveMatrix( i, j ), rElasticVariables,
    								msIndexVoigt2D[i][0], msIndexVoigt2D[i][1], msIndexVoigt2D[j][0], msIndexVoigt2D[j][1]);
    	  }

      }

	  	
  }


  //**************COMPUTE ALGORITHMIC CONSTITUTIVE MATRIX PULL-BACK*********************
  //************************************************************************************

  void HyperElasticPlaneStrain2DLaw::CalculateConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
						       const Matrix & rInverseDeformationGradientF,
						       Matrix& rConstitutiveMatrix)
  {

    rConstitutiveMatrix.clear();

    static const unsigned int msIndexVoigt2D [6][2] = { {0, 0}, {1, 1}, {0, 1} };

    for(unsigned int i=0; i<3; i++)
      {
	for(unsigned int j=0; j<3; j++)
	  {
	    rConstitutiveMatrix( i, j ) = ConstitutiveComponent(rConstitutiveMatrix( i, j ), rElasticVariables, rInverseDeformationGradientF,
								msIndexVoigt2D[i][0], msIndexVoigt2D[i][1], msIndexVoigt2D[j][0], msIndexVoigt2D[j][1]);
	  }

      }

    	  
  }



} // Namespace Kratos
