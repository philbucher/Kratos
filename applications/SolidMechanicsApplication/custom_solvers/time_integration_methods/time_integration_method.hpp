//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_TIME_INTEGRATION_METHOD_H_INCLUDED)
#define  KRATOS_TIME_INTEGRATION_METHOD_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/process_info.h"
#include "includes/variables.h"
#include "includes/node.h"
#include "custom_utilities/process_info_extensions.hpp"

namespace Kratos
{
  ///@addtogroup SolidMechanicsApplication
  ///@{
  
  ///@name Kratos Globals
  ///@{
  
  ///@}
  ///@name Type Definitions
  ///@{
  
  ///@}
  ///@name  Enum's
  ///@{
  
  ///@}
  ///@name  Functions
  ///@{
  
  ///@}
  ///@name Kratos Classes
  ///@{
    
  /// Short class definition.
  /** Detail class definition.     
   * This class performs predict and update of dofs variables, their time derivatives and time integrals      
   */
  template<class TVariableType, class TValueType>
  class TimeIntegrationMethod
  {
  public:
 
    ///@name Type Definitions
    ///@{
    
    /// NodeType
    typedef Node<3> NodeType;

    /// KratosVariable or KratosVariableComponent
    typedef const TVariableType*     VariablePointer;

    typedef const TValueType*           ValuePointer;

    typedef void (TimeIntegrationMethod::*MethodPointer) (NodeType& rNode);
    
    KRATOS_CLASS_POINTER_DEFINITION( TimeIntegrationMethod );

    typedef typename TimeIntegrationMethod::Pointer   TimeIntegrationMethodPointer;

    ///@}
    ///@name Life Cycle
    ///@{

    
    /// Default Constructor.
    TimeIntegrationMethod()
    {
      mpVariable = nullptr;
      mpFirstDerivative = nullptr;
      mpSecondDerivative = nullptr;

      mpPrimaryVariable = nullptr;
      mpInputVariable = nullptr;
    }

    /// Constructor.
    TimeIntegrationMethod(const TVariableType& rVariable)
    {
      mpVariable = &rVariable;
      mpFirstDerivative = nullptr;
      mpSecondDerivative = nullptr;
      
      //default dof variable
      mpPrimaryVariable = &rVariable;

      this->SetPointerMethods();
    }

    /// Constructor.
    TimeIntegrationMethod(const TVariableType& rVariable, const TVariableType& rFirstDerivative, const TVariableType& rSecondDerivative)
    {
      mpVariable = &rVariable;
      mpFirstDerivative = &rFirstDerivative;
      mpSecondDerivative = &rSecondDerivative;

      //default dof variable
      mpPrimaryVariable = &rVariable;

      this->SetPointerMethods();
    }
    
    /// Constructor.
    TimeIntegrationMethod(const TVariableType& rVariable, const TVariableType& rFirstDerivative, const TVariableType& rSecondDerivative, const TVariableType& rPrimaryVariable)
    {
      mpVariable = &rVariable;
      mpFirstDerivative = &rFirstDerivative;
      mpSecondDerivative = &rSecondDerivative;

      if( HasVariableName(rPrimaryVariable.Name()) ){
        //default dof variable
        mpPrimaryVariable = &rPrimaryVariable;
      }
      else{
        KRATOS_ERROR << "The primary variable supplied: "<<rPrimaryVariable.Name()<<" is not any of the time integration variables" << std::endl;
      }
      
      this->SetPointerMethods();
    }
    
    /// Copy Constructor.
    TimeIntegrationMethod(TimeIntegrationMethod& rOther)
      :mpVariable(rOther.mpVariable)
      ,mpFirstDerivative(rOther.mpFirstDerivative)
      ,mpSecondDerivative(rOther.mpSecondDerivative)
      ,mpPrimaryVariable(rOther.mpPrimaryVariable)
      ,mpInputVariable(rOther.mpInputVariable)
    {
    }

    /// Clone
    virtual TimeIntegrationMethodPointer Clone()
    {
      return TimeIntegrationMethodPointer( new TimeIntegrationMethod(*this) );
    }

    /// Destructor.
    virtual ~TimeIntegrationMethod(){}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // set parameters (to call it once with the original input parameters)
    virtual void CalculateParameters(ProcessInfo& rCurrentProcessInfo)
    {
    
    }
    
    // set parameters (do not calculate parameters here, only read them)
    virtual void SetParameters(const ProcessInfo& rCurrentProcessInfo)
    {
    
    }

    virtual void SetProcessInfoParameters(ProcessInfo& rCurrentProcessInfo)
    {
    
    }
    
    // get parameters for variables (RHS)
    virtual double& GetFirstDerivativeKineticParameter(double& rParameter)
    {
      rParameter = 0.0;
      return rParameter;
    }

    virtual double& GetSecondDerivativeKineticParameter(double& rParameter)
    {
      rParameter = 0.0;
      return rParameter;
    }

    // get parameters for matrices (LHS)
    virtual double& GetFirstDerivativeInertialParameter(double& rParameter)
    {
      rParameter = 1.0;
      return rParameter;
    }

    virtual double& GetSecondDerivativeInertialParameter(double& rParameter)
    {
      rParameter = 1.0;
      return rParameter;
    }
   
    // set input variable (constrained or dof variable)
    void SetInputVariable(const TVariableType& rVariable)
    {
      mpInputVariable = &rVariable;

      this->SetPointerAssignMethod();
    }
        

    // get primary variable name
    std::string GetPrimaryVariableName()
    {
      return (*this->mpPrimaryVariable).Name();
    }

    // get primary variable name
    std::string GetVariableName()
    {
      return (*this->mpVariable).Name();
    }
    
    // check if the integration method has the variable
    bool HasVariableName(const std::string& rVariableName)
    {
      if( this->mpVariable != nullptr ){ 
        if( rVariableName == (*this->mpVariable).Name() ){
          return true;
        }
        else if( this->mpFirstDerivative != nullptr ){ 
          if( rVariableName == (*this->mpFirstDerivative).Name() ){
            return true;
          }
          else if( this->mpSecondDerivative != nullptr ){
            if( rVariableName == (*this->mpSecondDerivative).Name() ){
              return true;
            }
          }
        }
      }     
      return false;
    }

    
    // check if the integration method has the step variable (step variable)
    virtual bool HasStepVariable()
    {
      return false;
    }
    
    // set step variable (step variable)
    virtual void SetStepVariable(const TVariableType& rStepVariable)
    {
      KRATOS_ERROR << " Calling SetStepVariable from time integration base class " <<std::endl;
    }

    // assign
    virtual void Assign(NodeType& rNode)
    {
     KRATOS_TRY

     (this->*this->mpAssign)(rNode);
      
     KRATOS_CATCH( "" )
    }

    // predict
    virtual void Predict(NodeType& rNode)
    {
     KRATOS_TRY
     
     (this->*this->mpPredict)(rNode);
  
     KRATOS_CATCH( "" )
    }

    // update
    virtual void Update(NodeType& rNode)
    {
     KRATOS_TRY
       
     (this->*this->mpUpdate)(rNode);  

     KRATOS_CATCH( "" )
    }

    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * @return 0 all ok
     */
    virtual int Check( const ProcessInfo& rCurrentProcessInfo )
    {
      KRATOS_TRY

      // Check that all required variables have been registered
      if( mpVariable == nullptr ){
        KRATOS_ERROR << " time integration method Variable not set " <<std::endl;
      }
      else{
        KRATOS_CHECK_VARIABLE_KEY((*mpVariable));
      }

      if( mpPrimaryVariable == nullptr ){
        KRATOS_ERROR << " time integration method PrimaryVariable not set " <<std::endl;
      }
      else{
        KRATOS_CHECK_VARIABLE_KEY((*mpPrimaryVariable));
      }
      
      // if( mpInputVariable == nullptr ){
      //   KRATOS_ERROR << " time integration method InputVariable not set " <<std::endl;
      // }
      // else{
      //   KRATOS_CHECK_VARIABLE_KEY((*mpInputVariable));
      // }
      
      // if( mpFirstDerivative == nullptr ){
      //   KRATOS_ERROR << " time integration method FirstDerivative not set " <<std::endl;
      // }
      // else{
      //   KRATOS_CHECK_VARIABLE_KEY((*mpFirstDerivative));
      // }

      // if( mpSecondDerivative == nullptr ){
      //   KRATOS_ERROR << " time integration method SecondDerivative not set " <<std::endl;
      // }
      // else{
      //   KRATOS_CHECK_VARIABLE_KEY((*mpSecondDerivative));
      // }
      
      return 0;
      
      KRATOS_CATCH("")
    }
    
    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{


    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "TimeIntegrationMethod";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "TimeIntegrationMethod";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
      rOStream << "TimeIntegrationMethod Data";     
    }

    
    ///@}
    ///@name Friends
    ///@{


    ///@}
    
  protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    // method variables and derivatives
    
    VariablePointer mpVariable;

    VariablePointer mpFirstDerivative;

    VariablePointer mpSecondDerivative;

    // primary variable (calculated variable 'dof')
    
    VariablePointer mpPrimaryVariable;

    
    // input variable (imposed variable or calculated variable)
    
    VariablePointer mpInputVariable;


    // method pointer
    
    MethodPointer mpAssign;
    
    MethodPointer mpPredict;
    
    MethodPointer mpUpdate;
    
    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    // set methods from primary variable
    void SetPointerMethods()
    {
      if( this->mpPrimaryVariable != nullptr ){ 
        
        if( this->mpVariable != nullptr ){ 
          if( *this->mpPrimaryVariable == *this->mpVariable ){
            mpPredict = &TimeIntegrationMethod::PredictFromVariable;
            mpUpdate  = &TimeIntegrationMethod::UpdateFromVariable;
          }
          else if( this->mpFirstDerivative != nullptr ){ 
            if( *this->mpPrimaryVariable == *this->mpFirstDerivative ){
              mpPredict = &TimeIntegrationMethod::PredictFromFirstDerivative;
              mpUpdate  = &TimeIntegrationMethod::UpdateFromFirstDerivative;
            }
            else if( this->mpSecondDerivative != nullptr ){
              if( *this->mpPrimaryVariable == *this->mpSecondDerivative ){
                mpPredict = &TimeIntegrationMethod::PredictFromSecondDerivative;
                mpUpdate  = &TimeIntegrationMethod::UpdateFromSecondDerivative;
              }
            }
          }
        }
      }
      
    }

    // set methods from input variable
    void SetPointerAssignMethod()
    {
      if( this->mpInputVariable != nullptr ){ 
        
        if( this->mpVariable != nullptr ){ 
          if( *this->mpInputVariable == *this->mpVariable ){
            mpAssign  = &TimeIntegrationMethod::AssignFromVariable;
          }
          else if( this->mpFirstDerivative != nullptr ){ 
            if( *this->mpInputVariable == *this->mpFirstDerivative ){
              mpAssign  = &TimeIntegrationMethod::AssignFromFirstDerivative;
            }
            else if( this->mpSecondDerivative != nullptr ){
              if( *this->mpInputVariable == *this->mpSecondDerivative ){
                mpAssign  = &TimeIntegrationMethod::AssignFromSecondDerivative;
              }
            }
          }
        }
      }
      
    }
    
    virtual void AssignFromVariable(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling predict from variable from time integration base class " <<std::endl;
    }
    
    virtual void AssignFromFirstDerivative(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling predict from first derivative from time integration base class " <<std::endl;
    }

    virtual void AssignFromSecondDerivative(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling predict from second derivative from time integration base class " <<std::endl;
    }

    virtual void AssignVariable(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling predict variable from time integration base class " <<std::endl;
    }

    virtual void AssignFirstDerivative(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling predict first derivative from time integration base class " <<std::endl;
    }

    virtual void AssignSecondDerivative(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling predict second derivative from time integration base class " <<std::endl;
    }
    
    virtual void PredictFromVariable(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling predict from variable from time integration base class " <<std::endl;
    }
    
    virtual void PredictFromFirstDerivative(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling predict from first derivative from time integration base class " <<std::endl;
    }

    virtual void PredictFromSecondDerivative(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling predict from second derivative from time integration base class " <<std::endl;
    }

    virtual void PredictVariable(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling predict variable from time integration base class " <<std::endl;
    }

    virtual void PredictFirstDerivative(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling predict first derivative from time integration base class " <<std::endl;
    }

    virtual void PredictSecondDerivative(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling predict second derivative from time integration base class " <<std::endl;
    }

    virtual void UpdateFromVariable(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling update from variable from time integration base class " <<std::endl;
    }
    
    virtual void UpdateFromFirstDerivative(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling update from first derivative from time integration base class " <<std::endl;
    }

    virtual void UpdateFromSecondDerivative(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling update from second derivative from time integration base class " <<std::endl;
    }
    
    virtual void UpdateVariable(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling update variable from time integration base class " <<std::endl;
    }

    virtual void UpdateFirstDerivative(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling update first derivative from time integration base class " <<std::endl;
    }

    virtual void UpdateSecondDerivative(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling update second derivative from time integration base class " <<std::endl;
    }


    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{
  
    ///@}

  private:

    ///@name Static Member Variables
    ///@{
  
    ///@}
    ///@name Member Variables
    ///@{
  
    ///@}
    ///@name Private Operators
    ///@{
  
    ///@}
    ///@name Private Operations
    ///@{
  
    ///@}
    ///@name Private  Access
    ///@{
  
    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
      rSerializer.save("Variable", mpVariable->Name());
      rSerializer.save("FirstDerivative", mpFirstDerivative->Name());
      rSerializer.save("SecondDerivative", mpSecondDerivative->Name());
      rSerializer.save("PrimaryVariable", mpPrimaryVariable->Name());
      rSerializer.save("InputVariable", mpInputVariable->Name());
    };

    virtual void load(Serializer& rSerializer)
    {
      std::string Name;
      rSerializer.load("Variable", Name);
      mpVariable = static_cast<VariablePointer>(KratosComponents<VariableData>::pGet(Name));
      rSerializer.load("FirstDerivative", Name);
      mpFirstDerivative = static_cast<VariablePointer>(KratosComponents<VariableData>::pGet(Name));
      rSerializer.load("SecondDerivative", Name);
      mpSecondDerivative = static_cast<VariablePointer>(KratosComponents<VariableData>::pGet(Name));
      rSerializer.load("PrimaryVariable", Name);
      mpPrimaryVariable = static_cast<VariablePointer>(KratosComponents<VariableData>::pGet(Name));
      rSerializer.load("InputVariable", Name);
      mpInputVariable = static_cast<VariablePointer>(KratosComponents<VariableData>::pGet(Name));

      this->SetPointerMethods();
      this->SetPointerAssignMethod();
    };
    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{
  
    ///@} 
  }; // Class TimeIntegrationMethod
  
  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  template<class TVariableType, class TValueType>
  inline std::istream & operator >> (std::istream & rIStream, TimeIntegrationMethod<TVariableType,TValueType>& rThis)
  {
    return rIStream;
  }

  template<class TVariableType, class TValueType>
  inline std::ostream & operator << (std::ostream & rOStream, const TimeIntegrationMethod<TVariableType,TValueType>& rThis)
  {
    return rOStream << rThis.Info();
  }

  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_TIME_INTEGRATION_METHOD_H_INCLUDED defined
