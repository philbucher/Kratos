//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                    
//

#if !defined(KRATOS_KERNEL_H_INCLUDED )
#define  KRATOS_KERNEL_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <unordered_map>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_application.h"


namespace Kratos
{

///@name Kratos Classes
///@{

/// Kernel is in charge of synchronization the whole system of Kratos itself and its appliction.
/** Kernel is the first component of the Kratos to be created and then used to plug the application into Kratos.
    Kernel takes the list of variables defined in the Variables file in Kratos and then by adding each application
    synchronizes the variables of this application with its variables and add the new ones to the Kratos.
    After adding all applications its time to initialize the Kratos to assign variables key to the list of all variables
    in Kratos and all added applications. Finally the initialized variables with keys are synchronized in each
    application in time of calling InitializeApplication method for each of them.
    The sequence of using Kernel is as follow:
1. Creating the Kernel using its default constructor
2. Adding applications to Kernel using AddApplication method
3. Initializing the Kernel using Initialize method
4. Initializing the applications using InitializeApplication method

    It is very important to perform all this step exactly in the same order as described above.

    @see AddApplication
    @see Initialize
    @see InitializeApplication
    @see KratosApplication
*/
class KRATOS_API(KRATOS_CORE) Kernel
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Kernel
    KRATOS_CLASS_POINTER_DEFINITION(Kernel);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** The default constructor creates a list of registerd variables in variables.cpp
        by calling the RegisterVariables method of application class.

        @see KratosApplication

    */
    Kernel();

    /// Copy constructor.
    /** This constructor is empty
    */
    Kernel(Kernel const& rOther) {}


    /// Destructor.
    virtual ~Kernel() {}


    ///@}
    ///@name Operations
    ///@{

    /// Pluging an application into Kratos.
    /** This method first call the register method of the new application in order to create the
        components list of the application and then syncronizes the lists of its components with Kratos ones.
        The synchronized lists are
      - Variables
      - Elements
      - Conditions

    @param NewApplication The application to be added and synchronized
    */
    void AddApplication(KratosApplication::Pointer pNewApplication)
    {
        //typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > array_1d_component_type;

        pNewApplication->Register();
        mApplicationsList[pNewApplication->Name()] = pNewApplication;
		
		/*
        KratosComponents<VariableData>::GetComponents().insert(NewApplication.GetVariables().begin(),NewApplication.GetVariables().end());
		
        KratosComponents<Variable<int> >::GetComponents().insert(NewApplication.GetComponents(Variable<int>("NONE")).begin(),
                NewApplication.GetComponents(Variable<int>("NONE")).end());
        KratosComponents<Variable<unsigned int> >::GetComponents().insert(NewApplication.GetComponents(Variable<unsigned int>("NONE")).begin(),
               NewApplication.GetComponents(Variable<unsigned int>("NONE")).end());
        KratosComponents<Variable<double> >::GetComponents().insert(NewApplication.GetComponents(Variable<double>("NONE")).begin(),
                NewApplication.GetComponents(Variable<double>("NONE")).end());
        KratosComponents<Variable<array_1d<double, 3> > >::GetComponents().insert(NewApplication.GetComponents(Variable<array_1d<double, 3> >("NONE")).begin(),
                NewApplication.GetComponents(Variable<array_1d<double, 3> >("NONE")).end());
        KratosComponents<Variable<Vector> >::GetComponents().insert(NewApplication.GetComponents(Variable<Vector>("NONE")).begin(),
                NewApplication.GetComponents(Variable<Vector>("NONE")).end());
        KratosComponents<Variable<Matrix> >::GetComponents().insert(NewApplication.GetComponents(Variable<Matrix>("NONE")).begin(),
                NewApplication.GetComponents(Variable<Matrix>("NONE")).end());
        
		Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3> > temp_adaptor(DISPLACEMENT, 0); // the displacement is not important, only an array_1d variable is needed!
 KratosApplication() = delete;
        KratosComponents<array_1d_component_type>::GetComponents().insert(NewApplication.GetComponents(array_1d_component_type("NONE", temp_adaptor)).begin(),				  
				NewApplication.GetComponents(array_1d_component_type("NONE", temp_adaptor)).end());

        KratosComponents<Element>::GetComponents().insert(NewApplication.GetElements().begin(),
                NewApplication.GetElements().end());
        KratosComponents<Condition>::GetComponents().insert(NewApplication.GetConditions().begin(),
                NewApplication.GetConditions().end());

//        KratosComponents<Variable<double> >::GetComponents().insert(NewApplication.GetComponents(Variable<double>("NONE")).begin(),
//                NewApplication.GetComponents(Variable<double>("NONE")).end());


        Serializer::GetRegisteredObjects().insert(NewApplication.GetRegisteredObjects().begin(), NewApplication.GetRegisteredObjects().end());
        Serializer::GetRegisteredObjectsName().insert(NewApplication.GetRegisteredObjectsName().begin(), NewApplication.GetRegisteredObjectsName().end());
		*/
    }

    /// Assign sequential key to the registered variables.
    /** This method assigns a sequential key to all registerd variables in kratos and all added applications.
        It is very important to call this function after adding ALL necessary applications using AddApplication
        methods before calling this function. Otherwise it leads to uninitialized variables with key 0!

        @see AddApplication
        @see InitializeApplication
    */
    void Initialize()
    {
        unsigned int j = 0;
        for(KratosComponents<VariableData>::ComponentsContainerType::iterator i = KratosComponents<VariableData>::GetComponents().begin() ;
                i != KratosComponents<VariableData>::GetComponents().end() ; i++)
            //const_cast<VariableData&>(i->second.get()).SetKey(++j);
            i->second->SetKey(++j);
    }

    /// Initializes and synchronizes the list of variables, elements and conditions in each application.
    /** This method gives the application the list of all variables, elements and condition which is registered
        by kratos and all other added applications.


        @see AddApplication
        @see Initialize
    */
    void InitializeApplication(KratosApplication& NewApplication)
    {
        /*
		NewApplication.SetComponents(KratosComponents<VariableData>::GetComponents());
        NewApplication.SetComponents(KratosComponents<Element>::GetComponents());
        NewApplication.SetComponents(KratosComponents<Condition>::GetComponents());

        NewApplication.GetRegisteredObjects().insert(Serializer::GetRegisteredObjects().begin(), Serializer::GetRegisteredObjects().end());
        NewApplication.GetRegisteredObjectsName().insert(Serializer::GetRegisteredObjectsName().begin(), Serializer::GetRegisteredObjectsName().end());
		*/
    }


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;

    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
    std::unordered_map< std::string, KratosApplication::Pointer > mApplicationsList;
//     KratosApplication mKratosApplication;


    ///@}
    ///@name Private Operations
    ///@{

    void RegisterVariables();

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    Kernel& operator=(Kernel const& rOther);


    ///@}

}; // Class Kernel

///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  Kernel& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Kernel& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_KERNEL_H_INCLUDED  defined 


