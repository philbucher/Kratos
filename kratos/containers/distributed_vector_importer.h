//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
#if !defined(KRATOS_DISTRIBUTED_VECTOR_IMPORTER_H_INCLUDED )
#define  KRATOS_DISTRIBUTED_VECTOR_IMPORTER_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "containers/distributed_system_vector.h"
#include "utilities/parallel_utilities.h"
#include "utilities/atomic_utilities.h"

namespace Kratos
{
///@addtogroup ApplicationNameApplication
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

/// Provides a DistributedVectorImporter which implements FEM assemble capabilities
template<class TDataType=double, class TIndexType=std::size_t>
class DistributedVectorImporter
{
public:
    ///@name Type Definitions
    ///@{
    typedef TIndexType IndexType;
    typedef int MpiIndexType;

    /// Pointer definition of DistributedVectorImporter
    KRATOS_CLASS_POINTER_DEFINITION(DistributedVectorImporter);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    template< class TGlobalIndicesVectorType>
    DistributedVectorImporter(
        const DataCommunicator& rComm,
        const TGlobalIndicesVectorType& rGlobalIndices,
        const DistributedNumbering<IndexType>& rNumbering)
        :
        mrComm(rComm)
    {
        mpNumbering = Kratos::make_unique< DistributedNumbering<IndexType> >(rNumbering);
        mImportedDataSize = rGlobalIndices.size();
        std::unordered_map<int, std::vector<IndexType>> to_recv_by_color; //do not need to store this

        for(unsigned int local_i=0; local_i<rGlobalIndices.size(); ++local_i)
        {
            IndexType global_i = rGlobalIndices[local_i];
            MpiIndexType owner_rank = mpNumbering->OwnerRank(global_i);
            IndexType remote_local_i = mpNumbering->RemoteLocalId(global_i, owner_rank);

            mlocal_i_by_color[owner_rank].push_back(local_i);
            to_recv_by_color[owner_rank].push_back(remote_local_i);
        }

        mIdOfLocallyOwnedTerms = mlocal_i_by_color[GetComm().Rank()];
        mLocallyOwnedIds = to_recv_by_color[GetComm().Rank()];

        //compute communication plan
        std::vector<MpiIndexType> send_list;
        for(const auto& item : to_recv_by_color)
        {
            MpiIndexType cpu_id = item.first;
            if(cpu_id != GetComm().Rank())
                send_list.push_back(cpu_id);
        }
        mvector_comm_colors = MPIColoringUtilities::ComputeCommunicationScheduling(send_list, rComm);

        //ensure that we have lists for all colors;
        for(auto color : mvector_comm_colors)
        {
            if(color >= 0) //-1 would imply no communication
            {
                mlocal_i_by_color[color]; //note that here we touch the entry and we create it if not existing
                to_recv_by_color[color];
            }
        }

        //communicate the remote_local_id so that the other node knows what to send
        for(auto color : mvector_comm_colors)
        {
            if(color >= 0) //-1 would imply no communication
            {
                //NOTE: this can be made nonblocking 
                mto_send_by_color[color] = rComm.SendRecv(to_recv_by_color[color], color, color); //TODO, we know all the sizes, we shall use that!
            }
        }
    }

    ///this function returns a local array containing the values identified by the rGlobalIndices list passed in the constructor
    DenseVector<TDataType> ImportData(const DistributedSystemVector<TDataType, TIndexType>& data_vector) const
    {
        DenseVector<TDataType> ImportedData(mImportedDataSize);

        std::vector<TDataType> send_buffer;
        std::vector<TDataType> recv_buffer;
        for(auto color : mvector_comm_colors)
        {
            if(color >= 0) //-1 would imply no communication
            {
                const auto& local_ids = mlocal_i_by_color.find(color)->second;
                const auto& to_send_ids = mto_send_by_color.find(color)->second; 

                send_buffer.resize(to_send_ids.size());
                recv_buffer.resize(local_ids.size());

                for(IndexType i=0; i<to_send_ids.size(); ++i)
                    send_buffer[i] = data_vector(to_send_ids[i]);

                //NOTE: this can be made nonblocking
                GetComm().SendRecv(send_buffer, color, 0, recv_buffer, color, 0); //TODO, we know all the sizes, we shall use that!

                //write the recv_data onto the output
                for(IndexType i=0; i<recv_buffer.size(); ++i)
                    ImportedData(local_ids[i]) = recv_buffer[i];
            }
        }
        //treat local data
        for(IndexType i=0; i<mLocallyOwnedIds.size(); ++i )

            ImportedData(mIdOfLocallyOwnedTerms[i]) = data_vector( mLocallyOwnedIds[i] );

        return ImportedData;

    }

    /// Destructor.
    virtual ~DistributedVectorImporter(){}

    ///@}
    ///@name Operators
    ///@{
    const DataCommunicator& GetComm() const{
        return mrComm;
    }

    
    ///@}
    ///@name Operations
    ///@{

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
    buffer << "DistributedVectorImporter" ;
    return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "DistributedVectorImporter";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

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


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


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
    const DataCommunicator& mrComm;
    typename DistributedNumbering<IndexType>::UniquePointer mpNumbering;


    IndexType mImportedDataSize;
    std::unordered_map<int, std::vector<IndexType>> mto_send_by_color;
    std::unordered_map<MpiIndexType, std::vector<IndexType>> mlocal_i_by_color;
    std::vector<IndexType> mLocallyOwnedIds;
    std::vector<IndexType> mIdOfLocallyOwnedTerms;
    std::vector<int> mvector_comm_colors;

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
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    DistributedVectorImporter& operator=(DistributedVectorImporter const& rOther){}

    /// Copy constructor.
    DistributedVectorImporter(DistributedVectorImporter const& rOther){}

    ///@}

}; // Class DistributedVectorImporter

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TDataType, class TIndexType>
inline std::istream& operator >> (std::istream& rIStream,
                DistributedVectorImporter<TDataType,TIndexType>& rThis)
                {
                    return rIStream;
                }

/// output stream function
template<class TDataType, class TIndexType>
inline std::ostream& operator << (std::ostream& rOStream,
                const DistributedVectorImporter<TDataType,TIndexType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_DISTRIBUTED_VECTOR_IMPORTER_H_INCLUDED  defined


