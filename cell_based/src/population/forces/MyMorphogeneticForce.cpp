/*

Copyright (c) 2005-2018, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "MyMorphogeneticForce.hpp"
#include "MyNoPBCToroidal2dVertexMesh.hpp"

template<unsigned DIM>
MyMorphogeneticForce<DIM>::MyMorphogeneticForce()
   : AbstractForce<DIM>(),
     mIfEquilibrateForAWhile(false),
     mTimeForEquilibrium(0.0)
{
}

template<unsigned DIM>
MyMorphogeneticForce<DIM>::~MyMorphogeneticForce()
{
}

template<unsigned DIM>
void MyMorphogeneticForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{    
    // Throw an exception message if not using a VertexBasedCellPopulation
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("MorphogeneticForce is to be used with a VertexBasedCellPopulation only");
    }

    double t_now = SimulationTime::Instance()->GetTime();
    if ( !(mIfEquilibrateForAWhile && t_now<=mTimeForEquilibrium) )
    {
        // Define some helper variables
        VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
        unsigned num_nodes = p_cell_population->GetNumNodes();
        unsigned num_top_boundary_nodes = 0;
        unsigned num_bottom_boundary_nodes = 0;

        double t_now = SimulationTime::Instance()->GetTime();
        if (this->mIfEquilibrateForAWhile && t_now>this->mTimeForEquilibrium)
        {
            // count the nodes of top and bottom boundaries
            for (unsigned node_index=0; node_index<num_nodes; node_index++)
            {
                Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);
                
                if (p_this_node->IsBoundaryNode())
                {
                    double y_coord= p_this_node->rGetLocation()[1];

                    if (y_coord > this->mCenterYCoordination)
                    {
                        num_top_boundary_nodes += 1;
                    }
                    else
                    {
                        num_bottom_boundary_nodes +=1;
                    }
                }
            }
        }
        
        std::cout<<"time="<<t_now<<std::endl;
        std::cout<<"num_top_boundary_nodes="<<num_top_boundary_nodes<<std::endl;
        std::cout<<"num_bottom_boundary_nodes="<<num_bottom_boundary_nodes<<std::endl;

        // Iterate over vertices in the cell population and distribute morphogenetic force on boundary nodes evenly
        for (unsigned node_index=0; node_index<num_nodes; node_index++)
        {
            Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);
            c_vector<double, DIM> force_on_node = zero_vector<double>(DIM);
            
            if (this->mAddPullingForceEvenlyOnNodesOfLeadingCell)
            {
                if (p_this_node->IsBoundaryNode())
                {
                    double y_coord= p_this_node->rGetLocation()[1];
                    if (y_coord > this->mCenterYCoordination)
                    {
                        force_on_node[1] += mPullingForceOnLeadingCell/num_top_boundary_nodes;
                    }
                    else
                    {
                        force_on_node[1] -= mPullingForceOnLeadingCell/num_bottom_boundary_nodes;
                    }
                }
            }

            // tmp: possilbe output for testing morphogenetic force:
            bool output_while_testing_morphogenetic_force = false;
            if (output_while_testing_morphogenetic_force)
            {
                VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(* p_this_node->rGetContainingElementIndices().begin());
                if (p_element->GetIsLeadingCell() && p_element->GetGroupNumber()>0)
                {
                    double t_now = SimulationTime::Instance()->GetTime();
                    if (p_this_node->rGetContainingElementIndices().size()==1)
                    {
                        std::cout << "time: " << t_now << std::endl;
                        std::cout << "node of Leading Cell" << std::endl
                                << " node index=" << p_this_node->GetIndex() << std::endl
                                << " node location: x=" << p_this_node->rGetLocation()[0] << ", y=" << p_this_node->rGetLocation()[1] << std::endl 
                                << " total morphogenetic force =" << mPullingForceOnLeadingCell << std::endl;
                    }
                    else if (p_this_node->rGetContainingElementIndices().size()==2)
                    {
                        std::cout << "time: " << t_now << std::endl;
                        std::cout << "node of Leading Cell" << std::endl
                                << " node index=" << p_this_node->GetIndex() << std::endl
                                << " node location: x=" << p_this_node->rGetLocation()[0] << ", y=" << p_this_node->rGetLocation()[1] << std::endl 
                                << " total morphogenetic force =" << mPullingForceOnLeadingCell << std::endl;
                    }
                }
            }

            p_cell_population->GetNode(node_index)->AddAppliedForceContribution(force_on_node);

        }
        // end of 'Iterate over nodes(vertices) in the cell population'
    }
}

template<unsigned DIM>
void MyMorphogeneticForce<DIM>::SetAddPullingForceEvenlyOnNodesOfLeadingCell(bool addPullingForceEvenlyOnNodesOfLeadingCell)
{
    mAddPullingForceEvenlyOnNodesOfLeadingCell = addPullingForceEvenlyOnNodesOfLeadingCell;
}

template<unsigned DIM>
void MyMorphogeneticForce<DIM>::SetPullingForceOnLeadingCell(double pullingForceOnLeadingCell)
{
    mPullingForceOnLeadingCell = pullingForceOnLeadingCell;
}

template<unsigned DIM>
void MyMorphogeneticForce<DIM>::SetCenterYCoordination(double centerYCoordination)
{
    mCenterYCoordination = centerYCoordination;
}

template<unsigned DIM>
void MyMorphogeneticForce<DIM>::SetIfEquilibrateForAWhile(bool ifEquilibrateForAWhile)
{
    mIfEquilibrateForAWhile = ifEquilibrateForAWhile;
}

template<unsigned DIM>
void MyMorphogeneticForce<DIM>::SetEndTimeForEquilibrium(double timeForEquilibrium)
{
    mTimeForEquilibrium = timeForEquilibrium;
}

template<unsigned DIM>
void MyMorphogeneticForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MorphogeneticForce>" << mPullingForceOnLeadingCell << "</MorphogeneticForce>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class MyMorphogeneticForce<1>;
template class MyMorphogeneticForce<2>;
template class MyMorphogeneticForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MyMorphogeneticForce)
