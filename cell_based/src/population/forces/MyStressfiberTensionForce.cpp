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

#include "MyStressfiberTensionForce.hpp"
#include "MyToroidal2dVertexMesh.hpp"
#include <cassert>

template<unsigned DIM>
MyStressfiberTensionForce<DIM>::MyStressfiberTensionForce()
   : AbstractForce<DIM>(),
     mIfEquilibrateForAWhile(false),
     mTimeForEquilibrium(0.0),
     mFlag(0)
{
}

template<unsigned DIM>
MyStressfiberTensionForce<DIM>::~MyStressfiberTensionForce()
{
}

template<unsigned DIM>
void MyStressfiberTensionForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{    
    // Throw an exception message if not using a VertexBasedCellPopulation
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("MyStressfiberTensionForce is to be used with a VertexBasedCellPopulation only");
    }

    assert(DIM>1);

    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    double t_now = SimulationTime::Instance()->GetTime();

    if (mIfEquilibrateForAWhile && (t_now>=mTimeForEquilibrium) && (mFlag==0) )
    {
        unsigned global_index = 0;
        for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
            elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
            ++elem_iter)
        {
            unsigned elem_index = elem_iter->GetIndex();
            VertexElement<DIM, DIM>* p_element = p_cell_population->rGetMesh().GetElement(elem_index);
            std::cout<<"there are "<< p_element->GetNumStressfibers() << " stress fibers in element "<<p_element->GetIndex()<<std::endl;
            std::cout<<"creating stress fibers..."<<std::endl;

            // std::vector<VertexElement<1,2>*> stressfibers;

            std::vector<Node<DIM>*> stressfiber_nodes;
            unsigned num_nodes_elem = p_element->GetNumNodes();
            unsigned lowest_local_index = 0;
            unsigned highest_local_index = 0;
            double y_lowest = p_element->GetNodeLocation(lowest_local_index)[1];
            double y_highest = p_element->GetNodeLocation(highest_local_index)[1];
            for (unsigned local_index=0; local_index<num_nodes_elem; local_index++)
            {
                double y_coord = p_element->GetNodeLocation(local_index)[1];
                if (y_coord<y_lowest)
                {
                    lowest_local_index = local_index;
                    y_lowest = y_coord;
                }
                if (y_coord>y_highest)
                {
                    highest_local_index = local_index;
                    y_highest = y_coord;
                }

                // unsigned node_global_index = p_element->GetNodeGlobalIndex(local_index);
                // Node<DIM>* p_this_node = p_cell_population->GetNode(node_global_index);
                // double y_coord= p_this_node->rGetLocation()[1];
            }
            unsigned next_node_local_index_of_lowest_node = (lowest_local_index+1)%num_nodes_elem;
            unsigned previous_node_local_index_of_highest_node = (num_nodes_elem+highest_local_index-1)%num_nodes_elem;

            stressfiber_nodes.push_back(p_element->GetNode(lowest_local_index));
            stressfiber_nodes.push_back(p_element->GetNode(next_node_local_index_of_lowest_node));
            stressfiber_nodes.push_back(p_element->GetNode(previous_node_local_index_of_highest_node));
            stressfiber_nodes.push_back(p_element->GetNode(highest_local_index));

            c_vector<double,2> end_points_ratio = zero_vector<double>(2);
            end_points_ratio[0] = 0.5;
            end_points_ratio[1] = 0.5;

            VertexElement<DIM-1,DIM>* p_stressfiber = new VertexElement<DIM-1,DIM>(global_index, stressfiber_nodes, end_points_ratio);
            // stressfibers.push_back(p_stressfiber);

            p_element->AddStressfiber(p_stressfiber);
            std::cout<< p_element->GetNumStressfibers() << " stress fibers are created in element "<<p_element->GetIndex()<<std::endl;
            
            VertexElement<DIM-1, DIM>* p_get_stressfiber = p_element->GetStressfiber(0);
            Node<DIM>* p_node0 = p_get_stressfiber->GetStressfiberNode(0);
            Node<DIM>* p_node1 = p_get_stressfiber->GetStressfiberNode(1);
            Node<DIM>* p_node2 = p_get_stressfiber->GetStressfiberNode(2);
            Node<DIM>* p_node3 = p_get_stressfiber->GetStressfiberNode(3);
            c_vector<double,2> ratio = p_get_stressfiber->GetStressfiberEndpointsratio();

            std::cout<< "Node 0: Index-" << p_node0->GetIndex() << ", Location:(" << p_node0->rGetLocation()[0]<<", "<<p_node0->rGetLocation()[1]<< ")" << std::endl;
            std::cout<< "Node 1: Index-" << p_node1->GetIndex() << ", Location:(" << p_node1->rGetLocation()[0]<<", "<<p_node1->rGetLocation()[1]<< ")" << std::endl;
            std::cout<< "Node 2: Index-" << p_node2->GetIndex() << ", Location:(" << p_node2->rGetLocation()[0]<<", "<<p_node2->rGetLocation()[1]<< ")" << std::endl;
            std::cout<< "Node 3: Index-" << p_node3->GetIndex() << ", Location:(" << p_node3->rGetLocation()[0]<<", "<<p_node3->rGetLocation()[1]<< ")" << std::endl;
            std::cout<< "Ratio1= " << ratio[0] << std::endl;
            std::cout<< "Ratio2= " << ratio[1] << std::endl;

            std::cout<<"deleting stress fibers..."<<std::endl;
            p_element->DeleteStressfiber(0);
            std::cout<<"there are "<< p_element->GetNumStressfibers() << " stress fibers in element "<<p_element->GetIndex()<<std::endl;
            std::cout<<"\n"<<std::endl;

            global_index++;
        }

        mFlag++;
    }

    if ( !(mIfEquilibrateForAWhile && t_now<=mTimeForEquilibrium) )
    {
        // Define some helper variables
        VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
        unsigned num_nodes = p_cell_population->GetNumNodes();

        // std::cout<<"time="<<t_now<<std::endl;

        // Iterate over vertices in the cell population and distribute morphogenetic force on boundary nodes evenly
        for (unsigned node_index=0; node_index<num_nodes; node_index++)
        {
            // Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);
            c_vector<double, DIM> force_on_node = zero_vector<double>(DIM);

            p_cell_population->GetNode(node_index)->AddAppliedForceContribution(force_on_node);
        }
        // end of 'Iterate over nodes(vertices) in the cell population'
    }
}

template<unsigned DIM>
void MyStressfiberTensionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<StressfiberTensionForce>" << " " << "</StressfiberTensionForce>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

template<unsigned DIM>
void MyStressfiberTensionForce<DIM>::SetIfEquilibrateForAWhile(bool ifEquilibrateForAWhile)
{
    mIfEquilibrateForAWhile = ifEquilibrateForAWhile;
}

template<unsigned DIM>
void MyStressfiberTensionForce<DIM>::SetEndTimeForEquilibrium(double timeForEquilibrium)
{
    mTimeForEquilibrium = timeForEquilibrium;
}

template<unsigned DIM>
void MyStressfiberTensionForce<DIM>::SetFlagForStressfiberCreation(unsigned flag)
{
    mFlag = flag;
}

// Explicit instantiation
template class MyStressfiberTensionForce<1>;
template class MyStressfiberTensionForce<2>;
template class MyStressfiberTensionForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MyStressfiberTensionForce)
