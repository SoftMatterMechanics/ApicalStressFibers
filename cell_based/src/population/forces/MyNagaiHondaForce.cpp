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

#include "MyNagaiHondaForce.hpp"
#include "MyToroidal2dVertexMesh.hpp"

template<unsigned DIM>
MyNagaiHondaForce<DIM>::MyNagaiHondaForce()
   : AbstractForce<DIM>(),
     mNagaiHondaDeformationEnergyParameter(1.0),
     mNagaiHondaMembraneSurfaceEnergyParameter(0.1),
     mNagaiHondaCellCellAdhesionEnergyParameter(0.0),
     mNagaiHondaCellBoundaryAdhesionEnergyParameter(0.0),
     mUseFixedTargetArea(true),
     mIfUseFaceElementToGetAdhesionParameter(false)
{
}

template<unsigned DIM>
MyNagaiHondaForce<DIM>::~MyNagaiHondaForce()
{
}

template<unsigned DIM>
void MyNagaiHondaForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{    
    // Throw an exception message if not using a VertexBasedCellPopulation
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("NagaiHondaForce is to be used with a VertexBasedCellPopulation only");
    }

    double t_now = SimulationTime::Instance()->GetTime();

    // Define some helper variables
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    unsigned num_nodes = p_cell_population->GetNumNodes();
    unsigned num_elements = p_cell_population->GetNumElements();

    // Begin by computing the area and perimeter of each element in the mesh, to avoid having to do this multiple times
    std::vector<double> element_areas(num_elements);
    std::vector<double> element_perimeters(num_elements);
    std::vector<double> target_areas(num_elements);
    // std::vector<double> element_perimeter_elasticity(num_elements);
    for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
         elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        element_areas[elem_index] = p_cell_population->rGetMesh().GetVolumeOfElement(elem_index);
        element_perimeters[elem_index] = p_cell_population->rGetMesh().GetSurfaceAreaOfElement(elem_index);
        try
        {
            // If we haven't specified a growth modifier, there won't be any target areas in the CellData array and CellData
            // will throw an exception that it doesn't have "target area" entries.  We add this piece of code to give a more
            // understandable message. There is a slight chance that the exception is thrown although the error is not about the
            // target areas.
            if (mUseFixedTargetArea)
            {
                target_areas[elem_index] = mFixedTargetArea;
            }
            else
            {
                target_areas[elem_index] = p_cell_population->GetCellUsingLocationIndex(elem_index)->GetCellData()->GetItem("target area");
                // std::cout << "elem_index = " << elem_index << ", target_area=" << target_areas[elem_index] << std::endl;
            }    
        }
        catch (Exception&)
        {
            EXCEPTION("You need to add an AbstractTargetAreaModifier to the simulation in order to use NagaiHondaForce");
        }

        // CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(elem_index);
        // element_perimeter_elasticity[elem_index] = p_cell->GetCellData()->GetItem("perimeter_elasticity");
    }

    // // original non-homogeneous SSA method using leading length:
    // double y_coord_leading_edge = 0.0;
    // for (unsigned node_index=0; node_index<num_nodes; node_index++)
    // {
    //     Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);
    //     if (p_this_node->rGetLocation()[1] > y_coord_leading_edge)
    //         y_coord_leading_edge= p_this_node->rGetLocation()[1];
    // }

    // // For My Detach Pattern Method, using group number:
    // unsigned number_of_groups = p_cell_population->rGetMesh().GetNumberOfGroups();
    // std::vector<double> leading_tops_of_groups(number_of_groups);
    // if (mUseMyDetachPatternMethod)
    // {
    //     for (unsigned i=0; i<number_of_groups; i++)
    //         leading_tops_of_groups[i] = p_cell_population->rGetMesh().GetLeadingTopOfTheGroup(i);
    // }

    // unsigned node_index_for_leading_top_of_the_group0 = 0;
    // if (mAddPullingForceOnNodeIndividually)
    //     node_index_for_leading_top_of_the_group0 = p_cell_population->rGetMesh().GetNodeIndexForLeadingTopOfTheGroup(0);

    // unsigned num_nodes_leading_cell_top = 0;
    // if (mAddPullingForceEvenlyOnNodesOfLeadingCell)
    //     num_nodes_leading_cell_top = p_cell_population->rGetMesh().GetNumNodesOfLeadingCellTopOfGroup(0);

    // Iterate over vertices in the cell population
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);

        /*
         * The force on this Node is given by the gradient of the total free
         * energy of the CellPopulation, evaluated at the position of the vertex. This
         * free energy is the sum of the free energies of all CellPtrs in
         * the cell population. The free energy of each CellPtr is comprised of three
         * parts - a cell deformation energy, a membrane surface tension energy
         * and an adhesion energy.
         *
         * Note that since the movement of this Node only affects the free energy
         * of the CellPtrs containing it, we can just consider the contributions
         * to the free energy gradient from each of these CellPtrs.
         */
        c_vector<double, DIM> deformation_contribution = zero_vector<double>(DIM);
        c_vector<double, DIM> membrane_surface_tension_contribution = zero_vector<double>(DIM);
        c_vector<double, DIM> adhesion_contribution = zero_vector<double>(DIM);

        // Find the indices of the elements owned by this node
        std::set<unsigned> containing_elem_indices = p_cell_population->GetNode(node_index)->rGetContainingElementIndices();

        // Iterate over these elements
        for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
             iter != containing_elem_indices.end();
             ++iter)
        {
            // Get this element, its index and its number of nodes
            VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
            unsigned elem_index = p_element->GetIndex();
            unsigned num_nodes_elem = p_element->GetNumNodes();

            // Find the local index of this node in this element
            unsigned local_index = p_element->GetNodeLocalIndex(node_index);

            // Add the force contribution from this cell's deformation energy (note the minus sign)
            c_vector<double, DIM> element_area_gradient = p_cell_population->rGetMesh().GetAreaGradientOfElementAtNode(p_element, local_index);
            // double area_stiffness = GetNagaiHondaDeformationEnergyParameter()*pow(abs(element_areas[elem_index]-target_areas[elem_index])/target_areas[elem_index],2.0/2.0);            
            // double area_stiffness = GetNagaiHondaDeformationEnergyParameter()*pow(target_areas[elem_index],2.0/2.0); 
            double area_stiffness = GetNagaiHondaDeformationEnergyParameter();
            deformation_contribution -= area_stiffness*(element_areas[elem_index] - target_areas[elem_index])*element_area_gradient;
            
            // std::cout<<"Ka["<<elem_index<<"]="<<area_stiffness<<std::endl;

            // Get the previous and next nodes in this element
            unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
            Node<DIM>* p_previous_node = p_element->GetNode(previous_node_local_index);

            unsigned next_node_local_index = (local_index+1)%num_nodes_elem;
            Node<DIM>* p_next_node = p_element->GetNode(next_node_local_index);

            // Compute the adhesion parameter for each of these edges
            // 0.5 means the edge is shared by two cells
            double previous_edge_adhesion_parameter = GetAdhesionParameter(p_previous_node, p_this_node, *p_cell_population);
            double next_edge_adhesion_parameter = GetAdhesionParameter(p_this_node, p_next_node, *p_cell_population);

            // Compute the gradient of each these edges, computed at the present node
            c_vector<double, DIM> previous_edge_gradient = -p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, previous_node_local_index);
            c_vector<double, DIM> next_edge_gradient = p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, local_index);

            // Add the force contribution from cell-cell and cell-boundary adhesion (note the minus sign)
            adhesion_contribution -= previous_edge_adhesion_parameter*previous_edge_gradient + next_edge_adhesion_parameter*next_edge_gradient;
            
            // Add the force contribution from this cell's membrane surface tension (note the minus sign)
            c_vector<double, DIM> element_perimeter_gradient;
            element_perimeter_gradient = previous_edge_gradient + next_edge_gradient;
            double cell_target_perimeter = 0.0;
            // membrane_surface_tension_contribution -= GetNagaiHondaMembraneSurfaceEnergyParameter()*(element_perimeters[elem_index] - cell_target_perimeter)*element_perimeter_gradient;

            // double perimeter_stiffness = GetNagaiHondaMembraneSurfaceEnergyParameter()*pow(mFixedTargetArea/target_areas[elem_index],1.0/2.0);
            // double perimeter_stiffness = GetNagaiHondaMembraneSurfaceEnergyParameter();
            // double target_element_perimeter_elasticity = GetNagaiHondaMembraneSurfaceEnergyParameter()*pow(mFixedTargetPerimeter/element_perimeters[elem_index],0.5);
            // double perimeter_stiffness = element_perimeter_elasticity[elem_index] + (target_element_perimeter_elasticity-element_perimeter_elasticity[elem_index])*0.01;
            // double perimeter_stiffness = target_element_perimeter_elasticity;

            double perimeter_stiffness = 0.0;
            CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(elem_index);
            if (t_now < this->mTimeForEquilibrium)
            {
                perimeter_stiffness = GetNagaiHondaMembraneSurfaceEnergyParameter()*pow(mFixedTargetPerimeter/element_perimeters[elem_index],1.0);
                p_cell->GetCellData()->SetItem("perimeter_elasticity", perimeter_stiffness);
            }
            else
            {
                perimeter_stiffness = p_cell->GetCellData()->GetItem("perimeter_elasticity");
            }
            membrane_surface_tension_contribution -= perimeter_stiffness*(element_perimeters[elem_index] - cell_target_perimeter)*element_perimeter_gradient;

            // CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(elem_index);
            // p_cell->GetCellData()->SetItem("perimeter_elasticity", perimeter_stiffness);
            p_cell->GetCellData()->SetItem("target_shape_index", -previous_edge_adhesion_parameter/(perimeter_stiffness*sqrt(target_areas[elem_index])));

            // my changes
            // double  rest_edge_length = sqrt(target_areas[elem_index]/(3*sqrt(3)/2));
            // double  rest_edge_length = 0;
            // double  previous_edge_length = p_cell_population->rGetMesh().GetDistanceBetweenNodes(p_element->GetNodeGlobalIndex(previous_node_local_index), node_index);
            // double  next_edge_length = p_cell_population->rGetMesh().GetDistanceBetweenNodes(node_index, p_element->GetNodeGlobalIndex(next_node_local_index));
            // membrane_surface_tension_contribution -= GetNagaiHondaMembraneSurfaceEnergyParameter()*((previous_edge_length - rest_edge_length)*previous_edge_gradient + (next_edge_length - rest_edge_length)*next_edge_gradient);
            // membrane_surface_tension_contribution -= GetNagaiHondaMembraneSurfaceEnergyParameter()*(previous_edge_gradient + next_edge_gradient);

        }// end of 'Iterate over these elements'

        c_vector<double, DIM> force_on_node = deformation_contribution + membrane_surface_tension_contribution + adhesion_contribution;

        // // if (mAddPullingForceOnNodeIndividually && p_this_node->GetIndex()==node_index_for_leading_top_of_the_group0)
        // //     force_on_node[1] += mPullingForceOnLeadingCell;
        // if (mAddPullingForceEvenlyOnNodesOfLeadingCell)
        // {
        //     bool node_belongs_to_leading_cell = false;
        //     for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
        //         iter != containing_elem_indices.end();
        //         ++iter)
        //     {
        //         if (p_cell_population->GetElement(*iter)->GetIsLeadingCellTop())
        //         {
        //             node_belongs_to_leading_cell = true;
        //             break;
        //         }
        //     }

        //     if (node_belongs_to_leading_cell)
        //     {
        //         assert(num_nodes_leading_cell_top!=0);
        //         double t_now = SimulationTime::Instance()->GetTime();
        //         if ( !(mIfEquilibrateForAWhile && t_now<=mTimeForEquilibrium) )
        //             force_on_node[1] += mPullingForceOnLeadingCell/num_nodes_leading_cell_top;
        //     }
        // }

        // // tmp: possilbe output for testing morphogenetic force:
        // bool output_while_testing_morphogenetic_force = false;
        // if (output_while_testing_morphogenetic_force)
        // {
        //     VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(* p_this_node->rGetContainingElementIndices().begin());
        //     if (p_element->GetIsLeadingCell() && p_element->GetGroupNumber()>0)
        //     {
        //         double t_now = SimulationTime::Instance()->GetTime();
        //         if (p_this_node->rGetContainingElementIndices().size()==1)
        //         {
        //             std::cout << "time: " << t_now << std::endl;
        //             std::cout << "node of Leading Cell" << std::endl
        //                     << " node index=" << p_this_node->GetIndex() << std::endl
        //                     << " node location: x=" << p_this_node->rGetLocation()[0] << ", y=" << p_this_node->rGetLocation()[1] << std::endl 
        //                     << " total morphogenetic force =" << mPullingForceOnLeadingCell << std::endl;
        //         }
        //         else if (p_this_node->rGetContainingElementIndices().size()==2)
        //         {
        //             std::cout << "time: " << t_now << std::endl;
        //             std::cout << "node of Leading Cell" << std::endl
        //                     << " node index=" << p_this_node->GetIndex() << std::endl
        //                     << " node location: x=" << p_this_node->rGetLocation()[0] << ", y=" << p_this_node->rGetLocation()[1] << std::endl 
        //                     << " total morphogenetic force =" << mPullingForceOnLeadingCell << std::endl;
        //         }
        //     }
        // }

        p_cell_population->GetNode(node_index)->AddAppliedForceContribution(force_on_node);

    }// end of 'Iterate over nodes(vertices) in the cell population'

}


template<unsigned DIM>
double MyNagaiHondaForce<DIM>::GetAdhesionParameter(Node<DIM>* pNodeA, Node<DIM>* pNodeB, VertexBasedCellPopulation<DIM>& rVertexCellPopulation)
{
    double adhesion_parameter = 0.0;

    // Find the indices of the elements owned by each node
    std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();

    // Find common elements
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_nodeA.begin(),
                          elements_containing_nodeA.end(),
                          elements_containing_nodeB.begin(),
                          elements_containing_nodeB.end(),
                          std::inserter(shared_elements, shared_elements.begin()));

    // Check that the nodes have a common edge
    assert(!shared_elements.empty());

    // If the edge corresponds to a single element, then the cell is on the boundary
    if (shared_elements.size() == 1)
    {
        adhesion_parameter = GetNagaiHondaCellBoundaryAdhesionEnergyParameter();
    //    return adhesion_parameter;
    }
    // // my changes
    // if (mIfUseFaceElementToGetAdhesionParameter)
    // {
    //     unsigned element_global_index = *shared_elements.begin();
    //     VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rVertexCellPopulation);
    //     VertexElement<DIM, DIM>* p_element = p_cell_population->rGetMesh().GetElement(element_global_index);
    //     if ((p_element->GetNodeLocalIndex(pNodeB->GetIndex()) - p_element->GetNodeLocalIndex(pNodeA->GetIndex())+p_element->GetNumNodes())%p_element->GetNumNodes()==1)
    //     {
    //         VertexElement<DIM-1, DIM>* pFace = p_element->GetFace(p_element->GetFaceLocalIndexUsingStartAndEndNodeGlobalIndex(pNodeA->GetIndex(), pNodeB->GetIndex()));
    //         adhesion_parameter =GetNagaiHondaCellCellAdhesionEnergyParameter() * pFace->GetUnifiedCellCellAdhesionEnergyParameter();
    //     }
    //     else
    //     {
    //         if ( (p_element->GetNodeLocalIndex(pNodeA->GetIndex()) - p_element->GetNodeLocalIndex(pNodeB->GetIndex())+p_element->GetNumNodes())%p_element->GetNumNodes()!=1)
    //         {
    //             std::cout << std::endl << "ERROR: Method MyNagaiHondaForce::GetAdhesionParameter";
    //             std::cout << std::endl << "Not Reachable!" << std::endl;
    //         }
    //         assert((p_element->GetNodeLocalIndex(pNodeA->GetIndex()) - p_element->GetNodeLocalIndex(pNodeB->GetIndex())+p_element->GetNumNodes())%p_element->GetNumNodes()==1);
    //         VertexElement<DIM-1, DIM>* pFace = p_element->GetFace(p_element->GetFaceLocalIndexUsingStartAndEndNodeGlobalIndex(pNodeB->GetIndex(), pNodeA->GetIndex()));
    //         adhesion_parameter =GetNagaiHondaCellCellAdhesionEnergyParameter() * pFace->GetUnifiedCellCellAdhesionEnergyParameter();    
    //     }
    // }
    else
        adhesion_parameter = GetNagaiHondaCellCellAdhesionEnergyParameter();

    return adhesion_parameter;
}

template<unsigned DIM>
double MyNagaiHondaForce<DIM>::GetNagaiHondaDeformationEnergyParameter()
{
    return mNagaiHondaDeformationEnergyParameter;
}

template<unsigned DIM>
double MyNagaiHondaForce<DIM>::GetNagaiHondaMembraneSurfaceEnergyParameter()
{
    return mNagaiHondaMembraneSurfaceEnergyParameter;
}

template<unsigned DIM>
double MyNagaiHondaForce<DIM>::GetNagaiHondaCellCellAdhesionEnergyParameter()
{
    return mNagaiHondaCellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
double MyNagaiHondaForce<DIM>::GetNagaiHondaCellBoundaryAdhesionEnergyParameter()
{
    return mNagaiHondaCellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void MyNagaiHondaForce<DIM>::SetNagaiHondaDeformationEnergyParameter(double deformationEnergyParameter)
{
    mNagaiHondaDeformationEnergyParameter = deformationEnergyParameter;
}

template<unsigned DIM>
void MyNagaiHondaForce<DIM>::SetNagaiHondaMembraneSurfaceEnergyParameter(double membraneSurfaceEnergyParameter)
{
    mNagaiHondaMembraneSurfaceEnergyParameter = membraneSurfaceEnergyParameter;
}

template<unsigned DIM>
void MyNagaiHondaForce<DIM>::SetNagaiHondaCellCellAdhesionEnergyParameter(double cellCellAdhesionEnergyParameter)
{
    mNagaiHondaCellCellAdhesionEnergyParameter = cellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
void MyNagaiHondaForce<DIM>::SetNagaiHondaCellBoundaryAdhesionEnergyParameter(double cellBoundaryAdhesionEnergyParameter)
{
    mNagaiHondaCellBoundaryAdhesionEnergyParameter = cellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void MyNagaiHondaForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<NagaiHondaDeformationEnergyParameter>" << mNagaiHondaDeformationEnergyParameter << "</NagaiHondaDeformationEnergyParameter>\n";
    *rParamsFile << "\t\t\t<NagaiHondaMembraneSurfaceEnergyParameter>" << mNagaiHondaMembraneSurfaceEnergyParameter << "</NagaiHondaMembraneSurfaceEnergyParameter>\n";
    *rParamsFile << "\t\t\t<NagaiHondaCellCellAdhesionEnergyParameter>" << mNagaiHondaCellCellAdhesionEnergyParameter << "</NagaiHondaCellCellAdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<NagaiHondaCellBoundaryAdhesionEnergyParameter>" << mNagaiHondaCellBoundaryAdhesionEnergyParameter << "</NagaiHondaCellBoundaryAdhesionEnergyParameter>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class MyNagaiHondaForce<1>;
template class MyNagaiHondaForce<2>;
template class MyNagaiHondaForce<3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MyNagaiHondaForce)
