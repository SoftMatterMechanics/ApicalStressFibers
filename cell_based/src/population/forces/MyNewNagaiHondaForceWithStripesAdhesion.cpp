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

#include "MyNewNagaiHondaForceWithStripesAdhesion.hpp"
#include "MyXToroidal2dVertexMesh.hpp"

template<unsigned DIM>
MyNewNagaiHondaForceWithStripesAdhesion<DIM>::MyNewNagaiHondaForceWithStripesAdhesion()
   : AbstractForce<DIM>(),
     mNagaiHondaDeformationEnergyParameter(1.0),
     mNagaiHondaMembraneSurfaceEnergyParameter(0.1),
     mNagaiHondaCellCellAdhesionEnergyParameter(0.0),
     mNagaiHondaCellBoundaryAdhesionEnergyParameter(0.0),
     mUseFixedTargetArea(true),
     mCaseNumberOfMembraneSurfaceEnergyForm(0),
     mIfUseFaceElementToGetAdhesionParameter(false),
     mOutputInformationForNagaiHondaForce(false)
{
}

template<unsigned DIM>
MyNewNagaiHondaForceWithStripesAdhesion<DIM>::~MyNewNagaiHondaForceWithStripesAdhesion()
{
}

template<unsigned DIM>
void MyNewNagaiHondaForceWithStripesAdhesion<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{    
    // Throw an exception message if not using a VertexBasedCellPopulation
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("NagaiHondaForce is to be used with a VertexBasedCellPopulation only");
    }

    // Define some helper variables
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    unsigned num_nodes = p_cell_population->GetNumNodes();
    unsigned num_elements = p_cell_population->GetNumElements();

    // my changes
    if (num_elements<20)
        return;

    // Begin by computing the area and perimeter of each element in the mesh, to avoid having to do this multiple times
    std::vector<double> element_areas(num_elements);
    std::vector<double> element_perimeters(num_elements);
    std::vector<double> target_areas(num_elements);
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
                target_areas[elem_index] = p_cell_population->GetCellUsingLocationIndex(elem_index)->GetCellData()->GetItem("target area");

        }
        catch (Exception&)
        {
            EXCEPTION("You need to add an AbstractTargetAreaModifier to the simulation in order to use NagaiHondaForce");
        }

    }

    // original non-homogeneous SSA method using leading length:
    double y_coord_leading_edge = 0.0;
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);
        if (p_this_node->rGetLocation()[1] > y_coord_leading_edge)
            y_coord_leading_edge= p_this_node->rGetLocation()[1];
    }

    // For My Detach Pattern Method, using group number:
    unsigned number_of_groups = p_cell_population->rGetMesh().GetNumberOfGroups();
    std::vector<double> leading_tops_of_groups(number_of_groups);
    if (mUseMyDetachPatternMethod)
    {
        for (unsigned i=0; i<number_of_groups; i++)
            leading_tops_of_groups[i] = p_cell_population->rGetMesh().GetLeadingTopOfTheGroup(i);
    }

    unsigned node_index_for_leading_top_of_the_group0 = 0;
    if (mAddPullingForceOnNodeIndividually)
        node_index_for_leading_top_of_the_group0 = p_cell_population->rGetMesh().GetNodeIndexForLeadingTopOfTheGroup(0);

    unsigned num_nodes_leading_cell_top = 0;
    if (mAddPullingForceEvenlyOnNodesOfLeadingCell)
        num_nodes_leading_cell_top = p_cell_population->rGetMesh().GetNumNodesOfLeadingCellTopOfGroup(0);

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
        c_vector<double, DIM> area_adhesion_contribution = zero_vector<double>(DIM);
        bool node_has_SSA = false;
        c_vector<double, DIM> averaged_inner_strip_substrate_adhesion_contribution = zero_vector<double>(DIM);
        c_vector<double, DIM> reservoir_substrate_adhesion_contribution = zero_vector<double>(DIM);

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
            deformation_contribution -= GetNagaiHondaDeformationEnergyParameter()*(element_areas[elem_index] - target_areas[elem_index])*element_area_gradient;
                 
            // Get the previous and next nodes in this element
            unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
            Node<DIM>* p_previous_node = p_element->GetNode(previous_node_local_index);

            unsigned next_node_local_index = (local_index+1)%num_nodes_elem;
            Node<DIM>* p_next_node = p_element->GetNode(next_node_local_index);

            // Compute the adhesion parameter for each of these edges
            // changes to be made:
            double previous_edge_adhesion_parameter = GetAdhesionParameter(p_previous_node, p_this_node, *p_cell_population);
            double next_edge_adhesion_parameter = GetAdhesionParameter(p_this_node, p_next_node, *p_cell_population);

            // Compute the gradient of each these edges, computed at the present node
            c_vector<double, DIM> previous_edge_gradient = -p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, previous_node_local_index);
            c_vector<double, DIM> next_edge_gradient = p_cell_population->rGetMesh().GetNextEdgeGradientOfElementAtNode(p_element, local_index);

            // Add the force contribution from cell-cell and cell-boundary adhesion (note the minus sign)
            adhesion_contribution -= previous_edge_adhesion_parameter*previous_edge_gradient + next_edge_adhesion_parameter*next_edge_gradient;
            
            // Add the force contribution from this cell's membrane surface tension (note the minus sign)
            c_vector<double, DIM> element_perimeter_gradient = previous_edge_gradient + next_edge_gradient;

            double cell_target_perimeter = 0.0;            
            unsigned case_number_of_membrane_surface_energy_form = this->GetCaseNumberOfMembraneSurfaceEnergyForm();
            
            // case 0:typical *membrane surface tension*
            if (case_number_of_membrane_surface_energy_form ==0)
            {
                membrane_surface_tension_contribution -= GetNagaiHondaMembraneSurfaceEnergyParameter()*1.0*(element_perimeters[elem_index] - cell_target_perimeter)*element_perimeter_gradient;
            }
            
            // case 1: MA feedback incorporated in *membrane surface tension*
            if (case_number_of_membrane_surface_energy_form ==1)
            {
                // to do later.
                // myosin_activity = p_cell->GetMyosinActivity();
                // membrane_surface_tension_contribution -= GetNagaiHondaMembraneSurfaceEnergyParameter()*myosin_activity*(element_perimeters[elem_index] - cell_target_perimeter)*element_perimeter_gradient;
            }
            
            // case 2: EMA feedback incorporated in *membrane surface tension*
            if (case_number_of_membrane_surface_energy_form ==2)// New energy form 1
            {
                double sum1 = 0.0;
                for (unsigned i = 0; i < p_element->GetNumNodes(); i++)
                {
                    unsigned this_node_global_index = p_element->GetNodeGlobalIndex(i);
                    unsigned next_node_global_index = p_element->GetNodeGlobalIndex((i + 1) % p_element->GetNumNodes());
                    double edge_length = p_cell_population->rGetMesh().GetDistanceBetweenNodes(this_node_global_index, next_node_global_index);
                    unsigned face_local_index = p_element->GetFaceLocalIndexUsingStartAndEndNodeGlobalIndex(this_node_global_index,next_node_global_index);
                    sum1 += sqrt(p_element->GetFace(face_local_index)->GetUnifiedEdgeMyosinActivty())*edge_length;
                }
                
                unsigned this_node_global_index = p_element->GetNodeGlobalIndex(local_index);
                unsigned next_node_global_index = p_element->GetNodeGlobalIndex((local_index+1)%p_element->GetNumNodes());
                unsigned previous_node_global_index = p_element->GetNodeGlobalIndex((local_index-1+p_element->GetNumNodes())%p_element->GetNumNodes());

                unsigned previous_face_local_index = p_element->GetFaceLocalIndexUsingStartAndEndNodeGlobalIndex(previous_node_global_index ,this_node_global_index);
                unsigned next_face_local_index = p_element->GetFaceLocalIndexUsingStartAndEndNodeGlobalIndex(this_node_global_index,next_node_global_index);
                c_vector<double, DIM> sum2 = sqrt(p_element->GetFace(previous_face_local_index)->GetUnifiedEdgeMyosinActivty())
                        *previous_edge_gradient + sqrt(p_element->GetFace(next_face_local_index)->GetUnifiedEdgeMyosinActivty())*next_edge_gradient;
                membrane_surface_tension_contribution -= GetNagaiHondaMembraneSurfaceEnergyParameter()*sum1*sum2;                
            }

            // case 3: EMA feedback incorporated in *membrane surface tension*
            if (case_number_of_membrane_surface_energy_form == 3)// New energy form 2
            {
                unsigned this_node_global_index = p_element->GetNodeGlobalIndex(local_index);
                unsigned next_node_global_index = p_element->GetNodeGlobalIndex((local_index+1)%p_element->GetNumNodes());
                unsigned previous_node_global_index = p_element->GetNodeGlobalIndex((local_index-1+p_element->GetNumNodes())%p_element->GetNumNodes());

                unsigned previous_face_local_index = p_element->GetFaceLocalIndexUsingStartAndEndNodeGlobalIndex(previous_node_global_index ,this_node_global_index);
                unsigned next_face_local_index = p_element->GetFaceLocalIndexUsingStartAndEndNodeGlobalIndex(this_node_global_index,next_node_global_index);

                double unified_edge_myosin_activity_of_previous_edge = p_element->GetFace(previous_face_local_index)->GetUnifiedEdgeMyosinActivty();
                double unified_edge_myosin_activity_of_next_edge= p_element->GetFace(next_face_local_index)->GetUnifiedEdgeMyosinActivty();

                double previous_edge_length = p_cell_population->rGetMesh().GetDistanceBetweenNodes(previous_node_global_index, this_node_global_index);
                double next_edge_length = p_cell_population->rGetMesh().GetDistanceBetweenNodes(this_node_global_index, next_node_global_index);

                c_vector<double, DIM> sum = unified_edge_myosin_activity_of_previous_edge*previous_edge_length*previous_edge_gradient \
                    + unified_edge_myosin_activity_of_next_edge*next_edge_length*next_edge_gradient;
                membrane_surface_tension_contribution -= 6*GetNagaiHondaMembraneSurfaceEnergyParameter()*sum;
            }

            /*----------------------------------------------------------------------------------------------------------*/
            /*---------------------------------Start of substrate adhesion contribution---------------------------------*/
            if (mIfConsiderSubstrateAdhesion)
            {
                assert( !((mIfUseNewSSADistributionRule||mUseMyDetachPatternMethod)&&mIfSubstrateAdhesionIsHomogeneous) );
                c_vector<double, DIM> substrate_adhesion_area_gradient = zero_vector<double>(DIM);
                c_vector<double, DIM> weighted_substrate_adhesion_area_gradient = zero_vector<double>(DIM);
                // Parameters
                double stripe_width = this->mStripWidth;
                double stripe_distance = this->mStripDistance;
                double strip_start_x_location = this->mStripStartXLocation;
                double strip_start_y_location = this->mStripStartYLocation;
                bool if_substrate_adhesion_is_homogeneous = this->mIfSubstrateAdhesionIsHomogeneous;
                double substrate_adhesion_leading_top_length = this->mSubstrateAdhesionLeadingTopLength;

                double small_change = mSmallChangeForAreaCalculation;
                double preferred_sample_dis = small_change/5.0;

                c_vector<double, DIM> previous_node_location = p_previous_node->rGetLocation();
                c_vector<double, DIM> this_node_location = p_this_node->rGetLocation();
                c_vector<double, DIM> next_node_location = p_next_node->rGetLocation();
                // consider periodicity: modify the node location!
                if (dynamic_cast<MyXToroidal2dVertexMesh*>(&p_cell_population->rGetMesh()) != nullptr)
                {
                    bool triangle_straddles_left_right_boundary = false;

                    c_vector<double, 2> vector1 = this_node_location - previous_node_location;
                    if (fabs(vector1[0]) > 0.5*mWidth)
                        triangle_straddles_left_right_boundary = true;
                    c_vector<double, 2> vector2 = next_node_location - this_node_location;
                    if (fabs(vector2[0]) > 0.5*mWidth)
                        triangle_straddles_left_right_boundary = true;
                    if (triangle_straddles_left_right_boundary)
                    {
                        if (previous_node_location[0] < mCenterOfWidth)
                            previous_node_location[0] += mWidth;
                        if (this_node_location[0] < mCenterOfWidth)
                            this_node_location[0] += mWidth;
                        if (next_node_location[0] < mCenterOfWidth)
                            next_node_location[0] += mWidth;
                    }
                }

                c_vector<double, DIM> centroid_triangle = 1/3.0*(previous_node_location+this_node_location+next_node_location);
                int stripe_num = (int) round((centroid_triangle[0] - strip_start_x_location)/stripe_distance);

                double stripe_location = strip_start_x_location+ stripe_distance*stripe_num;
                double stripe_left = stripe_location - stripe_width/2;
                double stripe_right = stripe_location + stripe_width/2;
                double left_stripe_right = strip_start_x_location+ stripe_distance*(stripe_num-1) + stripe_width/2;
                double right_stripe_left = strip_start_x_location+ stripe_distance*(stripe_num+1) - stripe_width/2;
                double sample_area_bottom = 0.0;
                double sample_area_top = 0.0;
                double sample_area_left = 0.0;
                double sample_area_right = 0.0;

                double expanded_triangle_box_bottom = -small_change+std::min(std::min(previous_node_location[1],this_node_location[1]),next_node_location[1]);
                double expanded_triangle_box_top = small_change+std::max(std::max(previous_node_location[1],this_node_location[1]),next_node_location[1]);
                double expanded_triangle_box_left = -small_change+std::min(std::min(previous_node_location[0],this_node_location[0]),next_node_location[0]);
                double expanded_triangle_box_right = small_change+std::max(std::max(previous_node_location[0],this_node_location[0]),next_node_location[0]);
                
                bool has_substrate_adhesion_area = true;

                if (expanded_triangle_box_top<strip_start_y_location)
                    has_substrate_adhesion_area = false;
                else
                {
                    sample_area_top = expanded_triangle_box_top;
                    sample_area_bottom = expanded_triangle_box_bottom>strip_start_y_location ? expanded_triangle_box_bottom : strip_start_y_location;
                }

                // codes follows only suit for the case that the triangle laps over only one strip!
                if (expanded_triangle_box_right<stripe_left)
                    has_substrate_adhesion_area = false;
                else if (expanded_triangle_box_left<stripe_left)
                {
                    sample_area_left = stripe_left;
                    if (expanded_triangle_box_right<stripe_right)
                        sample_area_right = expanded_triangle_box_right;
                    else
                        sample_area_right = stripe_right;
                }
                else if (expanded_triangle_box_right<stripe_right)
                {
                    sample_area_left = expanded_triangle_box_left;
                    sample_area_right = expanded_triangle_box_right;
                }
                else if (expanded_triangle_box_left<stripe_right)
                {
                    sample_area_left = expanded_triangle_box_left;
                    sample_area_right = stripe_right;
                }
                else
                    has_substrate_adhesion_area = false;

                // particular case:
                bool if_in_a_particular_case1 = false;
                bool if_in_a_particular_case2 = false;
                bool if_in_a_particular_case = false;
                double strip_interval_left = 0.0;
                double strip_interval_right = 0.0;
                if (expanded_triangle_box_top>=strip_start_y_location && expanded_triangle_box_right>stripe_left && expanded_triangle_box_left<left_stripe_right)
                    if_in_a_particular_case1 = true;
                if (expanded_triangle_box_top>=strip_start_y_location && expanded_triangle_box_left<stripe_right && expanded_triangle_box_right>right_stripe_left)
                    if_in_a_particular_case2 = true;
                if (if_in_a_particular_case1||if_in_a_particular_case2)
                {
                    if_in_a_particular_case = true;
                    has_substrate_adhesion_area = true;
                    sample_area_left = expanded_triangle_box_left;
                    sample_area_right = expanded_triangle_box_right;
                    if (if_in_a_particular_case1)
                    {
                        strip_interval_left = left_stripe_right;
                        strip_interval_right = stripe_left;
                    }
                    if (if_in_a_particular_case2)
                    {
                        assert(!if_in_a_particular_case1);
                        strip_interval_left = stripe_right;
                        strip_interval_right = right_stripe_left;
                    }

                }

                unsigned num_across = 0;
                unsigned num_up = 0;
                unsigned sample_num = 0;
                double sample_area = 0.0;

                if (has_substrate_adhesion_area)
                {
                    num_across = round((sample_area_right - sample_area_left) / preferred_sample_dis);
                    num_up = round((sample_area_top - sample_area_bottom) / preferred_sample_dis);
                    sample_num = num_across * num_up;
                    sample_area = (sample_area_top - sample_area_bottom)*(sample_area_right - sample_area_left);
                }

                bool if_node_a_inner_node = false;
                if (containing_elem_indices.size()==3)
                    if_node_a_inner_node = true;
                
                // Calculate (weighted_)substrate_adhesion_area_gradient!
                if (!if_node_a_inner_node && has_substrate_adhesion_area && (num_across>0) && (num_up>0))
                {
                    node_has_SSA = true;
                    c_vector<c_vector<double, DIM>, 3> points;
                    points[0]=previous_node_location;
                    points[1]=this_node_location;
                    points[2]=next_node_location;
                    c_vector<double, DIM> vec1 = this_node_location-previous_node_location;
                    c_vector<double, DIM> vec2 = next_node_location-this_node_location;                    

                    unsigned group_number_of_node = p_element->GetGroupNumber();
                    double leading_top_of_the_group = leading_tops_of_groups[group_number_of_node];

                    unsigned case_number = 0;
                    double SSA = 0.0; // for case 0
                    double SSA_above_this_node = 0.0; // for case 1
                    double SSA_below_this_node = 0.0; // for case 1
                    if (mIfUseNewSSADistributionRule)
                    {
                        if (p_this_node->rGetContainingElementIndices().size()==1)
                        {
                            case_number = 0;
                            CellPtr p_cell = p_cell_population->GetCellUsingLocationIndex(p_element->GetIndex());
                            if (p_element->GetIsLeadingCell() || p_element->GetIsJustReAttached())
                            {
                                SSA = mBasicSSA + (mSSAForMatureLamellipodium - mBasicSSA)*p_element->GetLamellipodiumStrength();
                                if (mSmallSSAAtFirst && SimulationTime::Instance()->GetTime()<mInitialTimeForSmallSSA)
                                {
                                    SSA = mBasicSSA + (mSmallSSAForInitialTime - mBasicSSA)*p_element->GetLamellipodiumStrength();
                                }
                                if(mKeepMovingForward)
                                {
                                    if (p_element->GetIsLeadingCellBottom())
                                    {
                                        if (p_this_node->rGetLocation()[1]>mSlowlyMovingForwardAfterThisHeight)
                                            SSA = std::max(mBasicSSA, SSA-0.5);
                                        else
                                            SSA = std::max(mBasicSSA, SSA-mSSABottomDecrease);
                                    }
                                }
                            }
                            else
                            {
                                assert(!p_element->GetIsLeadingCell() && !p_element->GetIsJustReAttached());
                                SSA = mBasicSSA;
                            }
                        }
                        else if (p_this_node->rGetContainingElementIndices().size()==2)
                        {
                            case_number = 1;
                            VertexElement<DIM, DIM>* p_element_1 = p_cell_population->GetElement(* p_this_node->rGetContainingElementIndices().begin());
                            CellPtr p_cell_1 = p_cell_population->GetCellUsingLocationIndex(p_element_1->GetIndex());
                            VertexElement<DIM, DIM>* p_element_2 = p_cell_population->GetElement(* p_this_node->rGetContainingElementIndices().begin()++);
                            CellPtr p_cell_2 = p_cell_population->GetCellUsingLocationIndex(p_element_2->GetIndex());
                            double centroid_y_cell_1 = p_cell_population->rGetMesh().GetCentroidOfElement(p_element_1->GetIndex())[1];
                            double centroid_y_cell_2 = p_cell_population->rGetMesh().GetCentroidOfElement(p_element_2->GetIndex())[1];
                            double SSA_cell_1 = 0.0;
                            double SSA_cell_2 = 0.0;
                            if (p_element_1->GetIsLeadingCell() || p_element_1->GetIsJustReAttached())
                            {
                                SSA_cell_1 = mBasicSSA + (mSSAForMatureLamellipodium - mBasicSSA)*p_element_1->GetLamellipodiumStrength();
                                if (mSmallSSAAtFirst && SimulationTime::Instance()->GetTime()<mInitialTimeForSmallSSA)
                                {
                                    SSA_cell_1 = mBasicSSA + (mSmallSSAForInitialTime - mBasicSSA)*p_element_1->GetLamellipodiumStrength();
                                }

                            }
                            else
                            {
                                assert(!p_element_1->GetIsLeadingCell() && !p_element_1->GetIsJustReAttached());
                                SSA_cell_1 = mBasicSSA;
                            }
                            if (p_element_2->GetIsLeadingCell() || p_element_2->GetIsJustReAttached())
                            {
                                SSA_cell_2 = mBasicSSA + (mSSAForMatureLamellipodium - mBasicSSA)*p_element_2->GetLamellipodiumStrength();
                                if (mSmallSSAAtFirst && SimulationTime::Instance()->GetTime()<mInitialTimeForSmallSSA)
                                {
                                    SSA_cell_2 = mBasicSSA + (mSmallSSAForInitialTime - mBasicSSA)*p_element_2->GetLamellipodiumStrength();
                                }

                            }
                            else
                            {
                                assert(!p_element_2->GetIsLeadingCell() && !p_element_2->GetIsJustReAttached());
                                SSA_cell_2 = mBasicSSA;
                            }

                            if (centroid_y_cell_1 > centroid_y_cell_2)
                            {
                                SSA_above_this_node = SSA_cell_1;
                                SSA_below_this_node = SSA_cell_2;
                            }
                            else
                            {
                                SSA_above_this_node = SSA_cell_2;
                                SSA_below_this_node = SSA_cell_1;
                            }
                        }
                    }

                    // Calculate initial adhesive area
                    double adhesive_sample_num = 0.0;
                    double weighted_adhesive_sample_num = 0.0;
                    for (unsigned i = 0; i <sample_num; i++)
                    {
                        double x_coord = sample_area_left + (sample_area_right - sample_area_left)/num_across*(i%num_across+0.5);
                        double y_coord = sample_area_bottom + (sample_area_top - sample_area_bottom)/num_up*(i/num_across+0.5);

                        c_vector<double, DIM> point;
                        point[0] = x_coord;
                        point[1] = y_coord;
                        c_vector<bool, 3> point_at_left_of_vector;
                        for (unsigned j = 0; j < 3; j++)
                        {
                            c_vector<double, DIM> vec1 = point - points[j];
                            c_vector<double, DIM> vec2 = points[(j+1)%3] - points[j];

                            if ( (vec1[0]*vec2[1]-vec1[1]*vec2[0]) > 0.0)
                                point_at_left_of_vector[j] = false;
                            else
                                point_at_left_of_vector[j] = true;
                        }
                        if (point_at_left_of_vector[0]==true && point_at_left_of_vector[1]==true && point_at_left_of_vector[2]==true)
                        {
                            if (!if_in_a_particular_case)
                                adhesive_sample_num += 1.0;
                            else // in a particular case!
                            {
                                if (point[0]<strip_interval_left || point[0]>strip_interval_right)
                                    adhesive_sample_num += 1.0;
                            }
                            
                            if (mIfUseNewSSADistributionRule)
                            {
                                assert(!mUseMyDetachPatternMethod);
                                if (case_number==0)
                                    weighted_adhesive_sample_num += 1.0*SSA;
                                else
                                {
                                    assert(case_number==1);
                                    if (y_coord >= p_this_node->rGetLocation()[1])
                                        weighted_adhesive_sample_num += 1.0*SSA_above_this_node;
                                    else
                                        weighted_adhesive_sample_num += 1.0*SSA_below_this_node;
                                }
                            }
                            else if (mUseMyDetachPatternMethod)
                            {
                                if (y_coord > (leading_top_of_the_group - substrate_adhesion_leading_top_length))
                                {
                                    weighted_adhesive_sample_num += 1.0*mSSAForMatureLamellipodium;
                                }
                                else
                                    weighted_adhesive_sample_num += 1.0*mBasicSSA;
                            }
                            else
                            {
                                if (y_coord > (y_coord_leading_edge - substrate_adhesion_leading_top_length))
                                {
                                    weighted_adhesive_sample_num += 1.0*mSSAForMatureLamellipodium;
                                }
                                else
                                    weighted_adhesive_sample_num += 1.0*mBasicSSA;
                            }
                            
                        }
                        else if (point_at_left_of_vector[0]==false && point_at_left_of_vector[1]==false && point_at_left_of_vector[2]==false)
                        {
                            if (!if_in_a_particular_case)
                                adhesive_sample_num += -1.0;
                            else // in a particular case!
                            {
                                if (point[0]<strip_interval_left || point[0]>strip_interval_right)
                                    adhesive_sample_num += -1.0;
                            }
                            

                            if (mIfUseNewSSADistributionRule)
                            {
                                assert(!mUseMyDetachPatternMethod);
                                if (case_number==0)
                                    weighted_adhesive_sample_num += -1.0*SSA;
                                else
                                {
                                    assert(case_number==1);
                                    if (y_coord >= p_this_node->rGetLocation()[1])
                                        weighted_adhesive_sample_num += -1.0*SSA_above_this_node;
                                    else
                                        weighted_adhesive_sample_num += -1.0*SSA_below_this_node;
                                }
                            }
                            else if (mUseMyDetachPatternMethod)
                            {
                                if (y_coord > (leading_top_of_the_group - substrate_adhesion_leading_top_length))
                                {
                                    weighted_adhesive_sample_num += -1.0*mSSAForMatureLamellipodium;
                                }
                                else
                                    weighted_adhesive_sample_num += -1.0*mBasicSSA;
                            }
                            else
                            {
                                if (y_coord > (y_coord_leading_edge - substrate_adhesion_leading_top_length))
                                {
                                    weighted_adhesive_sample_num += -1.0*mSSAForMatureLamellipodium;
                                }
                                else
                                    weighted_adhesive_sample_num += -1.0*mBasicSSA;
                            }

                        }
                    }
                    double substrate_adhesion_area = adhesive_sample_num/double(sample_num) * sample_area;
                    double weighted_substrate_adhesion_area = weighted_adhesive_sample_num/double(sample_num) * sample_area;

                    // Calculate adhesive area *changed* with small displacement of the node along the x or y axis
                    for (unsigned j = 0; j<2; j++)
                    {
                        c_vector<double, DIM> unit_vector_in_small_change_direction;
                        unit_vector_in_small_change_direction[0] = cos(j*M_PI/2);
                        unit_vector_in_small_change_direction[1] = sin(j*M_PI/2);
                        // my new changes for cosistent movement of nodes at edges of the strip.
                        if (mConsiderConsistencyForSSA)
                        {
                            if (points[1][0]-stripe_right>-small_change && points[1][0]-stripe_right<small_change)
                                unit_vector_in_small_change_direction[0] *= -1.0;
                        }

                        c_vector<c_vector<double, DIM>, 3> new_points = points;
                        new_points[1] += small_change*unit_vector_in_small_change_direction;

                        double adhesive_sample_num = 0.0;
                        double weighted_adhesive_sample_num = 0.0;

                        for (unsigned i = 0; i <sample_num; i++)
                        {
                            //double x_coord = sample_area_left + (sample_area_right - sample_area_left)*rand()/double(RAND_MAX);
                            //double y_coord = sample_area_bottom + (sample_area_top - sample_area_bottom)*rand()/double(RAND_MAX);
                            double x_coord = sample_area_left + (sample_area_right - sample_area_left)/num_across*(i%num_across+0.5);
                            double y_coord = sample_area_bottom + (sample_area_top - sample_area_bottom)/num_up*(i/num_across+0.5);

                            c_vector<double, DIM> point;
                            point[0] = x_coord;
                            point[1] = y_coord;
                            // bool point_is_outside = false;
                            c_vector<bool, 3> point_at_left_of_vector;
                            for (unsigned j = 0; j < 3; j++)
                            {
                                c_vector<double, DIM> vec1 = point - new_points[j];
                                c_vector<double, DIM> vec2 = new_points[(j+1)%3] - new_points[j];

                                if ( (vec1[0]*vec2[1]-vec1[1]*vec2[0]) > 0.0)
                                    point_at_left_of_vector[j] = false;
                                else
                                    point_at_left_of_vector[j] = true;
                            }
                            if (point_at_left_of_vector[0]==true && point_at_left_of_vector[1]==true && point_at_left_of_vector[2]==true)
                            {
                                if (!if_in_a_particular_case)
                                    adhesive_sample_num += 1.0;
                                else // in a particular case!
                                {
                                    if (point[0]<strip_interval_left || point[0]>strip_interval_right)
                                        adhesive_sample_num += 1.0;
                                }
                            
                                if (mIfUseNewSSADistributionRule)
                                {
                                    assert(!mUseMyDetachPatternMethod);
                                    if (case_number==0)
                                        weighted_adhesive_sample_num += 1.0*SSA;
                                    else
                                    {
                                        assert(case_number==1);
                                        if (y_coord >= p_this_node->rGetLocation()[1])
                                            weighted_adhesive_sample_num += 1.0*SSA_above_this_node;
                                        else
                                            weighted_adhesive_sample_num += 1.0*SSA_below_this_node;
                                    }
                                }
                                else if (mUseMyDetachPatternMethod)
                                {
                                    if (y_coord > (leading_top_of_the_group - substrate_adhesion_leading_top_length))
                                    {
                                        weighted_adhesive_sample_num += 1.0*mSSAForMatureLamellipodium;
                                    }
                                    else
                                        weighted_adhesive_sample_num += 1.0*mBasicSSA;
                                }
                                else
                                {
                                    if (y_coord > (y_coord_leading_edge - substrate_adhesion_leading_top_length))
                                    {
                                        weighted_adhesive_sample_num += 1.0*mSSAForMatureLamellipodium;
                                    }
                                    else
                                        weighted_adhesive_sample_num += 1.0*mBasicSSA;
                                }

                            }
                            else if (point_at_left_of_vector[0]==false && point_at_left_of_vector[1]==false && point_at_left_of_vector[2]==false)
                            {
                                if (!if_in_a_particular_case)
                                    adhesive_sample_num += -1.0;
                                else // in a particular case!
                                {
                                    if (point[0]<strip_interval_left || point[0]>strip_interval_right)
                                        adhesive_sample_num += -1.0;
                                }
                                
                                if (mIfUseNewSSADistributionRule)
                                {
                                    assert(!mUseMyDetachPatternMethod);
                                    if (case_number==0)
                                        weighted_adhesive_sample_num += -1.0*SSA;
                                    else
                                    {
                                        assert(case_number==1);
                                        if (y_coord >= p_this_node->rGetLocation()[1])
                                            weighted_adhesive_sample_num += -1.0*SSA_above_this_node;
                                        else
                                            weighted_adhesive_sample_num += -1.0*SSA_below_this_node;
                                    }
                                }
                                else if (mUseMyDetachPatternMethod)
                                {
                                    if (y_coord > (leading_top_of_the_group - substrate_adhesion_leading_top_length))
                                    {
                                        weighted_adhesive_sample_num += -1.0*mSSAForMatureLamellipodium;
                                    }
                                    else
                                        weighted_adhesive_sample_num += -1.0*mBasicSSA;
                                }
                                else
                                {
                                    if (y_coord > (y_coord_leading_edge - substrate_adhesion_leading_top_length))
                                    {
                                        weighted_adhesive_sample_num += -1.0*mSSAForMatureLamellipodium;
                                    }
                                    else
                                        weighted_adhesive_sample_num += -1.0*mBasicSSA;
                                }

                            }
                        }
                        double substrate_adhesion_area_new = adhesive_sample_num/double(sample_num) * sample_area;
                        double weighted_substrate_adhesion_area_new = weighted_adhesive_sample_num/double(sample_num) * sample_area;
                        
                        substrate_adhesion_area_gradient += (substrate_adhesion_area_new - substrate_adhesion_area)/small_change*unit_vector_in_small_change_direction;
                        weighted_substrate_adhesion_area_gradient += (weighted_substrate_adhesion_area_new - weighted_substrate_adhesion_area)/small_change*unit_vector_in_small_change_direction;
                            
                    }// end of calculation of gradient after considering small displacement of the node along the x or y axis

                }// end of statement 'if (has_substrate_adhesion_area)'

                if (if_substrate_adhesion_is_homogeneous == true)
                    area_adhesion_contribution -= mHomogeneousSubstrateAdhesionParameter*substrate_adhesion_area_gradient;
                else if (!mSSAStrengthenedOnlyInYDirection)
                    area_adhesion_contribution -= weighted_substrate_adhesion_area_gradient;
                else
                {
                    area_adhesion_contribution[0] -= mBasicSSA*substrate_adhesion_area_gradient[0];
                    area_adhesion_contribution[1] -= weighted_substrate_adhesion_area_gradient[1];
                }
                
                // If consider reservior substrate adhesion:
                if (mIfConsiderReservoirSubstrateAdhesion)
                {
                    // Attention
                    c_vector<double, DIM> reservoir_substrate_adhesion_area_gradient = zero_vector<double>(DIM);
                    // Parameters
                    double reservoir_end_y_location = this->mStripStartYLocation;

                    double small_change = mSmallChangeForAreaCalculation;
                    double preferred_sample_dis = small_change/5.0;

                    double sample_area_bottom = 0.0;
                    double sample_area_top = 0.0;
                    double sample_area_left = 0.0;
                    double sample_area_right = 0.0;

                    c_vector<double, DIM> previous_node_location = p_previous_node->rGetLocation();
                    c_vector<double, DIM> this_node_location = p_this_node->rGetLocation();
                    c_vector<double, DIM> next_node_location = p_next_node->rGetLocation();
                    // consider periodicity: modify the node location!
                    if (dynamic_cast<MyXToroidal2dVertexMesh*>(&p_cell_population->rGetMesh()) != nullptr)
                    {
                        bool triangle_straddles_left_right_boundary = false;

                        c_vector<double, 2> vector1 = this_node_location - previous_node_location;
                        if (fabs(vector1[0]) > 0.5*mWidth)
                            triangle_straddles_left_right_boundary = true;
                        c_vector<double, 2> vector2 = next_node_location - this_node_location;
                        if (fabs(vector2[0]) > 0.5*mWidth)
                            triangle_straddles_left_right_boundary = true;
                        if (triangle_straddles_left_right_boundary)
                        {
                            if (previous_node_location[0] < mCenterOfWidth)
                                previous_node_location[0] += mWidth;
                            if (this_node_location[0] < mCenterOfWidth)
                                this_node_location[0] += mWidth;
                            if (next_node_location[0] < mCenterOfWidth)
                                next_node_location[0] += mWidth;
                        }
                    }
                    
                    double expanded_triangle_box_bottom = -small_change+std::min(std::min(previous_node_location[1],this_node_location[1]),next_node_location[1]);
                    double expanded_triangle_box_top = small_change+std::max(std::max(previous_node_location[1],this_node_location[1]),next_node_location[1]);
                    double expanded_triangle_box_left = -small_change+std::min(std::min(previous_node_location[0],this_node_location[0]),next_node_location[0]);
                    double expanded_triangle_box_right = small_change+std::max(std::max(previous_node_location[0],this_node_location[0]),next_node_location[0]);
                    
                    // bool if_in_the_inner_of_reservoir = false;
                    // {
                    //     double top_y_location_of_these_nodes = 0.0;
                    //     if (containing_elem_indices.size() == 3)
                    //     {
                    //         for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
                    //             iter != containing_elem_indices.end();
                    //             ++iter)
                    //         {
                    //             VertexElement<DIM, DIM>* p_element = p_cell_population->GetElement(*iter);
                    //             unsigned num_nodes_elem = p_element->GetNumNodes();
                    //             unsigned local_index = p_element->GetNodeLocalIndex(node_index);
                    //             unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
                    //             Node<DIM>* p_previous_node = p_element->GetNode(previous_node_local_index);
                    //             unsigned next_node_local_index = (local_index+1)%num_nodes_elem;
                    //             Node<DIM>* p_next_node = p_element->GetNode(next_node_local_index);
                    //             c_vector<double, DIM> previous_node_location = p_previous_node->rGetLocation();
                    //             c_vector<double, DIM> this_node_location = p_this_node->rGetLocation();
                    //             c_vector<double, DIM> next_node_location = p_next_node->rGetLocation();
                                
                    //             double expanded_triangle_box_top = small_change+std::max(std::max(previous_node_location[1],this_node_location[1]),next_node_location[1]);
                    //             top_y_location_of_these_nodes = std::max(top_y_location_of_these_nodes, expanded_triangle_box_top);
                    //         }
                    //         if (top_y_location_of_these_nodes < reservoir_end_y_location)
                    //             if_in_the_inner_of_reservoir = true;
                    //     }
                    // }

                    bool has_substrate_adhesion_area = true;
                    if (expanded_triangle_box_bottom>reservoir_end_y_location)
                        has_substrate_adhesion_area = false;
                    else if (expanded_triangle_box_top <= 0.0) // note: we get errors here previously.
                        has_substrate_adhesion_area = false;
                    else                    
                    {
                        sample_area_top = expanded_triangle_box_top < reservoir_end_y_location ? expanded_triangle_box_top : reservoir_end_y_location;
                        sample_area_bottom = std::max(0.0, expanded_triangle_box_bottom);
                    }
                    sample_area_left = expanded_triangle_box_left;
                    sample_area_right = expanded_triangle_box_right;

                    unsigned num_across = 0;
                    unsigned num_up = 0;
                    unsigned sample_num = 0;
                    double sample_area = 0.0;

                    if (has_substrate_adhesion_area)
                    {
                        num_across = (unsigned)round((sample_area_right - sample_area_left) / preferred_sample_dis);
                        num_up = (unsigned)round((sample_area_top - sample_area_bottom) / preferred_sample_dis);
                        sample_num = num_across * num_up;
                        sample_area = (sample_area_top - sample_area_bottom)*(sample_area_right - sample_area_left);
                    }

                    bool if_node_a_inner_node = false;
                    if (containing_elem_indices.size()==3)
                        if_node_a_inner_node = true;

                    bool if_node_contact_to_reservoir_bottom = false;
                    if (fabs(p_this_node->rGetLocation()[1]) < 0.01/sqrt(M_PI/mFixedTargetArea))
                        if_node_contact_to_reservoir_bottom = true;

                    bool node_is_at_mesh_bottom = false;
                    if (containing_elem_indices.size()<3 && p_this_node->rGetLocation()[1]<3.0/4.0*strip_start_y_location)
                        node_is_at_mesh_bottom = true;
                    bool node_is_at_mesh_top = false;
                    if (containing_elem_indices.size()<3 && p_this_node->rGetLocation()[1]>3.0/4.0*strip_start_y_location)
                        node_is_at_mesh_top = true;

                    // Calculate reservoir_substrate_adhesion_area_gradient!
                    if (!if_node_a_inner_node && !(mIfIgnoreReservoirSubstrateAdhesionAtBottom&&node_is_at_mesh_bottom)
                                && !(mIfIgnoreReservoirSubstrateAdhesionAtTop&&node_is_at_mesh_top) 
                                && !if_node_contact_to_reservoir_bottom && has_substrate_adhesion_area && (num_across>0) && (num_up>0))
                    {
                        c_vector<c_vector<double, DIM>, 3> points;
                        points[0]=previous_node_location;
                        points[1]=this_node_location;
                        points[2]=next_node_location;
                        c_vector<double, DIM> vec1 = this_node_location-previous_node_location;
                        c_vector<double, DIM> vec2 = next_node_location-this_node_location;                    

                        // Calculate initial adhesive area
                        double adhesive_sample_num = 0.0;
                        for (unsigned i = 0; i <sample_num; i++)
                        {
                            double x_coord = sample_area_left + (sample_area_right - sample_area_left)/num_across*(i%num_across+0.5);
                            double y_coord = sample_area_bottom + (sample_area_top - sample_area_bottom)/num_up*(i/num_across+0.5);

                            c_vector<double, DIM> point;
                            point[0] = x_coord;
                            point[1] = y_coord;
                            c_vector<bool, 3> point_at_left_of_vector;
                            for (unsigned j = 0; j < 3; j++)
                            {
                                c_vector<double, DIM> vec1 = point - points[j];
                                c_vector<double, DIM> vec2 = points[(j+1)%3] - points[j];

                                if ( (vec1[0]*vec2[1]-vec1[1]*vec2[0]) > 0.0)
                                    point_at_left_of_vector[j] = false;
                                else
                                    point_at_left_of_vector[j] = true;
                            }
                            if (point_at_left_of_vector[0]==true && point_at_left_of_vector[1]==true && point_at_left_of_vector[2]==true)
                                adhesive_sample_num += 1.0;
                            else if (point_at_left_of_vector[0]==false && point_at_left_of_vector[1]==false && point_at_left_of_vector[2]==false)
                                adhesive_sample_num += -1.0;
                        }
                        double substrate_adhesion_area = adhesive_sample_num/double(sample_num) * sample_area;

                        // Calculate adhesive area *changed* with small displacement of the node along the x or y axis
                        for (unsigned j = 0; j<2; j++)
                        {
                            c_vector<double, DIM> unit_vector_in_small_change_direction;
                            unit_vector_in_small_change_direction[0] = cos(j*M_PI/2);
                            unit_vector_in_small_change_direction[1] = sin(j*M_PI/2);

                            c_vector<c_vector<double, DIM>, 3> new_points = points;
                            new_points[1] += small_change*unit_vector_in_small_change_direction;

                            double adhesive_sample_num = 0.0;

                            for (unsigned i = 0; i <sample_num; i++)
                            {
                                double x_coord = sample_area_left + (sample_area_right - sample_area_left)/num_across*(i%num_across+0.5);
                                double y_coord = sample_area_bottom + (sample_area_top - sample_area_bottom)/num_up*(i/num_across+0.5);

                                c_vector<double, DIM> point;
                                point[0] = x_coord;
                                point[1] = y_coord;
                                c_vector<bool, 3> point_at_left_of_vector;
                                for (unsigned j = 0; j < 3; j++)
                                {
                                    c_vector<double, DIM> vec1 = point - new_points[j];
                                    c_vector<double, DIM> vec2 = new_points[(j+1)%3] - new_points[j];

                                    if ( (vec1[0]*vec2[1]-vec1[1]*vec2[0]) > 0.0)
                                        point_at_left_of_vector[j] = false;
                                    else
                                        point_at_left_of_vector[j] = true;
                                }
                                if (point_at_left_of_vector[0]==true && point_at_left_of_vector[1]==true && point_at_left_of_vector[2]==true)
                                    adhesive_sample_num += 1.0;
                                else if (point_at_left_of_vector[0]==false && point_at_left_of_vector[1]==false && point_at_left_of_vector[2]==false)
                                    adhesive_sample_num += -1.0;
                            }
                            double substrate_adhesion_area_new = adhesive_sample_num/double(sample_num) * sample_area;
                            reservoir_substrate_adhesion_area_gradient += (substrate_adhesion_area_new - substrate_adhesion_area)/small_change*unit_vector_in_small_change_direction;
                        }// end of calculate adhesive area *changed* with small displacement of the node along the x or y axis
                        
                    }// end of statement 'if (has_substrate_adhesion_area)'

                    reservoir_substrate_adhesion_contribution -= mReservoirSubstrateAdhesionParameter*reservoir_substrate_adhesion_area_gradient;
                    
                }// end of If consider reservior substrate adhesion

            } //end of "if (mIfConsiderSubstrateAdhesion)"

            /*---------------------------------End of substrate adhesion contribution---------------------------------*/
            /*--------------------------------------------------------------------------------------------------------*/

        }// end of 'Iterate over these elements'

        c_vector<double, DIM> force_on_node = deformation_contribution + membrane_surface_tension_contribution + adhesion_contribution +area_adhesion_contribution +reservoir_substrate_adhesion_contribution;

        if (mAddPullingForceOnNodeIndividually && p_this_node->GetIndex()==node_index_for_leading_top_of_the_group0)
            force_on_node[1] += mPullingForceOnLeadingCell;
        else if (mAddPullingForceEvenlyOnNodesOfLeadingCell)
        {
            bool node_belongs_to_leading_cell = false;
            for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
                iter != containing_elem_indices.end();
                ++iter)
            {
                if (p_cell_population->GetElement(*iter)->GetIsLeadingCellTop())
                {
                    node_belongs_to_leading_cell = true;
                    break;
                }
            }

            if (node_belongs_to_leading_cell)
            {
                assert(num_nodes_leading_cell_top!=0);
                double t_now = SimulationTime::Instance()->GetTime();
                if ( !(mIfEquilibrateForAWhile && t_now<=mTimeForEquilibrium) )
                    force_on_node[1] += mPullingForceOnLeadingCell/num_nodes_leading_cell_top;
            }
        }

        // tmp: possilbe output for testing new SSA:
        bool output_while_testing_new_SSA = false;
        if (output_while_testing_new_SSA)
        {
            VertexElement<DIM, DIM>* pElement = p_cell_population->GetElement(* p_this_node->rGetContainingElementIndices().begin());
            if (pElement->GetIsLeadingCell() && pElement->GetGroupNumber()>0)
            {
                double t_now = SimulationTime::Instance()->GetTime();
                if (p_this_node->rGetContainingElementIndices().size()==1)
                {
                    std::cout << "Time: " << t_now << std::endl;
                    std::cout << "Node of LeadingCell, node index=" << p_this_node->GetIndex() 
                            << ", SSA case number=0, node location: x=" << p_this_node->rGetLocation()[0] << ", y=" << p_this_node->rGetLocation()[1]
                            << std::endl << "Lamellipodium strength=" << pElement->GetLamellipodiumStrength()
                            << std::endl << "SSA: x=" << area_adhesion_contribution[0] << ", y=" << area_adhesion_contribution[1] << std::endl;
                }
                else if (p_this_node->rGetContainingElementIndices().size()==2)
                {
                    std::cout << "Time: " << t_now << std::endl;
                    std::cout << "Node of LeadingCell, node index=" << p_this_node->GetIndex() 
                            << ", SSA case number=1, node location: x=" << p_this_node->rGetLocation()[0] << ", y=" << p_this_node->rGetLocation()[1]
                            << std::endl << "Lamellipodium strength1=" << pElement->GetLamellipodiumStrength() 
                            << ", Lamellipodium strength2=" << p_cell_population->GetElement(* p_this_node->rGetContainingElementIndices().begin()++)->GetLamellipodiumStrength()
                            << std::endl << "SSA: x=" << area_adhesion_contribution[0] << ", y=" << area_adhesion_contribution[1] << std::endl;
                }
            }
        }

        // to do: only typical nodes needs to be considered this carefully!
        if (mConsiderConsistencyForSSA && node_has_SSA)
        {
            c_vector<double, DIM> forces_except_SSA = deformation_contribution + membrane_surface_tension_contribution + adhesion_contribution +reservoir_substrate_adhesion_contribution;
            double forces_x_except_SSA = forces_except_SSA[0];
            double adaptive_dt = SimulationTime::Instance()->GetTimeStep();
            double node_x_location = p_this_node->rGetLocation()[0];
            double strip_right = mStripStartXLocation+mStripWidth/2.0;
            double strip_left = mStripStartXLocation-mStripWidth/2.0;
            double overall_dx =0.0;
            double overall_force_x = 0.0;

            double inner_strip_small_change = mSmallChangeForAreaCalculation;
            if (node_x_location>strip_right-mSmallChangeForAreaCalculation && node_x_location<strip_right+mSmallChangeForAreaCalculation)
                inner_strip_small_change = std::min(mSmallChangeForAreaCalculation-(node_x_location-strip_right), mSmallChangeForAreaCalculation);
            else if (node_x_location>mStripStartXLocation-mStripWidth/2.0-mSmallChangeForAreaCalculation && node_x_location<mStripStartXLocation-mStripWidth/2.0+mSmallChangeForAreaCalculation)
                inner_strip_small_change = std::min(mSmallChangeForAreaCalculation+(node_x_location-strip_left), mSmallChangeForAreaCalculation);
            averaged_inner_strip_substrate_adhesion_contribution[0] = area_adhesion_contribution[0] * mSmallChangeForAreaCalculation/inner_strip_small_change;
            averaged_inner_strip_substrate_adhesion_contribution[1] = area_adhesion_contribution[1];
            double SSA_x_inner_strip = averaged_inner_strip_substrate_adhesion_contribution[0];
            if (!(fabs(SSA_x_inner_strip) < 1000.0/pow(M_PI/mFixedTargetArea,1.5)) && inner_strip_small_change/mSmallChangeForAreaCalculation>0.1)
            {
                std::cout << "SSA_x_inner_strip is too large! SSA_x_inner_strip=" << SSA_x_inner_strip << std::endl;
                std::cout << "area_adhesion_contribution[0]=" << area_adhesion_contribution[0] << ", inner_strip_small_change/mSmallChange=" << inner_strip_small_change/mSmallChangeForAreaCalculation << std::endl;
            }

            if (inner_strip_small_change>0.2*mSmallChangeForAreaCalculation)
            {
                // if at around strip right!
                if (node_x_location>strip_right-mSmallChangeForAreaCalculation && node_x_location<strip_right+mSmallChangeForAreaCalculation)
                {
                    if (node_x_location>=strip_right)
                    {
                        if (forces_x_except_SSA>=0.0)
                            force_on_node[0] = forces_x_except_SSA;
                        else if (node_x_location+forces_x_except_SSA*adaptive_dt >= strip_right)
                            force_on_node[0] = forces_x_except_SSA;
                        else if (forces_x_except_SSA+SSA_x_inner_strip >= 0.0)
                        {
                            overall_dx = strip_right - node_x_location;
                            overall_force_x = overall_dx/adaptive_dt;
                            force_on_node[0] = overall_force_x;
                        }
                        else
                        {
                            assert(forces_x_except_SSA+SSA_x_inner_strip < 0.0);
                            assert(adaptive_dt-(strip_right-node_x_location)/forces_x_except_SSA >= 0.0);
                            overall_dx = strip_right - node_x_location +(adaptive_dt-(strip_right-node_x_location)/forces_x_except_SSA)*(forces_x_except_SSA+SSA_x_inner_strip);
                            overall_force_x = overall_dx/adaptive_dt;
                            force_on_node[0] = overall_force_x;
                        }
                    }
                    else
                    {
                        assert(node_x_location<strip_right);
                        if (forces_x_except_SSA<0.0)
                        {
                            if (forces_x_except_SSA+SSA_x_inner_strip < 0.0)
                                force_on_node[0] = forces_x_except_SSA + SSA_x_inner_strip;
                            else
                            {
                                assert(forces_x_except_SSA+SSA_x_inner_strip >= 0.0);
                                overall_dx = std::min((forces_x_except_SSA+SSA_x_inner_strip)*adaptive_dt, strip_right-node_x_location);
                                overall_force_x = overall_dx/adaptive_dt;
                                force_on_node[0] = overall_force_x;
                            }
                        }
                        else
                        {
                            assert(forces_x_except_SSA>=0.0);
                            if (forces_x_except_SSA+SSA_x_inner_strip <=0.0)
                                force_on_node[0] = forces_x_except_SSA +SSA_x_inner_strip;
                            else if ((forces_x_except_SSA+SSA_x_inner_strip)*adaptive_dt <= strip_right-node_x_location)
                                force_on_node[0] = forces_x_except_SSA +SSA_x_inner_strip;
                            else
                            {
                                assert(adaptive_dt-(strip_right-node_x_location)/(forces_x_except_SSA+SSA_x_inner_strip)>=0.0);
                                overall_dx = strip_right - node_x_location + forces_x_except_SSA*(adaptive_dt-(strip_right-node_x_location)/(forces_x_except_SSA+SSA_x_inner_strip));
                                overall_force_x = overall_dx/adaptive_dt;
                                force_on_node[0] = overall_force_x;
                            }
                        }
                        
                    }
                }
                // if at around strip left!            
                else if (node_x_location>strip_left-mSmallChangeForAreaCalculation && node_x_location<strip_left+mSmallChangeForAreaCalculation)
                {
                    if (node_x_location<=strip_left)
                    {
                        if (forces_x_except_SSA<=0.0)
                            force_on_node[0] = forces_x_except_SSA;
                        else if (node_x_location+forces_x_except_SSA*adaptive_dt <= strip_left)
                            force_on_node[0] = forces_x_except_SSA;
                        else if (forces_x_except_SSA+SSA_x_inner_strip <= 0.0)
                        {
                            overall_dx = strip_left - node_x_location;
                            overall_force_x = overall_dx/adaptive_dt;
                            force_on_node[0] = overall_force_x;
                        }
                        else
                        {
                            assert(forces_x_except_SSA+SSA_x_inner_strip > 0.0);
                            assert(adaptive_dt-(strip_left-node_x_location)/forces_x_except_SSA >= 0.0);
                            overall_dx = strip_left - node_x_location +(adaptive_dt-(strip_left-node_x_location)/forces_x_except_SSA)*(forces_x_except_SSA+SSA_x_inner_strip);
                            overall_force_x = overall_dx/adaptive_dt;
                            force_on_node[0] = overall_force_x;
                        }
                    }
                    else
                    {
                        assert(node_x_location>strip_left);
                        if (forces_x_except_SSA>0.0)
                        {
                            if (forces_x_except_SSA+SSA_x_inner_strip > 0.0)
                                force_on_node[0] = forces_x_except_SSA + SSA_x_inner_strip;
                            else
                            {
                                assert(forces_x_except_SSA+SSA_x_inner_strip <= 0.0);
                                overall_dx = std::max((forces_x_except_SSA+SSA_x_inner_strip)*adaptive_dt, strip_left-node_x_location);
                                overall_force_x = overall_dx/adaptive_dt;
                                force_on_node[0] = overall_force_x;
                            }
                        }
                        else
                        {
                            assert(forces_x_except_SSA<=0.0);
                            if (forces_x_except_SSA+SSA_x_inner_strip >=0.0)
                                force_on_node[0] = forces_x_except_SSA +SSA_x_inner_strip;
                            else if ((forces_x_except_SSA+SSA_x_inner_strip)*adaptive_dt >= strip_left-node_x_location)
                                force_on_node[0] = forces_x_except_SSA +SSA_x_inner_strip;
                            else
                            {
                                assert(adaptive_dt-(strip_left-node_x_location)/(forces_x_except_SSA+SSA_x_inner_strip)>=0.0);
                                overall_dx = strip_left - node_x_location + forces_x_except_SSA*(adaptive_dt-(strip_left-node_x_location)/(forces_x_except_SSA+SSA_x_inner_strip));
                                overall_force_x = overall_dx/adaptive_dt;
                                force_on_node[0] = overall_force_x;
                            }
                        }
                        
                    }
                }
                else
                    force_on_node[0] = forces_x_except_SSA + area_adhesion_contribution[0];
            }
            
        }
        
        p_cell_population->GetNode(node_index)->AddAppliedForceContribution(force_on_node);

        // tmp
        if (mOutputInformationForNagaiHondaForce)
        {
            if (norm_2(force_on_node)>2)
            {
                std::cout << "Weird! Force is too large! Node Index: " << p_this_node->GetIndex() << std::endl;
                std::cout << "Node location: " << p_this_node->rGetLocation()[0] << ", " << p_this_node->rGetLocation()[1] << std::endl;
                std::cout << "X direction:" << std::endl;
                std::cout << "force on node=" << force_on_node[0]
                        << " deformation_contribution=" << deformation_contribution[0]
                        << " membrane_surface_tension_contribution=" << membrane_surface_tension_contribution[0]
                        << " adhesion_contribution=" << adhesion_contribution[0]
                        << " area_adhesion_contribution=" << area_adhesion_contribution[0]
                        << " reservoir_substrate_adhesion_contribution=" << reservoir_substrate_adhesion_contribution[0] << std::endl;
                std::cout << "Y direction:" << std::endl;
                std::cout << "force on node=" << force_on_node[1]
                        << " deformation_contribution=" << deformation_contribution[1]
                        << " membrane_surface_tension_contribution=" << membrane_surface_tension_contribution[1]
                        << " adhesion_contribution=" << adhesion_contribution[1]
                        << " area_adhesion_contribution=" << area_adhesion_contribution[1]
                        << " reservoir_substrate_adhesion_contribution=" << reservoir_substrate_adhesion_contribution[1] << std::endl;
            }
        }

    }// end of 'Iterate over nodes(vertices) in the cell population'

}

template<unsigned DIM>
double MyNewNagaiHondaForceWithStripesAdhesion<DIM>::GetAdhesionParameter(Node<DIM>* pNodeA, Node<DIM>* pNodeB, VertexBasedCellPopulation<DIM>& rVertexCellPopulation)
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
        return adhesion_parameter;
    }

    // my changes
    if (mIfUseFaceElementToGetAdhesionParameter)
    {
        unsigned element_global_index = *shared_elements.begin();
        VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rVertexCellPopulation);
        VertexElement<DIM, DIM>* p_element = p_cell_population->rGetMesh().GetElement(element_global_index);
        if ((p_element->GetNodeLocalIndex(pNodeB->GetIndex()) - p_element->GetNodeLocalIndex(pNodeA->GetIndex())+p_element->GetNumNodes())%p_element->GetNumNodes()==1)
        {
            VertexElement<DIM-1, DIM>* pFace = p_element->GetFace(p_element->GetFaceLocalIndexUsingStartAndEndNodeGlobalIndex(pNodeA->GetIndex(), pNodeB->GetIndex()));
            adhesion_parameter =GetNagaiHondaCellCellAdhesionEnergyParameter() * pFace->GetUnifiedCellCellAdhesionEnergyParameter();
        }
        else
        {
            if ( (p_element->GetNodeLocalIndex(pNodeA->GetIndex()) - p_element->GetNodeLocalIndex(pNodeB->GetIndex())+p_element->GetNumNodes())%p_element->GetNumNodes()!=1)
            {
                std::cout << std::endl << "ERROR: Method MyNewNagaiHondaForceWithStripesAdhesion::GetAdhesionParameter";
                std::cout << std::endl << "Not Reachable!" << std::endl;
            }
            assert((p_element->GetNodeLocalIndex(pNodeA->GetIndex()) - p_element->GetNodeLocalIndex(pNodeB->GetIndex())+p_element->GetNumNodes())%p_element->GetNumNodes()==1);
            VertexElement<DIM-1, DIM>* pFace = p_element->GetFace(p_element->GetFaceLocalIndexUsingStartAndEndNodeGlobalIndex(pNodeB->GetIndex(), pNodeA->GetIndex()));
            adhesion_parameter =GetNagaiHondaCellCellAdhesionEnergyParameter() * pFace->GetUnifiedCellCellAdhesionEnergyParameter();    
        }
    }
    else
        adhesion_parameter =GetNagaiHondaCellCellAdhesionEnergyParameter();

    return adhesion_parameter;
}

template<unsigned DIM>
double MyNewNagaiHondaForceWithStripesAdhesion<DIM>::GetNagaiHondaDeformationEnergyParameter()
{
    return mNagaiHondaDeformationEnergyParameter;
}

template<unsigned DIM>
double MyNewNagaiHondaForceWithStripesAdhesion<DIM>::GetNagaiHondaMembraneSurfaceEnergyParameter()
{
    return mNagaiHondaMembraneSurfaceEnergyParameter;
}

template<unsigned DIM>
double MyNewNagaiHondaForceWithStripesAdhesion<DIM>::GetNagaiHondaCellCellAdhesionEnergyParameter()
{
    return mNagaiHondaCellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
double MyNewNagaiHondaForceWithStripesAdhesion<DIM>::GetNagaiHondaCellBoundaryAdhesionEnergyParameter()
{
    return mNagaiHondaCellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void MyNewNagaiHondaForceWithStripesAdhesion<DIM>::SetNagaiHondaDeformationEnergyParameter(double deformationEnergyParameter)
{
    mNagaiHondaDeformationEnergyParameter = deformationEnergyParameter;
}

template<unsigned DIM>
void MyNewNagaiHondaForceWithStripesAdhesion<DIM>::SetNagaiHondaMembraneSurfaceEnergyParameter(double membraneSurfaceEnergyParameter)
{
    mNagaiHondaMembraneSurfaceEnergyParameter = membraneSurfaceEnergyParameter;
}

template<unsigned DIM>
void MyNewNagaiHondaForceWithStripesAdhesion<DIM>::SetNagaiHondaCellCellAdhesionEnergyParameter(double cellCellAdhesionEnergyParameter)
{
    mNagaiHondaCellCellAdhesionEnergyParameter = cellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
void MyNewNagaiHondaForceWithStripesAdhesion<DIM>::SetNagaiHondaCellBoundaryAdhesionEnergyParameter(double cellBoundaryAdhesionEnergyParameter)
{
    mNagaiHondaCellBoundaryAdhesionEnergyParameter = cellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void MyNewNagaiHondaForceWithStripesAdhesion<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<NagaiHondaDeformationEnergyParameter>" << mNagaiHondaDeformationEnergyParameter << "</NagaiHondaDeformationEnergyParameter>\n";
    *rParamsFile << "\t\t\t<NagaiHondaMembraneSurfaceEnergyParameter>" << mNagaiHondaMembraneSurfaceEnergyParameter << "</NagaiHondaMembraneSurfaceEnergyParameter>\n";
    *rParamsFile << "\t\t\t<NagaiHondaCellCellAdhesionEnergyParameter>" << mNagaiHondaCellCellAdhesionEnergyParameter << "</NagaiHondaCellCellAdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<NagaiHondaCellBoundaryAdhesionEnergyParameter>" << mNagaiHondaCellBoundaryAdhesionEnergyParameter << "</NagaiHondaCellBoundaryAdhesionEnergyParameter>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class MyNewNagaiHondaForceWithStripesAdhesion<1>;
template class MyNewNagaiHondaForceWithStripesAdhesion<2>;
template class MyNewNagaiHondaForceWithStripesAdhesion<3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MyNewNagaiHondaForceWithStripesAdhesion)
