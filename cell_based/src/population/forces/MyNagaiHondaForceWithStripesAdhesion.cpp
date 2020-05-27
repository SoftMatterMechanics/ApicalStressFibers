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

#include "MyNagaiHondaForceWithStripesAdhesion.hpp"

template<unsigned DIM>
MyNagaiHondaForceWithStripesAdhesion<DIM>::MyNagaiHondaForceWithStripesAdhesion()
   : AbstractForce<DIM>(),
     mNagaiHondaDeformationEnergyParameter(1.0),
     mNagaiHondaMembraneSurfaceEnergyParameter(0.1),
     mNagaiHondaCellCellAdhesionEnergyParameter(0.0),
     mNagaiHondaCellBoundaryAdhesionEnergyParameter(0.0)
{
}

template<unsigned DIM>
MyNagaiHondaForceWithStripesAdhesion<DIM>::~MyNagaiHondaForceWithStripesAdhesion()
{
}

template<unsigned DIM>
void MyNagaiHondaForceWithStripesAdhesion<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    bool output_information =false;
    
    // Throw an exception message if not using a VertexBasedCellPopulation
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("NagaiHondaForce is to be used with a VertexBasedCellPopulation only");
    }

    // Define some helper variables
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    unsigned num_nodes = p_cell_population->GetNumNodes();
    unsigned num_elements = p_cell_population->GetNumElements();

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
            target_areas[elem_index] = p_cell_population->GetCellUsingLocationIndex(elem_index)->GetCellData()->GetItem("target area");
            if (this->mFixedTargetArea != 0.0)
            {
                target_areas[elem_index] = mFixedTargetArea;
            }
        }
        catch (Exception&)
        {
            EXCEPTION("You need to add an AbstractTargetAreaModifier to the simulation in order to use NagaiHondaForce");
        }

    }

    double y_coord_leading_edge = 0.0;
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        Node<DIM>* p_this_node = p_cell_population->GetNode(node_index);
        if (p_this_node->rGetLocation()[1] > y_coord_leading_edge)
            y_coord_leading_edge= p_this_node->rGetLocation()[1];
    }


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
            // default shape index of a circle:
            double cell_target_perimeter = 2*sqrt(M_PI*target_areas[elem_index]);
            if (this->mFixedTargetPerimeter !=0.0)   
            {
                cell_target_perimeter = this->mFixedTargetPerimeter;
            }
            if (this->mFixedTargetShapeIndex !=0.0)       
            {
                cell_target_perimeter = this->mFixedTargetShapeIndex * sqrt(target_areas[elem_index]);
            }
            
            CellPtr p_cell = p_cell_population->GetCellUsingLocationIndex(elem_index);

            double current_time = SimulationTime::Instance()->GetTime();

            // double myosin_activity = 0.0;

            // typical *membrane surface tension*
            if (current_time < (10000.0 + 1e-10))
            {
                // My changes: compute the stress tensor and 
                // iterate over all the faces(edges) of this element to get the EMA and C-C adh_para.
                membrane_surface_tension_contribution -= GetNagaiHondaMembraneSurfaceEnergyParameter()*1.0*(element_perimeters[elem_index] - cell_target_perimeter)*element_perimeter_gradient;
            }
            // MA feedback incorporated in *membrane surface tension*
            else if (current_time < (5000.0 + 1e-10))
            {
                // myosin_activity = p_cell->GetMyosinActivity();
                // membrane_surface_tension_contribution -= GetNagaiHondaMembraneSurfaceEnergyParameter()*myosin_activity*(element_perimeters[elem_index] - cell_target_perimeter)*element_perimeter_gradient;
            }
            // EMA feedback incorporated in *membrane surface tension*
            else
            {
                // double sum1 = 0.0;
                // c_vector<double, DIM> sum2 = zero_vector<double>(DIM);;

                // unsigned this_node_index = p_element->GetNodeGlobalIndex(0);
                // for (unsigned i = 0; i < num_nodes_elem; i++)
                // {
                //     unsigned next_node_index = p_element->GetNodeGlobalIndex((i + 1) % num_nodes_elem);

                //     double edge_length = p_cell_population->rGetMesh().GetDistanceBetweenNodes(this_node_index, next_node_index);
                //     sum1 += sqrt(p_element->GetEdgeMyosinActivity(i))*edge_length;
                //     this_node_index = next_node_index;
                // }

                // sum2 += sqrt(p_element->GetEdgeMyosinActivity(previous_node_local_index))*previous_edge_gradient;
                // sum2 += sqrt(p_element->GetEdgeMyosinActivity(local_index))*next_edge_gradient;

                // membrane_surface_tension_contribution -= GetNagaiHondaMembraneSurfaceEnergyParameter()*sum1*sum2;
            }

            /*--------------------------------------------------------------------------------------------------------------*/
            /*--------------------------------------------------------------------------------------------------------------*/
            /*---------------------------------Cell Area adhesion contribution---------------------------------*/

            double stripes_adhesion_start_time = -100.0;
            if (current_time < (stripes_adhesion_start_time + 1e-10))
            {
                area_adhesion_contribution += zero_vector<double>(DIM);
            }
            else
            {
                // Parameters
                double stripe_width = this->mStripWidth;
                double stripe_distance = this->mStripDistance;
                double strip_start_x_location = this->mStripStartXLocation;
                double strip_start_y_location = this->mStripStartYLocation;
                double substrate_adhesion_parameter = this->mSubstrateAdhesionParameter;
                bool if_substrate_adhesion_parameter_change = this->mIfSubstrateAdhesionParameterChange;
                double substrate_adhesion_leading_top_length = 1.0;

                bool use_fine_mesh =this->mUseFineMesh;
                double small_change = 0.02;
                double preferred_sample_dis = 0.004;
                if (!use_fine_mesh)
                {
                    small_change = 0.05;
                    preferred_sample_dis = 0.01;
                }

                c_vector<double, DIM> centroid = p_cell_population->rGetMesh().GetCentroidOfElement(elem_index);
                int stripe_num = round((centroid[0] - strip_start_x_location)/stripe_distance);

                double stripe_location = strip_start_x_location+ stripe_distance*stripe_num;
                double stripe_left = stripe_location - stripe_width/2;
                double stripe_right = stripe_location + stripe_width/2;

                double sample_area_bottom = 0.0;
                double sample_area_top = 0.0;
                double sample_area_left = 0.0;
                double sample_area_right = 0.0;

                c_vector<double, DIM> previous_node_location = p_previous_node->rGetLocation();
                c_vector<double, DIM> this_node_location = p_this_node->rGetLocation();
                c_vector<double, DIM> next_node_location = p_next_node->rGetLocation();
                
                double expanded_triangle_box_bottom = -small_change+std::min(std::min(previous_node_location[1],this_node_location[1]),next_node_location[1]);
                double expanded_triangle_box_top = small_change+std::max(std::max(previous_node_location[1],this_node_location[1]),next_node_location[1]);
                double expanded_triangle_box_left = -small_change+std::min(std::min(previous_node_location[0],this_node_location[0]),next_node_location[0]);
                double expanded_triangle_box_right = small_change+std::max(std::max(previous_node_location[0],this_node_location[0]),next_node_location[0]);
                
                bool has_substrate_adhesion_area = true;
                if (fabs(expanded_triangle_box_left-expanded_triangle_box_right) > stripe_distance/2)
                has_substrate_adhesion_area = false;

                if (expanded_triangle_box_top<strip_start_y_location)
                    has_substrate_adhesion_area = false;
                else
                {
                    sample_area_top = expanded_triangle_box_top;
                    sample_area_bottom = expanded_triangle_box_bottom>strip_start_y_location ? expanded_triangle_box_bottom : strip_start_y_location;
                }

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
                // Attention
                c_vector<double, DIM> substrate_adhesion_area_gradient = zero_vector<double>(DIM);
                c_vector<double, DIM> weighted_substrate_adhesion_area_gradient = zero_vector<double>(DIM);
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

                if (has_substrate_adhesion_area && (num_across>0) && (num_up>0))
                {
                    c_vector<c_vector<double, DIM>, 3> points;
                    points[0]=previous_node_location;
                    points[1]=this_node_location;
                    points[2]=next_node_location;
                    c_vector<double, DIM> vec1 = this_node_location-previous_node_location;
                    c_vector<double, DIM> vec2 = next_node_location-this_node_location;                    
                    //check
                    // if ((vec1[0]*vec2[1]-vec1[1]*vec2[0])< 0)
                    // std::cout << std::endl<< "Note:not counter clockwise!" << std::endl;

                    // Calculate initial adhesive area
                    unsigned adhesive_sample_num = 0;
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
                        bool point_is_outside = false;
                        for (unsigned j = 0; j < 3; j++)
                        {
                            c_vector<double, DIM> vec1 = point - points[j];
                            c_vector<double, DIM> vec2 = points[(j+1)%3] - points[j];

                            if ( (vec1[0]*vec2[1]-vec1[1]*vec2[0]) > 0.0)
                            {
                                point_is_outside = true;
                                break;
                            }
                        }
                        if (point_is_outside == false)
                        {
                            adhesive_sample_num += 1;
                            if (y_coord > (y_coord_leading_edge - substrate_adhesion_leading_top_length))
                            {
                                weighted_adhesive_sample_num += this-> mSubstrateAdhesionParameterTopAreaMagnifies;
                            }
                            else
                                weighted_adhesive_sample_num += 1.0;
                        }
                    }
                    double substrate_adhesion_area = adhesive_sample_num/double(sample_num) * sample_area;
                    double weighted_substrate_adhesion_area = weighted_adhesive_sample_num/double(sample_num) * sample_area;

                    // Calculate adhesive area *changed* with small displacement of the node along the x or y axis
                    for (unsigned j = 0; j<2; j++)
                    {
                        c_vector<double, DIM> small_vector;
                        small_vector[0] = small_change*cos(j*M_PI/2);
                        small_vector[1] = small_change*sin(j*M_PI/2);

                        c_vector<c_vector<double, DIM>, 3> new_points = points;
                        new_points[1] += small_vector;

                        unsigned adhesive_sample_num = 0;
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
                            bool point_is_outside = false;
                            for (unsigned j = 0; j < 3; j++)
                            {
                                c_vector<double, DIM> vec1 = point - new_points[j];
                                c_vector<double, DIM> vec2 = new_points[(j+1)%3] - new_points[j];

                                if ( (vec1[0]*vec2[1]-vec1[1]*vec2[0]) > 0.0)
                                {
                                    point_is_outside = true;
                                    break;
                                }
                            }
                            if (point_is_outside == false)
                            {
                                adhesive_sample_num += 1;
                                if (y_coord > (y_coord_leading_edge - substrate_adhesion_leading_top_length))
                                {
                                    weighted_adhesive_sample_num += this-> mSubstrateAdhesionParameterTopAreaMagnifies;
                                }
                                else
                                    weighted_adhesive_sample_num += 1.0;
                            }
                        }
                        double substrate_adhesion_area_new = adhesive_sample_num/double(sample_num) * sample_area;
                        double weighted_substrate_adhesion_area_new = weighted_adhesive_sample_num/double(sample_num) * sample_area;
                        substrate_adhesion_area_gradient += (substrate_adhesion_area_new - substrate_adhesion_area)/small_change*small_vector;
                        weighted_substrate_adhesion_area_gradient += (weighted_substrate_adhesion_area_new - weighted_substrate_adhesion_area)/small_change*small_vector;
                    }

                }// end of statement 'if (has_substrate_adhesion_area)'

                if (if_substrate_adhesion_parameter_change == false)
                    area_adhesion_contribution -= substrate_adhesion_parameter*substrate_adhesion_area_gradient;
                else
                    area_adhesion_contribution -= substrate_adhesion_parameter*weighted_substrate_adhesion_area_gradient;
            }
            if (output_information)
            {
                std::cout << "1.Element iteration, node_index: " << node_index << " ";
                std::cout << element_perimeters[elem_index] << " ";
                std::cout << element_areas[elem_index] << " ";
            }

        }// end of 'Iterate over these elements'

        c_vector<double, DIM> force_on_node = deformation_contribution + membrane_surface_tension_contribution + adhesion_contribution +area_adhesion_contribution;

        // ///=My changes: random force
        // double motility = 0.01;
        // // double Dr = 0.1;

        // // double U1 = (rand()%100)/double(100);
        // // double U2 = (rand()%100)/double(100);
        // // if (U1<1e-10)
        // //     U1 = 1e-10;
        // // double Z = sqrt((-2)*log(U1))*cos(2*M_PI*U2);
        // // double angle_change = Dr*Z;

        // double random_force_angle = rand();
        // //random_force_angle += angle_change;
        // //p_cell_population->GetNode(node_index)->SetRandomForceAngle(random_force_angle);
        // c_vector<double, DIM> random_force = zero_vector<double>(DIM);
        // //random_force_angle = M_PI/2;
        // random_force[0] = motility*cos(random_force_angle);
        // random_force[1] = motility*sin(random_force_angle);
        // //std::cout << random_force[0] << ' '<< motility << " " << random_force_angle << ' '<< cos(random_force_angle) << '\n';

        // //force_on_node = zero_vector<double>(DIM);
        // force_on_node += random_force;
        // //force_on_node += zero_vector<double>(DIM);

        /************************/
        /*std::cout << norm_2(deformation_contribution) << ' ' << norm_2(membrane_surface_tension_contribution) << ' ' << norm_2(adhesion_contribution) << ' ' << norm_2(area_adhesion_contribution) << ' ' << norm_2(random_force) << ' ' << area_adhesion_contribution[0] << '\n' ;*/
        // My test
        if (output_information)
        {
            std::cout << std::endl;
            std::cout << "2.Force analyze for node: " << node_index <<":";
            std::cout << deformation_contribution[0] << " " <<membrane_surface_tension_contribution[0] <<" "<<adhesion_contribution[0] <<" "<<area_adhesion_contribution[0]<<" "<< "Total force: "<<force_on_node[0]<<std::endl;
        }

        // int ii = 0;
        // for (std::set<unsigned>::iterator iter = containing_elem_indices.begin();
        //      iter != containing_elem_indices.end();
        //      ++iter)
        //      ii+=1;

        // c_vector<double, DIM> the_point_location = p_this_node->rGetLocation();
        // if (the_point_location[0]<3.9 && the_point_location[0]>0.1)
        // {
        //     if (the_point_location[1]>0.1)
        //     {
        //         if (ii==3 &&(abs(area_adhesion_contribution[0]) >1e-5 || abs(area_adhesion_contribution[1]) >1e-5))
        //         {
        //             if (output_information)
        //                 std::cout << "3.Error anlysation for node_index: " << node_index << "contaning elements: " << ii << " , location and force: " << the_point_location[0] << " " << the_point_location[1] << " " << area_adhesion_contribution[0]<< " "<<area_adhesion_contribution[1]<<std::endl ;
        //         }
        //     }
        // }

        p_cell_population->GetNode(node_index)->AddAppliedForceContribution(force_on_node);

    }// end of 'Iterate over vertices in the cell population'
}

template<unsigned DIM>
double MyNagaiHondaForceWithStripesAdhesion<DIM>::GetAdhesionParameter(Node<DIM>* pNodeA, Node<DIM>* pNodeB, VertexBasedCellPopulation<DIM>& rVertexCellPopulation)
{
    double adhesion_parameter = 0.0;
    unsigned global_index_nodeA = pNodeA->GetIndex();
    unsigned global_index_nodeB = pNodeB->GetIndex();

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

    // Change it later:
    unsigned element_global_index = *shared_elements.begin();//shared_elements[0]
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rVertexCellPopulation);
    VertexElement<DIM, DIM>* p_element = p_cell_population->rGetMesh().GetElement(element_global_index);
    // std::cout << "global_index_nodeA= " << global_index_nodeA << " global_index_nodeB" << global_index_nodeB <<std::endl;
    for (unsigned index=0; index < p_element->GetNumFaces(); index++)
    {
        VertexElement<DIM-1,DIM>* p_face = p_element->GetFace(index);
        
        // std::cout << "face local index= " << index << " p_face->GetNodeGlobalIndex(0)=" << p_face->GetNodeGlobalIndex(0) << " p_face->GetNodeGlobalIndex(1)=" << p_face->GetNodeGlobalIndex(1) << std::endl;
        if ((p_face->GetNodeGlobalIndex(0)==global_index_nodeA && p_face->GetNodeGlobalIndex(1)==global_index_nodeB)||(p_face->GetNodeGlobalIndex(0)==global_index_nodeB && p_face->GetNodeGlobalIndex(1)==global_index_nodeA))
        {
            adhesion_parameter = p_face->GetUnifiedCellCellAdhesionEnergyParameter();//not considering boundary edge
            break;
        }
        if (index == p_element->GetNumFaces()-1)
        {
            std::cout<<"/////ERR:method GetAdhesionParameter in My.Force.Adhesion gets error for finding relevant face(edge)" <<std::endl;
            std::cout<< "possible reason: T1swap" <<std::endl;
        }
    }
    // double adhesion_parameter = GetNagaiHondaCellCellAdhesionEnergyParameter();

    // // If the edge corresponds to a single element, then the cell is on the boundary
    // if (shared_elements.size() == 1)
    // {
    //     adhesion_parameter = GetNagaiHondaCellBoundaryAdhesionEnergyParameter();
    // }

    return adhesion_parameter;
}

template<unsigned DIM>
double MyNagaiHondaForceWithStripesAdhesion<DIM>::GetNagaiHondaDeformationEnergyParameter()
{
    return mNagaiHondaDeformationEnergyParameter;
}

template<unsigned DIM>
double MyNagaiHondaForceWithStripesAdhesion<DIM>::GetNagaiHondaMembraneSurfaceEnergyParameter()
{
    return mNagaiHondaMembraneSurfaceEnergyParameter;
}

template<unsigned DIM>
double MyNagaiHondaForceWithStripesAdhesion<DIM>::GetNagaiHondaCellCellAdhesionEnergyParameter()
{
    return mNagaiHondaCellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
double MyNagaiHondaForceWithStripesAdhesion<DIM>::GetNagaiHondaCellBoundaryAdhesionEnergyParameter()
{
    return mNagaiHondaCellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void MyNagaiHondaForceWithStripesAdhesion<DIM>::SetNagaiHondaDeformationEnergyParameter(double deformationEnergyParameter)
{
    mNagaiHondaDeformationEnergyParameter = deformationEnergyParameter;
}

template<unsigned DIM>
void MyNagaiHondaForceWithStripesAdhesion<DIM>::SetNagaiHondaMembraneSurfaceEnergyParameter(double membraneSurfaceEnergyParameter)
{
    mNagaiHondaMembraneSurfaceEnergyParameter = membraneSurfaceEnergyParameter;
}

template<unsigned DIM>
void MyNagaiHondaForceWithStripesAdhesion<DIM>::SetNagaiHondaCellCellAdhesionEnergyParameter(double cellCellAdhesionEnergyParameter)
{
    mNagaiHondaCellCellAdhesionEnergyParameter = cellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
void MyNagaiHondaForceWithStripesAdhesion<DIM>::SetNagaiHondaCellBoundaryAdhesionEnergyParameter(double cellBoundaryAdhesionEnergyParameter)
{
    mNagaiHondaCellBoundaryAdhesionEnergyParameter = cellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void MyNagaiHondaForceWithStripesAdhesion<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<NagaiHondaDeformationEnergyParameter>" << mNagaiHondaDeformationEnergyParameter << "</NagaiHondaDeformationEnergyParameter>\n";
    *rParamsFile << "\t\t\t<NagaiHondaMembraneSurfaceEnergyParameter>" << mNagaiHondaMembraneSurfaceEnergyParameter << "</NagaiHondaMembraneSurfaceEnergyParameter>\n";
    *rParamsFile << "\t\t\t<NagaiHondaCellCellAdhesionEnergyParameter>" << mNagaiHondaCellCellAdhesionEnergyParameter << "</NagaiHondaCellCellAdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<NagaiHondaCellBoundaryAdhesionEnergyParameter>" << mNagaiHondaCellBoundaryAdhesionEnergyParameter << "</NagaiHondaCellBoundaryAdhesionEnergyParameter>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class MyNagaiHondaForceWithStripesAdhesion<1>;
template class MyNagaiHondaForceWithStripesAdhesion<2>;
template class MyNagaiHondaForceWithStripesAdhesion<3>;


// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MyNagaiHondaForceWithStripesAdhesion)
