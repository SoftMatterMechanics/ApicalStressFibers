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
#include "MyNoPBCToroidal2dVertexMesh.hpp"
#include <cassert>
#include <algorithm>

template<unsigned DIM>
MyStressfiberTensionForce<DIM>::MyStressfiberTensionForce()
   : AbstractForce<DIM>(),
     mIfEquilibrateForAWhile(false),
     mStartTimeForStretching(0.0),
     mFlag(0),
     mNucleationPerimeterTension(0.0)
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
    double dt = SimulationTime::Instance()->GetTimeStep();

    std::string name_item;
    std::ostringstream oss;

    if ( !(mIfEquilibrateForAWhile && (t_now>=mStartTimeForStretching)) )
    {
        for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
            elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
            ++elem_iter)
        {
            unsigned elem_index = elem_iter->GetIndex();
            // VertexElement<DIM, DIM>* p_element = p_cell_population->rGetMesh().GetElement(elem_index);
            CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(elem_index);
            p_cell->GetCellData()->SetItem("nucleation_flag",mFlag);
            p_cell->GetCellData()->SetItem("peeling_flag",0);
        }
    }
    else
    {
        unsigned global_index = 0;

        OutputFileHandler output_file_handler(this->mOutputDirectory+"/", false);
        std::string output_file_for_angle_bisector;
        // output_file_for_angle_bisector = "anglebisector.dat";
        std::ostringstream file_string_angle_bisector;
        file_string_angle_bisector << "anglebisector" << std::to_string(mAreaSeed) << ".dat";
        output_file_for_angle_bisector = file_string_angle_bisector.str();
        mpAngleBisectorFile = output_file_handler.OpenOutputFile(output_file_for_angle_bisector, std::ios::app);
        if ( (t_now - 1.0*floor(t_now/1.0)) < dt )
        {
            *mpAngleBisectorFile << t_now << "\n";
        }

        std::string output_file_for_elem_info;
        // output_file_for_elem_info = "eleminfo.dat";
        std::ostringstream file_string_elem_info;
        file_string_elem_info << "eleminfo" << std::to_string(mAreaSeed) << ".dat";
        output_file_for_elem_info = file_string_elem_info.str();
        mpElementInfoFile = output_file_handler.OpenOutputFile(output_file_for_elem_info, std::ios::app);
        if ( (t_now - 1.0*floor(t_now/1.0)) < dt )
        {
            *mpElementInfoFile << t_now <<"\n";
        }

        std::string output_file_for_cell_stress;
        // output_file_for_elem_info = "cellstress.dat";
        std::ostringstream file_string_cell_stress;
        file_string_cell_stress << "cellstress" << std::to_string(mAreaSeed) << ".dat";
        output_file_for_cell_stress = file_string_cell_stress.str();
        mpCellStressFile = output_file_handler.OpenOutputFile(output_file_for_cell_stress, std::ios::app);
        if ( (t_now - 1.0*floor(t_now/1.0)) < dt )
        {
            *mpCellStressFile << t_now <<"\n";
        }

        for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
            elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
            ++elem_iter)
        {
            unsigned elem_index = elem_iter->GetIndex();
            VertexElement<DIM, DIM>* p_element = p_cell_population->rGetMesh().GetElement(elem_index);
            CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(elem_index);
            unsigned num_nodes_elem = p_element->GetNumNodes();
            unsigned num_stressfibers_elem = p_element->GetNumStressfibers();
            double cell_perimeter = p_cell_population->rGetMesh().GetSurfaceAreaOfElement(elem_index);
            double cell_target_perimeter = p_cell->GetCellData()->GetItem("cell_target_perimeter");
            double stressfiber_elasticity =  p_cell->GetCellData()->GetItem("perimeter_elasticity");
            // double stressfiber_elasticity =  this->mSfStiffness;
            double perimeter_tension = p_cell->GetCellData()->GetItem("cell_perimeter_tension");
            
            if (((t_now - 1.0*floor(t_now/1.0)) < dt) && (p_cell->GetCellData()->GetItem("is_boundary_cell")==0))
            {
                *mpElementInfoFile << elem_index << " ";
                *mpElementInfoFile << p_cell->GetCellData()->GetItem("target_area") << " ";
                *mpElementInfoFile << p_cell_population->rGetMesh().GetVolumeOfElement(elem_index) << " ";
                *mpElementInfoFile << p_cell->GetCellData()->GetItem("target_shape_index") << " ";
                *mpElementInfoFile << perimeter_tension << " ";
                *mpElementInfoFile << num_nodes_elem << " ";
                *mpElementInfoFile << p_cell->GetCellData()->GetItem("perimeter_elasticity") << " ";
            }

            double max_energy_release_rate = 0.0;
            std::vector<bool> IsPeeling(num_nodes_elem, false);
            // update the end points ratio and associated nodes of stress fibers contained by the element 
            // (special treatments will be adopted if T1 transition happens)
            if (num_stressfibers_elem>0)
            {
                std::vector<unsigned> delete_stressfiber_local_index;
                
                // find the stress fibers need to be deleted and update the ratio of normal ones.
                for (unsigned stressfiber_local_index=0; stressfiber_local_index<num_stressfibers_elem; stressfiber_local_index++)
                {
                    assert(stressfiber_local_index < p_element->GetNumStressfibers());
                    VertexElement<DIM-1, DIM>* p_stressfiber = p_element->GetStressfiber(stressfiber_local_index);
                    Node<DIM>* p_node0 = p_stressfiber->GetStressfiberNode(0);
                    Node<DIM>* p_node1 = p_stressfiber->GetStressfiberNode(1);
                    Node<DIM>* p_node2 = p_stressfiber->GetStressfiberNode(2);
                    Node<DIM>* p_node3 = p_stressfiber->GetStressfiberNode(3);
                    c_vector<double,2> ratio = p_stressfiber->GetStressfiberEndpointsratio();

                    unsigned containing0 = 0;
                    unsigned containing1 = 0;
                    unsigned containing3 = 0;
                    for (unsigned local_index=0; local_index<num_nodes_elem; local_index++)
                    {
                        if (p_node0->GetIndex() == p_element->GetNodeGlobalIndex(local_index))
                        {
                            containing0 = 1;
                        }
                        if (p_node1->GetIndex() == p_element->GetNodeGlobalIndex(local_index))
                        {
                            containing1 = 1;
                        }
                        if (p_node3->GetIndex() == p_element->GetNodeGlobalIndex(local_index))
                        {
                            containing3 = 1;
                        }
                    }

                    unsigned local_index_node0 = p_element->GetNodeLocalIndex(p_node0->GetIndex());
                    unsigned local_index_node1 = p_element->GetNodeLocalIndex(p_node1->GetIndex());
                    unsigned local_index_node2 = p_element->GetNodeLocalIndex(p_node2->GetIndex());
                    unsigned local_index_node3 = p_element->GetNodeLocalIndex(p_node3->GetIndex());
                    unsigned local_index_next_node0 = (local_index_node0+1)%num_nodes_elem;
                    unsigned local_index_next_node2 = (local_index_node2+1)%num_nodes_elem;
                    bool node0_neighbour_to_node1 = true;
                    bool node2_neighbour_to_node3 = true;
                    bool node1_overlap_with_node2 = true;
                    if (local_index_node1 != local_index_next_node0)
                    {
                        node0_neighbour_to_node1 = false;
                    }
                    if (local_index_node1 != local_index_node2)
                    {
                        node1_overlap_with_node2 = false;
                    }
                    if (local_index_node3 != local_index_next_node2)
                    {
                        node2_neighbour_to_node3 = false;
                    }

                    if (containing1 == 0)
                    {
                        delete_stressfiber_local_index.push_back(stressfiber_local_index);
                    }
                    else
                    {
                        if (!node1_overlap_with_node2)
                        {
                            std::cout<<"Error:nodes 1 and 2 of stress fiber " << stressfiber_local_index << " is not overlapped to each other in elem " << elem_index << std::endl;
                            EXCEPTION("nodes 1 and 2 are not overlapped!");
                        }

                        if (containing0 == 0 || (!node0_neighbour_to_node1))
                        {
                            unsigned local_index_previous_node1 = (num_nodes_elem+local_index_node1-1)%num_nodes_elem;
                            Node<DIM>* new_p_node0 = p_element->GetNode(local_index_previous_node1);
                            p_stressfiber->UpdateStressfiberNode(new_p_node0,0);
                        }
                        // else
                        // {
                        //     if ((ratio[0]>0.0) && (ratio[0]<1.0))
                        //     {
                        //         double length_old0_new1 = norm_2(p_node1->rGetLocation()-p_node0->GetOldLocation());
                        //         double length_old01 = norm_2(p_node1->GetOldLocation()-p_node0->GetOldLocation());
                        //         double length_new01 = norm_2(p_node1->rGetLocation()-p_node0->rGetLocation());
                        //         ratio[0] = 1 - (length_old0_new1 - ratio[0]*length_old01)/length_new01;
                        //         ratio[0] = ratio[0]<0? 0:(ratio[0]>1? 1:ratio[0]);
                        //     }
                        // }

                        if (containing3 == 0 || (!node2_neighbour_to_node3))
                        {
                            unsigned local_index_next_node2 = (local_index_node2+1)%num_nodes_elem;
                            Node<DIM>* new_p_node3 = p_element->GetNode(local_index_next_node2);
                            p_stressfiber->UpdateStressfiberNode(new_p_node3,3);
                        }
                        // else
                        // {
                        //     if ((ratio[1]>0.0) && (ratio[1]<1.0))
                        //     {
                        //         double length_old2_new3 = norm_2(p_node3->rGetLocation()-p_node2->GetOldLocation());
                        //         double length_old23 = norm_2(p_node3->GetOldLocation()-p_node2->GetOldLocation());
                        //         double length_new23 = norm_2(p_node3->rGetLocation()-p_node2->rGetLocation());
                        //         ratio[1] = 1 - (length_old2_new3 - ratio[1]*length_old23)/length_new23;
                        //         ratio[1] = ratio[1]<0? 0:(ratio[1]>1? 1:ratio[1]);
                        //     }
                        // }
                        // p_stressfiber->UpdateStressfiberEndpointsratio(ratio);
                    }
                }


                // nucleate stress fibers when the node1 is newly added.
                for (unsigned local_index=0; local_index<num_nodes_elem; local_index++)
                {
                    bool contained_by_stress_fiber = false;
                    for (unsigned stressfiber_local_index=0; stressfiber_local_index<num_stressfibers_elem; stressfiber_local_index++)
                    {
                        assert(stressfiber_local_index < p_element->GetNumStressfibers());
                        VertexElement<DIM-1, DIM>* p_stressfiber = p_element->GetStressfiber(stressfiber_local_index);
                        Node<DIM>* p_node1 = p_stressfiber->GetStressfiberNode(1);

                        if (p_node1->GetIndex() == p_element->GetNodeGlobalIndex(local_index))
                        {
                            contained_by_stress_fiber = true;
                        }
                    }

                    if (!contained_by_stress_fiber)
                    {
                        unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
                        unsigned next_node_local_index = (local_index+1)%num_nodes_elem;
                        std::vector<Node<DIM>*> stressfiber_nodes;

                        stressfiber_nodes.push_back(p_element->GetNode(previous_node_local_index));
                        stressfiber_nodes.push_back(p_element->GetNode(local_index));
                        stressfiber_nodes.push_back(p_element->GetNode(local_index));
                        stressfiber_nodes.push_back(p_element->GetNode(next_node_local_index));

                        c_vector<double, DIM> vec_CD = p_element->GetNodeLocation(previous_node_local_index) - p_element->GetNodeLocation(local_index);
                        c_vector<double, DIM> vec_CE = p_element->GetNodeLocation(next_node_local_index) - p_element->GetNodeLocation(local_index);
                        c_vector<double,2> end_points_ratio = zero_vector<double>(2);
                        double StretchingLengthOfNucleation = mRestLengthOfNucleation/cell_target_perimeter*cell_perimeter;
                        end_points_ratio[0] = std::max((norm_2(vec_CD)-StretchingLengthOfNucleation)/norm_2(vec_CD),0.0);
                        end_points_ratio[1] = std::min(StretchingLengthOfNucleation/norm_2(vec_CE),1.0);

                        VertexElement<DIM-1,DIM>* p_stressfiber = new VertexElement<DIM-1,DIM>(global_index, stressfiber_nodes, end_points_ratio, mRestLengthOfNucleation);
                        p_element->AddStressfiber(p_stressfiber);
                        
                        global_index++;

                        std::cout << "a new stress fiber is created in elem " << elem_index << "." << std::endl;
                    }
                }


                // delete the stress fibers whose node 1 is no longer contained by this element.
                unsigned num_stress_fibers_waiting_deletion = delete_stressfiber_local_index.size();
                if (num_stress_fibers_waiting_deletion > 0)
                {
                    for (unsigned index = 0; index < num_stress_fibers_waiting_deletion; index++)
                    {
                        p_element->DeleteStressfiber(delete_stressfiber_local_index[index]);
                        std::cout << "stress fiber " << index << " in elem " << elem_index << " is deleted." << std::endl;
                    }
                }

                
                // peeling the stress fiber
                num_stressfibers_elem = p_element->GetNumStressfibers();
                for (unsigned stressfiber_local_index=0; stressfiber_local_index<num_stressfibers_elem; stressfiber_local_index++)
                {
                    assert(stressfiber_local_index < p_element->GetNumStressfibers());
                    VertexElement<DIM-1, DIM>* p_stressfiber = p_element->GetStressfiber(stressfiber_local_index);
                    Node<DIM>* p_node0 = p_stressfiber->GetStressfiberNode(0);
                    Node<DIM>* p_node1 = p_stressfiber->GetStressfiberNode(1);
                    Node<DIM>* p_node2 = p_stressfiber->GetStressfiberNode(2);
                    Node<DIM>* p_node3 = p_stressfiber->GetStressfiberNode(3);
                    c_vector<double,2> ratio = p_stressfiber->GetStressfiberEndpointsratio();
                    double restlength = p_stressfiber->GetStressfiberRestLength();

                    // calculate the force of a stress fiber
                    c_vector<double, DIM> locationA = (1-ratio[0])*p_node0->rGetLocation() + ratio[0]*p_node1->rGetLocation();
                    c_vector<double, DIM> locationB = (1-ratio[1])*p_node2->rGetLocation() + ratio[1]*p_node3->rGetLocation();
                    c_vector<double, DIM> vec_AB = locationB - locationA;
                    c_vector<double, DIM> vec_01 = p_node1->rGetLocation() - p_node0->rGetLocation();
                    c_vector<double, DIM> vec_23 = p_node3->rGetLocation() - p_node2->rGetLocation();
                    double length_AB = norm_2(vec_AB);
                    double length_01 = norm_2(vec_01);
                    double length_23 = norm_2(vec_23);
                    // double peeling_length = (1-ratio[0])*length_01 + ratio[1]*length_23;
                    // double sf_tension = std::max(stressfiber_elasticity*(cell_perimeter*length_AB/peeling_length-cell_target_perimeter), 0.0);
                    // double sf_tension = std::max(stressfiber_elasticity*(length_AB-cell_target_perimeter*peeling_length/cell_perimeter), 0.0);
                    double sf_tension = std::max(stressfiber_elasticity*(length_AB - restlength), 0.0);
                    if (ratio[0]==1.0 || ratio[1]==0.0)
                    {
                        sf_tension = 0.0;
                    }
                    double alpha = acos((vec_01[0]*vec_AB[0]+vec_01[1]*vec_AB[1])/(norm_2(vec_01)*norm_2(vec_AB)));
                    double beta = acos((vec_23[0]*vec_AB[0]+vec_23[1]*vec_AB[1])/(norm_2(vec_23)*norm_2(vec_AB)));
                    max_energy_release_rate = std::max(sf_tension*(1-cos(alpha)), max_energy_release_rate);
                    max_energy_release_rate = std::max(sf_tension*(1-cos(beta)), max_energy_release_rate);

                    // if ((sf_tension*(1-cos(alpha))>mAdhesionEnergy) && (sf_tension>mC0))
                    // {
                    //     double peeling_rate_A = pow(1/mk*(sf_tension/mC0-1), mRatePower);
                    //     double peeling_length_step_A = peeling_rate_A*cos(alpha)*dt;
                    //     ratio[0] = ratio[0] - peeling_length_step_A/length_01;
                    //     ratio[0] = ratio[0]<0? 0:(ratio[0]>1? 1:ratio[0]);

                    //     IsPeeling[p_element->GetNodeLocalIndex(p_node1->GetIndex())] = true;
                    //     std::cout<<"peeling happens in elem "<< elem_index << "." << std::endl;
                    // }
                    // if ((sf_tension*(1-cos(beta))>mAdhesionEnergy) && (sf_tension>mC0))
                    // {
                    //     double peeling_rate_B = pow(1/mk*(sf_tension/mC0-1), mRatePower);
                    //     double peeling_length_step_B = peeling_rate_B*cos(beta)*dt;
                    //     ratio[1] = ratio[1] + peeling_length_step_B/length_23;
                    //     ratio[1] = ratio[1]<0? 0:(ratio[1]>1? 1:ratio[1]);

                    //     IsPeeling[p_element->GetNodeLocalIndex(p_node1->GetIndex())] = true;
                    //     std::cout<<"peeling happens in elem "<< elem_index << "." << std::endl;
                    // }

                    if (sf_tension*(1-cos(alpha))>mAdhesionEnergy && ratio[0]>0.0) 
                    {
                        double peeling_rate_A = pow(1/mk*(sf_tension*(1-cos(alpha))/mAdhesionEnergy-1), 1/mRatePower);
                        double peeling_length_step_A = peeling_rate_A*dt;
                        ratio[0] = ratio[0] - peeling_length_step_A/length_01;
                        ratio[0] = ratio[0]<0? 0:(ratio[0]>1? 1:ratio[0]);
                        restlength += peeling_length_step_A/cell_perimeter*cell_target_perimeter;

                        IsPeeling[p_element->GetNodeLocalIndex(p_node1->GetIndex())] = true;
                        p_cell->GetCellData()->SetItem("peeling_flag",1);
                        // std::cout<<"peeling happens in elem "<< elem_index << "." << std::endl;
                    }
                    if (sf_tension*(1-cos(beta))>mAdhesionEnergy && ratio[1]<1.0)
                    {
                        double peeling_rate_B = pow(1/mk*(sf_tension*(1-cos(beta))/mAdhesionEnergy-1), 1/mRatePower);
                        double peeling_length_step_B = peeling_rate_B*dt;
                        ratio[1] = ratio[1] + peeling_length_step_B/length_23;
                        ratio[1] = ratio[1]<0? 0:(ratio[1]>1? 1:ratio[1]);
                        restlength += peeling_length_step_B/cell_perimeter*cell_target_perimeter;

                        IsPeeling[p_element->GetNodeLocalIndex(p_node1->GetIndex())] = true;
                        p_cell->GetCellData()->SetItem("peeling_flag",1);
                        // std::cout<<"peeling happens in elem "<< elem_index << "." << std::endl;
                    }
                    p_stressfiber->UpdateStressfiberEndpointsratio(ratio);
                    p_stressfiber->UpdateStressfiberRestLength(restlength);
                }

            }
            p_cell->GetCellData()->SetItem("max_energy_release_rate", max_energy_release_rate);

            // record the info of opening angle and bisector orientation and create stress fibers if the tension is higher than a threshold
            for (unsigned local_index=0; local_index<num_nodes_elem; local_index++)
            {
                //  (next node)E---C(this node)
                //                 |
                //                 D(previous node)
                // record the angle and bisector orientation
                unsigned node_global_index = p_element->GetNodeGlobalIndex(local_index);
                unsigned next_node_local_index = (local_index+1)%num_nodes_elem;
                unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
                c_vector<double, DIM> vec_CD = p_element->GetNodeLocation(previous_node_local_index) - p_element->GetNodeLocation(local_index);
                c_vector<double, DIM> vec_CE = p_element->GetNodeLocation(next_node_local_index) - p_element->GetNodeLocation(local_index);
                double angle = acos((vec_CD[0]*vec_CE[0]+vec_CD[1]*vec_CE[1])/(norm_2(vec_CD)*norm_2(vec_CE)))*180/M_PI;
                c_vector<double, DIM> bisector = vec_CD/norm_2(vec_CD) + vec_CE/norm_2(vec_CE);
                double bisector_orientation = asin(bisector[0]/norm_2(bisector))*180/M_PI;

                oss.str("");
                oss << "node_" << local_index;
                name_item = oss.str();
                p_cell->GetCellData()->SetItem(name_item, node_global_index);

                oss.str("");
                oss << "angle_" << local_index;
                name_item = oss.str();
                p_cell->GetCellData()->SetItem(name_item, angle);

                oss.str("");
                oss << "bisector_" << local_index;
                name_item = oss.str();
                p_cell->GetCellData()->SetItem(name_item, bisector_orientation);

                if (((t_now - 1.0*floor(t_now/1.0)) < dt) && (p_cell->GetCellData()->GetItem("is_boundary_cell")==0))
                {
                    *mpAngleBisectorFile << elem_index  << " ";
                    *mpAngleBisectorFile << local_index  << " ";
                    *mpAngleBisectorFile << angle  << " ";
                    *mpAngleBisectorFile << bisector_orientation  << " ";
                    *mpAngleBisectorFile << "\n";
                }

                // nucleate stress fibers if the perimeter tension of the element is higher than a threshold
                if ( (num_stressfibers_elem == 0) && (perimeter_tension>this->mNucleationPerimeterTension) )
                {
                    std::vector<Node<DIM>*> stressfiber_nodes;

                    stressfiber_nodes.push_back(p_element->GetNode(previous_node_local_index));
                    stressfiber_nodes.push_back(p_element->GetNode(local_index));
                    stressfiber_nodes.push_back(p_element->GetNode(local_index));
                    stressfiber_nodes.push_back(p_element->GetNode(next_node_local_index));

                    c_vector<double,2> end_points_ratio = zero_vector<double>(2);
                    double StretchingLengthOfNucleation = mRestLengthOfNucleation/cell_target_perimeter*cell_perimeter;
                    end_points_ratio[0] = std::max((norm_2(vec_CD)-StretchingLengthOfNucleation)/norm_2(vec_CD),0.0);
                    end_points_ratio[1] = std::min(StretchingLengthOfNucleation/norm_2(vec_CE),1.0);

                    VertexElement<DIM-1,DIM>* p_stressfiber = new VertexElement<DIM-1,DIM>(global_index, stressfiber_nodes, end_points_ratio, mRestLengthOfNucleation);
                    p_element->AddStressfiber(p_stressfiber);
                    // std::cout<< p_element->GetNumStressfibers() << " stress fibers are created in element "<<p_element->GetIndex()<<std::endl;
                    
                    p_cell->GetCellData()->SetItem("nucleation_flag",mFlag+1);

                    global_index++;
                }
            }

            // calculate the stress fibers contribution
            num_stressfibers_elem = p_element->GetNumStressfibers();
            if (num_stressfibers_elem>0)
            {
                //     3---B---2
                //         |   |
                //         |   |
                //     0---A---1
                for (unsigned stressfiber_local_index=0; stressfiber_local_index<num_stressfibers_elem; stressfiber_local_index++)
                {
                    assert(stressfiber_local_index < p_element->GetNumStressfibers());
                    VertexElement<DIM-1, DIM>* p_stressfiber = p_element->GetStressfiber(stressfiber_local_index);
                    Node<DIM>* p_node0 = p_stressfiber->GetStressfiberNode(0);
                    Node<DIM>* p_node1 = p_stressfiber->GetStressfiberNode(1);
                    Node<DIM>* p_node2 = p_stressfiber->GetStressfiberNode(2);
                    Node<DIM>* p_node3 = p_stressfiber->GetStressfiberNode(3);
                    c_vector<double,2> ratio = p_stressfiber->GetStressfiberEndpointsratio();
                    double restlength = p_stressfiber->GetStressfiberRestLength();

                    c_vector<double, DIM> locationA = (1-ratio[0])*p_node0->rGetLocation() + ratio[0]*p_node1->rGetLocation();
                    c_vector<double, DIM> locationB = (1-ratio[1])*p_node2->rGetLocation() + ratio[1]*p_node3->rGetLocation();
                    c_vector<double, DIM> vec_AB = locationB - locationA;
                    // c_vector<double, DIM> vec_01 = p_node1->rGetLocation() - p_node0->rGetLocation();
                    // c_vector<double, DIM> vec_23 = p_node3->rGetLocation() - p_node2->rGetLocation();
                    double length_AB = norm_2(vec_AB);
                    // double length_01 = norm_2(vec_01);
                    // double length_23 = norm_2(vec_23);

                    // double peeling_length = (1-ratio[0])*length_01 + ratio[1]*length_23;
                    // double sf_tension = std::max(stressfiber_elasticity*(cell_perimeter*length_AB/peeling_length-cell_target_perimeter), 0.0);
                    // double sf_tension = std::max(stressfiber_elasticity*(length_AB-cell_target_perimeter*peeling_length/cell_perimeter), 0.0);
                    double sf_tension = std::max(stressfiber_elasticity*(length_AB - restlength), 0.0);
                    if (ratio[0]==1.0 || ratio[1]==0.0)
                    {
                        sf_tension = 0.0;
                    }

                    c_vector<double, DIM> sf_contribution_0 = sf_tension*(1-ratio[0])*vec_AB/length_AB;
                    c_vector<double, DIM> sf_contribution_1 = sf_tension*ratio[0]*vec_AB/length_AB;
                    c_vector<double, DIM> sf_contribution_2 = -sf_tension*(1-ratio[1])*vec_AB/length_AB;
                    c_vector<double, DIM> sf_contribution_3 = -sf_tension*ratio[1]*vec_AB/length_AB;
                    p_cell_population->GetNode(p_node0->GetIndex())->AddAppliedForceContribution(sf_contribution_0);
                    p_cell_population->GetNode(p_node1->GetIndex())->AddAppliedForceContribution(sf_contribution_1);
                    p_cell_population->GetNode(p_node2->GetIndex())->AddAppliedForceContribution(sf_contribution_2);
                    p_cell_population->GetNode(p_node3->GetIndex())->AddAppliedForceContribution(sf_contribution_3);

                    OutputFileHandler output_file_handler(this->mOutputDirectory+"/", false);
                    std::string output_file_for_peeling_bisector_orientation;
                    // output_file_for_peeling_bisector_orientation = "peelingbisectororientation.dat";
                    std::ostringstream file_string_peeling_bisector_orientation;
                    file_string_peeling_bisector_orientation << "peelingbisectororientation" << std::to_string(mAreaSeed) << ".dat";
                    output_file_for_peeling_bisector_orientation = file_string_peeling_bisector_orientation.str();
                    mpPeelingBisectorOrientationFile = output_file_handler.OpenOutputFile(output_file_for_peeling_bisector_orientation, std::ios::app);

                    unsigned node1_local_index = p_element->GetNodeLocalIndex(p_node1->GetIndex());
                    if ( (IsPeeling[node1_local_index])&&((t_now - 1.0*floor(t_now/1.0)) < dt) && (p_cell->GetCellData()->GetItem("is_boundary_cell")==0))
                    {
                        *mpPeelingBisectorOrientationFile << t_now << " ";
                        *mpPeelingBisectorOrientationFile << elem_index  << " ";
                        *mpPeelingBisectorOrientationFile << stressfiber_local_index  << " ";
                        *mpPeelingBisectorOrientationFile << p_cell_population->rGetMesh().GetVolumeOfElement(elem_index)  << " ";

                        oss.str("");
                        oss << "angle_" << node1_local_index;
                        name_item = oss.str();
                        *mpPeelingBisectorOrientationFile << p_cell->GetCellData()->GetItem(name_item)  << " ";

                        oss.str("");
                        oss << "bisector_" << node1_local_index;
                        name_item = oss.str();
                        *mpPeelingBisectorOrientationFile << p_cell->GetCellData()->GetItem(name_item) << " ";

                        *mpPeelingBisectorOrientationFile << ratio[0] << " ";
                        *mpPeelingBisectorOrientationFile << ratio[1] << " ";
                        *mpPeelingBisectorOrientationFile << sf_tension << " ";
                        *mpPeelingBisectorOrientationFile << "\n";
                    }

                    // OutputFileHandler output_file_handler(this->mOutputDirectory+"/", false);
                    std::string output_file_for_stress_fiber_tension;
                    // output_file_for_stress_fiber_tension = "stressfibertension.dat";
                    std::ostringstream file_string_fiber_tension;
                    file_string_fiber_tension << "stressfibertension" << std::to_string(mAreaSeed) << ".dat";
                    output_file_for_stress_fiber_tension = file_string_fiber_tension.str();
                    mpStressFiberTensionFile = output_file_handler.OpenOutputFile(output_file_for_stress_fiber_tension, std::ios::app);

                    if (((t_now - 1.0*floor(t_now/1.0)) < dt) && (p_cell->GetCellData()->GetItem("is_boundary_cell")==0))
                    {
                        *mpStressFiberTensionFile << t_now << " ";
                        *mpStressFiberTensionFile << elem_index << " ";
                        *mpStressFiberTensionFile << stressfiber_local_index  << " ";

                        oss.str("");
                        oss << "angle_" << node1_local_index;
                        name_item = oss.str();
                        *mpStressFiberTensionFile << p_cell->GetCellData()->GetItem(name_item)  << " ";

                        oss.str("");
                        oss << "bisector_" << node1_local_index;
                        name_item = oss.str();
                        *mpStressFiberTensionFile << p_cell->GetCellData()->GetItem(name_item) << " ";

                        *mpStressFiberTensionFile << asin(vec_AB[0]/length_AB)*180/M_PI << " ";
                        *mpStressFiberTensionFile << length_AB << " ";
                        *mpStressFiberTensionFile << sf_tension << " ";
                        *mpStressFiberTensionFile << ratio[0] << " ";
                        *mpStressFiberTensionFile << ratio[1] << " ";
                        *mpStressFiberTensionFile << "\n";
                    }
                }
            }

            UpdateStressStateOfCell(rCellPopulation, p_cell);

        }
        if ( (t_now - 1.0*floor(t_now/1.0)) < dt )
        {
            *mpElementInfoFile <<"\n";
            *mpCellStressFile <<"\n";
        }
    }

    /*
    // stressfiber check for each element
    for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
        elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
        ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        VertexElement<DIM, DIM>* p_element = p_cell_population->rGetMesh().GetElement(elem_index);
        // CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(elem_index);
        // p_cell->GetCellData()->SetItem("nucleation_flag",mFlag);
        unsigned num_nodes_elem = p_element->GetNumNodes();
        unsigned num_stressfibers_elem = p_element->GetNumStressfibers();

        if (num_stressfibers_elem>0)
        {
            if (num_nodes_elem != num_stressfibers_elem)
            {
                std::cout << "No of nodes and stress fibers are not equal!" << std::endl;
                EXCEPTION("No of nodes shall be equal to that of stress fibers.");
            }
            //     3---B---2
            //         |   |
            //         |   |
            //     0---A---1
            for (unsigned stressfiber_local_index=0; stressfiber_local_index<num_stressfibers_elem; stressfiber_local_index++)
            {
                assert(stressfiber_local_index < p_element->GetNumStressfibers());
                VertexElement<DIM-1, DIM>* p_stressfiber = p_element->GetStressfiber(stressfiber_local_index);
                Node<DIM>* p_node0 = p_stressfiber->GetStressfiberNode(0);
                Node<DIM>* p_node1 = p_stressfiber->GetStressfiberNode(1);
                Node<DIM>* p_node2 = p_stressfiber->GetStressfiberNode(2);
                Node<DIM>* p_node3 = p_stressfiber->GetStressfiberNode(3);
                c_vector<double,2> ratio = p_stressfiber->GetStressfiberEndpointsratio();

                unsigned local_index_node0 = p_element->GetNodeLocalIndex(p_node0->GetIndex());
                unsigned local_index_node1 = p_element->GetNodeLocalIndex(p_node1->GetIndex());
                unsigned local_index_node2 = p_element->GetNodeLocalIndex(p_node2->GetIndex());
                unsigned local_index_node3 = p_element->GetNodeLocalIndex(p_node3->GetIndex());
                unsigned local_index_next_node0 = (local_index_node0+1)%num_nodes_elem;
                unsigned local_index_next_node2 = (local_index_node2+1)%num_nodes_elem;

                if (local_index_node1 != local_index_next_node0)
                {
                    std::cout<<"Error:nodes 0 and 1 of stress fiber " << stressfiber_local_index << " is not neighbor to each other in elem " << elem_index << std::endl;
                    EXCEPTION("node 0 is not neighbour to node 1.");
                }
                if (local_index_node1 != local_index_node2)
                {
                    std::cout<<"Error:nodes 1 and 2 of stress fiber " << stressfiber_local_index << " is not overlapped to each other in elem " << elem_index << std::endl;
                    EXCEPTION("node 1 is not overlapped with node 2.");
                }
                if (local_index_node3 != local_index_next_node2)
                {
                    std::cout<<"Error:nodes 2 and 3 of stress fiber " << stressfiber_local_index << " is not neighbor to each other in elem " << elem_index << std::endl;
                    EXCEPTION("node 3 is not neighbour to node 2.");
                }
                if (ratio[0]>1.0 || ratio[0]<0.0)
                {
                    std::cout<<"ratio[0] is not reasonable!"<<std::endl;
                    EXCEPTION("ratio 0 is not reasonable.");
                }         
                if (ratio[1]>1.0 || ratio[1]<0.0)
                {
                    std::cout<<"ratio[1] is not reasonable!"<<std::endl;
                    EXCEPTION("ratio 1 is not reasonable.");
                }         
            }
        }
    }
    */
}

template<unsigned DIM>
void MyStressfiberTensionForce<DIM>::UpdateStressStateOfCell(AbstractCellPopulation<DIM,DIM>& rCellPopulation, CellPtr pCell)
{   
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("MyStressfiberTensionForce is to be used with a VertexBasedCellPopulation only");
    }

    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    unsigned elem_index = rCellPopulation.GetLocationIndexUsingCell(pCell);
    VertexElement<DIM, DIM>* p_element = p_cell_population->rGetMesh().GetElement(elem_index);

    double cell_stress_XX = 0.0;
    double cell_stress_YY = 0.0;
    double cell_stress_XY = 0.0;
    double cell_stress_YX = 0.0;
    double shape_tensor_XX = 0.0;
    double shape_tensor_YY = 0.0;
    double shape_tensor_XY = 0.0;

    double cell_area = p_cell_population->rGetMesh().GetVolumeOfElement(elem_index);
    double target_area = pCell->GetCellData()->GetItem("target_area");
    double perimeter_elasticity =  pCell->GetCellData()->GetItem("perimeter_elasticity");
    double cell_perimeter = p_cell_population->rGetMesh().GetSurfaceAreaOfElement(elem_index);
    c_vector<double, DIM> centroid = p_cell_population->rGetMesh().GetCentroidOfElement(elem_index);
    
    unsigned local_index_lowest_vertex = 0;
    double y_coor_lowest_vertex = 10000.0;
    for (unsigned local_index =0; local_index < p_element->GetNumNodes(); local_index++)
    {
        Node<DIM>* pNode = p_element->GetNode(local_index);
        if (pNode->rGetLocation()[1]<y_coor_lowest_vertex)
        {
            local_index_lowest_vertex = local_index;
            y_coor_lowest_vertex = pNode->rGetLocation()[1];
        }
    }
    pCell->GetCellData()->SetItem("LocalIndexOfLowestVertex", local_index_lowest_vertex);

    for (unsigned local_index =0; local_index < p_element->GetNumNodes(); local_index++)
    {
        //  --B
        //    |
        // C--A
        Node<DIM>* pNodeC = p_element->GetNode((local_index-1+p_element->GetNumNodes())%p_element->GetNumNodes());
        Node<DIM>* pNodeA = p_element->GetNode(local_index);
        Node<DIM>* pNodeB = p_element->GetNode((local_index+1)%p_element->GetNumNodes());
        c_vector<double,DIM> locationC = pNodeC->rGetLocation();
        c_vector<double,DIM> locationA = pNodeA->rGetLocation();
        c_vector<double,DIM> locationB = pNodeB->rGetLocation();
        double l_CA = p_cell_population->rGetMesh().GetDistanceBetweenNodes(pNodeC->GetIndex(), pNodeA->GetIndex());
        double l_AB = p_cell_population->rGetMesh().GetDistanceBetweenNodes(pNodeA->GetIndex(), pNodeB->GetIndex());
        c_vector<double,DIM> vec_CA =  p_cell_population->rGetMesh().GetVectorFromAtoB(locationC, locationA);
        c_vector<double,DIM> vec_AB =  p_cell_population->rGetMesh().GetVectorFromAtoB(locationA, locationB);

        std::set<unsigned> elements_containing_nodeC = pNodeC->rGetContainingElementIndices();
        std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
        std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();
        // Find common elements
        std::set<unsigned> shared_elements_CA;
        std::set_intersection(elements_containing_nodeC.begin(),
                            elements_containing_nodeC.end(),
                            elements_containing_nodeA.begin(),
                            elements_containing_nodeA.end(),
                            std::inserter(shared_elements_CA, shared_elements_CA.begin()));
        // Check that the nodes have a common edge
        assert(!shared_elements_CA.empty());
        if (shared_elements_CA.size() >= 3)
        {
            std::cout<< std::endl << "Get error in MyStressfiberTensionForce::UpdateStressStateOfCell";
            std::cout<< std::endl << "Get shared elements more than 2";
        }
        // Find common elements
        std::set<unsigned> shared_elements_AB;
        std::set_intersection(elements_containing_nodeA.begin(),
                            elements_containing_nodeA.end(),
                            elements_containing_nodeB.begin(),
                            elements_containing_nodeB.end(),
                            std::inserter(shared_elements_AB, shared_elements_AB.begin()));
        // Check that the nodes have a common edge
        assert(!shared_elements_AB.empty());
        if (shared_elements_AB.size() >= 3)
        {
            std::cout<< std::endl << "Get error in MyStressfiberTensionForce::UpdateStressStateOfCell";
            std::cout<< std::endl << "Get shared elements more than 2";
        }

        // force contributions on the vertex from the element
        // force from the area term
        c_vector<double,DIM> element_area_gradient = p_cell_population->rGetMesh().GetAreaGradientOfElementAtNode(p_element, local_index);
        c_vector<double,DIM> F_Press = -(cell_area-target_area)*element_area_gradient;

        // Tension from perimeter term
        c_vector<double,DIM> perimeter_gradient_at_node = p_cell_population->rGetMesh().GetPerimeterGradientOfElementAtNode(p_element, local_index);
        c_vector<double, DIM> F_Tens = -perimeter_elasticity*cell_perimeter*perimeter_gradient_at_node;

        // Cell-cell adhesion
        double Lambda = mCellCellAdhesionEnergyParameter;
        c_vector<double,DIM> F_Adhe = -Lambda*(vec_CA/l_CA -vec_AB/l_AB);

        c_vector<double,DIM> Force = F_Press +F_Tens +F_Adhe;

        // next: stress calculation
        c_vector<double,DIM> node_location_ralate_to_centroid =  p_cell_population->rGetMesh().GetVectorFromAtoB(centroid, locationA);
        cell_stress_XX += -1.0/cell_area*node_location_ralate_to_centroid[0]*Force[0];
        cell_stress_YY += -1.0/cell_area*node_location_ralate_to_centroid[1]*Force[1];
        cell_stress_XY += -1.0/cell_area*node_location_ralate_to_centroid[0]*Force[1];
        cell_stress_YX += -1.0/cell_area*node_location_ralate_to_centroid[1]*Force[0];

        shape_tensor_XX += 1.0/p_element->GetNumNodes()*node_location_ralate_to_centroid[0]*node_location_ralate_to_centroid[0];
        shape_tensor_YY += 1.0/p_element->GetNumNodes()*node_location_ralate_to_centroid[1]*node_location_ralate_to_centroid[1];
        shape_tensor_XY += 1.0/p_element->GetNumNodes()*node_location_ralate_to_centroid[0]*node_location_ralate_to_centroid[1];

        if (false)
        {
            std::string name_item;
            std::ostringstream oss;

            oss.str("");
            oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                    << "_F_Press_x";
            name_item = oss.str();
            pCell->GetCellData()->SetItem(name_item, F_Press[0]);
            oss.str("");
            oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                    << "_F_Press_y";
            name_item = oss.str();
            pCell->GetCellData()->SetItem(name_item, F_Press[1]);

            oss.str("");
            oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                    << "_F_Tens_x";
            name_item = oss.str();
            pCell->GetCellData()->SetItem(name_item, F_Tens[0]);
            oss.str("");
            oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                    << "_F_Tens_y";
            name_item = oss.str();
            pCell->GetCellData()->SetItem(name_item, F_Tens[1]);

            oss.str("");
            oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                    << "_F_Adhe_x";
            name_item = oss.str();
            pCell->GetCellData()->SetItem(name_item, F_Adhe[0]);
            oss.str("");
            oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                    << "_F_Adhe_y";
            name_item = oss.str();
            pCell->GetCellData()->SetItem(name_item, F_Adhe[1]);

            oss.str("");
            oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                    << "_Fx";
            name_item = oss.str();
            pCell->GetCellData()->SetItem(name_item, Force[0]);
            oss.str("");
            oss << (local_index+p_element->GetNumNodes()-local_index_lowest_vertex)%(p_element->GetNumNodes()) 
                    << "_Fy";
            name_item = oss.str();
            pCell->GetCellData()->SetItem(name_item, Force[1]);
        }
    } // end of iteration of vertices of the element, for calculation of force contributions and stress summation.

    double cell_stress_delta = (cell_stress_XX - cell_stress_YY)*(cell_stress_XX - cell_stress_YY) + 4*cell_stress_XY*cell_stress_XY;
    double stress_1 = 0.5*(cell_stress_XX + cell_stress_YY + sqrt(cell_stress_delta));
    double stress_2 = 0.5*(cell_stress_XX + cell_stress_YY - sqrt(cell_stress_delta));

    c_vector<double, DIM> stress_eigen_vec_1 = zero_vector<double>(DIM);
    double principal_axis_of_stress = 0.0;
    if (stress_1==cell_stress_XX)
        principal_axis_of_stress = 0.0;
    else if (stress_1==cell_stress_YY)
        principal_axis_of_stress = M_PI/2;
    else
    {
        stress_eigen_vec_1[0] = 1.0;
        stress_eigen_vec_1[1] = (stress_1-cell_stress_XX)/cell_stress_XY;
        principal_axis_of_stress = atan(stress_eigen_vec_1[1]/stress_eigen_vec_1[0]);
    }

    double shape_delta = (shape_tensor_XX-shape_tensor_YY)*(shape_tensor_XX-shape_tensor_YY) +4*shape_tensor_XY*shape_tensor_XY;
    double shape_lambda1 = 0.5*(shape_tensor_XX+shape_tensor_YY+sqrt(shape_delta));
    double shape_lambda2 = 0.5*(shape_tensor_XX+shape_tensor_YY-sqrt(shape_delta));
    double aspect_ratio = sqrt(shape_lambda1/shape_lambda2);

    c_vector<double, DIM> shape_eigen_vec_1 = zero_vector<double>(DIM);
    double principal_axis_of_shape = 0.0;
    if (shape_lambda1==shape_tensor_XX)
        principal_axis_of_shape = 0.0;
    else if (shape_lambda1==shape_tensor_YY)
        principal_axis_of_shape = M_PI/2;
    else
    {
        shape_eigen_vec_1[0] = 1.0;
        shape_eigen_vec_1[1] = (shape_lambda1-shape_tensor_XX)/shape_tensor_XY;
        principal_axis_of_shape = atan(shape_eigen_vec_1[1]/shape_eigen_vec_1[0]);
    }

    double t_now = SimulationTime::Instance()->GetTime();
    double dt = SimulationTime::Instance()->GetTimeStep();
    if (pCell->GetCellData()->GetItem("is_boundary_cell")==0)
    {
        pCell->GetCellData()->SetItem("StressXX", cell_stress_XX);
        pCell->GetCellData()->SetItem("StressYY", cell_stress_YY);
        pCell->GetCellData()->SetItem("StressXY", cell_stress_XY);
        pCell->GetCellData()->SetItem("Stress1", stress_1);
        pCell->GetCellData()->SetItem("Stress2", stress_2);
        pCell->GetCellData()->SetItem("PrincipalAxisOfStress", principal_axis_of_stress);
        pCell->GetCellData()->SetItem("shape_tensor_xx", shape_tensor_XX);
        pCell->GetCellData()->SetItem("shape_tensor_yy", shape_tensor_YY);
        pCell->GetCellData()->SetItem("shape_tensor_xy", shape_tensor_XY);
        pCell->GetCellData()->SetItem("AspectRatio", aspect_ratio);
        pCell->GetCellData()->SetItem("PrincipalAxisOfShape", principal_axis_of_shape);

        if ((t_now - 1.0*floor(t_now/1.0)) < dt)
        {
            OutputFileHandler output_file_handler(this->mOutputDirectory+"/", false);
            std::string output_file_for_cell_stress;
            // output_file_for_elem_info = "cellstress.dat";
            std::ostringstream file_string_cell_stress;
            file_string_cell_stress << "cellstress" << std::to_string(mAreaSeed) << ".dat";
            output_file_for_cell_stress = file_string_cell_stress.str();
            mpCellStressFile = output_file_handler.OpenOutputFile(output_file_for_cell_stress, std::ios::app);
            *mpCellStressFile << elem_index <<" ";
            *mpCellStressFile << cell_area <<" ";
            *mpCellStressFile << cell_stress_XX <<" ";
            *mpCellStressFile << cell_stress_YY <<" ";
            *mpCellStressFile << cell_stress_XY <<" ";
            *mpCellStressFile << stress_1 <<" ";
            *mpCellStressFile << stress_2 <<" ";
            *mpCellStressFile << principal_axis_of_stress <<" ";
        }
    }
    else
    {
        pCell->GetCellData()->SetItem("StressXX", 0);
        pCell->GetCellData()->SetItem("StressYY", 0);
        pCell->GetCellData()->SetItem("StressXY", 0);
        pCell->GetCellData()->SetItem("Stress1", 0);
        pCell->GetCellData()->SetItem("Stress2", 0);
        pCell->GetCellData()->SetItem("PrincipalAxisOfStress", 0);
        pCell->GetCellData()->SetItem("shape_tensor_xx", 0);
        pCell->GetCellData()->SetItem("shape_tensor_yy", 0);
        pCell->GetCellData()->SetItem("shape_tensor_xy", 0);
        pCell->GetCellData()->SetItem("AspectRatio", 1);
        pCell->GetCellData()->SetItem("PrincipalAxisOfShape", 0);
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
void MyStressfiberTensionForce<DIM>::SetAreaSeed(unsigned areaSeed)
{
    mAreaSeed = areaSeed;
}

template<unsigned DIM>
void MyStressfiberTensionForce<DIM>::SetIfEquilibrateForAWhile(bool ifEquilibrateForAWhile)
{
    mIfEquilibrateForAWhile = ifEquilibrateForAWhile;
}

template<unsigned DIM>
void MyStressfiberTensionForce<DIM>::SetStartTimeForStretching(double startTimeForStretching)
{
    mStartTimeForStretching = startTimeForStretching;
}

template<unsigned DIM>
void MyStressfiberTensionForce<DIM>::SetFlagForStressfiberCreation(unsigned flag)
{
    mFlag = flag;
}

template<unsigned DIM>
void MyStressfiberTensionForce<DIM>::SetStressfiberStiffness(double sfStiffness)
{
    mSfStiffness = sfStiffness;
}

template<unsigned DIM>
void MyStressfiberTensionForce<DIM>::SetNucleationThresholdOfPerimeterTension(double nucleationPerimeterTension)
{
    mNucleationPerimeterTension = nucleationPerimeterTension;
}

template<unsigned DIM>
void MyStressfiberTensionForce<DIM>::SetHalfWidth(double halfWidth)
{
    mHalfWidth = halfWidth;
}

template<unsigned DIM>
void MyStressfiberTensionForce<DIM>::SetRestLengthOfNucleation(double restLengthOfNucleation)
{
    mRestLengthOfNucleation = restLengthOfNucleation;
}

template<unsigned DIM>
void MyStressfiberTensionForce<DIM>::SetPeelingParameters(double adhesionEnergy, double k, double C0, double ratePower)
{
    mAdhesionEnergy = adhesionEnergy;
    mk = k;
    mC0 = C0;
    mRatePower = ratePower;
}

template<unsigned DIM>
void MyStressfiberTensionForce<DIM>::SetNagaiHondaCellCellAdhesionEnergyParameter(double cellCellAdhesionEnergyParameter)
{
    mCellCellAdhesionEnergyParameter = cellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
void MyStressfiberTensionForce<DIM>::SetOutputDirectory(std::string outputDirectory)
{
    mOutputDirectory = outputDirectory;
}

// Explicit instantiation
template class MyStressfiberTensionForce<1>;
template class MyStressfiberTensionForce<2>;
template class MyStressfiberTensionForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MyStressfiberTensionForce)
