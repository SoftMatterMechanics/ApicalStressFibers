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

    std::string name_item;
    std::ostringstream oss;

    if ( !(mIfEquilibrateForAWhile && (t_now>=mTimeForEquilibrium)) )
    {
        for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
            elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
            ++elem_iter)
        {
            //  (next node)E---C(this node)
            //                 |
            //                 D(previous node)
            unsigned elem_index = elem_iter->GetIndex();
            VertexElement<DIM, DIM>* p_element = p_cell_population->rGetMesh().GetElement(elem_index);
            CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(elem_index);
            
            unsigned num_nodes_elem = p_element->GetNumNodes();
            for (unsigned local_index=0; local_index<num_nodes_elem; local_index++)
            {
                unsigned node_global_index = p_element->GetNodeGlobalIndex(local_index);

                unsigned next_node_local_index = (local_index+1)%num_nodes_elem;
                unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
                c_vector<double, DIM> vec_CD = p_element->GetNodeLocation(previous_node_local_index) - p_element->GetNodeLocation(local_index);
                c_vector<double, DIM> vec_CE = p_element->GetNodeLocation(next_node_local_index) - p_element->GetNodeLocation(local_index);
                double angle = acos((vec_CD[0]*vec_CE[0]+vec_CD[1]*vec_CE[1])/(norm_2(vec_CD)*norm_2(vec_CE)))*180/M_PI;

                oss.str("");
                oss << "node_" << local_index;
                name_item = oss.str();
                p_cell->GetCellData()->SetItem(name_item, node_global_index);

                oss.str("");
                oss << "angle_" << local_index;
                name_item = oss.str();
                p_cell->GetCellData()->SetItem(name_item, angle);

                p_cell->GetCellData()->SetItem("nucleation_flag",mFlag);
            }
        }

        // Define some helper variables
        // VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
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
    else
    {
        unsigned global_index = 0;
        for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
            elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
            ++elem_iter)
        {
            unsigned elem_index = elem_iter->GetIndex();
            VertexElement<DIM, DIM>* p_element = p_cell_population->rGetMesh().GetElement(elem_index);
            CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(elem_index);

            double max_angle_diff = -1000.0;
            unsigned max_angle_diff_node = 0;
            double nucleation_angle = 0.0;
            double bisector_orientation = 0.0;
            unsigned num_nodes_elem = p_element->GetNumNodes();
            for (unsigned local_index=0; local_index<num_nodes_elem; local_index++)
            {
                unsigned node_global_index = p_element->GetNodeGlobalIndex(local_index);

                unsigned next_node_local_index = (local_index+1)%num_nodes_elem;
                unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
                c_vector<double, DIM> vec_CD = p_element->GetNodeLocation(previous_node_local_index) - p_element->GetNodeLocation(local_index);
                c_vector<double, DIM> vec_CE = p_element->GetNodeLocation(next_node_local_index) - p_element->GetNodeLocation(local_index);
                double angle = acos((vec_CD[0]*vec_CE[0]+vec_CD[1]*vec_CE[1])/(norm_2(vec_CD)*norm_2(vec_CE)))*180/M_PI;

                if (angle>90)
                {
                    for (unsigned old_local_index=0; old_local_index<num_nodes_elem; old_local_index++)   // here num_nodes_elem shall be the old one, but how to get it?
                    {
                        oss.str("");
                        oss << "node_" << old_local_index;
                        name_item = oss.str();
                        unsigned old_node_global_index = p_cell->GetCellData()->GetItem(name_item);

                        if (node_global_index == old_node_global_index)
                        {
                            oss.str("");
                            oss << "angle_" << local_index;
                            name_item = oss.str();
                            double angle_diff = angle - p_cell->GetCellData()->GetItem(name_item);
                            if (angle_diff > max_angle_diff)
                            {
                                max_angle_diff = angle_diff;
                                max_angle_diff_node = local_index;
                                nucleation_angle = angle;
                                c_vector<double, DIM> bisector = vec_CD/norm_2(vec_CD) + vec_CE/norm_2(vec_CE);
                                bisector_orientation = asin(bisector[1]/norm_2(bisector))*180/M_PI;
                            }
                        }
                    }
                }
            }

            for (unsigned local_index=0; local_index<num_nodes_elem; local_index++)
            {
                unsigned node_global_index = p_element->GetNodeGlobalIndex(local_index);

                unsigned next_node_local_index = (local_index+1)%num_nodes_elem;
                unsigned previous_node_local_index = (num_nodes_elem+local_index-1)%num_nodes_elem;
                c_vector<double, DIM> vec_CD = p_element->GetNodeLocation(previous_node_local_index) - p_element->GetNodeLocation(local_index);
                c_vector<double, DIM> vec_CE = p_element->GetNodeLocation(next_node_local_index) - p_element->GetNodeLocation(local_index);
                double angle = acos((vec_CD[0]*vec_CE[0]+vec_CD[1]*vec_CE[1])/(norm_2(vec_CD)*norm_2(vec_CE)))*180/M_PI;

                oss.str("");
                oss << "node_" << local_index;
                name_item = oss.str();
                p_cell->GetCellData()->SetItem(name_item, node_global_index);

                oss.str("");
                oss << "angle_" << local_index;
                name_item = oss.str();
                p_cell->GetCellData()->SetItem(name_item, angle);
            }

            if ( (p_element->GetNumStressfibers()<1) && (max_angle_diff>0.0) )
            {
                std::cout<<"there are "<< p_element->GetNumStressfibers() << " stress fibers in element "<<p_element->GetIndex()<<std::endl;
                std::cout<<"creating stress fibers..."<<std::endl;

                std::vector<Node<DIM>*> stressfiber_nodes;
                unsigned num_nodes_elem = p_element->GetNumNodes();

                // unsigned node_global_index = p_element->GetNodeGlobalIndex(local_index);
                // Node<DIM>* p_this_node = p_cell_population->GetNode(node_global_index);
                // double y_coord= p_this_node->rGetLocation()[1];
                unsigned next_node_local_index_of_max_angle_diff_node = (max_angle_diff_node+1)%num_nodes_elem;
                unsigned previous_node_local_index_of_max_angle_diff_node = (num_nodes_elem+max_angle_diff_node-1)%num_nodes_elem;

                stressfiber_nodes.push_back(p_element->GetNode(previous_node_local_index_of_max_angle_diff_node));
                stressfiber_nodes.push_back(p_element->GetNode(max_angle_diff_node));
                stressfiber_nodes.push_back(p_element->GetNode(max_angle_diff_node));
                stressfiber_nodes.push_back(p_element->GetNode(next_node_local_index_of_max_angle_diff_node));

                c_vector<double,2> end_points_ratio = zero_vector<double>(2);
                end_points_ratio[0] = 0.5;
                end_points_ratio[1] = 0.5;

                VertexElement<DIM-1,DIM>* p_stressfiber = new VertexElement<DIM-1,DIM>(global_index, stressfiber_nodes, end_points_ratio);
                // stressfibers.push_back(p_stressfiber);

                p_element->AddStressfiber(p_stressfiber);
                std::cout<< p_element->GetNumStressfibers() << " stress fibers are created in element "<<p_element->GetIndex()<<std::endl;
                
                p_cell->GetCellData()->SetItem("nucleation_node", max_angle_diff_node);
                p_cell->GetCellData()->SetItem("bisector_orientation", bisector_orientation);

                p_cell->GetCellData()->SetItem("nucleation_flag",mFlag+1);

                global_index++;
            }

            if (p_element->GetNumStressfibers()>0)
            {
                //     3---B---2
                //         |   |
                //         |   |
                //     0---A---1
                VertexElement<DIM-1, DIM>* p_get_stressfiber = p_element->GetStressfiber(0);
                Node<DIM>* p_node0 = p_get_stressfiber->GetStressfiberNode(0);
                Node<DIM>* p_node1 = p_get_stressfiber->GetStressfiberNode(1);
                Node<DIM>* p_node2 = p_get_stressfiber->GetStressfiberNode(2);
                Node<DIM>* p_node3 = p_get_stressfiber->GetStressfiberNode(3);
                c_vector<double,2> ratio = p_get_stressfiber->GetStressfiberEndpointsratio();

                c_vector<double, DIM> locationA = (1-ratio[0])*p_node0->rGetLocation() + ratio[0]*p_node1->rGetLocation();
                c_vector<double, DIM> locationB = (1-ratio[1])*p_node2->rGetLocation() + ratio[1]*p_node3->rGetLocation();
                c_vector<double, DIM> vec_AB = locationB - locationA;
                double length_AB = sqrt(vec_AB[0]*vec_AB[0]+vec_AB[1]*vec_AB[1]);
                double sf_tension = this->mSfTension;
           
                for (unsigned local_index=0; local_index<num_nodes_elem; local_index++)
                {
                    unsigned node_index = p_element->GetNodeGlobalIndex(local_index);
                    c_vector<double, DIM> sf_contribution = zero_vector<double>(DIM);
                    if (length_AB>0.0)
                    {
                        double location_x = p_element->GetNodeLocation(local_index)[0];
                        double location_y = p_element->GetNodeLocation(local_index)[1];
                        if (location_x==p_node0->rGetLocation()[0] && location_y==p_node0->rGetLocation()[1])
                        {
                            sf_contribution = sf_tension*(1-ratio[0])*vec_AB/length_AB;
                        }
                        if (location_x==p_node1->rGetLocation()[0] && location_y==p_node1->rGetLocation()[1])
                        {
                            sf_contribution = sf_tension*ratio[0]*vec_AB/length_AB;
                        }
                        if (location_x==p_node2->rGetLocation()[0] && location_y==p_node2->rGetLocation()[1])
                        {
                            sf_contribution = -sf_tension*(1-ratio[1])*vec_AB/length_AB;
                        }
                        if (location_x==p_node3->rGetLocation()[0] && location_y==p_node3->rGetLocation()[1])
                        {
                            sf_contribution = -sf_tension*ratio[1]*vec_AB/length_AB;
                        }
                    }
                    p_cell_population->GetNode(node_index)->AddAppliedForceContribution(sf_contribution);
                }
            }

            if (p_cell->GetCellData()->GetItem("nucleation_flag")==1)
            {
                OutputFileHandler output_file_handler(this->mOutputDirectory+"/", false);
                std::string output_file_for_bisector_orientation;
                output_file_for_bisector_orientation = "bisectororientation.dat";
                mpBisectorOrientationFile = output_file_handler.OpenOutputFile(output_file_for_bisector_orientation, std::ios::app);

                double t_now = SimulationTime::Instance()->GetTime();
                if (this->mIfEquilibrateForAWhile && t_now>this->mTimeForEquilibrium)
                { 
                    *mpBisectorOrientationFile << t_now << " ";
                }
                *mpBisectorOrientationFile << elem_index  << " ";
                *mpBisectorOrientationFile << max_angle_diff_node  << " ";
                *mpBisectorOrientationFile << max_angle_diff  << " ";
                *mpBisectorOrientationFile << nucleation_angle  << " ";
                *mpBisectorOrientationFile << bisector_orientation << " ";
                *mpBisectorOrientationFile << "\n";

                p_cell->GetCellData()->SetItem("nucleation_flag", p_cell->GetCellData()->GetItem("nucleation_flag")+1);
            }
        }
    }
    
    // std::vector<VertexElement<1,2>*> stressfibers;

    // if (local_index>5)
    // {
    //     std::string name_item1;
    //     std::string name_item2;
    //     oss.str("");
    //     oss << "node_" << local_index;
    //     name_item1 = oss.str();
    //     oss.str("");
    //     oss << "angle_" << local_index;
    //     name_item2 = oss.str();
    //     std::cout<< elem_index << " " << local_index << " " << p_cell->GetCellData()->GetItem(name_item1) << " " << p_cell->GetCellData()->GetItem(name_item2) << std::endl; 
    // }

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

template<unsigned DIM>
void MyStressfiberTensionForce<DIM>::SetStressfiberTension(double sfTension)
{
    mSfTension = sfTension;
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
