/*

Copyright (c) 2005-2020, University of Oxford.
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

#include "ForwardEulerNumericalMethod.hpp"
#include "MutableVertexMesh.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ForwardEulerNumericalMethod<ELEMENT_DIM,SPACE_DIM>::ForwardEulerNumericalMethod()
    : AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ForwardEulerNumericalMethod<ELEMENT_DIM,SPACE_DIM>::~ForwardEulerNumericalMethod()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ForwardEulerNumericalMethod<ELEMENT_DIM,SPACE_DIM>::UpdateAllNodePositions(double dt)
{
    if (!this->mUseUpdateNodeLocation)
    {
        // Apply forces to each cell, and save a vector of net forces F
        std::vector<c_vector<double, SPACE_DIM> > forces = this->ComputeForcesIncludingDamping();

        unsigned index = 0;
        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
             node_iter != this->mpCellPopulation->rGetMesh().GetNodeIteratorEnd();
             ++node_iter, ++index)
        {
            // Get the current node location and calculate the new location according to the forward Euler method
            const c_vector<double, SPACE_DIM>& r_old_location = node_iter->rGetLocation();
            c_vector<double, SPACE_DIM> displacement = dt * forces[index];

            // In the vertex-based case, the displacement may be scaled if the cell rearrangement threshold is exceeded
            this->DetectStepSizeExceptions(node_iter->GetIndex(), displacement, dt);

            c_vector<double, SPACE_DIM> new_location = r_old_location + displacement;
            this->SafeNodePositionUpdate(node_iter->GetIndex(), new_location);
        }
    }
    else
    {
        /*
         * If this type of cell population does not support the new numerical methods, delegate
         * updating node positions to the population itself.
         *
         * This only applies to NodeBasedCellPopulationWithBuskeUpdates.
         */
        this->mpCellPopulation->UpdateNodeLocations(dt);
    }
}

// my changes
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double ForwardEulerNumericalMethod<ELEMENT_DIM,SPACE_DIM>::GetNewAdaptiveTimestepAndUpdateAllNodePositions(double dt)
{
    double new_adaptive_timestep = dt;
    double maximum_velocity = 0.0;
    if (!this->mUseUpdateNodeLocation)
    {
        // Apply forces to each cell, and save a vector of net forces F
        std::vector<c_vector<double, SPACE_DIM> > forces = this->ComputeForcesIncludingDamping();

        unsigned i = 0;
        unsigned node_global_index = 0;
        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
                node_iter != this->mpCellPopulation->rGetMesh().GetNodeIteratorEnd();
                ++node_iter, ++i)
        {
            double velocity_abs = norm_2(forces[i]);
            if (velocity_abs > maximum_velocity)
                node_global_index = node_iter->GetIndex();
            maximum_velocity = std::max(maximum_velocity, velocity_abs);
        }
        // new_adaptive_timestep = static_cast<MutableVertexMesh<ELEMENT_DIM, SPACE_DIM>*>(& this->mpCellPopulation->rGetMesh())->GetCellRearrangementThreshold()/2.0/maximum_velocity;
        new_adaptive_timestep = this->mMaxMovementPerTimestep/maximum_velocity;
        new_adaptive_timestep = std::min(dt, new_adaptive_timestep);
        Node<SPACE_DIM>* p_Node = this->mpCellPopulation->rGetMesh().GetNode(node_global_index);
        if (mOutputNumericalMethodInformation)
        {
            std::cout << "UpdateNodes: timesteps elapsed: " << SimulationTime::Instance()->GetTimeStepsElapsed() << std::endl;
            std::cout << "node_global_index=" << node_global_index << " node_location=" << p_Node->rGetLocation()[0] << ", " << p_Node->rGetLocation()[1]
                    << " maximum_velocity=" << maximum_velocity << " new_adaptive_timestep=" << new_adaptive_timestep << std::endl;
            std::cout << std::endl;
        }

        // unsigned index = 0;
        // for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
        //      node_iter != this->mpCellPopulation->rGetMesh().GetNodeIteratorEnd();
        //      ++node_iter, ++index)
        // {
        //     // Get the current node location and calculate the new location according to the forward Euler method
        //     const c_vector<double, SPACE_DIM>& r_old_location = node_iter->rGetLocation();
        //     c_vector<double, SPACE_DIM> displacement = new_adaptive_timestep * forces[index];

        //     c_vector<double, SPACE_DIM> new_location = r_old_location + displacement;
        //     this->SafeNodePositionUpdate(node_iter->GetIndex(), new_location);
        // }

        // my changes
        double t_now = SimulationTime::Instance()->GetTime();
        if ( !(mIfEquilibrateForAWhile && t_now<=mTimeForEquilibrium) )
        {
            // find the square that differentiate the top, bottom, left and right boundary points
            double x_far_left = 0.0;
            double x_far_right = 0.0;
            double box_left = 0.0;
            double box_right = 0.0;
            double box_bottom = this->mCenterYCoordination;
            double box_top = this->mCenterYCoordination;
            
            unsigned node_index = 0;
            for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
                node_iter != this->mpCellPopulation->rGetMesh().GetNodeIteratorEnd();
                ++node_iter, ++node_index)
            {
                if (node_iter->IsBoundaryNode())
                {
                    double x_coord= node_iter->rGetLocation()[0];
                    double y_coord = node_iter->rGetLocation()[1];
                    if (x_coord < x_far_left)
                    {
                        x_far_left = x_coord;
                        box_left = x_far_left;
                    }
                    if (x_coord > x_far_right)
                    {
                        x_far_right = x_coord;
                        box_right = x_far_right;
                    }
                    if (y_coord > box_top)
                    {
                        box_top = y_coord;
                    }
                    if (y_coord < box_bottom)
                    {
                        box_bottom = y_coord;
                    }
                }
            }
            double new_box_left = 3.0/4*x_far_left + 1.0/4*x_far_right;
            double new_box_right = 1.0/4*x_far_left + 3.0/4*x_far_right;

            if (!(new_box_left==box_left && new_box_right==box_right))
            {
                box_left = new_box_left;
                box_right = new_box_right;

                for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
                    node_iter != this->mpCellPopulation->rGetMesh().GetNodeIteratorEnd();
                    ++node_iter)
                {
                    if (node_iter->IsBoundaryNode())
                    {
                        double x_coord = node_iter->rGetLocation()[0];
                        double y_coord = node_iter->rGetLocation()[1];
                        if ( (x_coord > box_left) && (x_coord < box_right) )
                        {
                            if ( (y_coord < this->mCenterYCoordination) && (y_coord > box_bottom) )
                            {
                                box_bottom = y_coord;
                            }
                            if ( (y_coord > this->mCenterYCoordination) && (y_coord < box_top) )
                            {
                                box_top = y_coord;
                            }
                        }
                    }
                }

                box_left = x_far_left;
                box_right = x_far_right;
                for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
                    node_iter != this->mpCellPopulation->rGetMesh().GetNodeIteratorEnd();
                    ++node_iter)
                {
                    if (node_iter->IsBoundaryNode())
                    {
                        double x_coord = node_iter->rGetLocation()[0];
                        double y_coord = node_iter->rGetLocation()[1];
                        if ( (y_coord > box_bottom) && (y_coord < box_top) )
                        {
                            if ( (x_coord < 0) && (x_coord > box_left) )
                            {
                                box_left = x_coord;
                            }
                            if ( (x_coord > 0) && (x_coord < box_right) )
                            {
                                box_right = x_coord;
                            }
                        }
                    }
                }

                new_box_left = box_left;
                new_box_right = box_right;
            }

            // std::cout << "left: " << box_left << std::endl;
            // std::cout << "right: " << box_right << std::endl;
            // std::cout << "bottom: " << box_bottom << std::endl;
            // std::cout << "top: " << box_top << std::endl;

            // count the nodes of top and bottom boundaries
            unsigned num_top_boundary_nodes = 0;
            unsigned num_bottom_boundary_nodes = 0;
            unsigned num_left_boundary_nodes = 0;
            unsigned num_right_boundary_nodes = 0;
            double left_damping_force = 0.0;
            double right_damping_force = 0.0;
            double top_damping_force = 0.0;
            double bottom_damping_force = 0.0;
            unsigned index = 0;
            for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
                node_iter != this->mpCellPopulation->rGetMesh().GetNodeIteratorEnd();
                ++node_iter, ++index)
            {
                if (node_iter->IsBoundaryNode())
                {
                    double x_coord= node_iter->rGetLocation()[0];
                    double y_coord= node_iter->rGetLocation()[1];

                    // // vertical stretching
                    // if ((x_coord <= box_left) && (y_coord < box_top) && (y_coord >box_bottom))
                    // {
                    //     num_left_boundary_nodes += 1;
                    //     left_damping_force += forces[index][0]; 
                    // }
                    // if ((x_coord >= box_right) && (y_coord < box_top) && (y_coord >box_bottom))
                    // {
                    //     num_right_boundary_nodes += 1;
                    //     right_damping_force += forces[index][0];
                    // }
                    // if (y_coord >= box_top)
                    // {
                    //     num_top_boundary_nodes += 1;
                    //     top_damping_force += forces[index][1];
                    // }
                    // if (y_coord <= box_bottom)
                    // {
                    //     num_bottom_boundary_nodes += 1;
                    //     bottom_damping_force += forces[index][1];
                    // }

                    // horizontal stretching
                    if (x_coord <= box_left)
                    {
                        num_left_boundary_nodes += 1;
                        left_damping_force += forces[index][0]; 
                    }
                    if (x_coord >= box_right)
                    {
                        num_right_boundary_nodes += 1;
                        right_damping_force += forces[index][0];
                    }
                    // if ((y_coord >= box_top) && (x_coord > box_left) && (x_coord < box_right))
                    if (y_coord >= box_top)
                    {
                        num_top_boundary_nodes += 1;
                        top_damping_force += forces[index][1];
                    }
                    // if ((y_coord <= box_bottom) && (x_coord > box_left) && (x_coord < box_right))
                    if (y_coord <= box_bottom)
                    {
                        num_bottom_boundary_nodes += 1;
                        bottom_damping_force += forces[index][1];
                    }
                }
            }

            // std::cout << "left nodes: " << num_left_boundary_nodes << std::endl;
            // std::cout << "right nodes: " << num_right_boundary_nodes << std::endl;
            // std::cout << "bottom nodes: " << num_bottom_boundary_nodes << std::endl;
            // std::cout << "top nodes: " << num_top_boundary_nodes << std::endl;

            // unsigned nodes_num= index;
            // std::vector<c_vector<double, SPACE_DIM> > old_location(nodes_num, zero_vector<double>(SPACE_DIM));
            // std::vector<c_vector<double, SPACE_DIM> > new_location(nodes_num, zero_vector<double>(SPACE_DIM));
            index = 0;
            for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
                node_iter != this->mpCellPopulation->rGetMesh().GetNodeIteratorEnd();
                ++node_iter, ++index)
            {
                // Get the current node location and calculate the new location according to the forward Euler method
                // unsigned node_global_index = node_iter->GetIndex();
                const c_vector<double, SPACE_DIM>& r_old_location = node_iter->rGetLocation(); 
                // old_location[node_global_index] = r_old_location;

                node_iter->SetOldLocation(r_old_location);
                // std::cout << "old location(" << node_global_index << ")=" << old_location[node_global_index][0] << ", " << old_location[node_global_index][1]<<std::endl; 
                // std::cout << "old location(" << node_global_index << ")=" << node_iter->GetOldLocation()[0] << ", " << node_iter->GetOldLocation()[1]<<std::endl; 

                c_vector<double, SPACE_DIM> displacement = zero_vector<double>(SPACE_DIM);
                double x_coord= node_iter->rGetLocation()[0];
                double y_coord= node_iter->rGetLocation()[1];
                
                if (this->mIfEquilibrateForAWhile && t_now>this->mTimeForEquilibrium && node_iter->IsBoundaryNode() )
                {
                    // // vertical stretching
                    // double HorizontalMorphogeneticForce1 = this->mHorizontalMorphogeneticForceGrowthRate*(t_now - this->mTimeForEquilibrium);
                    // double HorizontalMorphogeneticForce2 = this->mHorizontalMorphogeneticForceGrowthRate*(this->mRealEquilibriumTime - this->mTimeForEquilibrium);
                    // double VerticalMorphogeneticForce1 = this->mHorizontalMorphogeneticForceGrowthRate*(t_now - this->mTimeForEquilibrium);
                    // double VerticalMorphogeneticForce2 = this->mHorizontalMorphogeneticForceGrowthRate*(this->mRealEquilibriumTime - this->mTimeForEquilibrium)
                    //                                     + this->mVerticalMorphogeneticForceGrowthRate*(t_now - this->mRealEquilibriumTime);
                    
                    // double HorizontalMorphogeneticForce = (t_now < this->mRealEquilibriumTime? HorizontalMorphogeneticForce1:HorizontalMorphogeneticForce2);
                    // double VerticalMorphogeneticForce = (t_now < this->mRealEquilibriumTime? VerticalMorphogeneticForce1:VerticalMorphogeneticForce2);
                    
                    // if ((x_coord <= box_left) && (y_coord < box_top) && (y_coord > box_bottom))
                    // {
                    //     displacement[0] = new_adaptive_timestep * (-HorizontalMorphogeneticForce + left_damping_force)/num_left_boundary_nodes;
                    //     displacement[1] = new_adaptive_timestep * forces[index][1]; 
                    // }
                    // if ((x_coord >= box_right) && (y_coord < box_top) && (y_coord > box_bottom))
                    // {
                    //     displacement[0] = new_adaptive_timestep * (HorizontalMorphogeneticForce + right_damping_force)/num_right_boundary_nodes;
                    //     displacement[1] = new_adaptive_timestep * forces[index][1]; 
                    // }
                    // if (y_coord >= box_top)
                    // {
                    //     displacement[0] = new_adaptive_timestep * forces[index][0];
                    //     displacement[1] = new_adaptive_timestep * (VerticalMorphogeneticForce + top_damping_force)/num_top_boundary_nodes;
                    // }
                    // if (y_coord <= box_bottom) 
                    // {
                    //     displacement[0] = new_adaptive_timestep * forces[index][0];
                    //     displacement[1] = new_adaptive_timestep * (-VerticalMorphogeneticForce + bottom_damping_force)/num_bottom_boundary_nodes;
                    // }

                    // horizontal stretching
                    double HorizontalMorphogeneticForce1 = this->mHorizontalMorphogeneticForceGrowthRate*(t_now - this->mTimeForEquilibrium);
                    double HorizontalMorphogeneticForce2 = this->mHorizontalMorphogeneticForceGrowthRate*(this->mRealEquilibriumTime - this->mTimeForEquilibrium)
                                                            + this->mHorizontalMorphogeneticForceGrowthRate*(t_now - this->mRealEquilibriumTime);
                    // double HorizontalMorphogeneticForce3 = this->mHorizontalMorphogeneticForceGrowthRate*(this->mRealEquilibriumTime - this->mTimeForEquilibrium)
                    //                                         + this->mHorizontalMorphogeneticForceGrowthRate*(this->mLoadingGrowthStopTime - this->mRealEquilibriumTime);
                                                            
                    double VerticalMorphogeneticForce1 = this->mHorizontalMorphogeneticForceGrowthRate*(t_now - this->mTimeForEquilibrium);
                    double VerticalMorphogeneticForce2 = this->mHorizontalMorphogeneticForceGrowthRate*(this->mRealEquilibriumTime - this->mTimeForEquilibrium);
                    
                    // double HorizontalMorphogeneticForce = (t_now < this->mRealEquilibriumTime? HorizontalMorphogeneticForce1:(t_now < this->mLoadingGrowthStopTime? HorizontalMorphogeneticForce2:HorizontalMorphogeneticForce3));
                    double HorizontalMorphogeneticForce = (t_now < this->mRealEquilibriumTime? HorizontalMorphogeneticForce1:HorizontalMorphogeneticForce2);
                    double VerticalMorphogeneticForce = (t_now < this->mRealEquilibriumTime? VerticalMorphogeneticForce1:VerticalMorphogeneticForce2);
                    
                    if (x_coord <= box_left)
                    {
                        displacement[0] = new_adaptive_timestep * (-HorizontalMorphogeneticForce + left_damping_force)/num_left_boundary_nodes;
                        displacement[1] = new_adaptive_timestep * forces[index][1]; 
                    }
                    if (x_coord >= box_right)
                    {
                        displacement[0] = new_adaptive_timestep * (HorizontalMorphogeneticForce + right_damping_force)/num_right_boundary_nodes;
                        displacement[1] = new_adaptive_timestep * forces[index][1]; 
                    }
                    // if ((y_coord >= box_top) && (x_coord > box_left) && (x_coord < box_right))
                    if (y_coord >= box_top)
                    {
                        displacement[0] = new_adaptive_timestep * forces[index][0];
                        displacement[1] = new_adaptive_timestep * (VerticalMorphogeneticForce + top_damping_force)/num_top_boundary_nodes;
                    }
                    // if ((y_coord <= box_bottom) && (x_coord > box_left) && (x_coord < box_right))
                    if (y_coord <= box_bottom)
                    {
                        displacement[0] = new_adaptive_timestep * forces[index][0];
                        displacement[1] = new_adaptive_timestep * (-VerticalMorphogeneticForce + bottom_damping_force)/num_bottom_boundary_nodes;
                    }
                }
                else
                {
                    displacement = new_adaptive_timestep * forces[index];
                }

                c_vector<double, SPACE_DIM> r_new_location = r_old_location + displacement;
                this->SafeNodePositionUpdate(node_iter->GetIndex(), r_new_location);

                // new_location[node_global_index] = r_new_location;
            }
        }
        else
        {
            unsigned index = 0;
            for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
                node_iter != this->mpCellPopulation->rGetMesh().GetNodeIteratorEnd();
                ++node_iter, ++index)
            {
                // Get the current node location and calculate the new location according to the forward Euler method
                const c_vector<double, SPACE_DIM>& r_old_location = node_iter->rGetLocation();
                c_vector<double, SPACE_DIM> displacement = zero_vector<double>(SPACE_DIM);

                displacement = new_adaptive_timestep * forces[index];

                c_vector<double, SPACE_DIM> new_location = r_old_location + displacement;
                this->SafeNodePositionUpdate(node_iter->GetIndex(), new_location);
            }
        }

        // OutputFileHandler output_file_handler(this->mOutputDirectory+"/", false);
        // std::string output_file_for_reaction_forces;
        // output_file_for_reaction_forces = "reactionforces.dat";
        // mpReactionForcesFile = output_file_handler.OpenOutputFile(output_file_for_reaction_forces, std::ios::app);

        // double t_now = SimulationTime::Instance()->GetTime();
        // if (this->mIfEquilibrateForAWhile && t_now>this->mTimeForEquilibrium)
        // { 
        //     // *mpReactionForcesFile << t_now << "\n";
        //     *mpReactionForcesFile << t_now << " ";
        // }

        return new_adaptive_timestep;
    }
    else
    {
        /*
         * If this type of cell population does not support the new numerical methods, delegate
         * updating node positions to the population itself.
         *
         * This only applies to NodeBasedCellPopulationWithBuskeUpdates.
         */
        this->mpCellPopulation->UpdateNodeLocations(new_adaptive_timestep);

        return new_adaptive_timestep;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ForwardEulerNumericalMethod<ELEMENT_DIM, SPACE_DIM>::OutputNumericalMethodParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractNumericalMethod<ELEMENT_DIM,SPACE_DIM>::OutputNumericalMethodParameters(rParamsFile);
}

// Explicit instantiation
template class ForwardEulerNumericalMethod<1,1>;
template class ForwardEulerNumericalMethod<1,2>;
template class ForwardEulerNumericalMethod<2,2>;
template class ForwardEulerNumericalMethod<1,3>;
template class ForwardEulerNumericalMethod<2,3>;
template class ForwardEulerNumericalMethod<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ForwardEulerNumericalMethod)
