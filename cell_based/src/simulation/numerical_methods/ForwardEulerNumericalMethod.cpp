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
        OutputFileHandler output_file_handler(this->mOutputDirectory+"/", false);
        std::string output_file_for_reaction_forces;
        output_file_for_reaction_forces = "reactionforces.dat";
        mpReactionForcesFile = output_file_handler.OpenOutputFile(output_file_for_reaction_forces, std::ios::app);

        double t_now = SimulationTime::Instance()->GetTime();
        if (this->mIfEquilibrateForAWhile && t_now>this->mTimeForEquilibrium)
        { 
            // *mpReactionForcesFile << t_now << "\n";
            *mpReactionForcesFile << t_now << " ";
        }

        unsigned index = 0;
        c_vector<double, SPACE_DIM> up_sum_reaction_force = zero_vector<double>(SPACE_DIM);
        c_vector<double, SPACE_DIM> bottom_sum_reaction_force = zero_vector<double>(SPACE_DIM);
        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
             node_iter != this->mpCellPopulation->rGetMesh().GetNodeIteratorEnd();
             ++node_iter, ++index)
        {
            // Get the current node location and calculate the new location according to the forward Euler method
            const c_vector<double, SPACE_DIM>& r_old_location = node_iter->rGetLocation();
            c_vector<double, SPACE_DIM> displacement = zero_vector<double>(SPACE_DIM);
            
            if (this->mIfEquilibrateForAWhile && t_now>this->mTimeForEquilibrium && node_iter->IsBoundaryNode())
            {   
                double y_coord= node_iter->rGetLocation()[1];
                double vertical_velocity_of_boundary_nodes = this->mBoundaryVelocity;
                c_vector<double, SPACE_DIM> reaction_force = zero_vector<double>(SPACE_DIM);

                if (y_coord > this->mCenterYCoordination)
                {
                    // displacement[0] = new_adaptive_timestep * forces[index][0];
                    displacement[0] = 0.0;
                    displacement[1] = new_adaptive_timestep * vertical_velocity_of_boundary_nodes;
                    reaction_force[0] = -forces[index][0];
                    reaction_force[1] = vertical_velocity_of_boundary_nodes-forces[index][1];
                    up_sum_reaction_force += reaction_force; 
                    // *mpReactionForcesFile << "up" << " ";
                }
                else
                {
                    // displacement[0] = new_adaptive_timestep * forces[index][0];
                    displacement[0] = 0.0;
                    displacement[1] = -new_adaptive_timestep * vertical_velocity_of_boundary_nodes;
                    reaction_force[0] = -forces[index][0];
                    reaction_force[1] = -vertical_velocity_of_boundary_nodes-forces[index][1];
                    bottom_sum_reaction_force += reaction_force;
                    // *mpReactionForcesFile << "bottom" << " ";
                }

                // *mpReactionForcesFile << index  << " ";
                // *mpReactionForcesFile << reaction_force[0] << " " << reaction_force[1] << " ";
                // *mpReactionForcesFile << "\n";
            }
            else
            {
                displacement = new_adaptive_timestep * forces[index];
            }

            c_vector<double, SPACE_DIM> new_location = r_old_location + displacement;
            this->SafeNodePositionUpdate(node_iter->GetIndex(), new_location);
        }

        if (this->mIfEquilibrateForAWhile && t_now>this->mTimeForEquilibrium)
        { 
            *mpReactionForcesFile << up_sum_reaction_force[0] << " " << up_sum_reaction_force[1] << " "<< bottom_sum_reaction_force[0] << " " << bottom_sum_reaction_force[1] << " ";
            *mpReactionForcesFile << "\n";
        }

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
