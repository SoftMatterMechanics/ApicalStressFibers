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

#ifndef TESTMYSTRIPESADHESIONSIMULATION_HPP_
#define TESTMYSTRIPESADHESIONSIMULATION_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "SimpleTargetAreaModifier.hpp"

#include "FakePetscSetup.hpp"

//#include "MyEdgeMyosinActivityForce.hpp"
//#include "MyEdgeMyosinActivityModifier.hpp"
//#include "MyCircleBoundaryCondition.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "ToroidalHoneycombVertexMeshGenerator.hpp"
#include "NagaiHondaForce.hpp"
#include "MyNagaiHondaForceWithStripesAdhesion.hpp"
//#include "MyTimeRelatedPlaneBoundaryCondition.hpp"

class TestMyStripesAdhesionSimulation : public AbstractCellBasedTestSuite
{
public:

    void TestOscillation()
    {
        double rest_area = sqrt(3)/2/pow(1.25,2);
        //double rest_area = sqrt(3)/2;
        double initial_area = rest_area;
        unsigned num_ele_cross = 25;
        unsigned num_ele_up = 8;
        HoneycombVertexMeshGenerator generator(num_ele_cross, num_ele_up, false, 0.01, 0.001, initial_area);    // Parameters are: cells across, cells up
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        //double length = (num_ele_cross+0.5)*sqrt(initial_area*2/sqrt(3));
        double length = (num_ele_cross)*sqrt(initial_area*2/sqrt(3));

        // double height = (num_ele_up*sqrt(3)/2+1/sqrt(3)/2)*sqrt(initial_area*2/sqrt(3));

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        double stop_proliferate_time = 0.0;
        cells_generator.GenerateBasicRandomWithStopProliferateTime(cells, p_mesh->GetNumElements(), stop_proliferate_time, p_transit_type);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        //***********Attention****************************************************
        cell_population.SetRestrictVertexMovementBoolean(false);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("MyStripesAdhesionSimulationFile/K=6.75_Ga=0.119_SubAdh=-59.5_p0=3.8_Wid0.44Dis5.29_FineMesh_dt=0.05_T=600");
        simulator.SetEndTime(600.0);
        simulator.SetDt(0.05);
        simulator.SetSamplingTimestepMultiple(20);

        MAKE_PTR(MyNagaiHondaForceWithStripesAdhesion<2>, p_force);
        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(SimpleTargetAreaModifier<2>, p_growth_modifier, ());
        simulator.AddSimulationModifier(p_growth_modifier);

        //MAKE_PTR(MyEdgeMyosinActivityModifier, p_my_modifier);
        //simulator.AddSimulationModifier(p_my_modifier);

        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc);

        normal(0) = -1.0;
        normal(1) = 0.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc2);

        point(0) = length;
        normal(0) = 1.0;
        normal(1) = 0.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc3);

        // point(0) = 0.0;
        // point(1) = height;
        // normal(0) = 0.0;
        // normal(1) = 1.0;
        // double release_boundary_time = 0.0 - 1e-10;
        // MAKE_PTR_ARGS(MyTimeRelatedPlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal, release_boundary_time));
        // simulator.AddCellPopulationBoundaryCondition(p_bc4);

        simulator.Solve();
    }

};

#endif /* TESTMYSTRIPESADHESIONSIMULATION_HPP_ */
