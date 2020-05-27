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

#ifndef TESTMYXTOROIDALSTRIPESADHESIONSIMULATION_HPP_
#define TESTMYXTOROIDALSTRIPESADHESIONSIMULATION_HPP_

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
#include "MyXToroidalHoneycombVertexMeshGenerator.hpp"
#include "NagaiHondaForce.hpp"
#include "MyNagaiHondaForceWithStripesAdhesion.hpp"
//#include "MyTimeRelatedPlaneBoundaryCondition.hpp"

class TestMyXToroidalStripesAdhesionSimulation : public AbstractCellBasedTestSuite
{
public:

    void TestOscillation()
    {
        // double rest_area = sqrt(3)/2/pow(1.25,2);
        // //double rest_area = sqrt(3)/2;
        // double initial_area = rest_area;
        unsigned num_ele_cross = 12;
        unsigned num_ele_up = 8;
        double initial_area = M_PI;
        double cell_rearrangement_threshold = 0.02;
        // HoneycombVertexMeshGenerator generator(num_ele_cross, num_ele_up, false, cell_rearrangement_threshold, 0.001, initial_area);    // Parameters are: cells across, cells up
        // MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        MyXToroidalHoneycombVertexMeshGenerator generator(num_ele_cross, num_ele_up, initial_area, cell_rearrangement_threshold, 0.001);
        MyXToroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();

        //double length = (num_ele_cross+0.5)*sqrt(initial_area*2/sqrt(3));
        //double length = (num_ele_cross)*sqrt(initial_area*2/sqrt(3));

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

        // Force Parameters:
        MAKE_PTR(MyNagaiHondaForceWithStripesAdhesion<2>, p_force);
        double nagai_honda_deformation_energy_parameter = 1.0;
        double nagai_honda_membrane_surface_energy_parameter = 0.1;
        double nagai_honda_cell_cell_adhesion_energy_parameter =0.0;
        double nagai_honda_cell_boundary_adhesion_energy_parameter =0.0;
        double substrate_adhesion_parameter_top_area = -50.0;
        double substrate_adhesion_parameter = -1;
        double substrate_adhesion_parameter_top_area_magnifies = substrate_adhesion_parameter_top_area/substrate_adhesion_parameter;
        double leading_top_length = 4.0;

        double strip_width = sqrt(M_PI/(sqrt(3)/2))/2; //1.05;
        double strip_distance = 6*sqrt(M_PI/(sqrt(3)/2));
        double strip_start_x_location = -3*sqrt(M_PI/(sqrt(3)/2));
        double strip_start_y_location = 1.5*1/sqrt(3)*sqrt(initial_area/(sqrt(3)/2))*num_ele_up;
        double fixed_target_area = M_PI;
        double fixed_target_shape_index = 4.0; //sqrt(8)*pow(3.0,0.25);// regular hexagon
        bool use_fine_mesh = false;
        bool if_substrate_adhesion_parameter_change = true;

        p_force->SetNagaiHondaDeformationEnergyParameter(nagai_honda_deformation_energy_parameter);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(nagai_honda_membrane_surface_energy_parameter);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(nagai_honda_cell_cell_adhesion_energy_parameter);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(nagai_honda_cell_boundary_adhesion_energy_parameter);
        p_force->SetSubstrateAdhesionParameter(substrate_adhesion_parameter);
        p_force->SetStripWidth(strip_width);
        p_force->SetStripDistance(strip_distance);
        p_force->SetStripStartXLocation(strip_start_x_location);
        p_force->SetStripStartYLocation(strip_start_y_location);
        p_force->SetFixedTargetArea(fixed_target_area);
        p_force->SetFixedTargetShapeIndex(fixed_target_shape_index);
        p_force->SetUseFineMesh(use_fine_mesh);
        p_force->SetIfSubstrateAdhesionParameterChange(if_substrate_adhesion_parameter_change);
        //p_force->SetSubstrateAdhesionParameterChangePerUnitLength(substrate_adhesion_parameter_change_per_unit_length);
        p_force->SetSubstrateAdhesionParameterTopAreaMagnifies(substrate_adhesion_parameter_top_area_magnifies);
        p_force->SetSubstrateAdhesionLeadingTopLength(leading_top_length);
        simulator.AddForce(p_force);

        double dt = 0.1;
        double sampling_timestep_multiple = round(1/dt);
        std::string out_put_directory = "MyXToroidalStripesAdhesionSimulation/T1SwapSuccessForBoundary_NewFormOfSubAdh";
        out_put_directory += "_SubAdhLeadTopLen" + std::to_string(leading_top_length);
        out_put_directory += "_Dt=" + std::to_string(dt);
        out_put_directory += "_UseFineMesh=" + std::to_string(use_fine_mesh);
        out_put_directory += "_RearThresh=" + std::to_string(cell_rearrangement_threshold);
        out_put_directory += "_StripStartXLoc=" +std::to_string(strip_start_x_location);
        out_put_directory += "/_Ga="+std::to_string(nagai_honda_membrane_surface_energy_parameter);
        out_put_directory += "_ShapeIndex=" + std::to_string(fixed_target_shape_index);
        out_put_directory += "_SubAdhParaTopArea=" + std::to_string(substrate_adhesion_parameter_top_area);
        out_put_directory += "_SubAdhPara=" + std::to_string(substrate_adhesion_parameter);
        simulator.SetOutputDirectory(out_put_directory);
        std::cout << "OutputDirectory: " << out_put_directory << std::endl;
        // simulator.SetEndTime(800.0);
        simulator.SetEndTime(800.0);
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(sampling_timestep_multiple);




        MAKE_PTR_ARGS(SimpleTargetAreaModifier<2>, p_growth_modifier, ());
        simulator.AddSimulationModifier(p_growth_modifier);

        //MAKE_PTR(MyEdgeMyosinActivityModifier, p_my_modifier);
        //simulator.AddSimulationModifier(p_my_modifier);

        // c_vector<double,2> point = zero_vector<double>(2);
        // c_vector<double,2> normal = zero_vector<double>(2);
        // normal(1) = -1.0;
        // MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc, (&cell_population, point, normal));
        // simulator.AddCellPopulationBoundaryCondition(p_bc);

        // normal(0) = -1.0;
        // normal(1) = 0.0;
        // MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal));
        // simulator.AddCellPopulationBoundaryCondition(p_bc2);

        // point(0) = length;
        // normal(0) = 1.0;
        // normal(1) = 0.0;
        // MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal));
        // simulator.AddCellPopulationBoundaryCondition(p_bc3);

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

#endif /* TESTMYXTOROIDALSTRIPESADHESIONSIMULATION_HPP_ */
