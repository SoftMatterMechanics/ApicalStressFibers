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

#ifndef TESTMYFLUIDSTATE_HPP_
#define TESTMYFLUIDSTATE_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "FaceValueAndStressStateModifier.hpp"
#include "PolarityModifier.hpp"

#include "FakePetscSetup.hpp"

#include "PlaneBoundaryCondition.hpp"
#include "ToroidalHoneycombVertexMeshGenerator.hpp"
#include "MyXToroidalHoneycombVertexMeshGenerator.hpp"
#include "NagaiHondaForce.hpp"
#include "MyNagaiHondaForceWithStripesAdhesion.hpp"
#include "DiffusionForce.hpp"

class TestMyFluidState : public AbstractCellBasedTestSuite
{
public:

    void TestOscillation()
    {
        /*------------------------------START: Mesh Structure------------------------------*/
        bool use_two_periodic = false;
        bool if_use_compression = false;
        bool use_face_information = false;

        unsigned num_ele_cross = 10;
        if (use_two_periodic)
        {
          num_ele_cross *= 2;
        }
        unsigned num_ele_up = 10;
        if(if_use_compression)
        {
          num_ele_up*= 2;
        }
        double initial_area = M_PI;
        double cell_rearrangement_threshold = 0.05;
        ToroidalHoneycombVertexMeshGenerator generator(num_ele_cross, num_ele_up, cell_rearrangement_threshold, 0.001, initial_area);
        Toroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();
        p_mesh->SetUseFaceInformation(use_face_information);//
        /*------------------------------END: Mesh Structure------------------------------*/


        /*------------------------------START: Cells and Cell population--------------------------*/
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        double stop_proliferate_time = 0.0;
        cells_generator.GenerateBasicRandomWithStopProliferateTime(cells, p_mesh->GetNumElements(), stop_proliferate_time, p_transit_type);
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        bool restrict_vertex_movement = false;
        cell_population.SetRestrictVertexMovementBoolean(restrict_vertex_movement);//Attention!!!!!!!!!!!!
        /*------------------------------END: Cells and Cell population--------------------------*/


        /*--------------------------------START: TargetAreaModifier------------------------------*/
        OffLatticeSimulation<2> simulator(cell_population);
        MAKE_PTR_ARGS(SimpleTargetAreaModifier<2>, p_growth_modifier, ());
        p_growth_modifier->SetReferenceTargetArea(M_PI);
        p_growth_modifier->SetGrowthDuration(10.0);//
        simulator.AddSimulationModifier(p_growth_modifier);
        /*--------------------------------END: TargetAreaModifier------------------------------*/


        /*---------------------------------START: Add Numerical Method-----------------------------*/
        boost::shared_ptr<ForwardEulerNumericalMethod<2,2>> p_numerical_method(new ForwardEulerNumericalMethod<2,2>());
        p_numerical_method->SetCellPopulation(&cell_population);
        std::vector<boost::shared_ptr<AbstractForce<2, 2> > > force_collection = simulator.rGetForceCollection();
        //Note: Here the address of std::vector of force collections are different from that in simulator! 
        //But the corresponding shared_ptrs in the vectors should be the same.
        p_numerical_method->SetForceCollection(&force_collection);
        bool use_adaptive_timestep = true;
        p_numerical_method->SetUseAdaptiveTimestep(use_adaptive_timestep);//
        simulator.SetNumericalMethod(p_numerical_method);
        /*---------------------------------END: Add Numerical Method-----------------------------*/


        /*-----------------------------START: MyNagaiHondaForceWithStripesAdhesion---------------------*/
        MAKE_PTR(MyNagaiHondaForceWithStripesAdhesion<2>, p_force);
        // KA, A0
        double nagai_honda_deformation_energy_parameter = 1.0;
        bool use_fixed_target_area = true;
        double fixed_target_area = M_PI;
        double compress_ratio = 1.2;
        if (if_use_compression)
        {
          fixed_target_area *= compress_ratio;
        }

        // KP and LAMBDA
        double nagai_honda_membrane_surface_energy_parameter = 0.1;
        
        int case_number_of_membrane_surface_energy_form = 0;//2 for new membr_surf_energy form1, 3 for form2.

        double fixed_target_shape_index = 4.0;// {6/sqrt(6*sqrt(3)/4)}=3.72
        
        // Equation form: 0.5*Gamma*P^2, use Lambda instead of target perimeter
        double nagai_honda_cell_cell_adhesion_energy_parameter = -nagai_honda_membrane_surface_energy_parameter*(fixed_target_shape_index*sqrt(fixed_target_area));
        double nagai_honda_cell_boundary_adhesion_energy_parameter = nagai_honda_cell_cell_adhesion_energy_parameter;

        // Strips structure
        double strip_width = sqrt(M_PI/(sqrt(3)/2))/2; //1.05
        double strip_distance = 6*sqrt(M_PI/(sqrt(3)/2));//12.60
        double strip_start_x_location = 0.0;
        if (use_two_periodic)
          strip_start_x_location = -3*sqrt(M_PI/(sqrt(3)/2));
        double strip_start_y_location = 1e10+1.5*num_ele_up*1/sqrt(3)*sqrt(initial_area/(sqrt(3)/2));
        if (if_use_compression)
        {
          strip_start_y_location += 0.5*1/sqrt(3)*sqrt(initial_area/(sqrt(3)/2));
          strip_start_y_location += 0.5*(fixed_target_area/initial_area-1)*strip_start_y_location;
        }
        bool use_fine_mesh = false;

        // Substrate Adhesion Parameter
        bool if_substrate_adhesion_parameter_change = true;
        double leading_top_length = 4.0;
        double substrate_adhesion_parameter_top_area = -50.0;
        double substrate_adhesion_parameter = -1;
        double substrate_adhesion_parameter_top_area_magnifies = substrate_adhesion_parameter_top_area/substrate_adhesion_parameter;

        p_force->SetNagaiHondaDeformationEnergyParameter(nagai_honda_deformation_energy_parameter);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(nagai_honda_membrane_surface_energy_parameter);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(nagai_honda_cell_cell_adhesion_energy_parameter);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(nagai_honda_cell_boundary_adhesion_energy_parameter);
        p_force->SetUseFaceToGetAdhesionParameter(use_face_information);//***********
        p_force->SetSubstrateAdhesionParameter(substrate_adhesion_parameter);
        p_force->SetStripWidth(strip_width);
        p_force->SetStripDistance(strip_distance);
        p_force->SetStripStartXLocation(strip_start_x_location);
        p_force->SetStripStartYLocation(strip_start_y_location);
        p_force->SetUseFixedTargetArea(use_fixed_target_area);
        p_force->SetFixedTargetArea(fixed_target_area);
        p_force->SetFixedTargetShapeIndex(fixed_target_shape_index);
        p_force->SetUseFineMesh(use_fine_mesh);
        p_force->SetIfSubstrateAdhesionParameterChange(if_substrate_adhesion_parameter_change);
        //p_force->SetSubstrateAdhesionParameterChangePerUnitLength(substrate_adhesion_parameter_change_per_unit_length);
        p_force->SetSubstrateAdhesionParameterTopAreaMagnifies(substrate_adhesion_parameter_top_area_magnifies);
        p_force->SetSubstrateAdhesionLeadingTopLength(leading_top_length);
        p_force->SetCaseNumberOfMembraneSurfaceEnergyForm(case_number_of_membrane_surface_energy_form);
        p_force->SetOutputInformationForTestingCode(false);
        simulator.AddForce(p_force);
        /*-----------------------------END: MyNagaiHondaForceWithStripesAdhesion---------------------*/


        /*------------------------START: Feedback: FaceValueAndStressStateModifier: need modification---------------*/
        p_mesh->SetOutputInformationForTestingCode(false);
        MAKE_PTR_ARGS(FaceValueAndStressStateModifier<2>, p_face_value_and_stress_state_modifier, ());
        
        // Energy function parameter
        p_face_value_and_stress_state_modifier->SetNagaiHondaMembraneSurfaceEnergyParameter(nagai_honda_membrane_surface_energy_parameter);
        p_face_value_and_stress_state_modifier->SetFixedTargetPerimeter(fixed_target_shape_index*sqrt(fixed_target_area));
        
        // Feedback information
        double feedback_strength = 1e-2;
        double hill_coefficient = 2;
        double edge_length_at_rest = sqrt(M_PI/(6*sqrt(3)/4));
        p_face_value_and_stress_state_modifier->SetFeedbackStrength(feedback_strength);
        p_face_value_and_stress_state_modifier->SetHillCoefficient(hill_coefficient);
        p_face_value_and_stress_state_modifier->SetEdgeLengthAtRest(edge_length_at_rest);

        double feedback_strength_for_adhesion = 5*feedback_strength;
        double hill_coefficient_for_adhesion = hill_coefficient;
        p_face_value_and_stress_state_modifier->SetFeedbackStrengthForAdhesion(feedback_strength_for_adhesion);
        p_face_value_and_stress_state_modifier->SetHillCoefficientForAdhesion(hill_coefficient_for_adhesion);

        bool EMA_dont_decrease = true;// false for default
        bool CCA_dont_decrease = true;// true for default
        bool CCA_dont_increase = false;// true for default
        bool CCA_dont_inrease_until_shorter_than_a_threshold = true;
        double CCA_dont_increase_until_shorter_than_this_value = 0.1;
        p_face_value_and_stress_state_modifier->SetEMADontDecreaseWhenEdgeShrink_CCADontDecreaseWhenEdgeExpand_CCADontInreaseWhenEdgeShrink(EMA_dont_decrease,CCA_dont_decrease,CCA_dont_increase);
        p_face_value_and_stress_state_modifier->SetCCADontInreaseUntilShorterThanAThreshold(CCA_dont_inrease_until_shorter_than_a_threshold);
        p_face_value_and_stress_state_modifier->SetCCADontInreaseUntilShorterThanThisValue(CCA_dont_increase_until_shorter_than_this_value);
        p_face_value_and_stress_state_modifier->SetOutputInformationForTestingCode(false);
        if (use_face_information)
          simulator.AddSimulationModifier(p_face_value_and_stress_state_modifier);
        /*-----------------------END: Feedback: FaceValueAndStressStateModifier: need modification---------------*/


        /*------------------------------------START: Timestep---------------------------------------*/
        double dt = 0.1;
        double sampling_timestep_multiple = round(1/dt);
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(sampling_timestep_multiple);
        simulator.SetEndTime(800.0);
        /*------------------------------------END: Timestep---------------------------------------*/


        /*-----------------------------START: RandomForce and PolarityModifier--------------------------*/
        bool add_random_force = true;
        bool consider_polarity = true;//
        // Brownian diffusion
        double translational_diffusion_constant = 0.0;
        if (add_random_force)
        {
          MAKE_PTR(DiffusionForce<2>, p_force2);
          bool use_the_same_node_radius = true;
          double node_radius = 50*1e2;
          p_force2->SetConsiderPolarity(consider_polarity);
          p_force2->SetUseTheSameNodeRadius(use_the_same_node_radius);
          p_force2->SetTheSameNodeRadius(node_radius);
          translational_diffusion_constant = p_force2->GetDiffusionScalingConstant()/node_radius;
          simulator.AddForce(p_force2);
        }

        // Polarity force
        double rotational_diffusion_constant = 0.0;
        double polarity_magnitude = 0.0;
        if (consider_polarity)
        {
          MAKE_PTR_ARGS(PolarityModifier<2>, p_polarity_modifier, ());
          rotational_diffusion_constant = 0.01;
          double Dt = dt;
          polarity_magnitude = 0.1;
          double angle_for_initialization = M_PI;// 0 for all upwards; PI for all directions
          p_polarity_modifier->SetD(rotational_diffusion_constant);
          p_polarity_modifier->SetDt(Dt);
          p_polarity_modifier->SetPolarityMagnitude(polarity_magnitude);
          p_polarity_modifier->SetAngleForInitialization(angle_for_initialization);
          simulator.AddSimulationModifier(p_polarity_modifier);
        }
        /*-----------------------------END: RandomForce and PolarityModifier---------------------------*/


        /*--------------------------------START: Boundary condition-----------------------------*/
        if (if_use_compression)
        {
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc);

        c_vector<double,2> point2 = zero_vector<double>(2);
        c_vector<double,2> normal2 = zero_vector<double>(2);
        point2(1) = (1.5*num_ele_up+0.5)*1/sqrt(3)*sqrt(initial_area/(sqrt(3)/2));
        normal2(1) = 1.0;
        double stop_time = 5.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point2, normal2, stop_time));
        simulator.AddCellPopulationBoundaryCondition(p_bc2);
        }
        /*--------------------------------END: Boundary condition-----------------------------*/
        

        /*-------------------------------START: Output Directory---------------------------------*/
        std::ostringstream oss;
        std::string out_put_directory = "MyXToroidalStripesAdhesionSimulation/";

        //out_put_directory += "TestCodeAdhesionFeedback_ResetEdgeAfterT1Swap_19";
        //out_put_directory += "TestCompression2";
        out_put_directory += "TestFluidState4";

        if(if_use_compression)
        {
          oss << std::fixed << setprecision(2) << compress_ratio;
          out_put_directory += "_ComprRatio=" + oss.str();
        }

        if (add_random_force)
        {
          oss.str("");
          oss << std::scientific << setprecision(2) << translational_diffusion_constant;
          out_put_directory += "_TransDiffuConst=" + oss.str();
          if (consider_polarity)
          {
            out_put_directory += "_RotaDiffuConst=" + std::to_string(rotational_diffusion_constant);
            out_put_directory += "_PolarMagni=" + std::to_string(polarity_magnitude);
          }
        }

        oss.str("");
        oss << std::fixed << setprecision(3) << dt;
        out_put_directory += "/_Dt=" + oss.str();
        out_put_directory += "_CaseNumMSE=" + std::to_string(case_number_of_membrane_surface_energy_form);
        // oss.str("");
        // oss << std::scientific << setprecision(0) << feedback_strength;
        // out_put_directory += "_Alpha=2_FeedbackStren=" + oss.str();
        // oss.str("");
        // oss << std::fixed << setprecision(0) << hill_coefficient;
        // out_put_directory += "_HillCoeff=" + oss.str();
        // oss.str("");
        // oss << std::fixed << setprecision(3) << edge_length_at_rest;
        // out_put_directory += "_ELs=" + oss.str();
        // oss.str("");
        // oss << std::scientific << setprecision(0) << feedback_strength_for_adhesion;
        // out_put_directory += "_AdhFeedbStren=" + oss.str();    
        // oss.str("");
        // oss << std::fixed << setprecision(0) << hill_coefficient_for_adhesion;
        // out_put_directory += "_AdhHillCoeff=" + oss.str();
        // out_put_directory += "_EMADontDe=" + std::to_string(EMA_dont_decrease);
        // out_put_directory += "_CCADontDe=" + std::to_string(CCA_dont_decrease);
        // out_put_directory += "_CCADontIn=" + std::to_string(CCA_dont_increase);
        // if (CCA_dont_inrease_until_shorter_than_a_threshold)  
        // {
        //   oss.str("");
        //   oss << std::fixed << setprecision(1) << CCA_dont_increase_until_shorter_than_this_value;          
        //   out_put_directory += "_ACCInreaseThresh=" + oss.str();
        // }

        out_put_directory += "/_RestrVertexMove=" + std::to_string(restrict_vertex_movement); 
        out_put_directory += "_UseAdaptiveTimestep=" + std::to_string(use_adaptive_timestep);          
        // oss.str("");
        // oss << std::fixed << setprecision(1) << leading_top_length;
        // out_put_directory += "_SubAdhLeadTopLen" + oss.str();
        out_put_directory += "_UseFineMesh=" + std::to_string(use_fine_mesh);
        oss.str("");
        oss << std::fixed << setprecision(3) << cell_rearrangement_threshold;
        out_put_directory += "_RearThresh=" + oss.str();
        // oss.str("");
        // oss << std::fixed << setprecision(1) << strip_start_x_location;
        // out_put_directory += "_StripStartXLoc=" + oss.str();

        oss.str("");
        oss << std::fixed << setprecision(3) << nagai_honda_membrane_surface_energy_parameter;
        out_put_directory += "/_Ga=" + oss.str();
        oss.str("");
        oss << std::fixed << setprecision(1) << fixed_target_shape_index;
        out_put_directory += "_ShapeIndex=" + oss.str();
        // oss.str("");
        // oss << std::fixed << setprecision(1) << substrate_adhesion_parameter_top_area;
        // out_put_directory += "_SubAdhParaTopArea=" + oss.str();
        // oss.str("");
        // oss << std::fixed << setprecision(1) << substrate_adhesion_parameter;
        // out_put_directory += "_SubAdhPara=" + oss.str();

        simulator.SetOutputDirectory(out_put_directory);
        std::cout << std::endl << "OutputDirectory: " << out_put_directory << std::endl;
        /*-------------------------------END: Output Directory---------------------------------*/

        simulator.Solve();
    }

};

#endif /* TESTMYFLUIDSTATE_HPP_ */
