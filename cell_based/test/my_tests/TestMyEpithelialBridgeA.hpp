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

#ifndef TESTMYEPITHELIALBRIDGEA_HPP_
#define TESTMYEPITHELIALBRIDGEA_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "BernoulliTrialCellCycleModel.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "TargetAreaLinearGrowthModifier.hpp"
#include "FaceValueAndStressStateModifier.hpp"
#include "PolarityModifier.hpp"

#include "FakePetscSetup.hpp"

#include "PlaneBoundaryCondition.hpp"
#include "ToroidalHoneycombVertexMeshGenerator.hpp"
#include "MyXToroidalHoneycombVertexMeshGenerator.hpp"
#include "NagaiHondaForce.hpp"
#include "MyNagaiHondaForceWithStripesAdhesion.hpp"
#include "DiffusionForce.hpp"

// #include <ctime>

class TestMyEpithelialBridgeA : public AbstractCellBasedTestSuite
{
public:

    void TestStripSubstrateAdhesion()
    {
        /*-----------------------START: Frequently changed parameters-------------------------*/
        // PHASE DIAGRAM SEARCH:
        //                       OPTIONS:
        // double feedback_strength_for_myosin_activity = 0.0;
        // double set_substrate_adhesion_parameter_at_leading_top= -10.0;
        // double set_substrate_adhesion_parameter_below_leading_top = -1.0;
        // double set_reservoir_substrate_adhesion_parameter = -1.0;
        // double set_nagai_honda_membrane_surface_energy_parameter = 1.0;
        // double set_target_shape_index = 4.0; // 3.7224 for default
        // double set_polarity_magnitude = 0.5;

        double feedback_strength_for_myosin_activity = 0.0;
        double set_nagai_honda_membrane_surface_energy_parameter = 1.0;
        double set_target_shape_index = 4.0; // 3.7224 for default

        double set_substrate_adhesion_parameter_at_leading_top= -10.0;
        double set_substrate_adhesion_parameter_below_leading_top = -1.0;
        double set_reservoir_substrate_adhesion_parameter = -1.0;

        double set_polarity_magnitude = 0.5;

        /*-----------------------START: Frequently changed parameters-------------------------*/
        // timestep
        double set_dt = 0.1;
        bool apply_my_change_to_make_timestep_adaptive = true;
        bool consider_consistency_for_SSA = true;
        double max_movement_per_timestep = 0.05;
        bool restrict_vertex_movement = false;
        if (apply_my_change_to_make_timestep_adaptive)
          assert(restrict_vertex_movement==false);
        double small_change_for_area_calculation = 0.2;

        // cell rearrangement
        double set_cell_rearrangement_threshold = 0.05;

        // cell division
        bool run_with_birth =true;
        double time_for_one_division_of_cell_population = 25;
        double growth_rate_for_target_area_after_division = 0.1;
        bool use_my_division_rule_along_with_modifier = true;

        // structure
        bool if_use_larger_strip_distance = false;
        bool one_period_only = true;
        bool if_mesh_has_two_period = false;
        if (if_mesh_has_two_period)
          one_period_only = false;
        bool use_longer_mesh = false;
        
        // energy parameter
        // double set_nagai_honda_membrane_surface_energy_parameter = 0.1;
        bool set_use_fixed_target_area = false;
        // double set_target_shape_index = 4.0; // {6/sqrt(6*sqrt(3)/4)}=3.7224 for default
          // substrate adhesion
        bool if_consider_substrate_adhesion = true;
            // strip substrate adhesion
        bool if_substrate_adhesion_is_homogeneous = false;
              // if homogeneous
        double set_homogeneous_substrate_adhesion_parameter = -3.0;
              // if not homogeneous
        double set_substrate_adhesion_leading_top_length = 3.0;
        // double set_substrate_adhesion_parameter_at_leading_top= -3.0;
        // double set_substrate_adhesion_parameter_below_leading_top = -0.5;
                // detach pattern:
        bool use_my_detach_pattern_method = true;
        if (use_my_detach_pattern_method)
          assert(if_substrate_adhesion_is_homogeneous == false);
            // reservoir substrate adhesion
        bool if_consider_reservoir_substrate_adhesion = true;// true for default
        if (!if_consider_substrate_adhesion)
          assert(if_consider_reservoir_substrate_adhesion==false);
        bool if_ignore_reservoir_substrate_adhesion_at_top = false;// false for default
        bool if_ignore_reservoir_substrate_adhesion_at_bottom = true;// false for default
        // double set_reservoir_substrate_adhesion_parameter = -0.5;
        
        // feedback
        bool is_default_feedback_form = true;
        int case_number_of_membrane_surface_energy_form = 2; // 2, for default, 2 for new membr_surf_energy form1, 3 for form2.
        bool if_consider_feedback_of_face_values = true;
        bool if_apply_feedback_of_face_values_only_for_boundary_cells = true; // for testing fluid inside
        bool if_apply_feedback_of_face_values_only_for_top_boundary_cells = true; // for testing fluid inside
        bool if_consider_feedback_of_cell_cell_adhesion = true;
        bool EMA_dont_decrease = false;
        bool CCA_dont_decrease = false;
        bool CCA_increasing_has_a_threshold_of_edge_length = true;
        double CCA_increasing_threshold_of_edge_length_percentage = 0.5;
        // double feedback_strength_for_myosin_activity = 0;
        double hill_coefficient_for_myosin_activity = 8.0;
        double feedback_strength_for_adhesion = 0;
        double hill_coefficient_for_adhesion = 8.0;
          // for consistency
        if (feedback_strength_for_myosin_activity == 0.0)
          if_consider_feedback_of_face_values = false; // note: typically, we must have myosin feedback first.
        if (feedback_strength_for_adhesion == 0.0)
          if_consider_feedback_of_cell_cell_adhesion = false;
        if (if_consider_feedback_of_face_values)
          assert(case_number_of_membrane_surface_energy_form != 0);
        else
          case_number_of_membrane_surface_energy_form = 0;
        if (is_default_feedback_form)
          assert(!if_consider_feedback_of_cell_cell_adhesion && !EMA_dont_decrease);
        
        // random force and polarity
        bool add_random_force = true;
        bool consider_polarity = true;
        bool vanishing_motility_for_node_in_the_strip_interval = false;
        // double set_polarity_magnitude = 0.2;
        double set_rotational_diffusion_constant = 0.01;

        // output
        bool output_concise_swap_information_when_remesh = false;
        bool output_detailed_swap_information_when_remesh = false; //suggest "false" for concise output results
        bool output_information_for_nagai_honda_force = false;
        bool set_output_numerical_method_information = false;
        /*-----------------------END: Frequently changed parameters-------------------------*/


        /*------------------------------START: Mesh Structure------------------------------*/
        unsigned num_ele_cross = 6;
        if (if_use_larger_strip_distance)
          num_ele_cross = 10;
        if (if_mesh_has_two_period)
          num_ele_cross *= 2;

        unsigned num_ele_up = 8;
        if (use_longer_mesh)
          num_ele_up *=2;

        double initial_area = M_PI;
        double cell_rearrangement_threshold = set_cell_rearrangement_threshold;
        bool if_update_face_elements_in_mesh = if_consider_feedback_of_face_values;

        MyXToroidalHoneycombVertexMeshGenerator generator(num_ele_cross, num_ele_up, initial_area, cell_rearrangement_threshold, 0.001);
        MyXToroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();
        p_mesh->SetUpdateFaceElementsInMeshBoolean(if_update_face_elements_in_mesh);
        p_mesh->SetOutputConciseSwapInformationWhenRemesh(output_concise_swap_information_when_remesh);
        p_mesh->SetOutputDetailedSwapInformationWhenRemesh(output_detailed_swap_information_when_remesh);
        /*------------------------------END: Mesh Structure------------------------------*/


        /*------------------------------START: Cells and Cell population--------------------------*/
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<BernoulliTrialCellCycleModel, 2> cells_generator;
        bool use_bernoulli_trial_cell_cycle_model = true;
        cells_generator.SetUseBernoulliTrialCellCycleModel(use_bernoulli_trial_cell_cycle_model);
        double divide_probability_for_a_cell = 1.0/time_for_one_division_of_cell_population/(num_ele_cross*num_ele_up);
        double minimum_division_age = -0.01;
        cells_generator.SetDivideProbability(divide_probability_for_a_cell);
        cells_generator.SetMinimumDivisionAge(minimum_division_age);
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);
        
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetRestrictVertexMovementBoolean(restrict_vertex_movement);
        /*------------------------------END: Cells and Cell population--------------------------*/


        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetNoBirth(!run_with_birth);


        /*--------------------------------START: TargetAreaModifier------------------------------*/
        MAKE_PTR_ARGS(TargetAreaLinearGrowthModifier<2>, p_growth_modifier, ());
        p_growth_modifier->SetUseUseMyOwnRuleForUpdateTargetAreaOfCell(true);
        double reference_target_area_for_modifier = initial_area;
        p_growth_modifier->SetReferenceTargetArea(reference_target_area_for_modifier);
        p_growth_modifier->SetGrowthRate(growth_rate_for_target_area_after_division);
        p_growth_modifier->SetAgeToStartGrowing(0.0);
        simulator.AddSimulationModifier(p_growth_modifier);
        /*--------------------------------END: TargetAreaModifier------------------------------*/


        /*---------------------------------START: Add Numerical Method-----------------------------*/
        boost::shared_ptr<ForwardEulerNumericalMethod<2,2>> p_numerical_method(new ForwardEulerNumericalMethod<2,2>());
        p_numerical_method->SetCellPopulation(&cell_population);
        //Note: Here the address of std::vector of force collections are different from that in simulator!
        //But the corresponding shared_ptrs in the vectors should be the same.
        std::vector< boost::shared_ptr<AbstractForce<2, 2>> > force_collection = simulator.rGetForceCollection();
        p_numerical_method->SetForceCollection(&force_collection);
        p_numerical_method->SetMaxMovementPerTimestep(max_movement_per_timestep);
        p_numerical_method->SetOutputNumericalMethodInformation(set_output_numerical_method_information);
        simulator.SetNumericalMethod(p_numerical_method);
        /*---------------------------------END: Add Numerical Method-----------------------------*/


        /*-----------------------------START: MyNagaiHondaForceWithStripesAdhesion---------------------*/
        MAKE_PTR(MyNagaiHondaForceWithStripesAdhesion<2>, p_force);
        
        // KA, A0
        double nagai_honda_deformation_energy_parameter = 1;
        bool use_fixed_target_area = set_use_fixed_target_area;
        double fixed_target_area = initial_area;

        // KP 
        double nagai_honda_membrane_surface_energy_parameter = set_nagai_honda_membrane_surface_energy_parameter;
                
        // LAMBDA: Use equation form: 0.5*Gamma*P^2, target perimeter being 0!
        double target_shape_index = set_target_shape_index;
        double nagai_honda_cell_cell_adhesion_energy_parameter = -nagai_honda_membrane_surface_energy_parameter*(target_shape_index*sqrt(reference_target_area_for_modifier));
        if (use_fixed_target_area)
          nagai_honda_cell_cell_adhesion_energy_parameter = -nagai_honda_membrane_surface_energy_parameter*(target_shape_index*sqrt(fixed_target_area));
        double nagai_honda_cell_boundary_adhesion_energy_parameter = 0; //nagai_honda_cell_cell_adhesion_energy_parameter;
        bool if_use_face_element_to_get_adhesion_parameter = if_consider_feedback_of_face_values;

        // Strip Substrate Adhesion Parameter
        double homogeneous_substrate_adhesion_parameter = set_homogeneous_substrate_adhesion_parameter;
        double substrate_adhesion_leading_top_length = set_substrate_adhesion_leading_top_length;
        double substrate_adhesion_parameter_at_leading_top= set_substrate_adhesion_parameter_at_leading_top;
        double substrate_adhesion_parameter_below_leading_top = set_substrate_adhesion_parameter_below_leading_top;
        // Reservoir Substrate Adhesion Parameter
        double reservoir_substrate_adhesion_parameter = set_reservoir_substrate_adhesion_parameter;
        
        // Strips structure of substrate adhesion
        double strip_width = sqrt(initial_area/(sqrt(3)/2))/2; // =0.952~1.05
        double strip_distance = 6*sqrt(initial_area/(sqrt(3)/2)); // =11.428~12.60
        if(if_use_larger_strip_distance)
          strip_distance *= 200.0/120.0;
        double strip_start_x_location = 0.0;
        if (if_mesh_has_two_period)
          strip_start_x_location = -3*sqrt(initial_area/(sqrt(3)/2));
        double strip_start_y_location = 1.5*num_ele_up*1/sqrt(3)*sqrt(initial_area/(sqrt(3)/2));

        p_force->SetCaseNumberOfMembraneSurfaceEnergyForm(case_number_of_membrane_surface_energy_form);
        p_force->SetUseFaceElementToGetAdhesionParameterBoolean(if_use_face_element_to_get_adhesion_parameter);

        p_force->SetNagaiHondaDeformationEnergyParameter(nagai_honda_deformation_energy_parameter);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(nagai_honda_membrane_surface_energy_parameter);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(nagai_honda_cell_cell_adhesion_energy_parameter);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(nagai_honda_cell_boundary_adhesion_energy_parameter);        

        p_force->SetUseFixedTargetArea(use_fixed_target_area);
        p_force->SetFixedTargetArea(fixed_target_area);
        p_force->SetTargetShapeIndex(target_shape_index);
          // substrate adhesion info
        p_force->SetIfConsiderSubstrateAdhesion(if_consider_substrate_adhesion);
          // Strip Substrate Adhesion
        p_force->SetIfSubstrateAdhesionIsHomogeneous(if_substrate_adhesion_is_homogeneous);
        p_force->SetHomogeneousSubstrateAdhesionParameter(homogeneous_substrate_adhesion_parameter);
        p_force->SetSubstrateAdhesionLeadingTopLength(substrate_adhesion_leading_top_length);
        p_force->SetSubstrateAdhesionParameterAtLeadingTop(substrate_adhesion_parameter_at_leading_top);
        p_force->SetSubstrateAdhesionParameterBelowLeadingTop(substrate_adhesion_parameter_below_leading_top);
        p_force->SetConsiderConsistencyForSSA(consider_consistency_for_SSA);
        p_force->SetSmallChangeForAreaCalculation(small_change_for_area_calculation);
          // Reservoir Substrate Adhesion
        p_force->SetIfConsiderReservoirSubstrateAdhesion(if_consider_reservoir_substrate_adhesion);
        p_force->SetIfIgnoreReservoirSubstrateAdhesionAtTop(if_ignore_reservoir_substrate_adhesion_at_top);
        p_force->SetIfIgnoreReservoirSubstrateAdhesionAtBottom(if_ignore_reservoir_substrate_adhesion_at_bottom);
            // for consideration of periodicity
        p_force->SetCenterOfWidth(0.0);
        p_force->SetWidth(num_ele_cross*sqrt(initial_area/(sqrt(3)/2)));//check later!
        p_force->SetReservoirSubstrateAdhesionParameter(reservoir_substrate_adhesion_parameter);
          // strip info
        p_force->SetStripWidth(strip_width);
        p_force->SetStripDistance(strip_distance);
        p_force->SetStripStartXLocation(strip_start_x_location);
        p_force->SetStripStartYLocation(strip_start_y_location);
          // method for detachment behavior
        p_force->SetUseMyDetachPatternMethod(use_my_detach_pattern_method);

        p_force->SetOutputInformationForNagaiHondaForce(output_information_for_nagai_honda_force);
        simulator.AddForce(p_force);
        /*-----------------------------END: MyNagaiHondaForceWithStripesAdhesion---------------------*/


        /*-----------------------------START: RandomForce and PolarityModifier--------------------------*/
        // Brownian diffusion
        double translational_diffusion_constant = 0.0;
        if (add_random_force)
        {
          MAKE_PTR(DiffusionForce<2>, p_force2);
          bool use_the_same_node_radius = true;
          double node_radius = 50*1e2;
          p_force2->SetUseTheSameNodeRadius(use_the_same_node_radius);
          p_force2->SetTheSameNodeRadius(node_radius);
          p_force2->SetConsiderPolarity(consider_polarity);
          p_force2->SetVanishingMotilityForNodeInTheStripInterval(vanishing_motility_for_node_in_the_strip_interval);
          p_force2->SetOnePeriodOnly(one_period_only);
          p_force2->SetReservoirTop(strip_start_y_location);
          p_force2->SetStripStartLocation(strip_start_x_location);
          p_force2->SetStripWidth(strip_width);

          translational_diffusion_constant = p_force2->GetDiffusionScalingConstant()/node_radius;
          simulator.AddForce(p_force2);
        }

        // Polarity modifier
        double polarity_magnitude = 0.0;
        double rotational_diffusion_constant = 0.0;
        if (consider_polarity)
        {
          MAKE_PTR_ARGS(PolarityModifier<2>, p_polarity_modifier, ());
          polarity_magnitude = set_polarity_magnitude;
          rotational_diffusion_constant = set_rotational_diffusion_constant;
          p_polarity_modifier->SetPolarityMagnitude(polarity_magnitude);
          p_polarity_modifier->SetD(rotational_diffusion_constant);
          simulator.AddSimulationModifier(p_polarity_modifier);
        }
        /*-----------------------------END: RandomForce and PolarityModifier---------------------------*/


        /*------------------------START: !!!!!Feedback: FaceValueAndStressStateModifier: need modification---------------*/
        MAKE_PTR_ARGS(FaceValueAndStressStateModifier<2>, p_face_value_and_stress_state_modifier, ());
                
        p_face_value_and_stress_state_modifier->SetConsiderFeedbackOfFaceValues(if_consider_feedback_of_face_values);
        p_face_value_and_stress_state_modifier->SetConsiderFeedbackOfFaceValuesOnlyForBoundaryCells(if_apply_feedback_of_face_values_only_for_boundary_cells);
        p_face_value_and_stress_state_modifier->SetConsiderFeedbackOfFaceValuesOnlyForTopBoundaryCells(if_apply_feedback_of_face_values_only_for_top_boundary_cells);
        p_face_value_and_stress_state_modifier->SetConsiderFeedbackOfCellCellAdhesion(if_consider_feedback_of_cell_cell_adhesion);
        p_face_value_and_stress_state_modifier->SetEMADontDecrease_CCADontDecrease_HasAThreshold_Threshold(
            EMA_dont_decrease, CCA_dont_decrease, CCA_increasing_has_a_threshold_of_edge_length, CCA_increasing_threshold_of_edge_length_percentage);
        
        // Feedback information
        double edge_length_at_rest = sqrt(initial_area/(6*sqrt(3)/4)); // = 1.0996
        p_face_value_and_stress_state_modifier->SetEdgeLengthAtRest(edge_length_at_rest);
        p_face_value_and_stress_state_modifier->SetFeedbackStrengthForMyosinActivity(feedback_strength_for_myosin_activity);
        p_face_value_and_stress_state_modifier->SetHillCoefficientForMyosinActivity(hill_coefficient_for_myosin_activity);
        p_face_value_and_stress_state_modifier->SetFeedbackStrengthForAdhesion(feedback_strength_for_adhesion);
        p_face_value_and_stress_state_modifier->SetHillCoefficientForAdhesion(hill_coefficient_for_adhesion);

        // my stress state modifier
        p_face_value_and_stress_state_modifier->SetCalculateStressStateBoolean(true);
        p_face_value_and_stress_state_modifier->SetCaseNumberOfMembraneSurfaceEnergyForm(case_number_of_membrane_surface_energy_form);        
        p_face_value_and_stress_state_modifier->SetNagaiHondaMembraneSurfaceEnergyParameter(nagai_honda_membrane_surface_energy_parameter);
        p_face_value_and_stress_state_modifier->SetFixedTargetPerimeter(target_shape_index*sqrt(fixed_target_area));
        
        // my group number modifier
        p_face_value_and_stress_state_modifier->SetWriteGroupNumberToCell(use_my_detach_pattern_method);

        // my division state modifier
        p_face_value_and_stress_state_modifier->SetUseMyDivisionRuleAlongWithModifier(use_my_division_rule_along_with_modifier);
        p_face_value_and_stress_state_modifier->SetDivisionTime(time_for_one_division_of_cell_population);
        
        // outout information
        p_face_value_and_stress_state_modifier->SetOutputModifierInformationBoolean(false);

        simulator.AddSimulationModifier(p_face_value_and_stress_state_modifier);
        /*-----------------------END: !!!!!!Feedback: FaceValueAndStressStateModifier: need modification---------------*/


        /*------------------------------------START: Timestep---------------------------------------*/
        double dt = set_dt;
        double sampling_timestep_multiple = round(1/dt);

        simulator.SetApplyMyChangesToMakeTimestepAdaptive(apply_my_change_to_make_timestep_adaptive);
        simulator.SetDt(dt);

        simulator.SetApplySamplingTimeInsteadOfSamplingTimestep(apply_my_change_to_make_timestep_adaptive);
        simulator.SetSamplingTimestepMultiple(sampling_timestep_multiple);
        simulator.SetSamplingTime(1.0);

        simulator.SetEndTime(400.0);
        /*------------------------------------END: Timestep---------------------------------------*/


        /*--------------------------------START: Boundary condition-----------------------------*/
          c_vector<double,2> point = zero_vector<double>(2);
          c_vector<double,2> normal = zero_vector<double>(2);
          normal(1) = -1.0;
          double stop_time = DOUBLE_UNSET;
          MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc, (&cell_population, point, normal, stop_time));
          simulator.AddCellPopulationBoundaryCondition(p_bc);
        /*--------------------------------END: Boundary condition-----------------------------*/
        

        /*--------------------------START: Output Directory and Simulation Information File---------------------*/
        // output directory
        std::ostringstream oss;
        std::string output_directory = 
            "EpithelialBridgeSimulation/PHASE-DIAGRAM/Simulation Results Start From: 20-07-22/";

        oss.str("");
        oss << "_MyoFeStr=" << std::fixed << setprecision(2) << feedback_strength_for_myosin_activity
            << "_p0=" << std::fixed << setprecision(2) << target_shape_index
            << "_Ga=" << std::fixed << setprecision(2) << nagai_honda_membrane_surface_energy_parameter
            << "_SSATop=" << std::fixed << setprecision(1) << substrate_adhesion_parameter_at_leading_top
            << "_Bott=" << std::fixed << setprecision(1) << substrate_adhesion_parameter_below_leading_top
            << "_RSA=" << std::fixed << setprecision(1) << reservoir_substrate_adhesion_parameter
            << "_fp=" << ((polarity_magnitude>=0.01 || polarity_magnitude==0.0)? std::fixed : std::scientific) << setprecision(2) << polarity_magnitude;
        output_directory += oss.str();
        time_t raw_time = time(0);
        struct tm * now = localtime(& raw_time);
        oss.str("");
        oss << (now->tm_year + 1900 -2000) << '-' << (now->tm_mon + 1) << '-' <<  now->tm_mday 
            << ", " << now->tm_hour << ':' << now->tm_min << ':' << now->tm_sec;
        output_directory += "_timestamp=" + oss.str();

        std::string concise_output_directory = output_directory;
        
        // concise information written to directoory.
        oss.str("");
        oss << std::fixed << setprecision(3) << dt;
        output_directory += "/Dt=" + oss.str();
        oss.str("");
        oss << ((cell_rearrangement_threshold>=0.01)? std::fixed : std::scientific) << setprecision(2) << cell_rearrangement_threshold;
        output_directory += "_RearThr=" + oss.str();
        if (apply_my_change_to_make_timestep_adaptive)
        {
          oss.str("");
          oss << ((max_movement_per_timestep>=0.01)? std::fixed : std::scientific) << setprecision(2) << max_movement_per_timestep;
          output_directory += "_MaxMvDt=" + oss.str();
        }
        if (if_consider_substrate_adhesion)
        {
          oss.str("");
          oss << ((small_change_for_area_calculation>=0.01)? std::fixed : std::scientific) << setprecision(2) << small_change_for_area_calculation;
          output_directory += "_Dxy=" + oss.str();
        }
        if (!apply_my_change_to_make_timestep_adaptive)
          output_directory += "_MyAdaptDt=0";
        if (!restrict_vertex_movement)
          output_directory += "_ResMv=0";
        if (!consider_consistency_for_SSA)
          output_directory += "_ConsistMv=0";

        output_directory += "_|Divi=" + std::to_string(run_with_birth);
        if (if_consider_feedback_of_face_values && is_default_feedback_form)
          output_directory += "_HasDefaultFeedb";
        else
          output_directory += "_HasFeedb=" + std::to_string(if_consider_feedback_of_face_values);
        output_directory += "_HasRandF=" + std::to_string(add_random_force);
        output_directory += "_MSE=" + std::to_string(case_number_of_membrane_surface_energy_form);
        
        oss.str("");
        oss << std::fixed << setprecision(2) << nagai_honda_membrane_surface_energy_parameter;
        output_directory += "_|Ga=" + oss.str();
        oss.str("");
        oss << std::fixed << setprecision(2) << target_shape_index;
        output_directory += "_p0=" + oss.str();
        oss.str("");
        oss << ((fabs(nagai_honda_cell_cell_adhesion_energy_parameter)>=0.01)? std::fixed : std::scientific) << setprecision(2) << nagai_honda_cell_cell_adhesion_energy_parameter;
        output_directory += "_CCAdhe=" + oss.str();
        oss.str("");
        oss << ((nagai_honda_cell_boundary_adhesion_energy_parameter==0.0)? std::fixed : std::scientific) << setprecision(1) << nagai_honda_cell_boundary_adhesion_energy_parameter;
        output_directory += "_CBAdhe=" + oss.str();
        
        output_directory += "_|HasSubAdh=" + std::to_string(if_consider_substrate_adhesion);
        if (if_consider_substrate_adhesion)
        {
          if (if_substrate_adhesion_is_homogeneous)
            output_directory += "_SSAIsHomo";
          else
            output_directory += "_NonHomoSSA";       
          output_directory += "_HasRSA=" + std::to_string(if_consider_reservoir_substrate_adhesion);
        }
        
        // detailed information  written to directoory.
        output_directory += "/";
        if (run_with_birth)
        {
          output_directory += "CellDivi:";
          oss.str("");
          oss << std::fixed << setprecision(1) << time_for_one_division_of_cell_population;
          if (use_my_division_rule_along_with_modifier)
            output_directory += "FixedDiviDt=" + oss.str();
          else
            output_directory += "AveragedDiviDt=" + oss.str();
          oss.str("");
          oss << std::fixed << setprecision(2) << growth_rate_for_target_area_after_division;
          output_directory += "_GrRate=" + oss.str();
        }
        else
          output_directory += "NoCellDivi";

        if (if_consider_feedback_of_face_values)
        {
          // feedback form
          output_directory += "_|FeedbackForm:";
          if (is_default_feedback_form)
            output_directory += "Default";
          else
          {
            output_directory += "HasAdhFeedb=" + std::to_string(if_consider_feedback_of_cell_cell_adhesion);
            output_directory += "_EMACanDe=" + std::to_string(!EMA_dont_decrease);
            if (if_consider_feedback_of_cell_cell_adhesion)
            {
              output_directory += "_CCACanDe=" + std::to_string(!CCA_dont_decrease);
              if (CCA_increasing_has_a_threshold_of_edge_length)  
              {
                oss.str("");
                oss << std::fixed << setprecision(2) << CCA_increasing_threshold_of_edge_length_percentage;          
                output_directory += "_CCAInrThresh=" + oss.str();
              }
            }
          }
          // feedback parameters
          output_directory += "_|FeedbackPara:";
          oss.str("");
          oss << ((feedback_strength_for_myosin_activity<=100.0)? std::fixed : std::scientific) << setprecision(2) << feedback_strength_for_myosin_activity ;
          output_directory += "MyoFeStr=" + oss.str();
          oss.str("");
          oss << std::fixed << setprecision(1) << hill_coefficient_for_myosin_activity ;
          output_directory += "_MyoHill=" + oss.str();
          if (if_consider_feedback_of_cell_cell_adhesion)
          {
            oss.str("");
            oss << ((feedback_strength_for_adhesion<=100.0)? std::fixed : std::scientific) << setprecision(2) << feedback_strength_for_adhesion;
            output_directory += "_AdhFeStr=" + oss.str();
            oss.str("");
            oss << std::fixed << setprecision(1) << hill_coefficient_for_adhesion;
            output_directory += "_AdhHill=" + oss.str();
          }
        }

        if (add_random_force)
        {
          oss.str("");
          oss << std::scientific << setprecision(2) << translational_diffusion_constant;
          output_directory += "_|D_trans=" + oss.str();
          if (consider_polarity)
          {
            oss.str("");
            oss << ((polarity_magnitude>=0.01 || polarity_magnitude==0.0)? std::fixed : std::scientific) << setprecision(2) << polarity_magnitude;
            output_directory += "_fp=" + oss.str();
            oss.str("");
            oss << ((rotational_diffusion_constant>=0.01)? std::fixed : std::scientific) << setprecision(2) << rotational_diffusion_constant;
            output_directory += "_Dr=" + oss.str();
            if (vanishing_motility_for_node_in_the_strip_interval)
              output_directory += "_IntervalNodeMotility=0";
          }
        }   

        if(if_consider_substrate_adhesion)
        {
          if (if_substrate_adhesion_is_homogeneous)
          {
            oss.str("");
            oss << std::fixed << setprecision(1) << homogeneous_substrate_adhesion_parameter;
            output_directory += "_|HomoSSA=" + oss.str();
          }
          else
          {
            oss.str("");
            oss << std::fixed << setprecision(1) << substrate_adhesion_leading_top_length;
            output_directory += "_|SSA:LeadLeng=" + oss.str();
            oss.str("");
            oss << std::fixed << setprecision(1) << substrate_adhesion_parameter_at_leading_top;
            output_directory += "_Top=" + oss.str();
            oss.str("");
            oss << std::fixed << setprecision(1) << substrate_adhesion_parameter_below_leading_top;
            output_directory += "_Below=" + oss.str();
          }
          if(if_consider_reservoir_substrate_adhesion)
          {
            oss.str("");
            oss << std::fixed << setprecision(1) << reservoir_substrate_adhesion_parameter;
            output_directory += "_RSA=" + oss.str();
            if (if_ignore_reservoir_substrate_adhesion_at_top)
              output_directory += "_ConsiRSATop=0";
            if (if_ignore_reservoir_substrate_adhesion_at_bottom)
              output_directory += "_ConsiRSABott=0";

          }
        }

        bool omit_file_name_results_from_time_X = true;
        bool output_simulatin_information_to_file = true;
        simulator.SetOmitFileNameResultsFromTimeX(omit_file_name_results_from_time_X);
        simulator.SetOutputSimulationInformationToFile(output_simulatin_information_to_file);
        if (output_simulatin_information_to_file)
          simulator.InputSimulationInformation(output_directory);

        bool use_concise_output_directory = true;
        if (use_concise_output_directory)
        {
          simulator.SetOutputDirectory(concise_output_directory);
          std::cout << std::endl << "Concise output directory is set: " << concise_output_directory << std::endl;
        }
        else
        {
          simulator.SetOutputDirectory(output_directory);
          std::cout << std::endl << "Output directory is set: " << output_directory << std::endl;
        }
        /*--------------------------END: Output Directory and Simulation Information File---------------------*/

        simulator.Solve();
    }

};

#endif /* TESTMYEPITHELIALBRIDGEA_HPP_ */
