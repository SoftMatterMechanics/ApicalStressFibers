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

#ifndef TESTMYNEWEB_HPP_
#define TESTMYNEWEB_HPP_

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
#include "MyNewNagaiHondaForceWithStripesAdhesion.hpp"
#include "DiffusionForce.hpp"
#include "PlaneBasedCellKiller.hpp"

// #include <ctime>

class TestMyNewEB : public AbstractCellBasedTestSuite
{
public:

    void TestStripSubstrateAdhesion()
    {
        // assert(false);

        /* Fomulation change history: (A-PI)^2-1.0*S; =>(A'-1)^2-(1.0/PI)*S'; =>(A'-(1+1.0/PI))^2+...
        *                                                   ^                         ^
        *                                              reference_area         real_reference_area
        */
        double reference_area = 1.0;
        double real_reference_area = (M_PI+1.0)/M_PI;

        bool larger_strip_spacing_distance = false;
        double strip_dis_multiplied = 20.0/6.0;
        int move_mesh_right_for_N_periods = 0;
        double multiply_results_by = 1/sqrt(reference_area);

        double basic_SSA = +1.0/(M_PI/reference_area);

        double initial_area = reference_area;
        double strip_width = 5.5*sqrt(initial_area/(sqrt(3)/2)); // =0.952~1.05
        double strip_distance = 6*sqrt(initial_area/(sqrt(3)/2)); // =11.428~12.60
        
        // FOR PHASE DIAGRAM SEARCH:
        double target_shape_index = 4.5;//p0

        double pulling_force_on_leading_cell = 8/pow((M_PI/reference_area),1.5);// Fy

        bool if_consider_feedback_of_face_values = true;
        double feedback_strength_for_myosin_activity = 0.005/(M_PI/reference_area);//Fb

        double Gamma = 0.2/(M_PI/reference_area);//Ga

        bool run_with_birth =false;

        bool has_brownian_random_force = true;
        double translational_diffusion_constant = 0.0;
        bool has_polarity = true;
        double polarity_magnitude = 0.01;
        double rotational_diffusion_constant = 0.01/(M_PI/reference_area); //0.2*2.0*(M_PI/reference_area);

        double dt = 0.1*(M_PI/reference_area);
        double sampling_time = 1.0*(M_PI/reference_area);
        double end_time = 400.0*(M_PI/reference_area);
        double time_for_equilibrium = 0.0;

        bool restrict_vertex_movement = false; // cell population method
        bool apply_my_change_to_make_timestep_adaptive = true;
        double max_movement_per_timestep_using_my_adaptive_timestep_method = 0.05/sqrt((M_PI/reference_area));
        if (apply_my_change_to_make_timestep_adaptive)
          assert(restrict_vertex_movement==false && max_movement_per_timestep_using_my_adaptive_timestep_method!=0.0);

        
        /*-----------------------START: Frequently changed parameters-------------------------*/
        // Substrate adhesion:
        double SSA_for_mature_lamellipodium = -10.0/(M_PI/reference_area);
        
        // Strip substrate adhesion form:
/*----*/bool if_substrate_adhesion_is_homogeneous = true;
          // some basic mechanisms:
        bool classify_elements_with_group_numbers = true;
        bool mark_leading_cells = true;
        bool kill_cells_group_number_from1 = true;
        bool SSA_strengthened_only_in_y_direction = true;
          // homogeneous SSA case:
        bool add_pulling_force_on_node_individually = false;
/*----*/bool add_pulling_force_evenly_on_nodes_of_leading_cell = true;
          // non-homogeneous SSA case:
            // one leading top:
        double substrate_adhesion_leading_top_length = 3.0/sqrt((M_PI/reference_area));
            // several leading tops:
        bool use_my_detach_pattern_method = false;
            // distribution rule:
/*----*/bool use_new_SSA_distribution_rule = false;
        if (use_new_SSA_distribution_rule)
          assert(mark_leading_cells);
        double lamellipodium_maturation_rate = 1.0/(M_PI/reference_area);
        double lamellipodium_destruction_rate = 0.1/(M_PI/reference_area);
              // assistant methods for distribution rule:
        bool small_SSA_at_first = false;
        double initial_time_for_small_SSA = 5.0*(M_PI/reference_area);
        double small_SSA_for_initial_time = -6.0/(M_PI/reference_area);
        if (small_SSA_at_first)
          assert(use_new_SSA_distribution_rule && fabs(SSA_for_mature_lamellipodium) >= fabs(small_SSA_for_initial_time));
        bool keep_moving_forward = false;
        double SSA_bottom_decrease = 5.0/(M_PI/reference_area);
        double slowly_moving_forward_after_this_height = 25.0/sqrt((M_PI/reference_area));
        if (keep_moving_forward)
          assert(use_new_SSA_distribution_rule);

        // Cell division:
        double time_for_one_division_of_cell_population = 25*(M_PI/reference_area);
        double growth_rate_for_target_area_after_division = 0.1/(M_PI/reference_area);
        bool use_my_division_rule_along_with_modifier = true;

        // Strip substrate adhesion (SSA):
        assert( !(if_substrate_adhesion_is_homogeneous&&(use_my_detach_pattern_method||use_new_SSA_distribution_rule)) );
        assert( !(use_new_SSA_distribution_rule && use_my_detach_pattern_method) );
        /*-----------------------END: Frequently changed parameters-------------------------*/



        /*------------------------------START: Mesh Structure------------------------------*/
        unsigned num_ele_cross = 6;
        // bool larger_strip_spacing_distance = false;
        // double strip_dis_multiplied = 10.0/6.0;
        if (larger_strip_spacing_distance)
          num_ele_cross = (unsigned) round(num_ele_cross*strip_dis_multiplied);

        bool if_mesh_has_two_period = false;
        if (if_mesh_has_two_period)
          num_ele_cross *= 2;

        unsigned num_ele_up = 8;
        bool use_longer_mesh = false;
        if (use_longer_mesh)
          num_ele_up *=2;

        // double initial_area = reference_area;
        double cell_rearrangement_threshold = 0.05/sqrt((M_PI/reference_area));

        MyXToroidalHoneycombVertexMeshGenerator generator(num_ele_cross, num_ele_up, initial_area, cell_rearrangement_threshold, 0.001/sqrt((M_PI/reference_area)));
        MyXToroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();

        p_mesh->SetDistanceForT3SwapChecking(5.0/sqrt((M_PI/reference_area)));

        // int move_mesh_right_for_N_periods = 0;
        p_mesh->SetMoveMeshRightForNPeriods(move_mesh_right_for_N_periods);
        // double multiply_results_by = 1/sqrt(reference_area);
        p_mesh->SetMultiplyResultsBy(multiply_results_by);

        bool if_update_face_elements_in_mesh = if_consider_feedback_of_face_values;
        p_mesh->SetUpdateFaceElementsInMeshBoolean(if_update_face_elements_in_mesh);

        p_mesh->SetIfClassifyElementsWithGroupNumbers(classify_elements_with_group_numbers);
        p_mesh->SetMarkLeadingCells(mark_leading_cells);

        bool if_check_for_T4_swaps = false;
        bool perform_T4_swaps_only_when_cell_is_triangle = true;
        double typical_length_for_T4_swaps = 2.0/sqrt((M_PI/reference_area));
        double separation_ratio_with_rearr_thresh = 4.0;
        p_mesh->SetIfCheckForT4Swaps(if_check_for_T4_swaps);
        p_mesh->SetPerformT4SwapsOnlyWhenCellIsTriangle(perform_T4_swaps_only_when_cell_is_triangle);
        p_mesh->SetTypicalLengthForT4Swaps(typical_length_for_T4_swaps);
        p_mesh->SetSeparationRatioWithRearrThresh(separation_ratio_with_rearr_thresh);

        bool output_concise_swap_information_when_remesh = false;
        p_mesh->SetOutputConciseSwapInformationWhenRemesh(output_concise_swap_information_when_remesh);
        bool output_detailed_swap_information_when_remesh = false;
        p_mesh->SetOutputDetailedSwapInformationWhenRemesh(output_detailed_swap_information_when_remesh);
        /*------------------------------END: Mesh Structure------------------------------*/


        /*------------------------------START: Cells and Cell population--------------------------*/
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<BernoulliTrialCellCycleModel, 2> cells_generator;
        bool use_bernoulli_trial_cell_cycle_model = true;
        cells_generator.SetUseBernoulliTrialCellCycleModel(use_bernoulli_trial_cell_cycle_model);
        double divide_probability_for_a_cell = 1/(M_PI/reference_area)*1.0/time_for_one_division_of_cell_population/(num_ele_cross*num_ele_up);
        double minimum_division_age = -0.01;
        cells_generator.SetDivideProbability(divide_probability_for_a_cell);
        cells_generator.SetMinimumDivisionAge(minimum_division_age);
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);
        
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // bool restrict_vertex_movement = false;
        cell_population.SetRestrictVertexMovementBoolean(restrict_vertex_movement);
        /*------------------------------END: Cells and Cell population--------------------------*/


        OffLatticeSimulation<2> simulator(cell_population);
        // bool run_with_birth =true;
        simulator.SetNoBirth(!run_with_birth);


        /*--------------------------------START: TargetAreaModifier------------------------------*/
        MAKE_PTR_ARGS(TargetAreaLinearGrowthModifier<2>, p_growth_modifier, ());
        p_growth_modifier->SetUseUseMyOwnRuleForUpdateTargetAreaOfCell(true);
        p_growth_modifier->SetReferenceTargetArea(real_reference_area);
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
        
        // double max_movement_per_timestep_using_my_adaptive_timestep_method = 0.1;
        p_numerical_method->SetMaxMovementPerTimestep(max_movement_per_timestep_using_my_adaptive_timestep_method);

        bool output_numerical_method_information = false;
        p_numerical_method->SetOutputNumericalMethodInformation(output_numerical_method_information);
        simulator.SetNumericalMethod(p_numerical_method);
        /*---------------------------------END: Add Numerical Method-----------------------------*/


        /*-----------------------------START: MyNagaiHondaForceWithStripesAdhesion---------------------*/
        MAKE_PTR(MyNewNagaiHondaForceWithStripesAdhesion<2>, p_force);
        
        // Strips structure of substrate adhesion
        // double strip_width = sqrt(initial_area/(sqrt(3)/2))/2; // =0.952~1.05
        // double strip_distance = 6*sqrt(initial_area/(sqrt(3)/2)); // =11.428~12.60
        if(larger_strip_spacing_distance)
          strip_distance *= strip_dis_multiplied;
        double strip_start_x_location = 3.0*sqrt(initial_area/(sqrt(3)/2));
        if (if_mesh_has_two_period)
          strip_start_x_location = -3*sqrt(initial_area/(sqrt(3)/2));
        double strip_start_y_location = 1.5*num_ele_up*1/sqrt(3)*sqrt(initial_area/(sqrt(3)/2));

        int case_number_of_membrane_surface_energy_form = 2; // 2, for default, 2 for new membr_surf_energy form1, 3 for form2.
        p_force->SetCaseNumberOfMembraneSurfaceEnergyForm(case_number_of_membrane_surface_energy_form);
        bool if_use_face_element_to_get_adhesion_parameter = if_consider_feedback_of_face_values;
        p_force->SetUseFaceElementToGetAdhesionParameterBoolean(if_use_face_element_to_get_adhesion_parameter);

        double nagai_honda_deformation_energy_parameter = 1;
        p_force->SetNagaiHondaDeformationEnergyParameter(nagai_honda_deformation_energy_parameter);//KA
        // double Gamma = 0.2/(M_PI/reference_area);//Ga
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(Gamma);//KP
        double cell_cell_adhesion_parameter = -Gamma*(target_shape_index*sqrt(reference_area));
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(cell_cell_adhesion_parameter);
        double cell_boundary_adhesion_parameter = 0.0;
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(cell_boundary_adhesion_parameter);        

        bool use_fixed_target_area_without_modifier = false;
        p_force->SetUseFixedTargetArea(use_fixed_target_area_without_modifier); // used in the case where there is no target area modifier! (no division)
        p_force->SetFixedTargetArea(real_reference_area);
        // double target_shape_index = 0.0;
        p_force->SetTargetShapeIndex(target_shape_index);
          // substrate adhesion info
        bool if_consider_substrate_adhesion = true;
        p_force->SetIfConsiderSubstrateAdhesion(if_consider_substrate_adhesion);
          // Strip Substrate Adhesion
        //double basic_SSA = -1.0/(M_PI/reference_area);
        p_force->SetBasicSSA(basic_SSA);
        p_force->SetIfSubstrateAdhesionIsHomogeneous(if_substrate_adhesion_is_homogeneous);
        double homogeneous_substrate_adhesion_parameter = basic_SSA;
        p_force->SetHomogeneousSubstrateAdhesionParameter(homogeneous_substrate_adhesion_parameter);
        p_force->SetAddPullingForceOnNodeIndividually(add_pulling_force_on_node_individually);
        p_force->SetAddPullingForceEvenlyOnNodesOfLeadingCell(add_pulling_force_evenly_on_nodes_of_leading_cell);
        // double pulling_force_on_leading_cell = 11.0;
        p_force->SetPullingForceOnLeadingCell(pulling_force_on_leading_cell);
        p_force->SetSubstrateAdhesionLeadingTopLength(substrate_adhesion_leading_top_length);
        p_force->SetSSAForMatureLamellipodium(SSA_for_mature_lamellipodium);
        double small_change_for_area_calculation = 0.4/sqrt((M_PI/reference_area));
        p_force->SetSmallChangeForAreaCalculation(small_change_for_area_calculation);
            // for consideration of periodicity
        p_force->SetCenterOfWidth(0.0);
        p_force->SetWidth(num_ele_cross*sqrt(initial_area/(sqrt(3)/2)));//check later!
          // strip info
        p_force->SetStripWidth(strip_width);
        p_force->SetStripDistance(strip_distance);
        p_force->SetStripStartXLocation(strip_start_x_location);
        p_force->SetStripStartYLocation(strip_start_y_location);
          // method for detachment behavior
        p_force->SetUseMyDetachPatternMethod(use_my_detach_pattern_method);
          // new SSA distribution rule
        p_force->SetSSAStrengthenedOnlyInYDirection(SSA_strengthened_only_in_y_direction);
        p_force->SetIfUseNewSSADistributionRule(use_new_SSA_distribution_rule);
        p_force->SetSmallSSAAtFirst(small_SSA_at_first);
        p_force->SetInitialTimeForSmallSSA(initial_time_for_small_SSA);
        p_force->SetSmallSSAForInitialTime(small_SSA_for_initial_time);
        p_force->SetKeepMovingForward(keep_moving_forward);
        p_force->SetSSABottomDecrease(SSA_bottom_decrease);
        p_force->SetSlowlyMovingForwardAfterThisHeight(slowly_moving_forward_after_this_height);

        bool if_equilibrate_for_a_while = true;
        if (time_for_equilibrium == 0.0)
          if_equilibrate_for_a_while = false;
        p_force->SetIfEquilibrateForAWhile(if_equilibrate_for_a_while);
        p_force->SetTimeForEquilibrium(time_for_equilibrium);

        simulator.AddForce(p_force);
        /*-----------------------------END: MyNagaiHondaForceWithStripesAdhesion---------------------*/


        /*-----------------------------START: RandomForce and PolarityModifier--------------------------*/
        bool add_random_force = true;
        // for consistency:
        if (!has_brownian_random_force && polarity_magnitude==0.0)
          add_random_force = false;
        if (polarity_magnitude!=0.0)
          assert(add_random_force==true);
        if (polarity_magnitude==0.0)
          has_polarity = false;

        if (add_random_force)
        {
          MAKE_PTR(DiffusionForce<2>, p_force2);
          bool use_the_same_node_radius = true;
          double node_radius = 50*1e2*((M_PI/reference_area)*(M_PI/reference_area)); //50*1e2*((M_PI/reference_area)*(M_PI/reference_area));
          if (translational_diffusion_constant >0.0)
            node_radius = p_force2->GetDiffusionScalingConstant()/translational_diffusion_constant;
          translational_diffusion_constant = p_force2->GetDiffusionScalingConstant()/node_radius;
          p_force2->SetIsNoBrownianRandomForce(!has_brownian_random_force);
          p_force2->SetUseTheSameNodeRadius(use_the_same_node_radius);
          p_force2->SetTheSameNodeRadius(node_radius);
          p_force2->SetConsiderPolarity(has_polarity);

          p_force2->SetIfEquilibrateForAWhile(if_equilibrate_for_a_while);
          p_force2->SetTimeForEquilibrium(time_for_equilibrium);

          simulator.AddForce(p_force2);
        }

        // Polarity modifier
        if (has_polarity)
        {
          MAKE_PTR_ARGS(PolarityModifier<2>, p_polarity_modifier, ());
          p_polarity_modifier->SetPolarityMagnitude(polarity_magnitude);
          p_polarity_modifier->SetD(rotational_diffusion_constant);
          simulator.AddSimulationModifier(p_polarity_modifier);
        }
        /*-----------------------------END: RandomForce and PolarityModifier---------------------------*/


        /*------------------------START: !!!!!Feedback: FaceValueAndStressStateModifier: need modification---------------*/
        MAKE_PTR_ARGS(FaceValueAndStressStateModifier<2>, p_face_value_and_stress_state_modifier, ());

        bool is_default_feedback_form = true;
        bool if_apply_feedback_of_face_values_only_for_boundary_cells = true; // for testing fluid inside
        bool if_apply_feedback_of_face_values_only_for_top_boundary_cells = true; // for testing fluid inside
        bool if_consider_feedback_of_cell_cell_adhesion = true;
        bool EMA_dont_decrease = false;
        bool CCA_dont_decrease = false;
        bool CCA_increasing_has_a_threshold_of_edge_length = true;
        double CCA_increasing_threshold_of_edge_length_percentage = 0.5;
        // double feedback_strength_for_myosin_activity = 0.0;
        double hill_coefficient_for_myosin_activity = 8.0;
        double feedback_strength_for_adhesion = 0.0;
        double hill_coefficient_for_adhesion = 8.0;
          // for consistency of feedback form:
        if (if_consider_feedback_of_face_values)
          assert(case_number_of_membrane_surface_energy_form != 0);
        else
          case_number_of_membrane_surface_energy_form = 0;

        if (feedback_strength_for_myosin_activity == 0.0)
          if_consider_feedback_of_face_values = false; // note: typically, we must have myosin feedback first.
        if (feedback_strength_for_adhesion == 0.0)
          if_consider_feedback_of_cell_cell_adhesion = false;
        if (is_default_feedback_form)
          assert(!if_consider_feedback_of_cell_cell_adhesion && !EMA_dont_decrease);

        p_face_value_and_stress_state_modifier->SetConsiderFeedbackOfFaceValues(if_consider_feedback_of_face_values);
        p_face_value_and_stress_state_modifier->SetConsiderFeedbackOfFaceValuesOnlyForBoundaryCells(if_apply_feedback_of_face_values_only_for_boundary_cells);
        p_face_value_and_stress_state_modifier->SetConsiderFeedbackOfFaceValuesOnlyForTopBoundaryCells(if_apply_feedback_of_face_values_only_for_top_boundary_cells);
        p_face_value_and_stress_state_modifier->SetConsiderFeedbackOfCellCellAdhesion(if_consider_feedback_of_cell_cell_adhesion);
        p_face_value_and_stress_state_modifier->SetEMADontDecrease_CCADontDecrease_HasAThreshold_Threshold(
            EMA_dont_decrease, CCA_dont_decrease, CCA_increasing_has_a_threshold_of_edge_length, CCA_increasing_threshold_of_edge_length_percentage);
        
        // Feedback information
        double edge_length_at_rest = sqrt(initial_area/(6*sqrt(3)/4)); // = 1.0996
        p_face_value_and_stress_state_modifier->SetUseFixedTargetArea(use_fixed_target_area_without_modifier); // used in the case where there is no target area modifier! (no division)
        p_face_value_and_stress_state_modifier->SetFixedTargetArea(real_reference_area);
        p_face_value_and_stress_state_modifier->SetEdgeLengthAtRest(edge_length_at_rest);
        p_face_value_and_stress_state_modifier->SetFeedbackStrengthForMyosinActivity(feedback_strength_for_myosin_activity);
        p_face_value_and_stress_state_modifier->SetHillCoefficientForMyosinActivity(hill_coefficient_for_myosin_activity);
        p_face_value_and_stress_state_modifier->SetFeedbackStrengthForAdhesion(feedback_strength_for_adhesion);
        p_face_value_and_stress_state_modifier->SetHillCoefficientForAdhesion(hill_coefficient_for_adhesion);

        // my stress state modifier
        p_face_value_and_stress_state_modifier->SetCalculateStressStateBoolean(true);
        bool if_set_cell_data_of_detailed_force_contributions = false;
        p_face_value_and_stress_state_modifier->SetIfSetCellDataOfEachForceContributions(if_set_cell_data_of_detailed_force_contributions);
        p_face_value_and_stress_state_modifier->SetCaseNumberOfMembraneSurfaceEnergyForm(case_number_of_membrane_surface_energy_form);        
        p_face_value_and_stress_state_modifier->SetNagaiHondaMembraneSurfaceEnergyParameter(Gamma);
        p_face_value_and_stress_state_modifier->SetNagaiHondaCellCellAdhesionEnergyParameter(cell_cell_adhesion_parameter);
        p_face_value_and_stress_state_modifier->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(cell_boundary_adhesion_parameter);        
        p_face_value_and_stress_state_modifier->SetFixedTargetPerimeter(target_shape_index*sqrt(reference_area));//?
        
        // my group number modifier
        p_face_value_and_stress_state_modifier->SetWriteGroupNumberToCell(use_my_detach_pattern_method);

        // my division state modifier
        p_face_value_and_stress_state_modifier->SetUseMyDivisionRuleAlongWithModifier(use_my_division_rule_along_with_modifier);
        p_face_value_and_stress_state_modifier->SetDivisionTime(time_for_one_division_of_cell_population);

        // new SSA distribution rule
        p_face_value_and_stress_state_modifier->SetMarkLeadingCells(mark_leading_cells);
        p_face_value_and_stress_state_modifier->SetLamellipodiumMaturationRate(lamellipodium_maturation_rate);
        p_face_value_and_stress_state_modifier->SetLamellipodiumDestructionRate(lamellipodium_destruction_rate);
        // outout information
        p_face_value_and_stress_state_modifier->SetOutputModifierInformationBoolean(false);

        // myosin activity depression:
        bool has_myo_depression = false;
        double myosin_activity_depressing_rate = 0.05/(M_PI/reference_area);
        p_face_value_and_stress_state_modifier->SetHasMyosinActivityDepression(has_myo_depression);
        p_face_value_and_stress_state_modifier->SetMyosinActivityDepressedTime(400.0*(M_PI/reference_area));
        p_face_value_and_stress_state_modifier->SetMyosinActivityDepressingRate(myosin_activity_depressing_rate);

        simulator.AddSimulationModifier(p_face_value_and_stress_state_modifier);
        /*-----------------------END: !!!!!!Feedback: FaceValueAndStressStateModifier: need modification---------------*/


        /*------------------------------------START: CellKiller---------------------------------------*/
        c_vector<double, 2> point_vec = zero_vector<double>(2);
        double kill_height = 10000.0;
        point_vec[1] = kill_height;
        c_vector<double, 2> norm_vec = zero_vector<double>(2);
        norm_vec[1] = 1.0;
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, point_vec, norm_vec) );
        p_killer->SetKillCellsGroupNumberFrom1(kill_cells_group_number_from1);
        simulator.AddCellKiller(p_killer);
        /*------------------------------------END: CellKiller---------------------------------------*/


        /*------------------------------------START: Timestep---------------------------------------*/
        //bool apply_my_change_to_make_timestep_adaptive = true;
        simulator.SetApplyMyChangesToMakeTimestepAdaptive(apply_my_change_to_make_timestep_adaptive);
        simulator.SetApplySamplingTimeInsteadOfSamplingTimestep(apply_my_change_to_make_timestep_adaptive);

        //double dt = 0.1*(M_PI/reference_area);
        simulator.SetDt(dt);

        //double sampling_time = 1.0*(M_PI/reference_area);        
        simulator.SetSamplingTime(sampling_time);
        double sampling_timestep_multiple = (unsigned) round(sampling_time/dt);
        simulator.SetSamplingTimestepMultiple(sampling_timestep_multiple);

        //double end_time = 400.0*(M_PI/reference_area);
        simulator.SetEndTime(end_time);
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
        // Output directory:
        std::ostringstream oss;
        time_t raw_time = time(0);
        struct tm * now = localtime(& raw_time);

        std::string output_directory = 
            "EpithelialBridgeSimulation/PHASE-DIAGRAM/Simulation Results Start From: ";
        oss.str("");
        oss << (now->tm_year + 1900 -2000) << '-' << (now->tm_mon + 1) << '-' <<  now->tm_mday << '/';
        output_directory += oss.str();

        oss.str("");
        oss << "NewEB-IntervalResist-";

        if (if_substrate_adhesion_is_homogeneous)
        {
          if (add_pulling_force_on_node_individually)
            oss << "Fy_lead_node=" << std::fixed << setprecision(1) << pulling_force_on_leading_cell;
          if (add_pulling_force_evenly_on_nodes_of_leading_cell)
            oss << "Fy=" << std::fixed << setprecision(1) << pulling_force_on_leading_cell;
        }

        oss << "_Fb=" << std::fixed << setprecision(4) << feedback_strength_for_myosin_activity
            << "_p0=" << std::fixed << setprecision(2) << target_shape_index
            << "_Ga=" << ((Gamma>=0.01 || Gamma==0.0)? std::fixed : std::scientific) 
                << setprecision(2) << Gamma;
        if (add_random_force && has_brownian_random_force)
          oss << "_Dtr=" << std::scientific << setprecision(2) << translational_diffusion_constant;          
        else
          oss << "_NoBrown";
        if (add_random_force && polarity_magnitude!=0.0)
          oss << "_Fp=" << ((polarity_magnitude>=0.01)? std::fixed : std::scientific) << setprecision(2) << polarity_magnitude
              << "_Dr=" << ((rotational_diffusion_constant>=0.01)? std::fixed : std::scientific) << setprecision(2) << rotational_diffusion_constant;
        else
          oss << "_NoPolarity";

        if(if_check_for_T4_swaps)
        {
          oss << "_T4swaps=1";
          oss << "_ratio=" << std::fixed << setprecision(1) << separation_ratio_with_rearr_thresh
              << "_typical_leng=" << std::fixed << setprecision(1) << typical_length_for_T4_swaps;
        }
        if (if_consider_substrate_adhesion)
        {
          if (if_substrate_adhesion_is_homogeneous)
            oss << "_HomoSSA=" << std::fixed << setprecision(1) << homogeneous_substrate_adhesion_parameter;
          else
          {
            oss << "_NonHomoSSA:mature_lame=" << std::fixed << setprecision(1) << SSA_for_mature_lamellipodium
              << "_basic=" << std::fixed << setprecision(1) << basic_SSA;
            if (use_new_SSA_distribution_rule)
              oss << "_NewSSARule:matu_rate=" << std::fixed << setprecision(2) << lamellipodium_maturation_rate
                  << "_destru_rate=" << std::fixed << setprecision(2) << lamellipodium_destruction_rate;
            else if (use_my_detach_pattern_method)
              oss << "_Multiple_lead_tops=1";
          }
        }
        
        if (use_longer_mesh)
          oss << "_LongMesh";
        if (larger_strip_spacing_distance)
          oss << "_StripDis=" << std::fixed << setprecision(3) << strip_distance;
        if (!run_with_birth)
          oss << "_NoDivi";
        if (move_mesh_right_for_N_periods!=0)
          oss << "_MvRight=" << std::fixed << setprecision(0) << move_mesh_right_for_N_periods;
        if (has_myo_depression)
          oss << "_MyoDeprRate=" << std::scientific << setprecision(1)<< myosin_activity_depressing_rate;
        if (!if_apply_feedback_of_face_values_only_for_boundary_cells)
          oss << "_FbInside=1";
        output_directory += oss.str();

        oss.str("");
        oss << (now->tm_year + 1900 -2000) << '-' << (now->tm_mon + 1) << '-' <<  now->tm_mday 
            << ", " << now->tm_hour << ':' << now->tm_min << ':' << now->tm_sec;
        output_directory += "_Timestamp=" + oss.str();

        std::string concise_output_directory = output_directory;
        // Concise information has been written to directory.

        oss.str("");
        oss << std::fixed << setprecision(3) << dt;
        output_directory += "/Dt=" + oss.str();
        oss.str("");
        oss << ((cell_rearrangement_threshold>=0.01)? std::fixed : std::scientific) << setprecision(2) << cell_rearrangement_threshold;
        output_directory += "_RearThr=" + oss.str();
        if (!restrict_vertex_movement)
          output_directory += "_ResMv=0";
        if (apply_my_change_to_make_timestep_adaptive)
        {
          oss.str("");
          oss << ((max_movement_per_timestep_using_my_adaptive_timestep_method>=0.01)? std::fixed : std::scientific) << setprecision(2) << max_movement_per_timestep_using_my_adaptive_timestep_method;
          output_directory += "_AdaptDt:MaxMvDt=" + oss.str();
        }
        if (if_consider_substrate_adhesion)
        {
          oss.str("");
          oss << ((small_change_for_area_calculation>=0.01)? std::fixed : std::scientific) << setprecision(2) << small_change_for_area_calculation;
          output_directory += "_Dxy=" + oss.str();
        }

        output_directory += "_|Divi=" + std::to_string(run_with_birth);
        if (if_consider_feedback_of_face_values && is_default_feedback_form)
          output_directory += "_HasDefaultFeedb";
        else
          output_directory += "_HasFeedb=" + std::to_string(if_consider_feedback_of_face_values);
        output_directory += "_HasRandF=" + std::to_string(add_random_force);
        output_directory += "_MSE=" + std::to_string(case_number_of_membrane_surface_energy_form);
        
        oss.str("");
        oss << ((Gamma>=0.01 || Gamma==0.0)? std::fixed : std::scientific) 
                << setprecision(2) << Gamma;
        output_directory += "_|Ga=" + oss.str();
        oss.str("");
        oss << std::fixed << setprecision(2) << target_shape_index;
        output_directory += "_p0=" + oss.str();
        oss.str("");
        oss << ((fabs(cell_cell_adhesion_parameter)>=0.01)? std::fixed : std::scientific) << setprecision(2) << cell_cell_adhesion_parameter;
        output_directory += "_CCAdhe=" + oss.str();
        oss.str("");
        oss << ((fabs(cell_boundary_adhesion_parameter)>=0.01)? std::fixed : std::scientific) << setprecision(2) << cell_boundary_adhesion_parameter;
        output_directory += "_CBAdhe=" + oss.str();
        
        output_directory += "_|HasSubAdh=" + std::to_string(if_consider_substrate_adhesion);
        if (if_consider_substrate_adhesion)
        {
          if (if_substrate_adhesion_is_homogeneous)
            output_directory += "_SSAIsHomo";
          else
            output_directory += "_NonHomoSSA";       
        }
        // Detailed information has been written to directory.

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
          oss << ((feedback_strength_for_myosin_activity<=100.0)? std::fixed : std::scientific) << setprecision(4) << feedback_strength_for_myosin_activity ;
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
          if (has_brownian_random_force)
          {
            oss.str("");
            oss << std::scientific << setprecision(2) << translational_diffusion_constant;
            output_directory += "_|Dtr=" + oss.str();
          }
          if (has_polarity)
          {
            oss.str("");
            oss << ((polarity_magnitude>=0.01 || polarity_magnitude==0.0)? std::fixed : std::scientific) << setprecision(2) << polarity_magnitude;
            output_directory += "_fp=" + oss.str();
            oss.str("");
            oss << ((rotational_diffusion_constant>=0.01)? std::fixed : std::scientific) << setprecision(2) << rotational_diffusion_constant;
            output_directory += "_Dr=" + oss.str();
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
            oss << std::fixed << setprecision(1) << SSA_for_mature_lamellipodium;
            output_directory += "_mature_lamelli=" + oss.str();
            oss.str("");
            oss << std::fixed << setprecision(1) << basic_SSA;
            output_directory += "_basic=" + oss.str();
          }
        }

        // tmp
        if (small_SSA_at_first)
        {
          oss.str("");
          oss << "Small_SSA_mature_lame=" << std::fixed << setprecision(1) << small_SSA_for_initial_time 
              << "_for_initial_time=" << std::fixed << setprecision(1) << initial_time_for_small_SSA; 
          output_directory += oss.str();
        }
        // More detailed information has been written to directory.

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
          std::cout << std::endl << "Detailed output directory is set: " << output_directory << std::endl;
        }
        /*--------------------------END: Output Directory and Simulation Information File---------------------*/

        simulator.Solve();
    }

};

#endif /* TESTMYNEWEB_HPP_ */
