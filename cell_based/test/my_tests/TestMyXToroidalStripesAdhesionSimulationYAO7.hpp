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

#ifndef TESTMYXTOROIDALSTRIPESADHESIONSIMULATIONYAO7_HPP_
#define TESTMYXTOROIDALSTRIPESADHESIONSIMULATIONYAO7_HPP_

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

class TestMyXToroidalStripesAdhesionSimulationYAO7 : public AbstractCellBasedTestSuite
{
public:

    void TestOscillation()
    {
        /*-----------------------START: Frequently changed parameters-------------------------*/

        /* Temporarily changed important parameters:
         * 
         * dt=0.1; my_adaptive_dt=1; consider_consistency_for_SSA = 1; max_movement_per_timestep = 0.05;
         * set_cell_rearrangement_threshold = 0.05;
         * small_change_for_area_calculation = 0.2; // 0.15 default
         * 
         * considerSubAdh=1; StripHomo=0 (default); 0)(default)SSATop=-3; SSABottom=-0.5; SSALeadTopLeng=3.0; 1)SSA=-3;
         * considerResSubAdh=1; RSA=-0.5; IgnoreRSATop=0;
         * 
         * cell_division=1; time_one_division=25;
         * Ga=0.1;
         * 
         * use_my_detach_pattern_method = 1;
         * 
         * feeback=1; 1,8; 0,8; form=(0),0,1,1
         * p0=4.0;
         * HasRandomForce=1; fp=0.2; Dr=0.01;
         * 
        */

        // phase diagram:
        bool if_consider_feedback_of_face_values = false; // !feedback!
        double set_feedback_strength_for_myosin_activity = 1.0;
        if (if_consider_feedback_of_face_values == false)
          set_feedback_strength_for_myosin_activity = 0.0;
        double set_target_shape_index = 1.0; // {6/sqrt(6*sqrt(3)/4)}=3.72
        double set_polarity_magnitude = 0.1;

       
        // output directory
        std::ostringstream oss;
        std::string out_put_directory = "EpithelialBridgeSimulation/LookForMovementPatterns_StartFromJuly20/";
        time_t raw_time = time(0);   // get time now
        struct tm * now = localtime(& raw_time);
        oss.str("");
        oss << (now->tm_year + 1900 -2000) << '-' << (now->tm_mon + 1) << '-' <<  now->tm_mday 
            << ' ' << now->tm_hour << ':' << now->tm_min << ':' << now->tm_sec;
        out_put_directory += "__YAO7__StartTime=" + oss.str();
        out_put_directory += "__PHASE_DIAGRAM_";
        oss.str("");
        oss << std::scientific << setprecision(1) << set_feedback_strength_for_myosin_activity;
        out_put_directory += "_MyoFeStr=" + oss.str();
        oss.str("");
        oss << std::fixed << setprecision(2) << set_target_shape_index;
        out_put_directory += "_p0=" + oss.str();
        oss.str("");
        oss << std::scientific << setprecision(2) << set_polarity_magnitude;
        out_put_directory += "_fp=" + oss.str();
        
        // timestep
        double set_dt = 0.1; // 1.0/round(1.0/(0.01*120.0/50.0));
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
        double set_nagai_honda_membrane_surface_energy_parameter = 0.1;
        bool set_use_fixed_target_area = false;
        // double set_target_shape_index = 4.0; // {6/sqrt(6*sqrt(3)/4)}=3.72
          // substrate adhesion
        bool if_consider_substrate_adhesion = true;
            // strip substrate adhesion
        bool if_substrate_adhesion_is_homogeneous = false;
              // if homogeneous
        double set_homogeneous_substrate_adhesion_parameter = -3.0;//
              // if not homogeneous
        double set_substrate_adhesion_leading_top_length = 3.0;
        double set_substrate_adhesion_parameter_at_leading_top= -3.0;
        double set_substrate_adhesion_parameter_below_leading_top = -0.5;
                // detach pattern:
        bool use_my_detach_pattern_method = true;
            // reservoir substrate adhesion
        bool if_consider_reservoir_substrate_adhesion = true;// true for default
        bool if_ignore_reservoir_substrate_adhesion_at_top = false;// false for default
        bool if_ignore_reservoir_substrate_adhesion_at_bottom = false;// false for default
        double set_reservoir_substrate_adhesion_parameter = -0.5;//
        
        // feedback
        // bool if_consider_feedback_of_face_values = true; // !feedback!
        bool if_apply_feedback_of_face_values_only_for_boundary_cells = true; // for testing fluid inside
        bool if_apply_feedback_of_face_values_only_for_top_boundary_cells = true;
        bool if_update_unified_cell_cell_adhesion_of_face = false;//false for default feedback without AdhesionFeedback!
        int case_number_of_membrane_surface_energy_form = 2; // 2 for new membr_surf_energy form1, 3 for form2.
        if (if_consider_feedback_of_face_values)
          assert(case_number_of_membrane_surface_energy_form != 0);
        else
          case_number_of_membrane_surface_energy_form = 0;
        bool set_EMA_dont_decrease = false;// false for default feedback without AdhesionFeedback!
        bool set_CCA_dont_increase = true;// true for default feedback without AdhesionFeedback!
        bool set_CCA_dont_decrease = true;// true for default feedback without AdhesionFeedback!
        bool set_CCA_dont_inrease_until_shorter_than_a_threshold = true;
        double set_CCA_dont_increase_until_shorter_than_this_length = 0.2;
        // double set_feedback_strength_for_myosin_activity = 1;
        double set_hill_coefficient_for_myosin_activity = 8.0;
        double set_feedback_strength_for_adhesion = 0;
        double set_hill_coefficient_for_adhesion = 8.0;
        
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
        // CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
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
        std::vector< boost::shared_ptr<AbstractForce<2, 2>> > force_collection = simulator.rGetForceCollection();
        //Note: Here the address of std::vector of force collections are different from that in simulator!
        //But the corresponding shared_ptrs in the vectors should be the same.
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
          double Dt = set_dt;
          double angle_for_initialization = M_PI;// 0 for all upwards; PI for all directions
          p_polarity_modifier->SetPolarityMagnitude(polarity_magnitude);
          p_polarity_modifier->SetD(rotational_diffusion_constant);
          // tmp: delete this method later!
          p_polarity_modifier->SetDt(Dt);
          p_polarity_modifier->SetAngleForInitialization(angle_for_initialization);
          simulator.AddSimulationModifier(p_polarity_modifier);
        }
        /*-----------------------------END: RandomForce and PolarityModifier---------------------------*/


        /*------------------------START: !!!!!Feedback: FaceValueAndStressStateModifier: need modification---------------*/
        MAKE_PTR_ARGS(FaceValueAndStressStateModifier<2>, p_face_value_and_stress_state_modifier, ());
        
        // Energy function parameter
        p_face_value_and_stress_state_modifier->SetNagaiHondaMembraneSurfaceEnergyParameter(nagai_honda_membrane_surface_energy_parameter);
        p_face_value_and_stress_state_modifier->SetFixedTargetPerimeter(target_shape_index*sqrt(fixed_target_area));
        
        // Feedback information
        double edge_length_at_rest = sqrt(initial_area/(6*sqrt(3)/4)); // = 1.0996
        double feedback_strength_for_myosin_activity = set_feedback_strength_for_myosin_activity;
        double hill_coefficient_for_myosin_activity = set_hill_coefficient_for_myosin_activity;
        p_face_value_and_stress_state_modifier->SetEdgeLengthAtRest(edge_length_at_rest);
        p_face_value_and_stress_state_modifier->SetFeedbackStrengthForMyosinActivity(feedback_strength_for_myosin_activity);
        p_face_value_and_stress_state_modifier->SetHillCoefficientForMyosinActivity(hill_coefficient_for_myosin_activity);

        double feedback_strength_for_adhesion = set_feedback_strength_for_adhesion;
        double hill_coefficient_for_adhesion = set_hill_coefficient_for_adhesion;
        p_face_value_and_stress_state_modifier->SetFeedbackStrengthForAdhesion(feedback_strength_for_adhesion);
        p_face_value_and_stress_state_modifier->SetHillCoefficientForAdhesion(hill_coefficient_for_adhesion);

        bool EMA_dont_decrease = set_EMA_dont_decrease;
        bool CCA_dont_increase = set_CCA_dont_increase;
        bool CCA_dont_inrease_until_shorter_than_a_threshold = set_CCA_dont_inrease_until_shorter_than_a_threshold;
        bool CCA_dont_decrease = set_CCA_dont_decrease;
        double CCA_dont_increase_until_shorter_than_this_value = set_CCA_dont_increase_until_shorter_than_this_length;

        p_face_value_and_stress_state_modifier->SetConsiderFeedbackOfFaceValues(if_consider_feedback_of_face_values);
        p_face_value_and_stress_state_modifier->SetConsiderFeedbackOfFaceValuesOnlyForBoundaryCells(if_apply_feedback_of_face_values_only_for_boundary_cells);
        p_face_value_and_stress_state_modifier->SetConsiderFeedbackOfFaceValuesOnlyForTopBoundaryCells(if_apply_feedback_of_face_values_only_for_top_boundary_cells);
        p_face_value_and_stress_state_modifier->SetUpdateUnifiedCellCellAdhesionOfFace(if_update_unified_cell_cell_adhesion_of_face);
        p_face_value_and_stress_state_modifier->SetEMADontDecreaseWhenEdgeShrink_CCADontDecreaseWhenEdgeExpand_CCADontInreaseWhenEdgeShrink(EMA_dont_decrease,CCA_dont_decrease,CCA_dont_increase);
        p_face_value_and_stress_state_modifier->SetCCADontInreaseUntilShorterThanAThreshold(CCA_dont_inrease_until_shorter_than_a_threshold);
        p_face_value_and_stress_state_modifier->SetCCADontInreaseUntilShorterThanThisValue(CCA_dont_increase_until_shorter_than_this_value);
        
        p_face_value_and_stress_state_modifier->SetCalculateStressStateBoolean(true);
        p_face_value_and_stress_state_modifier->SetCaseNumberOfMembraneSurfaceEnergyForm(case_number_of_membrane_surface_energy_form);        
        p_face_value_and_stress_state_modifier->SetOutputModifierInformationBoolean(false);

        p_face_value_and_stress_state_modifier->SetWriteGroupNumberToCell(use_my_detach_pattern_method);

        p_face_value_and_stress_state_modifier->SetUseMyDivisionRuleAlongWithModifier(use_my_division_rule_along_with_modifier);
        p_face_value_and_stress_state_modifier->SetDivisionTime(time_for_one_division_of_cell_population);
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
        

        /*-------------------------------START: Output Directory---------------------------------*/
        // concise information
        oss.str("");
        oss << std::fixed << setprecision(3) << dt;
        out_put_directory += "/Dt=" + oss.str();
        out_put_directory += "_MyAdaptDt=" + std::to_string(apply_my_change_to_make_timestep_adaptive);        
        out_put_directory += "_ResMv=" + std::to_string(restrict_vertex_movement);
        out_put_directory += "_ConsistMv=" + std::to_string(consider_consistency_for_SSA);
        oss.str("");
        oss << std::scientific << setprecision(1) << max_movement_per_timestep;
        out_put_directory += "_MaxMvDt=" + oss.str();
        if (if_consider_substrate_adhesion)
        {
          oss.str("");
          oss << std::scientific << setprecision(1) << small_change_for_area_calculation;
          out_put_directory += "_Dxy=" + oss.str();
        }
        oss.str("");
        oss << std::scientific << setprecision(1) << cell_rearrangement_threshold;
        out_put_directory += "_RearThr=" + oss.str();
        out_put_directory += "_|_Divi=" + std::to_string(run_with_birth);
        oss.str("");
        oss << std::fixed << setprecision(3) << nagai_honda_membrane_surface_energy_parameter;
        out_put_directory += "_|_Ga=" + oss.str();
        oss.str("");
        oss << std::fixed << setprecision(2) << target_shape_index;
        out_put_directory += "_p0=" + oss.str();
        oss.str("");
        oss << std::scientific << setprecision(1) << nagai_honda_cell_cell_adhesion_energy_parameter;
        out_put_directory += "_CCAdhe=" + oss.str();
        oss.str("");
        oss << std::scientific << setprecision(1) << nagai_honda_cell_boundary_adhesion_energy_parameter;
        out_put_directory += "_CBAdhe=" + oss.str();
        out_put_directory += "_|_HasSA=" + std::to_string(if_consider_substrate_adhesion);
        if (if_consider_substrate_adhesion)
          out_put_directory += "_HomoSSA=" + std::to_string(if_substrate_adhesion_is_homogeneous);
        out_put_directory += "_|_HasRSA=" + std::to_string(if_consider_reservoir_substrate_adhesion);
        out_put_directory += "_|_MSE=" + std::to_string(case_number_of_membrane_surface_energy_form);
        out_put_directory += "_|_AddFeedb=" + std::to_string(if_consider_feedback_of_face_values);
        out_put_directory += "_|_AddRandF=" + std::to_string(add_random_force);
        
        // detailed information:
        out_put_directory += "/||";
        if (run_with_birth)
        {
          oss.str("");
          oss << std::fixed << setprecision(1) << time_for_one_division_of_cell_population;
          out_put_directory += "_|TDivi=" + oss.str();
          oss.str("");
          oss << std::fixed << setprecision(2) << growth_rate_for_target_area_after_division;
          out_put_directory += "_GrRate=" + oss.str();
          if (use_my_division_rule_along_with_modifier)
            out_put_directory += "_DiviDtFixed=1";
        }
        if(if_consider_substrate_adhesion)
        {
          if (if_substrate_adhesion_is_homogeneous)
          {
            oss.str("");
            oss << std::fixed << setprecision(1) << homogeneous_substrate_adhesion_parameter;
            out_put_directory += "_|HomoSSA=" + oss.str();
          }
          else
          {
            oss.str("");
            oss << std::fixed << setprecision(1) << substrate_adhesion_leading_top_length;
            out_put_directory += "_|SSALeadLeng=" + oss.str();
            oss.str("");
            oss << std::fixed << setprecision(1) << substrate_adhesion_parameter_at_leading_top;
            out_put_directory += "_SSATop=" + oss.str();
            oss.str("");
            oss << std::fixed << setprecision(1) << substrate_adhesion_parameter_below_leading_top;
            out_put_directory += "_SSABelow=" + oss.str();
          }
        }
        if(if_consider_reservoir_substrate_adhesion)
        {
          // out_put_directory += "_IgnRSABottom=" + std::to_string(if_ignore_reservoir_substrate_adhesion_at_bottom);
          oss.str("");
          oss << std::fixed << setprecision(1) << reservoir_substrate_adhesion_parameter;
          out_put_directory += "_|_RSA=" + oss.str();
          out_put_directory += "_ConsiRSATop=" + std::to_string(!if_ignore_reservoir_substrate_adhesion_at_top);
        }

        if (if_consider_feedback_of_face_values)
        {
          oss.str("");
          oss << std::scientific << setprecision(1) << feedback_strength_for_myosin_activity ;
          out_put_directory += "_|MyoFeStr=" + oss.str();
          oss.str("");
          oss << std::fixed << setprecision(0) << hill_coefficient_for_myosin_activity ;
          out_put_directory += "_MyoHill=" + oss.str();
          if (if_update_unified_cell_cell_adhesion_of_face)
          {
            oss.str("");
            oss << std::scientific << setprecision(1) << feedback_strength_for_adhesion;
            out_put_directory += "_AdhFeStr=" + oss.str();
            oss.str("");
            oss << std::fixed << setprecision(0) << hill_coefficient_for_adhesion;
            out_put_directory += "_AdhHill=" + oss.str();
          }
          
          if (EMA_dont_decrease==false && if_update_unified_cell_cell_adhesion_of_face==false)
            out_put_directory += "_DefaultFe";
          else
          {
            out_put_directory += "_EMACanIn=1";
            out_put_directory += "_EMACanDe=" + std::to_string(!EMA_dont_decrease);
            out_put_directory += "_CCACanIn=" + std::to_string(!CCA_dont_increase);
            out_put_directory += "_CCACanDe=" + std::to_string(!CCA_dont_decrease);
            if (CCA_dont_inrease_until_shorter_than_a_threshold)  
            {
              oss.str("");
              oss << std::fixed << setprecision(1) << CCA_dont_increase_until_shorter_than_this_value;          
              out_put_directory += "_CCAInrThresh=" + oss.str();
            }
          }
        }
        if (add_random_force)
        {
          oss.str("");
          oss << std::scientific << setprecision(2) << translational_diffusion_constant;
          out_put_directory += "_|D=" + oss.str();
          if (consider_polarity)
          {
            oss.str("");
            oss << std::scientific << setprecision(2) << polarity_magnitude;
            out_put_directory += "_fp=" + oss.str();
            oss.str("");
            oss << std::scientific << setprecision(2) << rotational_diffusion_constant;
            out_put_directory += "_Dr=" + oss.str();
            out_put_directory += "_MotiInterval=" + std::to_string(!vanishing_motility_for_node_in_the_strip_interval);
          }
        }

        simulator.SetOutputDirectory(out_put_directory);
        std::cout << std::endl << "OutputDirectory: " << out_put_directory << std::endl;
        /*-------------------------------END: Output Directory---------------------------------*/

        simulator.Solve();
    }

};

#endif /* TESTMYXTOROIDALSTRIPESADHESIONSIMULATIONYAO7_HPP_ */
