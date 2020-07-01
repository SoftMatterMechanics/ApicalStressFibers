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

#ifndef TESTMYXTOROIDALSTRIPESADHESIONSIMULATIONRAW_HPP_
#define TESTMYXTOROIDALSTRIPESADHESIONSIMULATIONRAW_HPP_

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

class TestMyXToroidalStripesAdhesionSimulationRAW : public AbstractCellBasedTestSuite
{
public:

    void TestOscillation()
    {
        /*-----------------------START: Frequently changed parameters-------------------------*/
        // output directory
        // temporarily changed important parameters:
        // feeback=1, restrict_vertex_movement=1, use_adaptive_timestep=1, dt=0.025, set_cell_rearrangement_threshold=0.01;
        std::ostringstream oss;
        std::string out_put_directory = "EpithelialBridgeSimulation/LookForMovementPatterns_StartFromJuly01/";
        time_t raw_time = time(0);   // get time now
        struct tm * now = localtime(& raw_time);
        oss.str("");
        oss << (now->tm_year + 1900 -2000) << '-' << (now->tm_mon + 1) << '-' <<  now->tm_mday 
            << ' ' << now->tm_hour << ':' << now->tm_min << ':' << now->tm_sec;
        out_put_directory += "__RAW__SimulationTime=" + oss.str();
        
        // timestep
        double set_dt = 0.025;
        bool restrict_vertex_movement = true; //???????
        bool use_adaptive_timestep = true; //?????????

        // cell rearrangement
        double set_cell_rearrangement_threshold = 0.01;

        // cell division
        bool run_with_birth =true;
        bool use_bernoulli_trial_cell_cycle_model = true;
        double time_for_one_division_of_cell_population = 10;
        double growth_rate_for_target_area_after_division = 0.1;

        // structure
        bool if_use_larger_strip_distance = false;
        bool if_mesh_has_two_period = false;
        bool use_longer_mesh = false;
        bool if_consider_compression = false;
        double compress_ratio = 1.2;
        
        // energy parameter
        double set_nagai_honda_membrane_surface_energy_parameter = 0.1;
        bool set_use_fixed_target_area = false;
        double set_target_shape_index = 4.0; // {6/sqrt(6*sqrt(3)/4)}=3.72
          // substrate adhesion
        bool if_consider_substrate_adhesion = true;
        bool set_use_fine_mesh_for_calculating_substrate_adhesion = false;
            // strip substrate adhesion
        bool if_substrate_adhesion_is_homogeneous = false;
              // homogeneous
        double set_homogeneous_substrate_adhesion_parameter = -75;
              // not homogeneous
        double set_substrate_adhesion_leading_top_length = 1.0;
        double set_substrate_adhesion_parameter_at_leading_top= -50.0;
        double set_substrate_adhesion_parameter = -10;
            // reservoir substrate adhesion
        bool if_consider_reservoir_substrate_adhesion = true;// true for default
        bool if_ignore_reservoir_substrate_adhesion_at_bottom = false;// false for default
        double set_reservoir_substrate_adhesion_parameter = -10.0;
          //tmp
        bool if_set_strip_start_y_location_artificially = false;//false for default
        double set_strip_start_y_location_artificially = 1.5*8*1/sqrt(3)*sqrt(M_PI/(sqrt(3)/2));
        
        // feedback
        bool if_consider_feedback_of_face_values = true;
        bool if_apply_feedback_of_face_values_only_for_boundary_cells = true; // for testing fluid inside
        bool if_apply_feedback_of_face_values_only_for_top_boundary_cells = true;
        int case_number_of_membrane_surface_energy_form = 2; // 2 for new membr_surf_energy form1, 3 for form2.
        if (if_consider_feedback_of_face_values)
          assert(case_number_of_membrane_surface_energy_form != 0);
        else
          case_number_of_membrane_surface_energy_form = 0;
        bool set_EMA_dont_decrease = false;// false for default!
        bool set_CCA_dont_decrease = true;// true for default!
        bool set_CCA_dont_increase = true;// true for default!
        bool set_CCA_dont_inrease_until_shorter_than_a_threshold = true;
        double set_CCA_dont_increase_until_shorter_than_this_length = 0.2;
        double set_feedback_strength_for_myosin_activity = 1.0;
        double set_feedback_strength_for_adhesion = 0.0;
        
        // polarity
        bool add_random_force = true;
        bool consider_polarity = true;
        double set_polarity_magnitude = 0.2;
        double set_rotational_diffusion_constant = 0.01;

        // output
        bool output_concise_swap_information_when_remesh = true;
        bool output_detailed_swap_information_when_remesh = false; //suggest "false" for concise output results
        bool output_information_for_nagai_honda_force = false;
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
        if(if_consider_compression)
          num_ele_up *= 1;

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
        cells_generator.SetUseBernoulliTrialCellCycleModel(use_bernoulli_trial_cell_cycle_model);
        double divide_probability_for_a_cell = 1.0/time_for_one_division_of_cell_population/(num_ele_cross*num_ele_up);
        double minimum_division_age = -0.01;
        cells_generator.SetDivideProbability(divide_probability_for_a_cell);
        cells_generator.SetMinimumDivisionAge(minimum_division_age);
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);
        
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetRestrictVertexMovementBoolean(restrict_vertex_movement);//Attention!!!!!!!!!!!!
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
        p_numerical_method->SetUseAdaptiveTimestep(use_adaptive_timestep);
        simulator.SetNumericalMethod(p_numerical_method);
        /*---------------------------------END: Add Numerical Method-----------------------------*/


        /*-----------------------------START: MyNagaiHondaForceWithStripesAdhesion---------------------*/
        MAKE_PTR(MyNagaiHondaForceWithStripesAdhesion<2>, p_force);
        
        // KA, A0
        double nagai_honda_deformation_energy_parameter = 1;
        bool use_fixed_target_area = set_use_fixed_target_area;
        double fixed_target_area = initial_area;
        if (if_consider_compression)
          fixed_target_area *= compress_ratio;

        // KP 
        double nagai_honda_membrane_surface_energy_parameter = set_nagai_honda_membrane_surface_energy_parameter;
                
        // LAMBDA: Use equation form: 0.5*Gamma*P^2, target perimeter being 0!
        double target_shape_index = set_target_shape_index;
        double nagai_honda_cell_cell_adhesion_energy_parameter = -nagai_honda_membrane_surface_energy_parameter*(target_shape_index*sqrt(reference_target_area_for_modifier));
        if (use_fixed_target_area)
          nagai_honda_cell_cell_adhesion_energy_parameter = -nagai_honda_membrane_surface_energy_parameter*(target_shape_index*sqrt(fixed_target_area));
        double nagai_honda_cell_boundary_adhesion_energy_parameter = 0; //nagai_honda_cell_cell_adhesion_energy_parameter;
        bool if_use_face_element_to_get_adhesion_parameter = if_consider_feedback_of_face_values;

        // Substrate Adhesion Parameter
        double substrate_adhesion_leading_top_length = set_substrate_adhesion_leading_top_length;
        double substrate_adhesion_parameter_at_leading_top= set_substrate_adhesion_parameter_at_leading_top;
        double substrate_adhesion_parameter = set_substrate_adhesion_parameter;
        if (if_substrate_adhesion_is_homogeneous)
          substrate_adhesion_parameter = set_homogeneous_substrate_adhesion_parameter;
        // reservoir substrate adhesion
        double reservoir_substrate_adhesion_parameter = set_reservoir_substrate_adhesion_parameter;
        
        // Strips structure of substrate adhesion
        bool use_fine_mesh_for_calculating_substrate_adhesion = set_use_fine_mesh_for_calculating_substrate_adhesion;
        double strip_width = sqrt(initial_area/(sqrt(3)/2))/2; // =0.952~1.05
        double strip_distance = 6*sqrt(initial_area/(sqrt(3)/2)); // =11.428~12.60
        if(if_use_larger_strip_distance)
          strip_distance *= 200.0/120.0;
        double strip_start_x_location = 0.0;
        if (if_mesh_has_two_period)
          strip_start_x_location = -3*sqrt(initial_area/(sqrt(3)/2));
        double strip_start_y_location = 1.5*num_ele_up*1/sqrt(3)*sqrt(initial_area/(sqrt(3)/2));
        if (if_consider_compression)
        {
          strip_start_y_location += 0.5*1/sqrt(3)*sqrt(initial_area/(sqrt(3)/2));
          strip_start_y_location += 0.5*(fixed_target_area/initial_area-1)*strip_start_y_location;
        }
        //tmp
        if (if_set_strip_start_y_location_artificially)
          strip_start_y_location = set_strip_start_y_location_artificially;

        p_force->SetCaseNumberOfMembraneSurfaceEnergyForm(case_number_of_membrane_surface_energy_form);
        p_force->SetUseFaceElementToGetAdhesionParameterBoolean(if_use_face_element_to_get_adhesion_parameter);

        p_force->SetNagaiHondaDeformationEnergyParameter(nagai_honda_deformation_energy_parameter);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(nagai_honda_membrane_surface_energy_parameter);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(nagai_honda_cell_cell_adhesion_energy_parameter);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(nagai_honda_cell_boundary_adhesion_energy_parameter);        

        p_force->SetUseFixedTargetArea(use_fixed_target_area);
        p_force->SetFixedTargetArea(fixed_target_area);
        p_force->SetTargetShapeIndex(target_shape_index);//?????????????????????????????
          // substrate adhesion info
        p_force->SetIfConsiderSubstrateAdhesion(if_consider_substrate_adhesion);
        p_force->SetUseFineMesh(use_fine_mesh_for_calculating_substrate_adhesion);
        p_force->SetIfSubstrateAdhesionIsHomogeneous(if_substrate_adhesion_is_homogeneous);
        p_force->SetSubstrateAdhesionLeadingTopLength(substrate_adhesion_leading_top_length);
        p_force->SetSubstrateAdhesionParameterAtLeadingTop(substrate_adhesion_parameter_at_leading_top);
        p_force->SetSubstrateAdhesionParameter(substrate_adhesion_parameter);
        p_force->SetIfConsiderReservoirSubstrateAdhesion(if_consider_reservoir_substrate_adhesion);
        p_force->SetIfIgnoreReservoirSubstrateAdhesionAtBottom(if_ignore_reservoir_substrate_adhesion_at_bottom);
        p_force->SetCenterOfWidth(0.0);
        p_force->SetWidth(num_ele_cross*sqrt(initial_area/(sqrt(3)/2)));//check later!
        p_force->SetReservoirSubstrateAdhesionParameter(reservoir_substrate_adhesion_parameter);
          // strip info
        p_force->SetStripWidth(strip_width);
        p_force->SetStripDistance(strip_distance);
        p_force->SetStripStartXLocation(strip_start_x_location);
        p_force->SetStripStartYLocation(strip_start_y_location);

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
          p_force2->SetConsiderPolarity(consider_polarity);
          p_force2->SetUseTheSameNodeRadius(use_the_same_node_radius);
          p_force2->SetTheSameNodeRadius(node_radius);
          translational_diffusion_constant = p_force2->GetDiffusionScalingConstant()/node_radius;
          simulator.AddForce(p_force2);
        }

        // Polarity force
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
        double hill_coefficient_for_myosin_activity = 2;
        p_face_value_and_stress_state_modifier->SetEdgeLengthAtRest(edge_length_at_rest);
        p_face_value_and_stress_state_modifier->SetFeedbackStrengthForMyosinActivity(feedback_strength_for_myosin_activity);
        p_face_value_and_stress_state_modifier->SetHillCoefficientForMyosinActivity(hill_coefficient_for_myosin_activity);

        double feedback_strength_for_adhesion = set_feedback_strength_for_adhesion;
        double hill_coefficient_for_adhesion = 2;
        p_face_value_and_stress_state_modifier->SetFeedbackStrengthForAdhesion(feedback_strength_for_adhesion);
        p_face_value_and_stress_state_modifier->SetHillCoefficientForAdhesion(hill_coefficient_for_adhesion);

        bool EMA_dont_decrease = set_EMA_dont_decrease;
        bool CCA_dont_decrease = set_CCA_dont_decrease;
        bool CCA_dont_increase = set_CCA_dont_increase;
        bool CCA_dont_inrease_until_shorter_than_a_threshold = set_CCA_dont_inrease_until_shorter_than_a_threshold;
        double CCA_dont_increase_until_shorter_than_this_value = set_CCA_dont_increase_until_shorter_than_this_length;

        p_face_value_and_stress_state_modifier->SetConsiderFeedbackOfFaceValues(if_consider_feedback_of_face_values);
        p_face_value_and_stress_state_modifier->SetConsiderFeedbackOfFaceValuesOnlyForBoundaryCells(if_apply_feedback_of_face_values_only_for_boundary_cells);
        p_face_value_and_stress_state_modifier->SetConsiderFeedbackOfFaceValuesOnlyForTopBoundaryCells(if_apply_feedback_of_face_values_only_for_top_boundary_cells);
        p_face_value_and_stress_state_modifier->SetEMADontDecreaseWhenEdgeShrink_CCADontDecreaseWhenEdgeExpand_CCADontInreaseWhenEdgeShrink(EMA_dont_decrease,CCA_dont_decrease,CCA_dont_increase);
        p_face_value_and_stress_state_modifier->SetCCADontInreaseUntilShorterThanAThreshold(CCA_dont_inrease_until_shorter_than_a_threshold);
        p_face_value_and_stress_state_modifier->SetCCADontInreaseUntilShorterThanThisValue(CCA_dont_increase_until_shorter_than_this_value);
        
        p_face_value_and_stress_state_modifier->SetCalculateStressStateBoolean(true);
        p_face_value_and_stress_state_modifier->SetCaseNumberOfMembraneSurfaceEnergyForm(case_number_of_membrane_surface_energy_form);        
        p_face_value_and_stress_state_modifier->SetOutputModifierInformationBoolean(false);
        simulator.AddSimulationModifier(p_face_value_and_stress_state_modifier);
        /*-----------------------END: !!!!!!Feedback: FaceValueAndStressStateModifier: need modification---------------*/


        /*------------------------------------START: Timestep---------------------------------------*/
        double dt = set_dt;
        double sampling_timestep_multiple = round(1/dt);
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(sampling_timestep_multiple);
        simulator.SetEndTime(400.0);
        /*------------------------------------END: Timestep---------------------------------------*/


        /*--------------------------------START: Boundary condition-----------------------------*/
        if (if_consider_compression)
        {
          c_vector<double,2> point = zero_vector<double>(2);
          c_vector<double,2> normal = zero_vector<double>(2);
          normal(1) = -1.0;
          double stop_time = DOUBLE_UNSET;
          MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc, (&cell_population, point, normal, stop_time));
          simulator.AddCellPopulationBoundaryCondition(p_bc);
        }
        /*--------------------------------END: Boundary condition-----------------------------*/
        

        /*-------------------------------START: Output Directory---------------------------------*/
        // concise information
        oss.str("");
        oss << std::fixed << setprecision(3) << dt;
        out_put_directory += "_|_Dt=" + oss.str();
        out_put_directory += "_ResMove=" + std::to_string(restrict_vertex_movement);
        out_put_directory += "_AdaptDt=" + std::to_string(use_adaptive_timestep);
        oss.str("");
        oss << std::fixed << setprecision(3) << cell_rearrangement_threshold;
        out_put_directory += "_RearThr=" + oss.str();
        out_put_directory += "_|_CellDivi=" + std::to_string(run_with_birth);
        oss.str("");
        oss << std::fixed << setprecision(3) << nagai_honda_membrane_surface_energy_parameter;
        out_put_directory += "_|_Gamma=" + oss.str();
        oss.str("");
        oss << std::fixed << setprecision(1) << target_shape_index;
        out_put_directory += "_ShapeIndex=" + oss.str();
        oss.str("");
        oss << std::scientific << setprecision(2) << nagai_honda_cell_cell_adhesion_energy_parameter;
        out_put_directory += "_CCAdhe=" + oss.str();
        oss.str("");
        oss << std::scientific << setprecision(2) << nagai_honda_cell_boundary_adhesion_energy_parameter;
        out_put_directory += "_CBAdhe=" + oss.str();
        out_put_directory += "_|_ConsiSubAdh=" + std::to_string(if_consider_substrate_adhesion);
        if (if_consider_substrate_adhesion)
          out_put_directory += "_StripSubAdhIsHomo=" + std::to_string(if_substrate_adhesion_is_homogeneous);
        out_put_directory += "_|_ConsiRSA=" + std::to_string(if_consider_reservoir_substrate_adhesion);
        out_put_directory += "_|_MSE=" + std::to_string(case_number_of_membrane_surface_energy_form);
        out_put_directory += "_|_AddFeedback=" + std::to_string(if_consider_feedback_of_face_values);
        out_put_directory += "_|_AddRandForce=" + std::to_string(add_random_force);
        
        // detailed information:
        oss.str("");
        oss << std::fixed << setprecision(3) << dt;
        out_put_directory += "/Dt=" + oss.str();
        if (run_with_birth)
        {
          oss.str("");
          oss << std::fixed << setprecision(1) << time_for_one_division_of_cell_population;
          out_put_directory += "_|_TOneDivi=" + oss.str();
          oss.str("");
          oss << std::fixed << setprecision(2) << growth_rate_for_target_area_after_division;
          out_put_directory += "_GrRate=" + oss.str();
        }
        if(if_consider_substrate_adhesion)
        {
          if (!if_substrate_adhesion_is_homogeneous)
          {
            oss.str("");
            oss << std::fixed << setprecision(1) << substrate_adhesion_leading_top_length;
            out_put_directory += "_|_SSALeadTopLeng" + oss.str();
            oss.str("");
            oss << std::fixed << setprecision(1) << substrate_adhesion_parameter_at_leading_top;
            out_put_directory += "_SSATop=" + oss.str();
            oss.str("");
            oss << std::fixed << setprecision(1) << substrate_adhesion_parameter;
            out_put_directory += "_SSAEnd=" + oss.str();
          }
          else
          {
            oss.str("");
            oss << std::fixed << setprecision(1) << substrate_adhesion_parameter;
            out_put_directory += "_|_HomoSSA=" + oss.str();
          }
        }
        if(if_consider_reservoir_substrate_adhesion)
        {
          // out_put_directory += "_IgnRSABottom=" + std::to_string(if_ignore_reservoir_substrate_adhesion_at_bottom);
          oss.str("");
          oss << std::fixed << setprecision(1) << reservoir_substrate_adhesion_parameter;
          out_put_directory += "_|_RSA=" + oss.str();
        }

        // if(if_consider_compression)
        // {
        //   oss.str("");
        //   oss << std::fixed << setprecision(2) << compress_ratio;
        //   out_put_directory += "_ComprRatio=" + oss.str();
        // }
        // oss.str("");
        // oss << std::fixed << setprecision(1) << strip_start_x_location;
        // out_put_directory += "_StripXLoc=" + oss.str();

        if (if_consider_feedback_of_face_values)
        {
          oss.str("");
          oss << std::scientific << setprecision(0) << feedback_strength_for_myosin_activity ;
          out_put_directory += "_|_FeStrenMyo=" + oss.str();
          oss.str("");
          oss << std::scientific << setprecision(0) << feedback_strength_for_adhesion;
          out_put_directory += "_FeStrenAdh=" + oss.str();
          // oss.str("");
          // oss << std::fixed << setprecision(0) << hill_coefficient_for_myosin_activity ;
          // out_put_directory += "_MyoHillCoeff=" + oss.str();
          // oss.str("");
          // oss << std::fixed << setprecision(0) << hill_coefficient_for_adhesion;
          // out_put_directory += "_AdhHillCoeff=" + oss.str();
          out_put_directory += "_EMACanDe=" + std::to_string(!EMA_dont_decrease);
          out_put_directory += "_CCACanDe=" + std::to_string(!CCA_dont_decrease);
          out_put_directory += "_CCACanIn=" + std::to_string(!CCA_dont_increase);
          if (CCA_dont_inrease_until_shorter_than_a_threshold)  
          {
            oss.str("");
            oss << std::fixed << setprecision(1) << CCA_dont_increase_until_shorter_than_this_value;          
            out_put_directory += "_CCAInrThresh=" + oss.str();
          }
        }
        if (add_random_force)
        {
          oss.str("");
          oss << std::scientific << setprecision(2) << translational_diffusion_constant;
          out_put_directory += "_|_D=" + oss.str();
          if (consider_polarity)
          {
            oss.str("");
            oss << std::scientific << setprecision(2) << polarity_magnitude;
            out_put_directory += "_fp=" + oss.str();
            oss.str("");
            oss << std::scientific << setprecision(2) << rotational_diffusion_constant;
            out_put_directory += "_Dr=" + oss.str();
          }
        }

        simulator.SetOutputDirectory(out_put_directory);
        std::cout << std::endl << "OutputDirectory: " << out_put_directory << std::endl;
        /*-------------------------------END: Output Directory---------------------------------*/

        simulator.Solve();
    }

};

#endif /* TESTMYXTOROIDALSTRIPESADHESIONSIMULATIONRAW_HPP_ */
