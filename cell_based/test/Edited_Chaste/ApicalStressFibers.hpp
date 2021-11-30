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

/*

This program is customized for investigating the roles of apical stress fibers
during the Drosophila development, particularly at the stage of 18-26hAPF (hours 
after pupa formation), by Dr. Chao FANG at Department of Mechanical Engineering,
The University of Hong Kong.

*/

#ifndef APICALSTRESSFIBERS_HPP_
#define APICALSTRESSFIBERS_HPP_

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
#include "PlaneBasedCellKiller.hpp"

// #include <ctime>

class ApicalStressFibers : public AbstractCellBasedTestSuite
{
public:

    void TestApicalStressFibers()
    {
        /*------------------------------START: Basic Settings----------------------------*/
        double target_shape_index = 4.25;//p0
        double reference_area = M_PI;
        double initial_area = reference_area;

      /* Strip Structure & Cell Mesh */
        bool   strip_width_doubled_for_multiple_leading_cells = false;
        double strip_width_mutiple = 8.0;
      //  double strip_width_multiple_for_sliding = 15.0;
        bool   if_use_larger_strip_distance = false;
        double strip_dis_multiplier = 40.0/12.0; // strip center distance, ensure the (num_ele_cross) being even number
        bool   use_longer_mesh = false; // for (num_ele_up) mesh
        int    num_ele_up_multiplier = 4;
        int    move_mesh_right_for_N_periods = 0; // for display of multiple periods
        bool   one_strip_only_in_a_period = true;

        unsigned num_ele_cross = 40; // must be even number
        if  (if_use_larger_strip_distance)
            num_ele_cross = (unsigned)round(num_ele_cross*strip_dis_multiplier);

        unsigned num_ele_up = 40;
        if  (use_longer_mesh)
            num_ele_up *= num_ele_up_multiplier;

        double center_of_width = 0.0;       // change made by Chao
        double width = num_ele_cross*sqrt(initial_area/(sqrt(3)/2));   //width of reservoir, change made by Chao

        double strip_width = 20*sqrt(initial_area/(sqrt(3)/2)); // default =0.9523 (1/2 cell width)
        if  (strip_width_doubled_for_multiple_leading_cells)
            strip_width = strip_width*strip_width_mutiple;

        double strip_distance = width;

        double strip_start_x_location = 0.0;
        if  (strip_width_doubled_for_multiple_leading_cells)
            strip_start_x_location += 1/2*sqrt(initial_area/(sqrt(3)/2));
        double strip_start_y_location = 1.5*num_ele_up*1/sqrt(3)*sqrt(initial_area/(sqrt(3)/2));

        /* Energy equation form: 1/2*KA*(A-A0)^2 + 1/2*Gamma*P^2 + 1.0*Lambda*P - alpha*Sa. ( Lambda=-Ga*P0, p0=P0/sqrt(A0) ) */
      /* 1. Area expansion resistance */
        double nagai_honda_deformation_energy_parameter = 1; // KA
        bool   use_fixed_target_area_without_modifier = true; // A0:

      /* 2. Myosin activity */
        double nagai_honda_membrane_surface_energy_parameter = 0.2/(M_PI/reference_area);//Gamma
        double edge_length_at_rest = sqrt(initial_area/(6*sqrt(3)/4)); // = 1.0996

      /* 3. Cell-Cell adhesion */
        double cell_cell_adhesion_parameter = -nagai_honda_membrane_surface_energy_parameter*(target_shape_index*sqrt(reference_area));
        double cell_boundary_adhesion_parameter = 0.0; // cell-cell adhesion at boundary

      /* 4. Morphogenetic Force */
        // Note that pulling force is realized by different ways for epithelial bridge and vortex formation
        // epithelial bridge: pulling force is directly applied on the nodes of leading cells
        // vortex formation: pulling force is treated equivalently by larger substrate adhesion on strip 
        bool   multiple_leading_cells = true;
        // unsigned leading_cell_number = 5;
        unsigned leading_cell_number = (unsigned) round( strip_width/(sqrt(initial_area/(sqrt(3)/2))) );
        if (!multiple_leading_cells)
           leading_cell_number = 1;
        // double pulling_force_on_leading_cell = 3.0*leading_cell_number/pow((M_PI/reference_area),1.5);// Fy
        bool   add_pulling_force_evenly_on_nodes_of_leading_cell = true;
        double pulling_force_on_leading_cell = 3.0/pow((M_PI/reference_area),1.5);// Fy

      /* 5. Random force */
        bool   add_random_force = true;

        // Brownian random force
        bool   has_brownian_random_force = false;
        double translational_diffusion_constant = 0.0;
        double set_node_radius = 1.0/pow(1.0,2.0)*50*1e2*((M_PI/reference_area)*(M_PI/reference_area)); // effective radius of node in Einstein relation

        // Cell polarity
        bool   has_polarity = true;
        // double polarity_magnitude = 0.2;
        bool   seed_manually = true;
        // unsigned seed_for_initial_random_polarity = 1u;
        // seed_for_initial_random_polarity += 10;
        double rotational_diffusion_constant = 0.4/(M_PI/reference_area);

        if (polarity_magnitude==0.0)
           has_polarity = false;
        
        if ((!has_brownian_random_force)&&(!has_polarity))
           add_random_force = false;
        if (polarity_magnitude!=0.0)
           assert(add_random_force == true);

        double polarity_magnitude = 0.4;
        unsigned seed_for_initial_random_polarity = 3u;

      /* 6. Time */
        double end_time_for_equilibrium = 0.0;
        bool   if_equilibrate_for_a_while = true;
        if (end_time_for_equilibrium <= 0.0)
           if_equilibrate_for_a_while = false;
        double polarity_magnitude_equilibrium = 0.5;
        
        double dt = 0.05*(M_PI/reference_area); // Previously 0.025
        // double end_time = 800.0*(M_PI/reference_area);
        double max_movement_per_timestep = 0.05/sqrt((M_PI/reference_area)); // Previously 0.05

        bool   apply_my_change_to_make_timestep_adaptive = true;

        double end_time = 800.0*(M_PI/reference_area);

      /* 7. Cell rearrangement */
        // the minimum threshold distance for element T1 rearrangement 
        double cell_rearrangement_threshold = 0.05/sqrt((M_PI/reference_area)); 

      /* 8. Output & display */
        bool   output_concise_swap_information_when_remesh = false;
        bool   output_detailed_swap_information_when_remesh = false;
        bool   output_numerical_method_information = false;
        bool   output_information_for_nagai_honda_force = false;     
        bool   if_set_cell_data_of_detailed_force_contributions = false; // for output and display each component in Paraview

        /*------------------------------END: Basic Settings----------------------------*/


        /*------------------------------START: Mesh Structure------------------------------*/
        // Strips structure of substrate adhesion
        bool   if_update_face_elements_in_mesh = if_consider_feedback_of_cell_cell_adhesion;
        
        MyXToroidalHoneycombVertexMeshGenerator generator(num_ele_cross, num_ele_up, initial_area, cell_rearrangement_threshold, 0.001/sqrt((M_PI/reference_area)));
        MyXToroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();
        p_mesh->SetMultiplyResultsBy(multiply_results_by);
        p_mesh->SetDistanceForT3SwapChecking(5.0/sqrt((M_PI/reference_area)));
        p_mesh->SetMoveMeshRightForNPeriods(move_mesh_right_for_N_periods);        
        p_mesh->SetUpdateFaceElementsInMeshBoolean(if_update_face_elements_in_mesh);
        p_mesh->SetIfClassifyElementsWithGroupNumbers(classify_elements_with_group_numbers);
        p_mesh->SetMarkLeadingCells(mark_leading_cells);
        p_mesh->SetIfCheckForT4Swaps(if_check_for_T4_swaps);
        p_mesh->SetOutputConciseSwapInformationWhenRemesh(output_concise_swap_information_when_remesh);
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
        cell_population.SetRestrictVertexMovementBoolean(restrict_vertex_movement);
        /*------------------------------END: Cells and Cell population--------------------------*/


        OffLatticeSimulation<2> simulator(cell_population);
        // output cell velocity:
        bool my_output_cell_velocity = true;
        simulator.SetMyOutputCellVelocities(my_output_cell_velocity);
        if (my_output_cell_velocity && seed_manually)
          simulator.SetMySeed(seed_for_initial_random_polarity);
        bool output_cell_velocity = true;
        simulator.SetOutputCellVelocities(output_cell_velocity);
        simulator.SetNoBirth(!run_with_birth);

        /*--------------------------------START: TargetAreaModifier------------------------------*/
        MAKE_PTR_ARGS(TargetAreaLinearGrowthModifier<2>, p_growth_modifier, ());
        p_growth_modifier->SetUseUseMyOwnRuleForUpdateTargetAreaOfCell(true);
        p_growth_modifier->SetReferenceTargetArea(reference_area);
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
        p_numerical_method->SetOutputNumericalMethodInformation(output_numerical_method_information);
        simulator.SetNumericalMethod(p_numerical_method);
        /*---------------------------------END: Add Numerical Method-----------------------------*/


        /*-----------------------------START: MyNagaiHondaForceWithStripesAdhesion---------------------*/
        MAKE_PTR(MyNagaiHondaForceWithStripesAdhesion<2>, p_force);
        
//         int  case_number_of_membrane_surface_energy_form = 1; // 1 for default
// /*?*/   bool if_use_face_element_to_get_adhesion_parameter = if_consider_feedback_of_cell_cell_adhesion;

//         p_force->SetCaseNumberOfMembraneSurfaceEnergyForm(case_number_of_membrane_surface_energy_form);
// /*?*/   p_force->SetUseFaceElementToGetAdhesionParameterBoolean(if_use_face_element_to_get_adhesion_parameter);

        p_force->SetNagaiHondaDeformationEnergyParameter(nagai_honda_deformation_energy_parameter); // KA
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(nagai_honda_membrane_surface_energy_parameter); // Gamma
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(cell_cell_adhesion_parameter); // Lambda
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(cell_boundary_adhesion_parameter);        

        p_force->SetUseFixedTargetArea(use_fixed_target_area_without_modifier); // used in the case where there is no target area modifier! (no division)
        p_force->SetFixedTargetArea(reference_area);
        p_force->SetTargetShapeIndex(target_shape_index);
          // Strip Substrate Adhesion
        p_force->SetAddPullingForceEvenlyOnNodesOfLeadingCell(add_pulling_force_evenly_on_nodes_of_leading_cell);
        p_force->SetLeadingCellNumber(leading_cell_number);
        p_force->SetPullingForceOnLeadingCell(pulling_force_on_leading_cell);
        p_force->SetSubstrateAdhesionLeadingTopLength(substrate_adhesion_leading_top_length);
          // for consideration of periodicity
        p_force->SetCenterOfWidth(center_of_width);
        p_force->SetWidth(width);//check later!
        p_force->SetReservoirSubstrateAdhesionParameter(reservoir_substrate_adhesion_parameter);
          // strip info
        p_force->SetStripWidth(strip_width);
        p_force->SetStripDistance(strip_distance);
        p_force->SetStripStartXLocation(strip_start_x_location);
        p_force->SetStripStartYLocation(strip_start_y_location);

        p_force->SetIfEquilibrateForAWhile(if_equilibrate_for_a_while);
        p_force->SetEndTimeForEquilibrium(end_time_for_equilibrium);

        p_force->SetOutputInformationForNagaiHondaForce(output_information_for_nagai_honda_force);
        simulator.AddForce(p_force);
        /*-----------------------------END: MyNagaiHondaForceWithStripesAdhesion---------------------*/


        /*-----------------------------START: RandomForce and PolarityModifier--------------------------*/
        // Brownian diffusion
        if (add_random_force)
        {
          MAKE_PTR(DiffusionForce<2>, p_force2);
          bool   use_the_same_node_radius = true;
          double node_radius = 50*1e2*((M_PI/reference_area)*(M_PI/reference_area));
          if  (translational_diffusion_constant >0.0)
              node_radius = p_force2->GetDiffusionScalingConstant()/translational_diffusion_constant;
          else
              node_radius = set_node_radius; // 50*1e2*((M_PI/reference_area)*(M_PI/reference_area));
          translational_diffusion_constant = p_force2->GetDiffusionScalingConstant()/node_radius;

          p_force2->SetHasBrownianRandomForce(has_brownian_random_force);
          p_force2->SetUseTheSameNodeRadius(use_the_same_node_radius);
          p_force2->SetTheSameNodeRadius(node_radius);
          p_force2->SetHasPolarity(has_polarity);
          p_force2->SetOnePeriodOnly(one_strip_only_in_a_period);

          p_force2->SetIfEquilibrateForAWhile(if_equilibrate_for_a_while);
          p_force2->SetEndTimeForEquilibrium(end_time_for_equilibrium);

          simulator.AddForce(p_force2);
        }

        // Polarity modifier
        if (has_polarity)
        {
          MAKE_PTR_ARGS(PolarityModifier<2>, p_polarity_modifier, ());
          p_polarity_modifier->SetPolarityMagnitude(polarity_magnitude);
          p_polarity_modifier->SetPolarityMagnitudeEquilibrium(polarity_magnitude_equilibrium);
          p_polarity_modifier->SetSeedManually(seed_manually);
          p_polarity_modifier->SetSeedForInitialRandomPolarity(seed_for_initial_random_polarity);
          p_polarity_modifier->SetD(rotational_diffusion_constant);
          simulator.AddSimulationModifier(p_polarity_modifier);
        }
        /*-----------------------------END: RandomForce and PolarityModifier---------------------------*/


        /*------------------------START: Modifier: need modification---------------*/
        MAKE_PTR_ARGS(FaceValueAndStressStateModifier<2>, p_face_value_and_stress_state_modifier, ());
        
        // equilibrium related
        p_face_value_and_stress_state_modifier->SetIfEquilibrateForAWhile(if_equilibrate_for_a_while);
        p_face_value_and_stress_state_modifier->SetEndTimeForEquilibrium(end_time_for_equilibrium);

        simulator.AddSimulationModifier(p_face_value_and_stress_state_modifier);
        /*-----------------------END: Modifier: need modification---------------*/


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
        unsigned sampling_timestep_multiple = (unsigned) round(sampling_time/dt);

        simulator.SetApplyMyChangesToMakeTimestepAdaptive(apply_my_change_to_make_timestep_adaptive);
        simulator.SetDt(dt);

        simulator.SetApplySamplingTimeInsteadOfSamplingTimestep(apply_my_change_to_make_timestep_adaptive);
        simulator.SetSamplingTimestepMultiple(sampling_timestep_multiple);
        simulator.SetSamplingTime(sampling_time);

        simulator.SetEndTime(end_time);
        /*------------------------------------END: Timestep---------------------------------------*/


        /*--------------------------------START: Boundary condition-----------------------------*/
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(1) = -1.0;
        double stop_time = DOUBLE_UNSET;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc, (&cell_population, point, normal, stop_time));
        c_vector<double,2> point2 = zero_vector<double>(2);
        c_vector<double,2> normal2 = zero_vector<double>(2);
        point2(1) = strip_start_y_location;
        normal2(1) = 1.0;
        double stop_time2 = end_time_for_equilibrium;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point2, normal2, stop_time2));
        simulator.AddCellPopulationBoundaryCondition(p_bc);
        simulator.AddCellPopulationBoundaryCondition(p_bc2);
        /*--------------------------------END: Boundary condition-----------------------------*/
        

        /*--------------------------START: Output Directory and Simulation Information File---------------------*/
        // Output directory:
        std::ostringstream oss;
        time_t raw_time = time(0);
        struct tm * now = localtime(& raw_time);

        std::string output_directory = "PHASE-DIAGRAM/Date: ";
        oss.str("");
        oss << (now->tm_year + 1900 -2000) << '-' << (now->tm_mon + 1) << '-' <<  now->tm_mday << '/';
        output_directory += oss.str();

        oss.str("");
        oss << "W=" << std::fixed << setprecision(0) << round(strip_width/0.95) << '/';
        output_directory += oss.str();

        oss.str("");

        if (multiple_leading_cells)
        oss << "_Fy=" << std::fixed  << setprecision(0) << pulling_force_on_leading_cell;

        oss << "_Fp=" << ((polarity_magnitude>=0.01 || polarity_magnitude==0.0)? std::fixed : std::scientific) << setprecision(2) << polarity_magnitude;
        if (polarity_magnitude!=0.0)
        {
          oss << "_Dr=" << ((rotational_diffusion_constant>=0.01)? std::fixed : std::scientific) << setprecision(2) << rotational_diffusion_constant;
          if (seed_manually)
            oss << "_PSeed=" << seed_for_initial_random_polarity;
          else
            oss << "_PSeed=N";
        }
        if (add_random_force&&(has_brownian_random_force))
          oss << "_Brown:D=" << std::scientific << setprecision(2) << translational_diffusion_constant;
        else
          oss << "_Brown=0";

        oss << "_p0=" << std::fixed << setprecision(2) << target_shape_index;        

        oss << "_Dt=" << std::scientific << setprecision(1) << dt;
        if (apply_my_change_to_make_timestep_adaptive)
          oss << "_MaxMv=" << ((max_movement_per_timestep>=0.01)? std::fixed : std::scientific) << setprecision(3) << max_movement_per_timestep;
        oss << "_T1Thresh=" << ((cell_rearrangement_threshold>=0.01)? std::fixed : std::scientific) << setprecision(3) << cell_rearrangement_threshold;

        // oss << "_A0=" << std::fixed << setprecision(2) << reference_area;
        oss << "_Ga=" << ((nagai_honda_membrane_surface_energy_parameter>=0.01 || nagai_honda_membrane_surface_energy_parameter==0.0)? std::fixed : std::scientific) 
                << setprecision(2) << nagai_honda_membrane_surface_energy_parameter;

        if (num_ele_up!=8)
          oss << "_NumUp=" << num_ele_up;
        if (num_ele_cross!=6)
          oss << "_NumCr=" << num_ele_cross;
        if ( fabs( strip_width - 0.5*sqrt(3)*sqrt(initial_area/(6*sqrt(3)/4)) )>1e-10 )
          oss << "_SWid=" << std::fixed << setprecision(0) << strip_width;
        if ( fabs(strip_distance - 6*sqrt(3)*sqrt(initial_area/(6*sqrt(3)/4)))>1e-10 )
          oss << "_SDis=" << std::fixed << setprecision(0) << strip_distance; 
        if (move_mesh_right_for_N_periods!=0)
          oss << "_MvRight=" << std::fixed << setprecision(0) << move_mesh_right_for_N_periods;

        output_directory += oss.str();

        oss.str("");
        oss << (now->tm_year + 1900 -2000) << '-' << (now->tm_mon + 1) << '-' <<  now->tm_mday 
            << ", " << now->tm_hour << ':' << now->tm_min << ':' << now->tm_sec;
        output_directory += "_Timestamp=" + oss.str();

        std::string concise_output_directory = output_directory;
        
        // Concise information written to directoory.
        oss.str("");
        oss << std::fixed << setprecision(4) << dt;
        output_directory += "/Dt=" + oss.str();
        oss.str("");
        oss << ((cell_rearrangement_threshold>=0.01)? std::fixed : std::scientific) << setprecision(3) << cell_rearrangement_threshold;
        output_directory += "_T1Thresh=" + oss.str();
        if (apply_my_change_to_make_timestep_adaptive)
        {
          oss.str("");
          oss << ((max_movement_per_timestep>=0.01)? std::fixed : std::scientific) << setprecision(3) << max_movement_per_timestep;
          output_directory += "_MaxMvDt=" + oss.str();
        }
        if (!apply_my_change_to_make_timestep_adaptive)
          output_directory += "_MyAdaptDt=0";

        output_directory += "_HasRandF=" + std::to_string(add_random_force);
        
        oss.str("");
        oss << ((nagai_honda_membrane_surface_energy_parameter>=0.01 || nagai_honda_membrane_surface_energy_parameter==0.0)? std::fixed : std::scientific) 
                << setprecision(2) << nagai_honda_membrane_surface_energy_parameter;
        output_directory += "_|Ga=" + oss.str();
        oss.str("");
        oss << std::fixed << setprecision(2) << target_shape_index;
        output_directory += "_p0=" + oss.str();
        
        // Detailed information  written to directoory.
        output_directory += "_|Structure:";
        oss << "StripDis=" << std::fixed << setprecision(3) << strip_distance;
        oss << "_StripWid=" << std::fixed << setprecision(3) << strip_width;

        if (add_random_force)
        {
          if (has_brownian_random_force)
          {
            oss.str("");
            oss << std::scientific << setprecision(2) << translational_diffusion_constant;
            output_directory += "_|D_trans=" + oss.str();
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

        // tmp
        if (if_equilibrate_for_a_while)
        {
          oss.str("");
          oss << "End_time_equilib=" << std::fixed << setprecision(1) << end_time_for_equilibrium;
          output_directory += oss.str();
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

#endif /* APICALSTRESSFIBERS_HPP_ */
