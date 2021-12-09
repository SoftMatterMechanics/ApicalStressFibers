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
// #include "FaceValueAndStressStateModifier.hpp"
#include "PolarityModifier.hpp"

#include "FakePetscSetup.hpp"

#include "PlaneBoundaryCondition.hpp"
#include "ToroidalHoneycombVertexMeshGenerator.hpp"
#include "MyToroidalHoneycombVertexMeshGenerator.hpp"
// #include "NagaiHondaForce.hpp"
#include "MyNagaiHondaForce.hpp"
#include "MyMorphogeneticForce.hpp"
#include "DiffusionForce.hpp"
#include "PlaneBasedCellKiller.hpp"

// #include <ctime>

class ApicalStressFibers : public AbstractCellBasedTestSuite
{
public:

    void TestApicalStressFibers()
    {
      /*-------------------------START: Basic Settings-----------------------*/
      /* Energy equation form: 1/2*KA*(A-A0)^2 + 1/2*KL*(L-L0)^2 + 2*Gamma*L */

      // 1. Cell mesh and size
        unsigned num_ele_across = 8; // cell number along anterior-posterior, must be an even number
        unsigned num_ele_up = 8; // cell number along medial-lateral, must be an even number

        double target_shape_index = 4.0;
        double initial_area = M_PI; // set the initial areas to be uniform, A = 3*sqrt(3)/2*l^2
        
        double center_y_coordination = 0.5*(1.5*num_ele_up+0.5)*sqrt(initial_area/(3*sqrt(3)/2));

      // 2. Area elasticity
        double area_elastic_modulus = 1.0; // KA
        bool   use_fixed_target_area_without_modifier = false;
        double target_area = M_PI; // A0, target areas are not uniform
        unsigned random_seed_for_target_area = 1;
        double min_target_area = initial_area/5;
        double max_target_area = initial_area*9/5;

      // 3. Edge elasticity
        double edge_elastic_modulus = 0.04; // KL
        // double rest_edge_length = sqrt(target_area/(6*sqrt(3)/4));

      // 4. Cell-Cell adhesion & constant cortical contraction
        double cell_cell_adhesion_energy_density = -2*target_shape_index*edge_elastic_modulus; // Gamma, this parameter consists of cell-cell adhesion and cortical contraction 
        double cell_boundary_adhesion_energy_density = -2*target_shape_index*edge_elastic_modulus; // Gamma at boundary
        bool   if_use_face_element_to_get_adhesion_parameter = false;

      // 5. Random force
        bool   add_random_force = true;
        bool   has_brownian_random_force = false; // brownian random force is used in cell center model
        double translational_diffusion_constant = 0.0;
        double set_node_radius = 2.0; // effective cell diameter (in units of 10 microns)

        bool   has_polarity = true;
        bool   seed_manually = true;
        unsigned seed_for_initial_random_polarity = 1;
        double polarity_magnitude_equilibrium = 0.02;  // for before equilibrium
        double polarity_magnitude = 0.0;  // for after equilibrium
        double rotational_diffusion_constant = 0.2;

      // 6. Morphogenetic Force
        bool     multiple_leading_cells = true;
        bool     add_pulling_force_evenly_on_nodes_of_leading_cell = true;
        double   pulling_force_on_leading_cell = 0.0;// MF (morphogenetic force)

      // 7. Cell rearrangement threshold length for T1 & T2 transitions
        double cell_rearrangement_threshold = 0.05; 
        double t2_threshold = 0.001;

      // 8. Time
        bool   if_equilibrate_for_a_while = true;
        double time_for_equilibrium = 200.0;
        if (time_for_equilibrium <= 0.0)
           if_equilibrate_for_a_while = false;
        
        double dt = 0.01;
        double end_time = 300.0;
        double max_movement_per_timestep = 0.001; 
        bool   apply_adaptive_timestep = true;
        double sampling_time = 1.0;
        unsigned sampling_timestep_multiple = (unsigned) round(sampling_time/dt);

      // 9. Output & display
        bool   if_update_face_elements_in_mesh = true;
        bool   output_concise_swap_information_when_remesh = false;
        bool   output_detailed_swap_information_when_remesh = false;
        bool   output_numerical_method_information = false;
        // bool   output_information_for_nagai_honda_force = false;     
        // bool   if_set_cell_data_of_detailed_force_contributions = false; // for output and display each component in Paraview
        bool   use_concise_output_directory = true;
      /*------------------------------END: Basic Settings----------------------------*/


      /*-----------------------START: Generate cell monolayer mesh-------------------*/
        MyToroidalHoneycombVertexMeshGenerator generator(num_ele_across, num_ele_up, initial_area, cell_rearrangement_threshold, t2_threshold);
        MyToroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();

        p_mesh->SetDistanceForT3SwapChecking(5.0); 
        p_mesh->SetUpdateFaceElementsInMeshBoolean(if_update_face_elements_in_mesh);
      //  p_mesh->SetIfClassifyElementsWithGroupNumbers(classify_elements_with_group_numbers);
      //  p_mesh->SetMarkLeadingCells(mark_leading_cells);
        p_mesh->SetOutputConciseSwapInformationWhenRemesh(output_concise_swap_information_when_remesh);
        p_mesh->SetOutputDetailedSwapInformationWhenRemesh(output_detailed_swap_information_when_remesh);
      /*------------------------END: Generate cell monolayer mesh--------------------*/


      /*---------------------------START: Cells and Cell population--------------------------*/
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<BernoulliTrialCellCycleModel, 2> cells_generator;
        
        bool use_bernoulli_trial_cell_cycle_model = true;
        double period_of_cell_division = DOUBLE_UNSET;
        double divide_probability_for_a_cell = 1.0/period_of_cell_division/(num_ele_across*num_ele_up);
        double minimum_division_age = -0.01;
        cells_generator.SetUseBernoulliTrialCellCycleModel(use_bernoulli_trial_cell_cycle_model);        
        cells_generator.SetDivisionProbability(divide_probability_for_a_cell);
        cells_generator.SetMinimumDivisionAge(minimum_division_age);
        cells_generator.SetRandomSeedForTargetAreas(random_seed_for_target_area);
        cells_generator.SetLimitsOfTargetAreas(min_target_area, max_target_area);

        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);
        
        bool restrict_vertex_movement = false;
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetRestrictVertexMovementBoolean(restrict_vertex_movement);
      /*------------------------------END: Cells and Cell population--------------------------*/


      /*---------------------------------START: Simulator settings----------------------------*/
        OffLatticeSimulation<2> simulator(cell_population);

        // output cell velocity:
        bool my_output_cell_velocity = true;
        simulator.SetMyOutputCellVelocities(my_output_cell_velocity);
        if (my_output_cell_velocity && seed_manually)
          simulator.SetMySeed(seed_for_initial_random_polarity);
        
        bool output_cell_velocity = true;
        bool run_with_birth = false;
        simulator.SetOutputCellVelocities(output_cell_velocity);
        simulator.SetNoBirth(!run_with_birth);

        // Timestep
        simulator.SetApplyAdaptiveTimestep(apply_adaptive_timestep);
        simulator.SetDt(dt);
        simulator.SetApplySamplingTimeInsteadOfSamplingTimestep(apply_adaptive_timestep);
        simulator.SetSamplingTimestepMultiple(sampling_timestep_multiple);
        simulator.SetSamplingTime(sampling_time);
        simulator.SetEndTime(end_time);
        
        // Numerical Method
        boost::shared_ptr<ForwardEulerNumericalMethod<2,2>> p_numerical_method(new ForwardEulerNumericalMethod<2,2>());
        p_numerical_method->SetCellPopulation(&cell_population);
        //Note: Here the address of std::vector of force collections are different from that in simulator!
        //But the corresponding shared_ptrs in the vectors should be the same.
        std::vector< boost::shared_ptr<AbstractForce<2, 2>> > force_collection = simulator.rGetForceCollection();
        p_numerical_method->SetForceCollection(&force_collection);
        p_numerical_method->SetMaxMovementPerTimestep(max_movement_per_timestep);
        p_numerical_method->SetOutputNumericalMethodInformation(output_numerical_method_information);
        simulator.SetNumericalMethod(p_numerical_method);
      /*---------------------------------END: Simulator settings-----------------------------*/


      /*--------------------------------START: Modifier Settings-----------------------------*/
        // // TargetAreaModifier
        // MAKE_PTR_ARGS(TargetAreaLinearGrowthModifier<2>, p_growth_modifier, ());

        // double growth_rate_for_target_area_after_division = 0.0;
        // p_growth_modifier->SetUseUseMyOwnRuleForUpdateTargetAreaOfCell(true);
        // p_growth_modifier->SetReferenceTargetArea(initial_area);// to be determined
        // p_growth_modifier->SetGrowthRate(growth_rate_for_target_area_after_division);
        // p_growth_modifier->SetAgeToStartGrowing(0.0);
        // simulator.AddSimulationModifier(p_growth_modifier);

        // // Contractility and Cell-Cell Adhesion Modifier
        // MAKE_PTR_ARGS(FaceValueAndStressStateModifier<2>, p_face_value_and_stress_state_modifier, ());
        
        // p_face_value_and_stress_state_modifier->SetIfEquilibrateForAWhile(if_equilibrate_for_a_while);
        // p_face_value_and_stress_state_modifier->SetEndTimeForEquilibrium(time_for_equilibrium);
        // simulator.AddSimulationModifier(p_face_value_and_stress_state_modifier);

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
      /*----------------------------------END: Modifier Settings------------------------------*/


      /*---------------------------------START: MyMorphogeneticForce-----------------------------*/
        MAKE_PTR(MyMorphogeneticForce<2>, p_morphogenetic_force);
        
        p_morphogenetic_force->SetAddPullingForceEvenlyOnNodesOfLeadingCell(add_pulling_force_evenly_on_nodes_of_leading_cell);
        p_morphogenetic_force->SetPullingForceOnLeadingCell(pulling_force_on_leading_cell);
        p_morphogenetic_force->SetCenterYCoordination(center_y_coordination);
        p_morphogenetic_force->SetIfEquilibrateForAWhile(if_equilibrate_for_a_while);
        p_morphogenetic_force->SetEndTimeForEquilibrium(time_for_equilibrium);

        simulator.AddForce(p_morphogenetic_force);
      /*-------------------------------------END: MyNagaiHondaForce------------------------------*/


      /*---------------------------------START: MyNagaiHondaForce-----------------------------*/
        MAKE_PTR(MyNagaiHondaForce<2>, p_nh_force);
        
        p_nh_force->SetNagaiHondaDeformationEnergyParameter(area_elastic_modulus); // KA
        p_nh_force->SetNagaiHondaMembraneSurfaceEnergyParameter(edge_elastic_modulus); // KL
        p_nh_force->SetNagaiHondaCellCellAdhesionEnergyParameter(cell_cell_adhesion_energy_density); // Gamma
        p_nh_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(cell_boundary_adhesion_energy_density); // Gamma at boundary

        p_nh_force->SetUseFixedTargetArea(use_fixed_target_area_without_modifier); // used in the case where there is no target area modifier! (no division)
        p_nh_force->SetFixedTargetArea(target_area); // to be determined
        p_nh_force->SetTargetShapeIndex(target_shape_index);
        p_nh_force->SetUseFaceElementToGetAdhesionParameterBoolean(if_use_face_element_to_get_adhesion_parameter);

        // // Morphogenetic force
        // p_force->SetAddPullingForceEvenlyOnNodesOfLeadingCell(add_pulling_force_evenly_on_nodes_of_leading_cell);
        // p_force->SetPullingForceOnLeadingCell(pulling_force_on_leading_cell);

        // p_force->SetIfEquilibrateForAWhile(if_equilibrate_for_a_while);
        // p_force->SetEndTimeForEquilibrium(time_for_equilibrium);

        // p_force->SetOutputInformationForNagaiHondaForce(output_information_for_nagai_honda_force);
        simulator.AddForce(p_nh_force);
      /*-------------------------------------END: MyNagaiHondaForce------------------------------*/


      /*-----------------------------------START: Brownian Random Force -------------------------*/
        // Brownian diffusion
        if (add_random_force)
        {
          MAKE_PTR(DiffusionForce<2>, p_random_force);
          bool   use_the_same_node_radius = true;
          double node_radius = 2.0; // effective cell diameter (in units 10 microns)
          if  (translational_diffusion_constant >0.0)
              node_radius = p_random_force->GetDiffusionScalingConstant()/translational_diffusion_constant;
          else
              node_radius = set_node_radius;
          translational_diffusion_constant = p_random_force->GetDiffusionScalingConstant()/node_radius;

          p_random_force->SetHasBrownianRandomForce(has_brownian_random_force);
          p_random_force->SetUseTheSameNodeRadius(use_the_same_node_radius);
          p_random_force->SetTheSameNodeRadius(node_radius);

          p_random_force->SetHasPolarity(has_polarity);

          p_random_force->SetIfEquilibrateForAWhile(if_equilibrate_for_a_while);
          p_random_force->SetEndTimeForEquilibrium(time_for_equilibrium);

          simulator.AddForce(p_random_force);
        }
      /*------------------------------------END: Brownian Random Force ---------------------------*/


      /*------------------------------------START: CellKiller---------------------------------------*/
        c_vector<double, 2> point_vec = zero_vector<double>(2);
        double kill_height = 10000.0;
        point_vec[1] = kill_height;
        c_vector<double, 2> norm_vec = zero_vector<double>(2);
        norm_vec[1] = 1.0;

        bool kill_cells_group_number_from1 = true;

        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, point_vec, norm_vec) );
        p_killer->SetKillCellsGroupNumberFrom1(kill_cells_group_number_from1);
        simulator.AddCellKiller(p_killer);
      /*------------------------------------END: CellKiller---------------------------------------*/


      /*----------------------------------START: Boundary condition-------------------------------*/
        c_vector<double,2> point1 = zero_vector<double>(2);
        c_vector<double,2> normal1 = zero_vector<double>(2);
        normal1(1) = -1.0;
        double stop_time1 = time_for_equilibrium;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point1, normal1, stop_time1));

        c_vector<double,2> point2 = zero_vector<double>(2);
        c_vector<double,2> normal2 = zero_vector<double>(2);
        point2(1) = (1.5*num_ele_up + 0.5)*sqrt(initial_area/(3*sqrt(3)/2));
        normal2(1) = 1.0;
        double stop_time2 = time_for_equilibrium;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point2, normal2, stop_time2));

        simulator.AddCellPopulationBoundaryCondition(p_bc1);
        simulator.AddCellPopulationBoundaryCondition(p_bc2);
      /*--------------------------------END: Boundary condition-----------------------------*/
    

      /*--------------------------START: Output Directory and Simulation Information File---------------------*/
        // Output directory:
        std::ostringstream oss;
        time_t raw_time = time(0);
        struct tm * now = localtime(& raw_time);

        std::string output_directory = "aSF/Date: ";
        oss.str("");
        oss << (now->tm_year + 1900 -2000) << '-' << (now->tm_mon + 1) << '-' <<  now->tm_mday << '/';
        output_directory += oss.str();

        oss.str("");
        oss << "Timestamp=" << now->tm_hour << ':' << now->tm_min << ':' << now->tm_sec;

        oss << "_NumUp=" << num_ele_up;
        oss << "_NumAc=" << num_ele_across;

        oss << "_KA=" << ((area_elastic_modulus>=0.01 || area_elastic_modulus==0.0)? std::fixed : std::scientific) 
                << setprecision(2) << area_elastic_modulus;
        oss << "_KL=" << ((edge_elastic_modulus>=0.01 || edge_elastic_modulus==0.0)? std::fixed : std::scientific) 
                << setprecision(2) << edge_elastic_modulus;
        oss << "_Gamma=" << ((fabs(cell_cell_adhesion_energy_density)>=0.01 || fabs(cell_cell_adhesion_energy_density)==0.0)? std::fixed : std::scientific) 
                << setprecision(2) << cell_cell_adhesion_energy_density;

        if (multiple_leading_cells)
        oss << "_MF=" << std::fixed  << setprecision(0) << pulling_force_on_leading_cell;

        if (if_equilibrate_for_a_while)
        {
          oss << "_eqtime=" << std::fixed << setprecision(1) << time_for_equilibrium;
        }

        oss << "_Fp=" << ((polarity_magnitude>=0.01 || polarity_magnitude==0.0)? std::fixed : std::scientific) << setprecision(2) << polarity_magnitude;
        if (polarity_magnitude!=0.0)
        {
          oss << "_Dr=" << ((rotational_diffusion_constant>=0.01)? std::fixed : std::scientific) << setprecision(2) << rotational_diffusion_constant;
          if (seed_manually)
            oss << "_PSeed=" << seed_for_initial_random_polarity;
          else
            oss << "_PSeed=NA";
        }      

        oss << "_Dt=" << std::scientific << setprecision(1) << dt;
        if (apply_adaptive_timestep)
          oss << "_MaxMv=" << ((max_movement_per_timestep>=0.01)? std::fixed : std::scientific) << setprecision(3) << max_movement_per_timestep;
        oss << "_T1Thresh=" << ((cell_rearrangement_threshold>=0.01)? std::fixed : std::scientific) << setprecision(3) << cell_rearrangement_threshold;

        output_directory += oss.str();
        std::string concise_output_directory = output_directory;

        bool omit_file_name_results_from_time_X = true;
        bool output_simulatin_information_to_file = true;
        simulator.SetOmitFileNameResultsFromTimeX(omit_file_name_results_from_time_X);
        simulator.SetOutputSimulationInformationToFile(output_simulatin_information_to_file);
        if (output_simulatin_information_to_file)
          simulator.InputSimulationInformation(output_directory);

        if (use_concise_output_directory)
        {
          simulator.SetOutputDirectory(concise_output_directory);
          std::cout << std::endl << "Concise output directory is set: " << concise_output_directory << std::endl;
        }
        else
        {
          std::cout << std::endl << "Output directory is not set. Please set up an output directory!" << std::endl;
          EXCEPTION("Output directory is not set.");
        }
      /*--------------------------END: Output Directory and Simulation Information File---------------------*/

        simulator.Solve();
    }

};

#endif /* APICALSTRESSFIBERS_HPP_ */
