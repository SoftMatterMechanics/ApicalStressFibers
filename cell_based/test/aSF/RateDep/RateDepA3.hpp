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
#include "MyNoPBCToroidalHoneycombVertexMeshGenerator.hpp"
// #include "MyToroidalHoneycombVertexMeshGenerator.hpp"
// #include "NagaiHondaForce.hpp"
#include "MyNagaiHondaForce.hpp"
// #include "MyMorphogeneticForce.hpp"
#include "MyStressfiberTensionForce.hpp"
#include "DiffusionForce.hpp"
#include "PlaneBasedCellKiller.hpp"

// #include <ctime>

class ApicalStressFibers : public AbstractCellBasedTestSuite
{
public:

    void TestApicalStressFibers()
    {
      /*-------------------------START: Basic Settings-----------------------*/
      /* Energy equation form: 1/2*Ka*(A-A0)^2 + 1/2*Kp*(P-P0)^2 + Gamma*L */

        double target_shape_index = 3.5;
        double cell_cell_adhesion_energy_density = -0.3;  // Gamma, this parameter consists of cell-cell adhesion and cortical contraction 
        double cell_boundary_adhesion_energy_density = -0.3;  // Gamma at boundary
        double polarity_magnitude_before_equilibrium = 0.04;  // for before equilibrium
        unsigned random_seed_for_target_area = 13;
        // 5. Stress fiber tension
        double sf_stiffness = 0.0;   // 0.0 means the stiffness of stress fibers is equal to the perimeter stiffness
        double nucleation_perimeter_tension = 0.0;  // threshold for stress fibers nucleation
        double rest_length_of_nucleation = 0.19;   // delta0 = 0.02;
        double adhesion_energy = 10;   // cortex-membrane adhesion energy
        double cyto_viscosity = 1;  // viscosity of cytoplasma

      // 1. Cell mesh and size
        unsigned num_ele_across = 16; // cell number along anterior-posterior, must be an even number
        unsigned num_ele_up = 16; // cell number along medial-lateral, must be an even number
        // double target_shape_index = 3.0;
        bool   seed_manually = true;
        // unsigned random_seed_for_target_area = 1;
        double min_target_area = 1.0;
        double max_target_area = 7.0;
        bool   use_fixed_target_area_without_modifier = false;
        double target_area = (min_target_area + max_target_area)/2; // A0, target areas are not uniform
        double initial_area = (min_target_area + max_target_area)/2; // set the initial areas to be uniform, A = 3*sqrt(3)/2*l^2
        double center_y_coordination = 0.5*(1.5*num_ele_up+0.5)*sqrt(initial_area/(3*sqrt(3)/2));
        double half_width = num_ele_across/2*sqrt(initial_area/(3*sqrt(3)/2))*sqrt(3)/2*2;

      // 2. Area elasticity
        double area_elastic_modulus = 1.0; // Ka

      // 4. Cell-cell adhesion & constant cortical contraction
      // this parameter value here is one half of its real value because the edge is shared by two cells
        // double cell_cell_adhesion_energy_density = -0.3;  // Gamma, this parameter consists of cell-cell adhesion and cortical contraction 
        // double cell_boundary_adhesion_energy_density = -0.3;  // Gamma at boundary
        bool   if_use_face_element_to_get_adhesion_parameter = false;

      // 5. Stress fiber tension


      // 6. morphogenetic force
        double horizontal_morphogenetic_force = 5.5;
        double vertical_morphogenetic_force = 2;
        double horizontal_morphogenetic_force_growth_rate = 0.001;
        double vertical_morphogenetic_force_growth_rate = 0.02;

      // 7. Random force
        bool   add_random_force = true;
        bool   has_brownian_random_force = false; // brownian random force is used in cell center model
        double translational_diffusion_constant = 0.0;
        double set_node_radius = 2.0; // effective cell diameter (in units of 10 microns)

        bool   has_polarity = true;
        unsigned seed_for_initial_random_polarity = 3;
        // double polarity_magnitude_before_equilibrium = 0.04;  // for before equilibrium
        double polarity_magnitude_after_equilibrium = 0.0;  // for after equilibrium
        double rotational_diffusion_constant = 0.5;

      // 8. Time
        bool   if_equilibrate_for_a_while = true;
        double time_for_rest = 0;
        double time_for_random_movement = 250.0;
        double time_for_relaxation = 150.0;
        double time_for_equilibrium = time_for_rest + time_for_random_movement + time_for_relaxation;
        if (time_for_equilibrium <= 0.0)
           if_equilibrate_for_a_while = false;

        double dt = 0.05;
        double real_equilibrium_time = time_for_equilibrium + vertical_morphogenetic_force/vertical_morphogenetic_force_growth_rate;
        double start_time_for_stretching = real_equilibrium_time;
        double end_time = real_equilibrium_time + (horizontal_morphogenetic_force - vertical_morphogenetic_force)/horizontal_morphogenetic_force_growth_rate; 
        double max_movement_per_timestep = 0.05; 
        bool   apply_adaptive_timestep = true;
        double sampling_time = 1.0;
        unsigned sampling_timestep_multiple = (unsigned) round(sampling_time/dt);
        
      // 9. Cell rearrangement threshold length for T1 & T2 transitions
        double cell_rearrangement_threshold = 0.01; 
        double t2_threshold = 0.001;
        double t3_threshold = 5.0;

      // 10. Output & display
        bool   if_update_face_elements_in_mesh = true;
        bool   output_concise_swap_information_when_remesh = false;
        bool   output_detailed_swap_information_when_remesh = false;
        bool   output_numerical_method_information = false;
        // bool   output_information_for_nagai_honda_force = false;     
        // bool   if_set_cell_data_of_detailed_force_contributions = false; // for output and display each component in Paraview
        bool   use_concise_output_directory = true;
      /*------------------------------END: Basic Settings----------------------------*/


      /*-----------------------START: Generate cell monolayer mesh-------------------*/
        MyNoPBCToroidalHoneycombVertexMeshGenerator generator(num_ele_across, num_ele_up, initial_area, cell_rearrangement_threshold, t2_threshold);
        MyNoPBCToroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();

        p_mesh->SetDistanceForT3SwapChecking(t3_threshold); 
        p_mesh->SetUpdateFaceElementsInMeshBoolean(if_update_face_elements_in_mesh);
        p_mesh->SetAreaSeed(random_seed_for_target_area);
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
        // cells_generator.SetPerimeterElasticityParameter(edge_elastic_modulus);

        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);
        
        bool restrict_vertex_movement = false;
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetRestrictVertexMovementBoolean(restrict_vertex_movement);
      /*------------------------------END: Cells and Cell population--------------------------*/


      /*---------------------------------START: Simulator settings----------------------------*/
        OffLatticeSimulation<2> simulator(cell_population);

        // output cell velocity:
        bool run_with_birth = false;
        bool output_cell_velocity = true;
        bool my_output_cell_velocity = true;
        bool output_cell_elongation = true;

        simulator.SetNoBirth(!run_with_birth);
        simulator.SetOutputCellVelocities(output_cell_velocity);
        simulator.SetMyOutputCellVelocities(my_output_cell_velocity);
        if (my_output_cell_velocity && seed_manually)
          simulator.SetAreaSeed(random_seed_for_target_area);
        simulator.SetOutputCellElongation(output_cell_elongation);

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
        p_numerical_method->SetCenterYCoordination(center_y_coordination);
        p_numerical_method->SetIfEquilibrateForAWhile(if_equilibrate_for_a_while);
        p_numerical_method->SetTimeForEquilibrium(time_for_equilibrium);
        p_numerical_method->SetRealEquilibriumTime(real_equilibrium_time);
        p_numerical_method->SetMorphogeneticForceGrowthRate(horizontal_morphogenetic_force_growth_rate, vertical_morphogenetic_force_growth_rate);
        simulator.SetNumericalMethod(p_numerical_method);
      /*---------------------------------END: Simulator settings-----------------------------*/


      /*--------------------------------START: Modifier Settings-----------------------------*/
        // Polarity modifier
        if (has_polarity)
        {
          MAKE_PTR_ARGS(PolarityModifier<2>, p_polarity_modifier, ());
          p_polarity_modifier->SetPolarityMagnitudeAfterEquilibrium(polarity_magnitude_after_equilibrium);
          p_polarity_modifier->SetPolarityMagnitudeBeforeEquilibrium(polarity_magnitude_before_equilibrium);
          p_polarity_modifier->SetSeedManually(seed_manually);
          p_polarity_modifier->SetSeedForInitialRandomPolarity(seed_for_initial_random_polarity);
          p_polarity_modifier->SetD(rotational_diffusion_constant);
          simulator.AddSimulationModifier(p_polarity_modifier);
        }
      /*----------------------------------END: Modifier Settings------------------------------*/


      // /*---------------------------------START: MyMorphogeneticForce-----------------------------*/
      //   MAKE_PTR(MyMorphogeneticForce<2>, p_morphogenetic_force);
        
      //   p_morphogenetic_force->SetAddPullingForceEvenlyOnNodesOfLeadingCell(add_pulling_force_evenly_on_nodes_of_leading_cell);
      //   p_morphogenetic_force->SetPullingForceOnLeadingCell(pulling_force_on_leading_cell);
      //   p_morphogenetic_force->SetCenterYCoordination(center_y_coordination);
      //   p_morphogenetic_force->SetIfEquilibrateForAWhile(if_equilibrate_for_a_while);
      //   p_morphogenetic_force->SetEndTimeForEquilibrium(time_for_equilibrium);

      //   simulator.AddForce(p_morphogenetic_force);
      // /*-------------------------------------END: MyNagaiHondaForce------------------------------*/


      /*---------------------------------START: MyNagaiHondaForce-----------------------------*/
        MAKE_PTR(MyNagaiHondaForce<2>, p_nh_force);
        
        p_nh_force->SetNagaiHondaDeformationEnergyParameter(area_elastic_modulus); // KA
        // p_nh_force->SetNagaiHondaMembraneSurfaceEnergyParameter(edge_elastic_modulus); // KL
        p_nh_force->SetNagaiHondaCellCellAdhesionEnergyParameter(cell_cell_adhesion_energy_density); // Gamma
        p_nh_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(cell_boundary_adhesion_energy_density); // Gamma at boundary

        p_nh_force->SetUseFixedTargetArea(use_fixed_target_area_without_modifier); // used in the case where there is no target area modifier! (no division)
        p_nh_force->SetFixedTargetArea(target_area); // to be determined
        p_nh_force->SetTargetShapeIndex(target_shape_index);
        // p_nh_force->SetFixedTargetPerimeter(target_perimeter);
        p_nh_force->SetTimeForRest(time_for_rest);
        p_nh_force->SetUseFaceElementToGetAdhesionParameterBoolean(if_use_face_element_to_get_adhesion_parameter);

        // p_force->SetOutputInformationForNagaiHondaForce(output_information_for_nagai_honda_force);
        simulator.AddForce(p_nh_force);
      /*-------------------------------------END: MyNagaiHondaForce------------------------------*/


      /*---------------------------------START: My Stressfiber Tension Force-----------------------------*/
        MAKE_PTR(MyStressfiberTensionForce<2>, p_sf_force);
        
        p_sf_force->SetAreaSeed(random_seed_for_target_area);
        p_sf_force->SetIfEquilibrateForAWhile(if_equilibrate_for_a_while);
        p_sf_force->SetStartTimeForStretching(start_time_for_stretching);
        p_sf_force->SetFlagForStressfiberCreation(0);
        p_sf_force->SetStressfiberStiffness(sf_stiffness);
        p_sf_force->SetNucleationThresholdOfPerimeterTension(nucleation_perimeter_tension);
        p_sf_force->SetRestLengthOfNucleation(rest_length_of_nucleation);
        p_sf_force->SetPeelingParameters(adhesion_energy, cyto_viscosity);
        p_sf_force->SetNagaiHondaCellCellAdhesionEnergyParameter(cell_cell_adhesion_energy_density);
        // p_sf_force->SetHalfWidth(half_width);

        // simulator.AddForce(p_sf_force);
      /*-------------------------------------END: My Stressfiber Tension Force------------------------------*/


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
          p_random_force->SetStartTimeForRandom(time_for_rest);
          p_random_force->SetEndTimeForRandom(time_for_rest + time_for_random_movement);

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
        double stop_time1 = time_for_rest + time_for_random_movement;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point1, normal1, stop_time1));

        c_vector<double,2> point2 = zero_vector<double>(2);
        c_vector<double,2> normal2 = zero_vector<double>(2);
        point2(1) = (1.5*num_ele_up + 0.5)*sqrt(initial_area/(3*sqrt(3)/2));
        normal2(1) = 1.0;
        double stop_time2 = time_for_rest + time_for_random_movement;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point2, normal2, stop_time2));

        c_vector<double,2> point3 = zero_vector<double>(2);
        c_vector<double,2> normal3 = zero_vector<double>(2);
        point3(0) = -half_width*1.0;
        normal3(0) = -1.0;
        double stop_time3 = time_for_rest + time_for_random_movement;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point3, normal3, stop_time3));

        c_vector<double,2> point4 = zero_vector<double>(2);
        c_vector<double,2> normal4 = zero_vector<double>(2);
        point4(0) = (half_width + sqrt(initial_area/(3*sqrt(3)/2))*sqrt(3)/2)*1.0;
        normal4(0) = 1.0;
        double stop_time4 = time_for_rest + time_for_random_movement;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point4, normal4, stop_time4));        

        simulator.AddCellPopulationBoundaryCondition(p_bc1);
        simulator.AddCellPopulationBoundaryCondition(p_bc2);
        simulator.AddCellPopulationBoundaryCondition(p_bc3);
        simulator.AddCellPopulationBoundaryCondition(p_bc4);        
      /*--------------------------------END: Boundary condition-----------------------------*/
    

      /*--------------------------START: Output Directory and Simulation Information File---------------------*/
        // Output directory:
        std::ostringstream oss;
        time_t raw_time = time(0);
        struct tm * now = localtime(& raw_time);

        std::string output_directory = "aSF/RateDep/";
        oss.str("");
        oss << "p0=" << std::fixed  << setprecision(1) << target_shape_index << ",";
        oss << "Gamma=" << ((fabs(cell_cell_adhesion_energy_density)>=0.01 || fabs(cell_cell_adhesion_energy_density)==0.0)? std::fixed : std::scientific) 
                << setprecision(1) << cell_cell_adhesion_energy_density << "/";
        output_directory += oss.str();
        
        oss.str("");
        oss << "Date=" << (now->tm_year + 1900 -2000) << '-' << (now->tm_mon + 1) << '-' <<  now->tm_mday;
        oss << "_Time=" << now->tm_hour << ':' << now->tm_min << ':' << now->tm_sec;
        // oss << "_NumUp=" << num_ele_up;
        // oss << "_NumAc=" << num_ele_across;
        // oss << "_Ka=" << ((area_elastic_modulus>=0.01 || area_elastic_modulus==0.0)? std::fixed : std::scientific) 
        //         << setprecision(2) << area_elastic_modulus;
        // oss << "_Kp=" << ((edge_elastic_modulus>=0.01 || edge_elastic_modulus==0.0)? std::fixed : std::scientific) 
        //         << setprecision(2) << edge_elastic_modulus;

        // oss << "_Ai=" << std::fixed  << setprecision(2) << initial_area;
        // oss << "_minA0=" << std::fixed  << setprecision(2) << min_target_area;
        // oss << "_maxA0=" << std::fixed  << setprecision(2) << max_target_area;
        oss << "_Aseed=" << random_seed_for_target_area;

        // oss << "_fp=" << std::fixed  << setprecision(3) << polarity_magnitude_before_equilibrium;
        // oss << "_Pseed=" << seed_for_initial_random_polarity;

        // oss << "_dt=" << std::fixed  << setprecision(3) << dt;
        // oss << "_maxd=" << std::fixed  << setprecision(3) << max_movement_per_timestep;
        // if (if_equilibrate_for_a_while)
        // {
        //    oss << "_eqtime=" << std::fixed << setprecision(1) << time_for_equilibrium;
        // }
        oss << "_realeqtime=" << std::fixed << setprecision(1) << real_equilibrium_time;
        oss << "_simtime=" << std::fixed << setprecision(1) << end_time; 
        
        oss << "_Fx=" << std::fixed << setprecision(1) << horizontal_morphogenetic_force;
        oss << "_Fy=" << std::fixed << setprecision(1) << vertical_morphogenetic_force;
        oss << "_vFx=" << std::fixed << setprecision(4) << horizontal_morphogenetic_force_growth_rate;
        oss << "_vFy=" << std::fixed << setprecision(2) << vertical_morphogenetic_force_growth_rate;

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
        p_mesh->SetOutputDirectory(concise_output_directory);

        p_sf_force->SetOutputDirectory(concise_output_directory);
        simulator.AddForce(p_sf_force);

        simulator.Solve();
    }

};

#endif /* APICALSTRESSFIBERS_HPP_ */
