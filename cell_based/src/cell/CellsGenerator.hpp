/*

Copyright (c) 2005-2020, University of Oxford.
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

#ifndef CELLSGENERATOR_HPP_
#define CELLSGENERATOR_HPP_

#include <vector>
#include "Cell.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "RandomNumberGenerator.hpp"
#include "BernoulliTrialCellCycleModel.hpp"

/**
 * A helper class for generating a vector of cells for a given mesh.
 *
 * It is templated over types of cell-cycle model and spatial dimension.
 */
template<class CELL_CYCLE_MODEL, unsigned DIM>
class CellsGenerator
{
private:

    bool mUseARandomDivisionRule;
    unsigned mDivideTimestepMultiple;
    bool mUseBernoulliTrialCellCycleModel;
    double mDivisionProbability;
    double mMinimumDivisionAge;
    unsigned mRandomSeed;
    double mMinTargetArea;
    double mMaxTargetArea;
    double mPerimeterElasticity;

public:

    // CellsGenerator();

    // ~CellsGenerator();

    /**
     * Fills a vector of cells with a specified cell-cycle model, to match
     * a given number of cells. Gives them birth times of 0 for node 0,
     * -1 for node 1, -2 for node 2 etc...
     *
     * @param rCells  An empty vector of cells to fill up.
     * @param numCells  The number of cells to generate.
     * @param locationIndices is used when a birth-time hint is needed for individual cell.
     *             Defaults to an empty vector -- otherwise must be of length numCells
     * @param pCellProliferativeType  the cell proliferative type to give each cell (defaults to StemCellProliferativeType)
     */
    void GenerateBasic(std::vector<CellPtr>& rCells,
                       unsigned numCells,
                       const std::vector<unsigned> locationIndices=std::vector<unsigned>(),
                       boost::shared_ptr<AbstractCellProperty> pCellProliferativeType=boost::shared_ptr<AbstractCellProperty>());

    /**
     * Fills a vector of cells with a specified cell-cycle model, to match
     * a given number of cells. Gives cells a random birth time drawn uniformly
     * from 0 to the AverageStemCellCycleTime.
     *
     * @param rCells  An empty vector of cells to fill up.
     * @param numCells  The number of cells to generate.
     * @param pCellProliferativeType  the cell proliferative type to give each cell (defaults to StemCellProliferativeType)
     */
    void GenerateBasicRandom(std::vector<CellPtr>& rCells,
                             unsigned numCells,
                             boost::shared_ptr<AbstractCellProperty> pCellProliferativeType=boost::shared_ptr<AbstractCellProperty>());

    /**
     * Fills a vector of cells with birth times to match a given vector of location indices.
     *
     * @param rCells  An empty vector of cells to fill up.
     * @param locationIndices  The indices of the cell population to assign real cells to.
     * @param pCellProliferativeType  the cell proliferative type to give each cell (defaults to StemCellProliferativeType)
     */
    void GenerateGivenLocationIndices(std::vector<CellPtr>& rCells,
                                      const std::vector<unsigned> locationIndices,
                                      boost::shared_ptr<AbstractCellProperty> pCellProliferativeType=boost::shared_ptr<AbstractCellProperty>());
    
    void SetUseBernoulliTrialCellCycleModel(bool useBernoulliTrialCellCycleModel)
    {
        mUseBernoulliTrialCellCycleModel = useBernoulliTrialCellCycleModel;
    }

    void SetDivisionProbability(double divisionProbability)
    {
        mDivisionProbability = divisionProbability;
    }

    void SetMinimumDivisionAge(double minimumDivisionAge)
    {
        mMinimumDivisionAge = minimumDivisionAge;
    }

    void SetRandomSeedForTargetAreas(unsigned randomSeed)
    {
        mRandomSeed = randomSeed;
    }

    void SetLimitsOfTargetAreas(double minTargetArea, double maxTargetArea)
    {
        mMinTargetArea = minTargetArea;
        mMaxTargetArea = maxTargetArea;
    }

    void SetPerimeterElasticityParameter(double perimeterElasticity)  
    {
        mPerimeterElasticity = perimeterElasticity;
    }                                      
                                                                        
};

template<class CELL_CYCLE_MODEL, unsigned DIM>
void CellsGenerator<CELL_CYCLE_MODEL,DIM>::GenerateBasic(std::vector<CellPtr>& rCells,
                                                         unsigned numCells,
                                                         const std::vector<unsigned> locationIndices,
                                                         boost::shared_ptr<AbstractCellProperty> pCellProliferativeType)
{
    rCells.clear();

    if (!locationIndices.empty())
    {
        // If location indices is given, then it needs to match the number of output cells
        if (numCells != locationIndices.size())
        {
            EXCEPTION("The size of the locationIndices vector must match the required number of output cells");
        }
    }
    rCells.reserve(numCells);

    // Create cells
    for (unsigned i=0; i<numCells; i++)
    {
        CELL_CYCLE_MODEL* p_cell_cycle_model = new CELL_CYCLE_MODEL;
        p_cell_cycle_model->SetDimension(DIM);

        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));

        if (!pCellProliferativeType)
        {
            p_cell->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
        }
        else
        {
            p_cell->SetCellProliferativeType(pCellProliferativeType);
        }

        double birth_time;
        if (!locationIndices.empty())
        {
            birth_time = 0.0 - locationIndices[i];
        }
        else
        {
            birth_time = 0.0 - i;
        }

        p_cell->SetBirthTime(birth_time);
        rCells.push_back(p_cell);
    }
}

template<class CELL_CYCLE_MODEL, unsigned DIM>
void CellsGenerator<CELL_CYCLE_MODEL,DIM>::GenerateBasicRandom(std::vector<CellPtr>& rCells,
                                                               unsigned numCells,
                                                               boost::shared_ptr<AbstractCellProperty> pCellProliferativeType)
{
    rCells.clear();

    rCells.reserve(numCells);

    RandomNumberGenerator::Instance()->Reseed(mRandomSeed);

    // Create cells
    for (unsigned i=0; i<numCells; i++)
    {
        CELL_CYCLE_MODEL* p_cell_cycle_model = new CELL_CYCLE_MODEL;
        p_cell_cycle_model->SetDimension(DIM);

        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));

        if (!pCellProliferativeType)
        {
            p_cell->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
        }
        else
        {
            p_cell->SetCellProliferativeType(pCellProliferativeType);
        }

        double birth_time = -p_cell_cycle_model->GetAverageStemCellCycleTime()*RandomNumberGenerator::Instance()->ranf();

        if (p_cell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
        {
            birth_time = -p_cell_cycle_model->GetAverageTransitCellCycleTime()*RandomNumberGenerator::Instance()->ranf();
        }

        p_cell->SetBirthTime(birth_time);

        // my changes
        if (mUseBernoulliTrialCellCycleModel)
        {
            if (dynamic_cast<BernoulliTrialCellCycleModel*>(p_cell->GetCellCycleModel())==nullptr)
                std::cout << std::endl << "ERROR: Please use BernoulliTrialCellCycleModel!" << std::endl;
            assert(dynamic_cast<BernoulliTrialCellCycleModel*>(p_cell->GetCellCycleModel())!=nullptr);
            static_cast<BernoulliTrialCellCycleModel*>(p_cell->GetCellCycleModel())->SetDivisionProbability(mDivisionProbability);
            static_cast<BernoulliTrialCellCycleModel*>(p_cell->GetCellCycleModel())->SetMinimumDivisionAge(mMinimumDivisionAge);
        }

        // my changes (to generate cells with random target areas in the range [min_target_area, max_target_area])
        // srand((unsigned)time(NULL));// if srand() used in main time loop, we may not need it here. Try later!
        double random_number = RandomNumberGenerator::Instance()->ranf();
        double cell_target_area = mMinTargetArea + (mMaxTargetArea - mMinTargetArea)*random_number;
        // double random_number = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
        // double cell_target_area = (mMinTargetArea+mMaxTargetArea)/2 + (mMaxTargetArea - mMinTargetArea)/2*1/3*random_number;
        // std::cout << "seed=" << mRandomSeed << ", random=" << random_number << std::endl;                                                                      
        p_cell->GetCellData()->SetItem("target_area", cell_target_area);

        // p_cell->GetCellData()->SetItem("perimeter_elasticity", mPerimeterElasticity);

        rCells.push_back(p_cell);
    }
}

template<class CELL_CYCLE_MODEL, unsigned DIM>
void CellsGenerator<CELL_CYCLE_MODEL,DIM>::GenerateGivenLocationIndices(std::vector<CellPtr>& rCells,
                                                                        const std::vector<unsigned> locationIndices,
                                                                        boost::shared_ptr<AbstractCellProperty> pCellProliferativeType)
{
    assert(!locationIndices.empty());

    unsigned num_cells = locationIndices.size();

    rCells.clear();
    rCells.reserve(num_cells);
    CellPropertyRegistry::Instance()->Clear();

    for (unsigned i=0; i<num_cells; i++)
    {
        CELL_CYCLE_MODEL* p_cell_cycle_model = new CELL_CYCLE_MODEL;
        p_cell_cycle_model->SetDimension(DIM);

        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

        CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));

        if (!pCellProliferativeType)
        {
            p_cell->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
        }
        else
        {
            p_cell->SetCellProliferativeType(pCellProliferativeType);
        }

        double birth_time = 0.0 - locationIndices[i];
        p_cell->SetBirthTime(birth_time);
        rCells.push_back(p_cell);
    }
}

#endif /* CELLSGENERATOR_HPP_ */
