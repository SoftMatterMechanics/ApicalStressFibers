/*

Copyright (c) 2005-2013, University of Oxford.
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

#ifndef TESTABSTRACTCONTRACTIONCELLFACTORY_HPP_
#define TESTABSTRACTCONTRACTIONCELLFACTORY_HPP_


#include <cxxtest/TestSuite.h>

#include "LabelBasedContractionCellFactory.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestAbstractContractionCellFactory : public CxxTest::TestSuite
{
public:
    void TestContractionCellFactory() throw (Exception)
    {
        {
            LabelBasedContractionCellFactory<2> factory(CONSTANT);

            ConstantActiveTension* p_model =
                dynamic_cast<ConstantActiveTension*> (factory.CreateContractionCellForElement(0u));
            TS_ASSERT(p_model);
            delete p_model;
        }

        {
            LabelBasedContractionCellFactory<2> factory(NONPHYSIOL1);

            NonPhysiologicalContractionModel* p_model =
                dynamic_cast<NonPhysiologicalContractionModel*> (factory.CreateContractionCellForElement(0u));
            TS_ASSERT(p_model);
            delete p_model;
        }

        {
            LabelBasedContractionCellFactory<2> factory(NASH2004);

            Nash2004ContractionModel* p_model =
                dynamic_cast<Nash2004ContractionModel*> (factory.CreateContractionCellForElement(0u));
            TS_ASSERT(p_model);
            delete p_model;
        }

        {
            LabelBasedContractionCellFactory<2> factory(KERCHOFFS2003);

            Kerchoffs2003ContractionModel* p_model =
                dynamic_cast<Kerchoffs2003ContractionModel*> (factory.CreateContractionCellForElement(0u));
            TS_ASSERT(p_model);
            delete p_model;
        }
    }
};



#endif /* TESTABSTRACTCONTRACTIONCELLFACTORY_HPP_ */
