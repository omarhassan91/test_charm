#ifndef OPENMM_MONTECARLOBAROSTATIMPL2_H_
#define OPENMM_MONTECARLOBAROSTATIMPL2_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "openmm/internal/ForceImpl.h"
#include "openmm/MonteCarloBarostat2.h"
#include "openmm/Kernel.h"
#include <sfmt/SFMT.h>
#include <string>

namespace OpenMM {

/**
 * This is the internal implementation of MonteCarloBarostat2.
 */
class MonteCarloBarostatImpl2 : public ForceImpl {
public:
    MonteCarloBarostatImpl2(const MonteCarloBarostat2& owner);
    void initialize(ContextImpl& context);
    virtual const MonteCarloBarostat2& getOwner() const {
        return owner;
    }
    void updateContextState(ContextImpl& context);
    double calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
        // This force doesn't apply forces to particles.
        return 0.0;
    }
    std::map<std::string, double> getDefaultParameters();
    std::vector<std::string> getKernelNames();
private:
    const MonteCarloBarostat2& owner;
    int step, numAttempted, numAccepted;
    int totalAccepted;
    /* 3 Dimensions of scales are needed for different baro states */
    Vec3 lengthScale;
    OpenMM_SFMT::SFMT random;
    Kernel kernel;
};

} // namespace OpenMM

#endif

