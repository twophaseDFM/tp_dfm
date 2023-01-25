// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief Utility function to compute the velocities at
 *        the cell centers for a given grid and solution vector.
 */
#ifndef DUMUX_COMPUTEVELOCITIES_HH
#define DUMUX_COMPUTEVELOCITIES_HH

#include <dune/common/indices.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/discretization/elementsolution.hh>


namespace Dumux {

//! Function to compute velocities in a single-phase context
template<class TypeTag, std::size_t id,
         class Assembler, class CouplingManager,
         class FVGridGeometry, class GridVariables,
         class SolutionVector, class VelocityVector>
void computeVelocities(Dune::index_constant<id> domainId,
                       const Assembler& assembler,
                       CouplingManager& couplingManager,
                       const FVGridGeometry& fvGridGeometry,
                       const GridVariables& gridVars,
                       const SolutionVector& x,
                       VelocityVector& velocities)
{
    velocities.resize(fvGridGeometry.gridView().size(0));
    for (const auto& element : elements(fvGridGeometry.gridView()))
    {
        auto fvGeometry = localView(fvGridGeometry);
        auto elemVolVars = localView(gridVars.curGridVolVars());
        auto elemFluxVarsCache = localView(gridVars.gridFluxVarsCache());

        couplingManager.bindCouplingContext(domainId, element, assembler);
        fvGeometry.bind(element);
        elemVolVars.bind(element, fvGeometry, x);

        if (FVGridGeometry::discMethod != DiscretizationMethod::box)
        {
            elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);
            GetPropType<TypeTag, Properties::VelocityOutput> velOutput(gridVars);
            velOutput.calculateVelocity(velocities, element, fvGeometry, elemVolVars, elemFluxVarsCache, 1);
        }
        else
        {
            double volume = 0.0;
            typename VelocityVector::value_type vel(0.0);
            const auto elemSol = elementSolution(element, x, fvGridGeometry);
            for (const auto& scv : scvs(fvGeometry))
            {
                auto velScv = evalGradients(element,
                                            element.geometry(),
                                            fvGridGeometry,
                                            elemSol,
                                            scv.center())[0];

                velScv.axpy(-1.0*elemVolVars[scv].density(0),
                            couplingManager.problem(domainId).spatialParams().gravity(scv.center()));

                velScv *= elemVolVars[scv].permeability();
                velScv *= elemVolVars[scv].mobility(1);
                velScv *= -1.0;
                velScv *= scv.volume();
                vel += velScv;
                volume += scv.volume();
            }
            vel /= volume;
            velocities[fvGridGeometry.elementMapper().index(element)] = vel;
        }
    }
}

} // end namespace Dumux

#endif
