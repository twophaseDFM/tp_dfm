// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup MultiDomain
 * \ingroup MultiDomainFacet
 * \ingroup TwoPTests
 * \brief The spatial parameters for the fracture sub-domain in the exercise
 *        on two-phase flow in fractured porous media.
 */
#ifndef DUMUX_COURSE_FRACTURESEXERCISE_FRACTURE_SPATIALPARAMS_HH
#define DUMUX_COURSE_FRACTURESEXERCISE_FRACTURE_SPATIALPARAMS_HH

#include <dumux/io/grid/griddata.hh>

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>
#include "fractureproblem.hh"
#include <dumux/porousmediumflow/2p/model.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup MultiDomainFacet
 * \ingroup TwoPTests
 * \brief The spatial params the two-phase facet coupling test
 */
template< class FVGridGeometry, class Scalar >
//template<class TypeTag>
class FractureSpatialParams
: public FVSpatialParams< FVGridGeometry, Scalar, FractureSpatialParams<FVGridGeometry, Scalar> >
//: public FVSpatialParams <GetPropType<TypeTag, Properties::GridGeometry>,
//  	  	  	  	  	  	  GetPropType<TypeTag, Properties::Scalar>,
//						  FractureSpatialParams<TypeTag>>
{
//	using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
//	using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ThisType = FractureSpatialParams< FVGridGeometry, Scalar >;
//	using ThisType = FractureSpatialParams<TypeTag>;
    using ParentType = FVSpatialParams< FVGridGeometry, Scalar, ThisType >;

    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using GridView = typename FVGridGeometry::GridView;
    using Grid = typename GridView::Grid;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PcKrSwCurve = FluidMatrix::VanGenuchtenDefault<Scalar>;

    static constexpr int dimWorld = GridView::dimensionworld;

//    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    //! export the type used for permeabilities
    using PermeabilityType = Scalar;

    //! the constructor
    FractureSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                          std::shared_ptr<const Dumux::GridData<Grid>> gridData,
                          const std::string& paramGroup)
    : ParentType(fvGridGeometry)
    , gridDataPtr_(gridData)
    , pcKrSwCurve_("Fracture.SpatialParams")
    , barrierPcKrSwCurve_("Fracture.SpatialParams.Barrier")
    {
        porosity_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Porosity");
        a1_ = (getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture1"));
        a2_ = (getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture2"));
        a3_ = (getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture3"));
        a4_ = (getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture4"));
        a5_ = (getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture5"));
        a6_ = (getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture6", 1e-3));
        kn_ = (getParamFromGroup<Scalar>(paramGroup, "Problem.Stiffness"));
        wte_ = (getParamFromGroup<Scalar>(paramGroup, "Problem.WettingPhaseThermalExpansionCoefficient"));
        nte_ = (getParamFromGroup<Scalar>(paramGroup, "Problem.NonWettingPhaseThermalExpansionCoefficient"));
        a_ = (getParamFromGroup<Scalar>(paramGroup, "Problem.a"));
        b_ = (getParamFromGroup<Scalar>(paramGroup, "Problem.b"));
    }

//    //! Function for defining the (intrinsic) permeability \f$[m^2]\f$.
//    template< class ElementSolution, class FluidState >
    template< class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
    	int pressureIdx = 0;
    	int saturationIdx = 1;
    	int temperatureIdx = 2;
        const auto& priVars = elemSol[scv.localDofIndex()];
        const auto fmi = fluidMatrixInteraction(element, scv, elemSol);
        const auto peff_ = priVars[saturationIdx] * (priVars[pressureIdx] + fmi.pc(1 - priVars[saturationIdx])) + (1 - priVars[saturationIdx]) * priVars[pressureIdx];

		const GlobalPosition& globalPos = scv.center();
		const auto domainHeight = 100.0 + 6000;
		Scalar densityW = 1000.0;
		const auto g = -9.81;
		const auto initialTemperature = 283.0 + (domainHeight - globalPos[dimWorld-1])*0.03;
		const auto initialPressure = 1e5 - (domainHeight - globalPos[dimWorld-1])*densityW*g;
		const auto ThermalExpan = (priVars[temperatureIdx] - initialTemperature) * (priVars[saturationIdx] * nte_ + (1 - priVars[saturationIdx]) * wte_);
		const auto deltaP = peff_ - initialPressure;
//        const auto deltaT_ = volVars.temperature() - initialTemperature;

        Scalar a1 = a1_ + a_ * deltaP/kn_ + b_ * ThermalExpan * a1_;
        Scalar a2 = a2_ + a_ * deltaP/kn_ + b_ * ThermalExpan * a2_;
        Scalar a3 = a3_ + a_ * deltaP/kn_ + b_ * ThermalExpan * a3_;
        Scalar a4 = a4_ + a_ * deltaP/kn_ + b_ * ThermalExpan * a4_;
        Scalar a5 = a5_ + a_ * deltaP/kn_ + b_ * ThermalExpan * a5_;
        Scalar a6 = a6_ + a_ * deltaP/kn_ + b_ * ThermalExpan * a6_;

        if (getElementDomainMarker(element) == 1)
        	return a1*a1/12;
        else if(getElementDomainMarker(element) == 2)
        	return a2*a2/12;
        else if(getElementDomainMarker(element) == 3)
        	return a3*a3/12;
        else if(getElementDomainMarker(element) == 4)
        	return a4*a4/12;
        else if(getElementDomainMarker(element) == 5)
        	return a5*a5/12;
        else
        	return a6*a6/12;
//    	return 1e-8;
    }

//    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
//    {
////    	if (globalPos[0] < 0 + eps_)
////    		return permeabilityWell_ ;
////    	else
//		return 1e-10;
//    }


    //! Return the porosity
    template< class ElementSolution >
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        return porosity_;
    }

    /*!
     * \brief Returns the fluid-matrix interaction law for the sub-control volume
     *
     * \param element The current finite element
     * \param scv The sub-control volume
     * \param elemSol The current element solution
     */
    template<class ElementSolution>
    auto fluidMatrixInteraction(const Element& element,
                                const SubControlVolume& scv,
                                const ElementSolution& elemSol) const
    {  
        return makeFluidMatrixInteraction(pcKrSwCurve_); 
    }

    //! Water is the wetting phase
    template< class FluidSystem >
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    {
        // we set water as the wetting phase here
        // which is phase0Idx in the H2oN2 fluid system
        return FluidSystem::phase0Idx;
    }

    Scalar temperature() const
    { return temperature_ ; }

    //! returns the domain marker for an element
    int getElementDomainMarker(const Element& element) const
    { return gridDataPtr_->getElementDomainMarker(element); }

    Scalar pressure(int phaseIdx) const
    {return pressure_ ;}

//    const CouplingManager& couplingManager() const
//    { return *couplingManagerPtr_; }

private:
    //! pointer to the grid data (contains domain markers)
    std::shared_ptr<const Dumux::GridData<Grid>> gridDataPtr_;
//    std::shared_ptr<CouplingManager> couplingManagerPtr_;

    Scalar porosity_;
    const PcKrSwCurve pcKrSwCurve_;
    const PcKrSwCurve barrierPcKrSwCurve_;
    Scalar pressure_;
    Scalar temperature_;
    Scalar a1_,a2_,a3_,a4_,a5_,a6_;
    Scalar kn_, nte_, wte_;
    Scalar a_, b_ ;
};

} // end namespace Dumux

#endif
