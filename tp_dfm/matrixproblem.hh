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
 * \brief The sub-problem for the matrix domain in the exercise on two-phase flow in fractured porous media.
 */
#ifndef DUMUX_COURSE_FRACTURESEXERCISE_MATRIX_PROBLEM_HH
#define DUMUX_COURSE_FRACTURESEXERCISE_MATRIX_PROBLEM_HH

// we need this in this test in order to define the domain
// id of the fracture problem (see function interiorBoundaryTypes())
#include <dune/common/indices.hh>

// include the base problem and properties we inherit from
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup MultiDomainFacet
 * \ingroup TwoPTests
 * \brief The sub-problem for the matrix domain in the exercise on two-phase flow in fractured porous media.
 */
template<class TypeTag>
class MatrixSubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using Scalar = typename GridVariables::Scalar;

    using FVGridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;

    static constexpr int dimWorld = GridView::dimensionworld;

    // some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum
    {
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        temperatureIdx = Indices::temperatureIdx,

        //! Equation indices
        contiCO2EqIdx = Indices::conti0EqIdx + FluidSystem::CO2Idx,
		contiH2OEqIdx = Indices::conti0EqIdx + FluidSystem::BrineIdx,
        energyEqIdx = Indices::energyEqIdx,

        //! Phase indices
        wPhaseIdx = FluidSystem::BrineIdx,
        nPhaseIdx = FluidSystem::CO2Idx,
    };

public:
    //! The constructor
    MatrixSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                     std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                     const std::string& paramGroup = "")
    : ParentType(fvGridGeometry, spatialParams, paramGroup)
//    , boundaryOverPressure_(getParamFromGroup<Scalar>(paramGroup, "Problem.BoundaryOverPressure"))
//    , boundarySaturation_(getParamFromGroup<Scalar>(paramGroup, "Problem.BoundarySaturation"))
    , InjectionRate_(getParamFromGroup<Scalar>(paramGroup, "Problem.InjectionRate"))
    , InjectionTemperature_(getParamFromGroup<Scalar>(paramGroup, "Problem.InjectionTemperature"))
    , InjectionPressure_(getParamFromGroup<Scalar>(paramGroup, "Problem.InjectionPressure"))
    , InjectionSaturation_(getParamFromGroup<Scalar>(paramGroup, "Problem.InjectionSaturation"))
    , IsInjectCO2_(getParamFromGroup<Scalar>(paramGroup, "Problem.IsInjectCO2"))
    , IsPureCO2_(getParamFromGroup<Scalar>(paramGroup, "Problem.IsPureCO2"))
    , IsNeumannInjection_(getParamFromGroup<Scalar>(paramGroup, "Problem.IsNeumannInjection"))
    , ProductionPressure_ (getParamFromGroup<Scalar>(paramGroup, "Problem.ProductionPressure"))
    , ProductionYMin_(getParamFromGroup<Scalar>(paramGroup, "Problem.ProductionYMin", 40))
    , ProductionYMax_(getParamFromGroup<Scalar>(paramGroup, "Problem.ProductionYMax", 60))
    , InjectionYMin_(getParamFromGroup<Scalar>(paramGroup, "Problem.InjectionYMin", 40))
    , InjectionYMax_(getParamFromGroup<Scalar>(paramGroup, "Problem.InjectionYMax", 60))
    , name_(getParamFromGroup<std::string>(paramGroup, "Problem.name"))
    {
        // initialize the fluid system, i.e. the tabulation
        // of water properties. Use the default p/T ranges.
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        FluidSystem::init();

        const std::string fileName = "1_matrix" + name_ + ".csv";
        fileOutput_.open(fileName,std::ios::app);

        /* unit of injection well is same as production well */
        fileOutput_ << "time(s),time(year),Sn,Temperature(K),FlowCO2(mole/s),MassCO2(mole),FlowH2O(mole/s),MassH2O(mole),HeatProductionRate(MW),"
        			   "Energy(MJ),InjectionCO2(mole),OverPressure(MPa),CPUtime,Sn_inj,T_inj,OPw_inj,OPnw_inj,S_CO2,FlowCO2_inj,MassCO2_inj,FlowH2O_inj,MassH2O_inj" << std::endl;
        fileOutput_.close();
    }

    template<class VTKWriter>
    void addVtkFields(VTKWriter& vtk)
    {
        vtk.addField(deltaT_, "deltaT");
        vtk.addField(overPressure_, "deltaP");
    }

    void updateVtkFields(const SolutionVector& curSol)
    {
    	int dofCodim = 0;
    	const auto& gridView = this->gridGeometry().gridView();
    	deltaT_.resize(gridView.size(dofCodim));
    	overPressure_.resize(gridView.size(dofCodim));

        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto elemSol = elementSolution(element, curSol, this->gridGeometry());
//
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);
//
            for (auto&& scv : scvs(fvGeometry))
            {
                const auto dofIdxGlobal = scv.dofIndex();
                const GlobalPosition& globalPos = scv.center();
//                VolumeVariables volVars;
                const auto& priVars = elemSol[scv.localDofIndex()];
                const auto initialValues = initialAtPos(globalPos);
//                volVars.update(elemSol, *this, element, scv);
                deltaT_[dofIdxGlobal] = priVars[temperatureIdx] - initialValues[temperatureIdx];
                overPressure_[dofIdxGlobal] = priVars[pressureIdx] - initialValues[pressureIdx];
            }
        }
    }

    //! Specifies the type of boundary condition at a given position
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;

        // we consider buoancy-driven upwards migration of nitrogen and set
        // Dirichlet BCs on the top and bottom boundary
        values.setAllNeumann();
//        if (globalPos[1] < 1e-6 || globalPos[1] > this->gridGeometry().bBoxMax()[1] - 1e-6)
		if (isProduction(globalPos))
			values.setAllDirichlet();

		if(!IsNeumannInjection_)
		{
			if(isInjection(globalPos))
				values.setAllDirichlet();
		}

        return values;
    }

    //! Specifies the type of interior boundary condition to be used on a sub-control volume face
    BoundaryTypes interiorBoundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();

        return values;
    }

    //! evaluates the Dirichlet boundary condition for a given position
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        // initialize values with the initial conditions
        auto values = initialAtPos(globalPos);

        if (isProduction(globalPos))
        	values[pressureIdx] += - ProductionPressure_;

    	if (isInjection(globalPos))
    	{
        	if (IsInjectCO2_)
			{
				values[pressureIdx] += InjectionPressure_;
				values[saturationIdx] = InjectionSaturation_;
				values[temperatureIdx] = InjectionTemperature_;
			}
        	else
        	{
				values[pressureIdx] += InjectionPressure_;
				values[saturationIdx] = 0;
				values[temperatureIdx] = InjectionTemperature_;
        	}
        }

        return values;
    }

    //! evaluate the initial conditions
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        // For the grid used here, the height of the domain is equal
        // to the maximum y-coordinate
        const auto domainHeight = this->gridGeometry().bBoxMax()[dimWorld-1] + 6000;

        // we assume a constant water density of 1000 for initial conditions!
        const auto& g = this->spatialParams().gravity(globalPos);
        PrimaryVariables values;
        Scalar densityW = 1000.0;
        values[pressureIdx] = 1e5 - (domainHeight - globalPos[dimWorld-1])*densityW*g[dimWorld-1];
        values[temperatureIdx] = 283.0 + (domainHeight - globalPos[dimWorld-1])*0.03;
        if (IsPureCO2_)
        	{values[saturationIdx] = 1.0 - eps_;}
        else
        	{values[saturationIdx] = 0.0;}

        return values;
    }

    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
    	NumEqVector values(0.0);

		const auto& globalPos = scvf.ipGlobal();

		if (IsNeumannInjection_)
		{
			if (globalPos[dimWorld-1] < 60 + eps_ && globalPos[dimWorld-1] > 40 - eps_ && globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_)
			{
				// compute enthalpy flux associated with this injection [(J/(kg*s)]

				// energy flux is mass flux times specific enthalpy
				// inject air. negative values mean injection
	            using FluidState = GetPropType<TypeTag, Properties::FluidState>;
	            FluidState fs;

	            const auto initialValues = initialAtPos(globalPos);
	            fs.setPressure(wPhaseIdx, initialValues[pressureIdx]);
	            fs.setPressure(nPhaseIdx, initialValues[pressureIdx]); // assume pressure equality here
	            fs.setTemperature(wPhaseIdx, InjectionTemperature_);
	            fs.setTemperature(nPhaseIdx, InjectionTemperature_);

	            if (IsInjectCO2_)
	            {
					values[contiCO2EqIdx] = -InjectionRate_ ; // kg/(s*m^2) flow rate * density
	        		values[energyEqIdx] = values[contiCO2EqIdx]*FluidSystem::enthalpy(fs, nPhaseIdx);
	            }
				else
				{
					values[contiH2OEqIdx] = -InjectionRate_ ; // kg/(s*m^2) flow rate * density
	        		values[energyEqIdx] = values[contiH2OEqIdx]*FluidSystem::enthalpy(fs, wPhaseIdx);
				}
			}
		}
//		if (globalPos[dimWorld-1] < 75 + eps_ && globalPos[dimWorld-1] > 25 - eps_ && globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_)
//		{
//			// compute enthalpy flux associated with this injection [(J/(kg*s)]
//			const auto& volVars = elemVolVars[scvf.insideScvIdx()];
//
//			values[contiCO2EqIdx] = InjectionRate_ * volVars.saturation(nPhaseIdx) ; // kg/(s*m^2)
//			values[contiH2OEqIdx] = InjectionRate_ * volVars.saturation(wPhaseIdx) ; // kg/(s*m^2)
//			values[energyEqIdx] = values[contiCO2EqIdx]*FluidSystem::enthalpy(volVars.fluidState(), nPhaseIdx) +
//								  values[contiH2OEqIdx]*FluidSystem::enthalpy(volVars.fluidState(), wPhaseIdx);
//			values[energyEqIdx] = values[conti0EqIdx]*IapwsH2O::liquidEnthalpy(volVars.temperature(),volVars.pressure());
//		}

		return values;
    }

    void calculateOutput(const GridVariables& gridVariables, const SolutionVector& sol, Scalar time)
    {
        NumEqVector Source(0.0);
        PrimaryVariables InitialValues;
		using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
		FluxVariables fluxVars;

        Scalar index_p = 0, index_pin = 0, index_pf = 0, index_if = 0;
        Scalar outputT = 0, outputSn = 0, outputOverPw = 0, outputOverPn = 0, outputEnw = 0, outputEw = 0;
        Scalar injectT = 0, injectSn = 0, injectOverPw = 0, injectOverPn = 0;
        Scalar flux_nw = 0, flux_w = 0, flux_inj = 0, flux_winj = 0, flux_nwinj = 0;
        Scalar Area_inj = 0;
        Scalar productionCO2 = 0, productionH2O = 0, injectionCO2 = 0, storageCO2 = 0, Evolume = 0;
        using std::abs;

        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
			auto fvGeometry = localView(this->gridGeometry());
			fvGeometry.bind(element);
			auto elemVolVars = localView(gridVariables.curGridVolVars());
			elemVolVars.bind(element, fvGeometry, sol);
	        auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());
	        elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

			for (auto&& scv : scvs(fvGeometry))
			{
				const auto& globalPos = scv.center();
				const auto& volVars = elemVolVars[scv];
				Source = this->source(element, fvGeometry, elemVolVars, scv);
				InitialValues = initialAtPos(globalPos);
				storageCO2 += volVars.saturation(nPhaseIdx) * scv.volume() * FluidSystem::CO2::gasDensity(volVars.temperature(),volVars.pressure(nPhaseIdx));
				Evolume = scv.volume();

				if (isProduction_area(globalPos))
				{
					index_p += 1;
					outputSn += volVars.saturation(nPhaseIdx);
					outputT += volVars.temperature();
					outputOverPw += (volVars.pressure(wPhaseIdx) - InitialValues[pressureIdx]);
					outputOverPn += (volVars.pressure(nPhaseIdx) - InitialValues[pressureIdx]);
				}

				if(isInjection_area(globalPos))
				{
					index_pin += 1;
					injectSn += volVars.saturation(nPhaseIdx);
					injectT += volVars.temperature();
					injectOverPw += (volVars.pressure(wPhaseIdx) - InitialValues[pressureIdx]);
					injectOverPn += (volVars.pressure(nPhaseIdx) - InitialValues[pressureIdx]);
				}
			}

/* ************** Calculate production Flux  ****************** */
			for (const auto& scvf : scvfs(fvGeometry))
			{
				if (scvf.boundary())
					continue;

				const auto& globalPos = scvf.ipGlobal();
		        const auto& volVars = elemVolVars[scvf.insideScvIdx()];
				const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
                const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
                const bool isOnProduction = insideScv.center()[0] < ProductionBoundary_ && outsideScv.center()[0] > ProductionBoundary_ && globalPos[1] > ProductionLowLimit_ - eps_;
				fluxVars.init(*this, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
				if (isOnProduction)
				{
					auto upwindTerm_nw = [] (const auto& volVars)
					{ return volVars.mobility(nPhaseIdx) * FluidSystem::CO2::gasDensity(volVars.temperature(),volVars.pressure(nPhaseIdx)); };
					auto upwindTerm_w = [] (const auto& volVars)
					{ return volVars.mobility(wPhaseIdx) * FluidSystem::Brine::liquidDensity(volVars.temperature(), volVars.pressure(wPhaseIdx)); };

					index_pf += 1;
					flux_nw += (fluxVars.advectiveFlux(nPhaseIdx, upwindTerm_nw) / FluidSystem::CO2::molarMass());
					flux_w += (fluxVars.advectiveFlux(wPhaseIdx, upwindTerm_w) / FluidSystem::Brine::molarMass());
					outputEnw += (fluxVars.advectiveFlux(nPhaseIdx, upwindTerm_nw) * (-FluidSystem::CO2::gasEnthalpy(InjectionTemperature_,volVars.pressure(nPhaseIdx)) + FluidSystem::CO2::gasEnthalpy(volVars.temperature(),volVars.pressure(nPhaseIdx))));
					outputEw += (fluxVars.advectiveFlux(wPhaseIdx, upwindTerm_w) * (-FluidSystem::Brine::liquidEnthalpy(InjectionTemperature_,volVars.pressure(wPhaseIdx)) + FluidSystem::Brine::liquidEnthalpy(volVars.temperature(),volVars.pressure(wPhaseIdx))));
				}

//				outputEnw = (flux_nw * FluidSystem::CO2::molarMass() * FluidSystem::CO2::gasEnthalpy(volVars.temperature(),volVars.pressure(nPhaseIdx)));
//				outputEw = (flux_w * FluidSystem::Brine::molarMass() * FluidSystem::Brine::gasEnthalpy(volVars.temperature(),volVars.pressure(wPhaseIdx)));
			}

		/* *************** Calculate Injection Flux ******************/
			for (const auto& scvf: scvfs(fvGeometry))
			{
				const auto& globalPos = scvf.ipGlobal();
				fluxVars.init(*this, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
				if (isInjection(globalPos))
				{
					if (scvf.boundary())
					 {
						flux_inj += (- neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf)[contiCO2EqIdx]*scvf.area());
						Area_inj += scvf.area();
						index_if += 1;

						/* Calculate the flow rate near injection well */
//						flux_winj += abs(fluxVars.advectiveFlux(wPhaseIdx, upwindTerm_w) / FluidSystem::Brine::molarMass());
//						flux_nwinj += abs(fluxVars.advectiveFlux(nPhaseIdx, upwindTerm_nw) / FluidSystem::CO2::molarMass());
					 }
				}
			}
        }

        outputSn /= index_p;
        outputT /= index_p;
        outputOverPw /= (index_p * 1e6);   // transfer to MPa
        outputOverPn /= (index_p * 1e6);
        injectSn /= index_pin;
        injectT /= index_pin;
        injectOverPw /= (index_pin * 1e6);
        injectOverPn /= (index_pin * 1e6);

        productionCO2 = flux_nw * timeStepSize_;       // mole
        productionH2O = flux_w * timeStepSize_;
        totalProductionCO2_ += productionCO2;
        totalProductionH2O_ += productionH2O;
        totalH2OinInjection_ += flux_winj * timeStepSize_;
        totalCO2inInjection_ += flux_nwinj * timeStepSize_;
        injectionCO2 = flux_inj * timeStepSize_;
        totalInjectionCO2_ += injectionCO2;
        Scalar HeatproductionRate = (outputEnw + outputEw)/1e6;
        outputEnw *= timeStepSize_ /1e6;
        outputEw *= timeStepSize_ /1e6;
        totalE_ += outputEnw + outputEw;

        std::cout << '\n' << " **** NAME **** " << '\n';
        std::cout << '\n' << name_ << '\n';

        std::cout << '\n' << " **** SATURATION **** " << '\n';
        std::cout << "outputSn = " << outputSn << ", storageCO2 = " << storageCO2 << ", scv.volume = " << Evolume << '\n' << '\n';
        std::cout << " **** TEMPERATURE AND PRESSURE **** " << '\n';
        std::cout << "outputT = " << outputT <<  " outputOverPw = " << outputOverPw << " outputOverPn = " << outputOverPn << '\n' << '\n';
//        std::cout << " **** MASSFLOW **** " << '\n' ;
//        std::cout << "CO2 produced in this time step = " << productionCO2 << " mole" << '\n' << "total CO2 produced sofar = " << totalProductionCO2_ << " mole" << '\n';
//        std::cout << "H2O produced in this time step = " << productionH2O << " mole" << '\n' << "total H2O produced sofar = " << totalProductionH2O_ << " mole" << '\n' << '\n';
//        std::cout << "CO2 injected in this time step = " << injectionCO2 << " mole" << '\n' << "total CO2 injected sofar = " << totalInjectionCO2_ << " mole" << '\n' << '\n';
//        std::cout << " **** ENERGY **** " << '\n' ;
//        std::cout << "Energy produced from CO2 in this time step = " << outputEnw << " MJ" << '\n';
//        std::cout << "Energy produced from H2O in this time step = " << outputEw << " MJ" << '\n';
//        std::cout << "Total Energy produced = " <<totalE_ << " MJ" << '\n' << '\n' << '\n';

//        fileOutput_ << "time time(a) Sn Temperature OverPressure MassCO2 MassH2O Energy InjectionCO2" << std::endl;

        const std::string fileName = "1_matrix" + name_ + ".csv";
//        const std::string fileName = "~/home/dejian/dumux/tp/build-cmake/tp/new/randomK/lowmesh/result/1" + name_ + ".csv";
        fileOutput_.open(fileName,std::ios::app);
//        fileOutput_ << time << " " << time/86400/365 << " " << outputSn << " " << outputT << " " << outputOverPw << " " << totalProductionCO2_ <<
//        			" " << totalProductionH2O_ << " " << outputE << " " << totalE_ << " " << totalInjectionCO2_ << std::endl;
        fileOutput_ << time << "," << time/86400/365 << "," << outputSn << "," << outputT  << "," << flux_nw << "," << totalProductionCO2_ << ","
        			<< flux_w << "," << totalProductionH2O_ << "," << HeatproductionRate << "," << totalE_ << "," << totalInjectionCO2_ << ","
					<< outputOverPw << "," << wallTime_ << "," << injectSn << "," << injectT << "," << injectOverPw << "," << injectOverPn << ","
					<< storageCO2 << "," << flux_winj << "," << totalH2OinInjection_ << "," << flux_nwinj << "," << totalCO2inInjection_ << ","
					<< std::endl;

        fileOutput_.close();
    }

    //! returns the temperature in \f$\mathrm{[K]}\f$ in the domain
    Scalar temperature() const
    { return temperature_; }

    //! sets the pointer to the coupling manager.
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManagerPtr_ = cm; }

    //! returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

    void setTime( Scalar time )
    {
        time_ = time;
    }

    /*!
     * \brief The time step size.
     *
     * This is used to calculate the source term.
     */
    void setTimeStepSize( Scalar timeStepSize )
     {
        timeStepSize_ = timeStepSize;
     }

    void setWallTime(Scalar wallTime)
    {
    	wallTime_ = wallTime;
    }

private:

    bool isProduction (const GlobalPosition globalPos) const
    {
    	return globalPos[1] < ProductionYMax_ + eps_ &&
			   globalPos[1] > ProductionYMin_ - eps_ &&
			   globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_ ;
    }

    bool isInjection (const GlobalPosition globalPos) const
    {
    	return globalPos[1] < InjectionYMax_ &&
    		   globalPos[1] > InjectionYMin_ &&
			   globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_;
    }

    bool isProduction_area (const GlobalPosition globalPos) const
    {
    	return globalPos[1] < ProductionYMax_ + eps_ &&
			   globalPos[1] > ProductionYMin_ - eps_ &&
			   globalPos[0] < ProductionXMax_ + eps_ &&
			   globalPos[0] > ProductionXMax_ - elementX_ - eps_;
    }

    bool isInjection_area (const GlobalPosition globalPos) const
    {
    	return globalPos[1] < InjectionYMax_ &&
    		   globalPos[1] > InjectionYMin_ &&
			   globalPos[0] < this->gridGeometry().bBoxMin()[0] + elementX_ + eps_ &&
			   globalPos[0] > this->gridGeometry().bBoxMin()[0] - eps_;
    }

    std::shared_ptr<CouplingManager> couplingManagerPtr_;

    Scalar boundaryOverPressure_;
    Scalar boundarySaturation_;
    Scalar InjectionTemperature_, InjectionPressure_, InjectionSaturation_;
    Scalar InjectionRate_;
    bool IsInjectCO2_;
    bool IsPureCO2_;
    bool IsNeumannInjection_;
    Scalar temperature_;
    Scalar saturation_;
    static constexpr Scalar eps_ = 1e-7;
    std::vector<Scalar> overPressure_, deltaT_;

    Scalar InjectionYMax_, InjectionYMin_, ProductionYMax_, ProductionYMin_, ProductionLowLimit_ = 40;
    Scalar ProductionPressure_;
    Scalar time_ = 0.0, timeStepSize_ = 0.0, wallTime_ = 0.0;
	Scalar elementX_ = 10, ProductionBoundary_ = 980, ProductionXMax_ = 980;

	Scalar totalProductionCO2_ = 0, totalProductionH2O_ = 0, totalInjectionCO2_ = 0, totalE_ = 0, totalCO2inInjection_ = 0, totalH2OinInjection_ = 0;

	std::string name_;
    std::ofstream fileOutput_;
};

} // end namespace Dumux

#endif
