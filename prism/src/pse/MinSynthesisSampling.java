//==============================================================================
//	
//	Copyright (c) 2014-
//	Authors:
//	* Andrej Tokarcik <andrejtokarcik@gmail.com> (Masaryk University)
//	
//------------------------------------------------------------------------------
//	
//	This file is part of PRISM.
//	
//	PRISM is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//	
//	PRISM is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//	
//	You should have received a copy of the GNU General Public License
//	along with PRISM; if not, write to the Free Software Foundation,
//	Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//	
//==============================================================================

package pse;

import java.util.Map.Entry;

import explicit.CTMC;
import explicit.CTMCModelChecker;
import prism.PrismException;
import simulator.SimulatorEngine;

public final class MinSynthesisSampling extends MinSynthesis {
	private explicit.ConstructModel constructModel;
	private double lastMinimalSampleProb = Double.POSITIVE_INFINITY;

	public MinSynthesisSampling(double probTolerance, int initState, SimulatorEngine simulatorEngine)
	{
		super(probTolerance, initState);
		constructModel = new explicit.ConstructModel(modelChecker, simulatorEngine);
	}

	/**
	 * Sampling approach to determining the maximal lower bound.
	 */
	@Override
	protected double getMinimalUpperBound(BoxRegionValues regionValues) throws PrismException
	{
		double minimalUpperBound = Double.POSITIVE_INFINITY;
		BoxRegion minimalUpperBoundRegion = null;
		for (Entry<BoxRegion, BoxRegionValues.StateValuesPair> entry : regionValues) {
			double currentUpperBound = (Double) entry.getValue().getMax().getValue(initState);
			if (currentUpperBound < minimalUpperBound) {
				minimalUpperBound = currentUpperBound;
				minimalUpperBoundRegion = entry.getKey();
			}
		}

		double minimalSampleProb = Double.POSITIVE_INFINITY;
		for (Point sample : minimalUpperBoundRegion.getPointSamples()) {
			CTMC ctmc = model.instantiate(sample, modelChecker.getModulesFile(), constructModel);
			CTMCModelChecker ctmcModelChecker = new CTMCModelChecker(modelChecker);
			double currentSampleProb = (Double) ctmcModelChecker.checkExpression(ctmc, propExpr).getValue(initState);
			if (currentSampleProb < minimalSampleProb)
				minimalSampleProb = currentSampleProb;
		}

		if (minimalSampleProb < lastMinimalSampleProb) {
			lastMinimalSampleProb = minimalSampleProb;
			return minimalSampleProb;
		}
		return lastMinimalSampleProb;
	}

	@Override
	public String toString()
	{
		return super.toString() + " (sampling)";
	}
}
