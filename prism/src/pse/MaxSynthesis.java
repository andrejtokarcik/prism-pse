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

import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;

import parser.ast.Expression;
import prism.PrismException;
import prism.PrismLog;

abstract class MaxSynthesis extends AbstractMinMaxSynthesis {
	// Solution structures
	protected List<BoxRegion> regionsMaximising = new LinkedList<BoxRegion>();
	protected List<BoxRegion> regionsNonmaximising = new LinkedList<BoxRegion>();
	
	public MaxSynthesis(Expression propExpr, double probTolerance, int initState) throws PrismException
	{
		super(propExpr, probTolerance, initState);
	}

	public abstract double getMaximalLowerBound(BoxRegionValues regionValues);
	
	@Override
	public void examineWholeComputation(BoxRegionValues regionValues) throws DecompositionNeeded
	{
		regionsMaximising.clear();
		regionsNonmaximising.clear();

		// NB: In the following, the term `bounds' refers to constraints on the probability
		// of the property's being satisfied in a given region.  This is not to be confused
		// with `bounds' in the sense of upper/lower values of parameter ranges characterising
		// the parameter regions/subspaces.

		// Determine the maximal lower bound
		double maximalLowerBound = getMaximalLowerBound(regionValues);

		// Determine the maximising regions
		for (Entry<BoxRegion, BoxRegionValues.StateValuesPair> entry : regionValues) {
			if ((Double) entry.getValue().getMax().getValue(initState) < maximalLowerBound)
				regionsNonmaximising.add(entry.getKey());
			else
				regionsMaximising.add(entry.getKey());
		}

		// Determine the deciding probability bounds
		BoxRegion regionToDecompose = null;
		double minimalLowerBoundOfMaximising = Double.POSITIVE_INFINITY;
		double maximalUpperBoundOfMaximising = Double.NEGATIVE_INFINITY;
		for (Entry<BoxRegion, BoxRegionValues.StateValuesPair> entry : regionValues) {
			if (!regionsMaximising.contains(entry.getKey()))
				continue;

			double currentLowerBound = (Double) entry.getValue().getMin().getValue(initState);
			if (currentLowerBound < minimalLowerBoundOfMaximising) {
				minimalLowerBoundOfMaximising = currentLowerBound;
			}

			double currentUpperBound = (Double) entry.getValue().getMax().getValue(initState);
			if (currentUpperBound > maximalUpperBoundOfMaximising) {
				regionToDecompose = entry.getKey();
				maximalUpperBoundOfMaximising = currentUpperBound;
			}
		}
		
		// Evaluate whether a decomposition is needed
		probDifference = maximalUpperBoundOfMaximising - minimalLowerBoundOfMaximising;
		if (probDifference > probTolerance) {
			// Decompose a maximising region with the maximal upper bound
			throw new DecompositionNeeded(regionToDecompose);
		}
	}

	@Override
	public void printSolution(PrismLog log)
	{
		log.println("\nSolution of the max synthesis problem for property " + propExpr + " using the naive approach:\n");

		log.print("Regions maximising the property satisfaction probability:");
		printRegions(log, regionsMaximising);
		log.print("Non-maximising regions:");
		printRegions(log, regionsNonmaximising);
		
		log.println("\nmax {upper prob bounds of maximising} - min {lower prob bounds of maximising} = " + probDifference);
		log.println("probability tolerance = " + probTolerance);
	}
}
