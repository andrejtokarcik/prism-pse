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
import parser.ast.ExpressionProb;
import parser.ast.RelOp;
import prism.PrismException;
import prism.PrismLog;

abstract class OptimisingSynthesis extends DecompositionProcedure
{
	// Synthesis parameters
	protected double probTolerance;
	protected String captionForOptimising;

	// Properties of the model being model-checked
	protected int initState;

	// Solution structures
	protected LabelledBoxRegions regionsOptimising;
	protected LabelledBoxRegions regionsNonoptimising;
	private double minimalLowerProbBoundOfOptimising;
	private double maximalUpperProbBoundOfOptimising;
	protected List<Double> demarcationProbBounds;

	public OptimisingSynthesis(double probTolerance, int initState)
	{
		this.probTolerance = probTolerance;
		this.initState = initState;
	}

	@Override
	public void initialiseModelChecking(PSEModelChecker modelChecker, PSEModel model, Expression propExpr) throws PrismException
	{
		super.initialiseModelChecking(modelChecker, model, propExpr);
		regionsOptimising = new LabelledBoxRegions();
		regionsNonoptimising = new LabelledBoxRegions();
		demarcationProbBounds = new LinkedList<Double>();
	}

	@Override
	protected void processPropertyExpression() throws PrismException
	{
		try {
			ExpressionProb probExpr = (ExpressionProb) propExpr;
			if (probExpr.getRelOp() != RelOp.EQ)
				throw new ClassCastException();
		} catch (ClassCastException e) {
			throw new PrismException("Min and max syntheses require a P operator of the form P=?");
		}
	}

	protected abstract void determineOptimisingRegions(BoxRegionValues regionValues) throws PrismException;

	@Override
	protected void verifyRegionValues(BoxRegionValues regionValues) throws DecompositionNeeded, PrismException
	{
		determineOptimisingRegions(regionValues);

		// Determine the deciding probability bounds
		minimalLowerProbBoundOfOptimising = Double.POSITIVE_INFINITY;
		maximalUpperProbBoundOfOptimising = Double.NEGATIVE_INFINITY;
		BoxRegion regionToDecomposeMin = null;
		BoxRegion regionToDecomposeMax = null;
		for (Entry<BoxRegion, BoxRegionValues.StateValuesPair> entry : regionValues) {
			if (!regionsOptimising.contains(entry.getKey()))
				continue;

			double currentLowerProbBound = (Double) entry.getValue().getMin().getValue(initState);
			if (currentLowerProbBound < minimalLowerProbBoundOfOptimising) {
				minimalLowerProbBoundOfOptimising = currentLowerProbBound;
				regionToDecomposeMin = entry.getKey();
			}

			double currentUpperProbBound = (Double) entry.getValue().getMax().getValue(initState);
			if (currentUpperProbBound > maximalUpperProbBoundOfOptimising) {
				maximalUpperProbBoundOfOptimising = currentUpperProbBound;
				regionToDecomposeMax = entry.getKey();
			}
		}

		// Evaluate whether a decomposition is needed
		if (maximalUpperProbBoundOfOptimising - minimalLowerProbBoundOfOptimising > probTolerance) {
			LabelledBoxRegions regionsToDecompose = new LabelledBoxRegions();
			regionsToDecompose.add(regionToDecomposeMin, "min lower prob bound");
			regionsToDecompose.add(regionToDecomposeMax, "max upper prob bound");
			throw new DecompositionNeeded("Probability tolerance was not satisfied, " +
					maximalUpperProbBoundOfOptimising + " - " + minimalLowerProbBoundOfOptimising + " > " + probTolerance,
					regionsToDecompose);
		}
	}

	@Override
	public void printSolution(PrismLog log)
	{
		printIntro(log);

		log.print("\nRegions " + captionForOptimising + " the property satisfaction probability");
		log.println(" (" + regionsOptimising.size() + "):");
		regionsOptimising.print(log);
		log.print("Non-" + captionForOptimising + " regions");
		log.println(" (" + regionsNonoptimising.size() + "):");
		regionsNonoptimising.print(log);

		log.println("\nMin lower prob bound of " + captionForOptimising + " regions = " + minimalLowerProbBoundOfOptimising);
		log.println("Max upper prob bound of " + captionForOptimising + " regions = " + maximalUpperProbBoundOfOptimising);
		log.println("Probability tolerance = " + probTolerance);
		
		log.println("\nDemarcation prob bounds in the order they were used to exclude regions:");
		log.println(demarcationProbBounds);
	}
}
