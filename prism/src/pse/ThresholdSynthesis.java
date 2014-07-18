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

import java.util.List;
import java.util.LinkedList;
import java.util.Map.Entry;

import parser.ast.Expression;
import parser.ast.ExpressionProb;
import parser.ast.RelOp;
import prism.PrismException;
import prism.PrismLog;

public final class ThresholdSynthesis extends DecompositionProcedure {
	// Synthesis parameters
	private boolean aboveIsTrue;
	private double threshold;
	private double volumeTolerance;

	// Properties of the model being model-checked
	private int initState;
	private double completeSpaceVolume;

	// Solution structures
	private List<BoxRegion> regionsBelowThreshold = new LinkedList<BoxRegion>();
	private List<BoxRegion> regionsAboveThreshold = new LinkedList<BoxRegion>();
	private List<BoxRegion> regionsUndecided = new LinkedList<BoxRegion>();
	private double undecidedVsComplete;

	public ThresholdSynthesis(double volumeTolerance, int initState, BoxRegion completeSpace) throws PrismException
	{
		this.volumeTolerance = volumeTolerance;
		this.initState = initState;
		this.completeSpaceVolume = completeSpace.getVolume();
	}

	@Override
	protected void processPropertyExpression() throws PrismException
	{
		try {
			ExpressionProb probExpr = (ExpressionProb) propExpr;
			RelOp relOp = probExpr.getRelOp();
			if (!relOp.isLowerBound() && !relOp.isUpperBound())
				throw new ClassCastException();
			aboveIsTrue = relOp.isLowerBound();
			Expression boundExpr = probExpr.getProb();
			if (boundExpr == null)
				throw new ClassCastException();
			threshold = (Double) boundExpr.evaluate(modelChecker.getConstantValues());
		} catch (ClassCastException e) {
			throw new PrismException("Threshold synthesis requires a P operator with a lower/upper bound");
		}
	}

	@Override
	public void examineWholeComputation(BoxRegionValues regionValues) throws DecompositionNeeded
	{
		// Determine the regions and compute the volume of undecided regions
		regionsBelowThreshold.clear();
		regionsAboveThreshold.clear();
		regionsUndecided.clear();
		BoxRegion regionToDecompose = null;
		double undecidedVolume = 0.0;
		double greatestVolume = Double.NEGATIVE_INFINITY;
		for (Entry<BoxRegion, BoxRegionValues.StateValuesPair> entry : regionValues) {
			if ((Double) entry.getValue().getMin().getValue(initState) >= threshold)
				regionsAboveThreshold.add(entry.getKey());
			else if ((Double) entry.getValue().getMax().getValue(initState) < threshold)
				regionsBelowThreshold.add(entry.getKey());
			else {
				regionsUndecided.add(entry.getKey());
				undecidedVolume += entry.getKey().getVolume();
				if (entry.getKey().getVolume() > greatestVolume)
					regionToDecompose = entry.getKey();
			}
		}

		// Evaluate whether a decomposition is needed
		undecidedVsComplete = undecidedVolume / completeSpaceVolume;
		if (undecidedVolume / completeSpaceVolume > volumeTolerance) {
			// Decompose a largest undecided region
			throw new DecompositionNeeded(regionToDecompose);
		}
	}

	@Override
	public void printSolution(PrismLog log)
	{
		super.printIntro(log);

		log.print("\nTrue regions:");
		printRegions(log, aboveIsTrue ? regionsAboveThreshold : regionsBelowThreshold);
		log.print("False regions:");
		printRegions(log, aboveIsTrue ? regionsBelowThreshold : regionsAboveThreshold);
		log.print("Undecided regions:");
		printRegions(log, regionsUndecided);

		log.println("\nvol(undecided) / vol(complete parameter space) = " + undecidedVsComplete);
		log.println("volume tolerance = " + volumeTolerance);
	}

	@Override
	public String toString()
	{
		return "threshold synthesis";
	}
}
