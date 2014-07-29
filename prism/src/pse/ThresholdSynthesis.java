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

public final class ThresholdSynthesis extends DecompositionProcedure
{
	// Synthesis parameters
	private boolean aboveIsTrue;
	private double threshold;
	private double volumeTolerance;

	// Properties of the model being model-checked
	private int initState;
	private double completeSpaceVolume;

	// Solution structures
	private List<BoxRegion> belowRegions;
	private List<BoxRegion> aboveRegions;
	private List<BoxRegion> undecidedRegions;
	private double undecidedVsComplete;

	public ThresholdSynthesis(double volumeTolerance, int initState, BoxRegion completeSpace) throws PrismException
	{
		this.volumeTolerance = volumeTolerance;
		this.initState = initState;
		this.completeSpaceVolume = completeSpace.getVolume();
	}

	@Override
	public void initialise(PSEModelChecker modelChecker, PSEModel model, Expression propExpr) throws PrismException
	{
		super.initialise(modelChecker, model, propExpr);
		belowRegions = new LinkedList<BoxRegion>();
		aboveRegions = new LinkedList<BoxRegion>();
		undecidedRegions = new LinkedList<BoxRegion>();
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
	public void verifyRegionValues(BoxRegionValues regionValues) throws DecompositionNeeded
	{
		// Determine the regions and compute the volume of undecided regions
		undecidedRegions.clear();
		BoxRegion regionToDecompose = null;
		double undecidedVolume = 0.0;
		double greatestVolume = Double.NEGATIVE_INFINITY;
		for (Entry<BoxRegion, BoxRegionValues.StateValuesPair> entry : regionValues) {
			if (aboveRegions.contains(entry.getKey()) || belowRegions.contains(entry.getKey()))
				continue;

			if ((Double) entry.getValue().getMin().getValue(initState) >= threshold)
				aboveRegions.add(entry.getKey());
			else if ((Double) entry.getValue().getMax().getValue(initState) < threshold)
				belowRegions.add(entry.getKey());
			else {
				undecidedRegions.add(entry.getKey());
				double currentVolume = entry.getKey().getVolume();
				undecidedVolume += currentVolume;
				if (currentVolume > greatestVolume) {
					greatestVolume = currentVolume;
					regionToDecompose = entry.getKey();
				}
			}
		}

		// Evaluate whether a decomposition is needed
		undecidedVsComplete = undecidedVolume / completeSpaceVolume;
		if (undecidedVolume / completeSpaceVolume > volumeTolerance) {
			throw new DecompositionNeeded(regionToDecompose, "largest volume of undecided");
		}
	}

	@Override
	public void printSolution(PrismLog log)
	{
		printIntro(log);

		List<BoxRegion> trueRegions = aboveIsTrue ? aboveRegions : belowRegions;
		List<BoxRegion> falseRegions = aboveIsTrue ? belowRegions : aboveRegions;

		log.println("\nTrue regions (" + trueRegions.size() + "):");
		printRegions(log, trueRegions);
		log.println("False regions (" + falseRegions.size() + "):");
		printRegions(log, falseRegions);
		log.println("Undecided regions (" + undecidedRegions.size() + "):");
		printRegions(log, undecidedRegions);

		log.println("\nvol(undecided) / vol(complete parameter space) = " + undecidedVsComplete);
		log.println("Volume tolerance = " + volumeTolerance);
	}

	@Override
	public String toString()
	{
		return "threshold synthesis";
	}
}
