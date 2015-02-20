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

import parser.ast.Expression;
import parser.ast.ExpressionProb;
import parser.ast.ExpressionReward;
import parser.ast.RelOp;
import prism.PrismException;
import prism.PrismLog;

/**
 * Decomposition procedure solving the threshold synthesis problem.
 */
public final class ThresholdSynthesis extends DecompositionProcedure
{
	// Synthesis parameters
	/** true iff the property formula is of the form P>=r[subformula] */
	private boolean aboveIsTrue;
	/** threshold for distinguishing "above" regions from those "below" */
	private double threshold;
	/** greatest acceptable value for {@link #undecidedVsComplete} */
	private double volumeTolerance;

	// Properties of the model to be checked
	/** model's initial state */
	private int initState;
	/** volume of the complete parameter space */
	private double completeSpaceVolume;

	// Solution structures
	/** regions marked as "below" {@code threshold} */
	private LabelledBoxRegions belowRegions;
	/** regions marked as "above" {@code threshold} */
	private LabelledBoxRegions aboveRegions;
	/** regions marked as yet "undecided" */
	private LabelledBoxRegions undecidedRegions;
	/** ratio of total volume of "undecided" regions to volume
	 *  of the complete parameter space */
	private double undecidedVsComplete;

	public ThresholdSynthesis(double volumeTolerance) throws PrismException
	{
		this.volumeTolerance = volumeTolerance;
	}

	@Override
	public void initialiseModelChecking(PSEModelChecker modelChecker, PSEModel model, Expression propExpr) throws PrismException
	{
		super.initialiseModelChecking(modelChecker, model, propExpr);
		initState = model.getFirstInitialState();
		completeSpaceVolume = model.getCompleteSpace().volume();
		belowRegions = new LabelledBoxRegions();
		aboveRegions = new LabelledBoxRegions();
		undecidedRegions = new LabelledBoxRegions();
	}

	@Override
	protected void processPropertyExpression() throws PrismException
	{
		try {
			RelOp relOp;
			Expression boundExpr;
			try {
				ExpressionProb probExpr = (ExpressionProb) propExpr;
				relOp = probExpr.getRelOp();
				boundExpr = probExpr.getProb();
			} catch (ClassCastException e) {
				ExpressionReward rewExpr = (ExpressionReward) propExpr;
				relOp = rewExpr.getRelOp();
				boundExpr = rewExpr.getReward();
			}
			if (!relOp.isLowerBound() && !relOp.isUpperBound())
				throw new ClassCastException();
			aboveIsTrue = relOp.isLowerBound();
			if (boundExpr == null)
				throw new ClassCastException();
			threshold = (Double) boundExpr.evaluate(modelChecker.getConstantValues());
		} catch (ClassCastException e) {
			throw new PrismException("Threshold synthesis requires a bounded P or R operator");
		}
	}

	@Override
	protected void verifyRegionValues(BoxRegionValues regionValues) throws DecompositionNeeded
	{
		// Determine the regions and compute the volume of undecided regions
		undecidedRegions.clear();
		BoxRegion regionToDecompose = null;
		double undecidedVolume = 0.0;
		double greatestVolume = Double.NEGATIVE_INFINITY;
		for (Entry<BoxRegion, BoxRegionValues.StateValuesPair> entry : regionValues) {
			if (aboveRegions.contains(entry.getKey()) || belowRegions.contains(entry.getKey())) {
				continue;
			}
			double lowerProbBound = (Double) entry.getValue().getMin().getValue(initState);
			double upperProbBound = (Double) entry.getValue().getMax().getValue(initState);
			if (lowerProbBound >= threshold) {
				aboveRegions.add(entry.getKey(), "lower prob bound = " + lowerProbBound + ", upper prob bound = " + upperProbBound);
			} else if (upperProbBound < threshold) {
				belowRegions.add(entry.getKey(), "lower prob bound = " + lowerProbBound + ", upper prob bound = " + upperProbBound);
			} else {
				undecidedRegions.add(entry.getKey(), "lower prob bound = " + lowerProbBound + ", upper prob bound = " + upperProbBound);
				double currentVolume = entry.getKey().volume();
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
			throw new DecompositionNeeded("The volume of undecided regions is too large,\n" +
					undecidedVolume + " / " + completeSpaceVolume + " > " + volumeTolerance,
					regionToDecompose, "largest undecided region");
		}
	}

	@Override
	public void printSolution(PrismLog log, boolean verbose)
	{
		printIntro(log);

		LabelledBoxRegions trueRegions = aboveIsTrue ? aboveRegions : belowRegions;
		LabelledBoxRegions falseRegions = aboveIsTrue ? belowRegions : aboveRegions;

		log.println("\nTrue regions (" + trueRegions.size() + "):");
		trueRegions.print(log);
		log.println("False regions (" + falseRegions.size() + "):");
		falseRegions.print(log);
		log.println("Undecided regions (" + undecidedRegions.size() + "):");
		undecidedRegions.print(log);

		log.println("\nvol(undecided) / vol(complete parameter space) = " + undecidedVsComplete);
		log.println("Volume tolerance = " + volumeTolerance);
	}

	@Override
	public String toString()
	{
		return "threshold synthesis";
	}
}
