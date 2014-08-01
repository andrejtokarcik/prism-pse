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

import java.util.Set;

import parser.ast.Expression;
import parser.ast.ExpressionFilter;
import prism.PrismException;
import prism.PrismLog;

abstract class DecompositionProcedure
{
	protected PSEModelChecker modelChecker;
	protected PSEModel model;
	protected Expression propExpr;

	@SuppressWarnings("serial")
	public static class DecompositionNeeded extends Exception
	{
		protected String reason;
		protected LabelledBoxRegions regionsToDecompose;
		protected BoxRegionValues examinedRegionValues;

		public DecompositionNeeded(String reason, LabelledBoxRegions regionsToDecompose)
		{
			this.reason = reason;
			this.regionsToDecompose = regionsToDecompose;
		}

		public DecompositionNeeded(String reason, BoxRegion region, String explanation)
		{
			this.reason = reason;
			this.regionsToDecompose = new LabelledBoxRegions(region, explanation);
		}

		public Set<BoxRegion> getRegionsToDecompose()
		{
			return regionsToDecompose.keySet();
		}

		public void printRegionsToDecompose(PrismLog log)
		{
			log.print("The following " + regionsToDecompose.size() + " regions are to be decomposed");
			log.println(" because " + reason + ":");
			regionsToDecompose.print(log);
		}

		public void setExaminedRegionValues(BoxRegionValues examinedRegionValues)
		{
			this.examinedRegionValues = examinedRegionValues;
		}

		public BoxRegionValues getExaminedRegionValues()
		{
			return examinedRegionValues;
		}
	}

	public void initialise(PSEModelChecker modelChecker, PSEModel model, Expression propExpr) throws PrismException
	{
		this.modelChecker = modelChecker;
		this.model = model;
		this.propExpr = propExpr;
		processPropertyExpression();
	}

	protected void processPropertyExpression() throws PrismException
	{
		// Wrap a filter round the property, if needed
		// (in order to extract the final result of model checking)
		propExpr = ExpressionFilter.addDefaultFilterIfNeeded(propExpr, model.getNumInitialStates() == 1);
	}

	public Expression getPropertyExpression()
	{
		return propExpr;
	}

	public void examineSingleIteration(BoxRegion region, double probsMin[], double probsMax[]) throws DecompositionNeeded, PrismException
	{
		verifySingleRegion(region, probsMin, probsMax);
	}

	protected void verifySingleRegion(BoxRegion region, double probsMin[], double probsMax[]) throws DecompositionNeeded, PrismException {}

	public void examineWholeComputation(BoxRegionValues regionValues) throws DecompositionNeeded, PrismException
	{
		try {
			verifyRegionValues(regionValues);
		} catch (DecompositionNeeded e) {
			e.setExaminedRegionValues(regionValues);
			throw e;
		}
	}

	protected void verifyRegionValues(BoxRegionValues regionValues) throws DecompositionNeeded, PrismException {}

	public void printSolution(PrismLog log)
	{
		// The default filter added above takes care of printing the solution
	}

	protected void printIntro(PrismLog log)
	{
		log.println("\nSolution of " + toString() + " for property " + propExpr + ":");
	}
}
