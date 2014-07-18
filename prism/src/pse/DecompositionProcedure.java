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

import parser.Values;
import parser.ast.Expression;
import parser.ast.ExpressionFilter;
import prism.PrismException;
import prism.PrismLog;

abstract class DecompositionProcedure {
	protected Expression propExpr;

	@SuppressWarnings("serial")
	public static class DecompositionNeeded extends Exception
	{
		protected BoxRegion region;

		public DecompositionNeeded(BoxRegion region)
		{
			this.region = region;
		}

		public BoxRegion getRegion()
		{
			return region;
		}
	}

	protected void processPropertyExpression(boolean singleInit, Values constantValues) throws PrismException
	{
		// Wrap a filter round the property, if needed
		// (in order to extract the final result of model checking)
		propExpr = ExpressionFilter.addDefaultFilterIfNeeded(propExpr, singleInit);
	}

	public Expression adjustPropertyExpression(Expression propExpr, boolean singleInit, Values constantValues) throws PrismException
	{
		this.propExpr = propExpr;
		processPropertyExpression(singleInit, constantValues);
		return this.propExpr;
	}

	public void examineSingleIteration(BoxRegion region, double probsMin[], double probsMax[]) throws DecompositionNeeded {}
	
	public void examineWholeComputation(BoxRegionValues regionValues) throws DecompositionNeeded {}

	public void printSolution(PrismLog log) {
		// The default filter added above takes care of printing the solution
	}

	protected void printIntro(PrismLog log) {
		log.println("\nSolution of " + toString() + " for property " + propExpr + ":");
	}

	protected void printRegions(PrismLog log, List<BoxRegion> regions) {
		if (regions.isEmpty()) {
			log.println("\n * [none]");
		} else {
			log.println();
			for(BoxRegion region : regions) {
				log.println(" * " + region);
		    }
		}
	}
}
