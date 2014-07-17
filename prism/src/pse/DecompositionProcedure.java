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

import parser.ast.Expression;
import parser.ast.ExpressionFilter;
import prism.PrismException;
import prism.PrismLog;

abstract class DecompositionProcedure {
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

	public Expression adjustPropertyExpression(Expression propExpr, PSEModel model) throws PrismException
	{
		// Wrap a filter round the property, if needed
		// (in order to extract the final result of model checking)
		return ExpressionFilter.addDefaultFilterIfNeeded(propExpr, model.getNumInitialStates() == 1);
	}
	
	public void examineSingleIteration(BoxRegion region, double probsMin[], double probsMax[]) throws DecompositionNeeded {}
	
	public void examineWholeComputation(BoxRegionValues regionValues) throws DecompositionNeeded {}
	
	public void printSolution(PrismLog log) {
		// The solution is printed when model-checking the default filter
	}
}
