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
import parser.ast.ExpressionProb;
import parser.ast.RelOp;
import prism.PrismException;

abstract class AbstractMinMaxSynthesis extends DecompositionProcedure {
	// Parameters of the synthesis procedures
	protected Expression propExpr;
	protected double probTolerance;

	// Properties of the model being model-checked
	protected int initState;

	// Solution structures
	protected double probDifference;

	public AbstractMinMaxSynthesis(Expression propExpr, double probTolerance, int initState) throws PrismException
	{
		this.propExpr = propExpr;
		processPropertyExpression(propExpr);
		this.probTolerance = probTolerance;
		this.initState = initState;
	}

	private void processPropertyExpression(Expression propExpr) throws PrismException
	{
		try {
			ExpressionProb probExpr = (ExpressionProb) propExpr;
			if (probExpr.getRelOp() != RelOp.EQ)
				throw new ClassCastException();
		} catch (ClassCastException e) {
			throw new PrismException("Min and max syntheses require a P operator of the form P=?");
		}
	}

	@Override
	public Expression adjustPropertyExpression(Expression propExpr, PSEModel model) throws PrismException
	{
		// Prevent even the default filter from being added
		return propExpr;
	}
}
