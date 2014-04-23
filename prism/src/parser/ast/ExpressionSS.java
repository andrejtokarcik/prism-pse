//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford, formerly University of Birmingham)
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

package parser.ast;

import parser.*;
import parser.visitor.*;
import prism.PrismLangException;

public class ExpressionSS extends Expression
{
	RelOp relOp = null;
	Expression prob = null;
	Expression expression = null;
	// Note: this "old-style" filter is just for display purposes
	// The parser creates an (invisible) new-style filter around this expression
	Filter filter = null;
	
	// Constructors
	
	public ExpressionSS()
	{
	}
	
	public ExpressionSS(Expression e, String r, Expression p)
	{
		expression = e;
		relOp = RelOp.parseSymbol(r);
		prob = p;
	}

	// Set methods
	
	public void setRelOp(RelOp relOp)
	{
		this.relOp = relOp;
	}

	public void setRelOp(String r)
	{
		relOp = RelOp.parseSymbol(r);
	}

	public void setProb(Expression p)
	{
		prob = p;
	}

	public void setExpression(Expression e)
	{
		expression = e;
	}
	
	public void setFilter(Filter f)
	{
		filter = f;
	}

	// Get methods
	
	public RelOp getRelOp()
	{
		return relOp;
	}
	
	public Expression getProb()
	{
		return prob;
	}

	public Expression getExpression()
	{
		return expression;
	}
	
	public Filter getFilter()
	{
		return filter;
	}

	// Methods required for Expression:
	
	/**
	 * Is this expression constant?
	 */
	public boolean isConstant()
	{
		return false;
	}

	@Override
	public boolean isProposition()
	{
		return false;
	}
	
	/**
	 * Evaluate this expression, return result.
	 * Note: assumes that type checking has been done already.
	 */
	public Object evaluate(EvaluateContext ec) throws PrismLangException
	{
		throw new PrismLangException("Cannot evaluate an S operator without a model");
	}

	/**
	  * Get "name" of the result of this expression (used for y-axis of any graphs plotted)
	  */
	public String getResultName()
	{
		return (prob == null) ? "Probability" : "Result";
	}

	@Override
	public boolean returnsSingleValue()
	{
		return false;
	}

	// Methods required for ASTElement:
	
	/**
	 * Visitor method.
	 */
	public Object accept(ASTVisitor v) throws PrismLangException
	{
		return v.visit(this);
	}
	
	/**
	 * Convert to string.
	 */
	public String toString()
	{
		String s = "";
		
		s += "S" + relOp;
		s += (prob==null) ? "?" : prob.toString();
		s += " [ " + expression;
		if (filter != null) s += " "+filter;
		s += " ]";
		
		return s;
	}

	/**
	 * Perform a deep copy.
	 */
	public Expression deepCopy()
	{
		ExpressionSS expr = new ExpressionSS();
		expr.setExpression(expression == null ? null : expression.deepCopy());
		expr.setRelOp(relOp);
		expr.setProb(prob == null ? null : prob.deepCopy());
		expr.setFilter(filter == null ? null : (Filter)filter.deepCopy());
		expr.setType(type);
		expr.setPosition(this);
		return expr;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result
				+ ((expression == null) ? 0 : expression.hashCode());
		result = prime * result + ((prob == null) ? 0 : prob.hashCode());
		result = prime * result + ((relOp == null) ? 0 : relOp.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		ExpressionSS other = (ExpressionSS) obj;
		if (expression == null) {
			if (other.expression != null)
				return false;
		} else if (!expression.equals(other.expression))
			return false;
		if (prob == null) {
			if (other.prob != null)
				return false;
		} else if (!prob.equals(other.prob))
			return false;
		if (relOp != other.relOp)
			return false;
		return true;
	}
}

//------------------------------------------------------------------------------
