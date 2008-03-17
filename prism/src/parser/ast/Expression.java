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

// Abstract class for PRISM language expressions

public abstract class Expression extends ASTElement
{
	// type constants
	public static final int INT = 1;
	public static final int DOUBLE = 2;
	public static final int BOOLEAN = 3;
	
	// and a function to get at their names
	public static String getTypeString(int i) {
		switch (i) {
			case INT: return "int";
			case DOUBLE: return "double";
			case BOOLEAN: return "bool";
			default: return "(unknown)";
		}
	}
	
	// another useful function telling you which types can be assigned to which others
	public static boolean canAssignTypes(int tl, int tr)
	{
		switch (tl) {
			// boolean can only be assigned boolean
			case Expression.BOOLEAN: return (tr == Expression.BOOLEAN);
			// int can only be assigned int
			case Expression.INT: return (tr == Expression.INT);
			// double can be assigned int or double
			case Expression.DOUBLE: return (tr == Expression.INT || tr == Expression.DOUBLE);
			// should never happen...
			default: return false;
		}
	}
	
	// Methods required for Expression (all subclasses should implement/override):
	
	/**
	 * Is this expression constant?
	 */
	public abstract boolean isConstant();

	/**
	 * Evaluate this expression, return result.
	 * Note: assumes that type checking has been done already.
	 */
	public abstract Object evaluate(Values constantValues, Values varValues) throws PrismLangException;
	
	/**
	  * Get "name" of the result of this expression (used for y-axis of any graphs plotted)
	  */
	public String getResultName()
	{
		return "Result";
	}

	 // Overrided version of deepCopy() from superclass ASTElement (to reduce casting).
	 
	/**
	 * Perform a deep copy.
	 */
	public abstract Expression deepCopy();
	
	// Utility methods:
	
	/**
	 * Check expression (property) for validity with respect to a particular model type
	 * (i.e. whether not it is a property that can be model checked for that model type).
	 */
	public void checkValid(int modelType) throws PrismLangException
	{
		CheckValid visitor = new CheckValid(modelType);
		accept(visitor);
	}
	
	// evaluate to an int
	// any typing issues cause an exception
	// [does nothing to the expression itself]
	public int evaluateInt(Values constantValues, Values varValues) throws PrismLangException
	{
		Object o;
		
		o = evaluate(constantValues, varValues);
		
		if (!(o instanceof Integer)) {
			throw new PrismLangException("Cannot evaluate to an integer", this);
		}
		
		return ((Integer)o).intValue();
	}

	// evaluate to a double
	// any typing issues cause an exception
	// [does nothing to the expression itself]
	public double evaluateDouble(Values constantValues, Values varValues) throws PrismLangException
	{
		Object o;
		
		o = evaluate(constantValues, varValues);
		
		if (o instanceof Boolean) {
			throw new PrismLangException("Cannot evaluate to a double", this);
		}
		if (o instanceof Integer) {
			return ((Integer)o).intValue();
		}
		if (o instanceof Double) {
			return ((Double)o).doubleValue();
		}
		
		throw new PrismLangException("Cannot evaluate to a double", this);
	}

	// evaluate to a boolean
	// any typing issues cause an exception
	// [does nothing to the expression itself]
	public boolean evaluateBoolean(Values constantValues, Values varValues) throws PrismLangException
	{
		Object o;
		
		o = evaluate(constantValues, varValues);
		
		if (!(o instanceof Boolean)) {
			throw new PrismLangException("Cannot evaluate to a boolean", this);
		}
		
		return ((Boolean)o).booleanValue();
	}
	
	// Static constructors for convenience
	public static Expression True() { return new ExpressionLiteral(Expression.BOOLEAN, true); }
	public static Expression False() { return new ExpressionLiteral(Expression.BOOLEAN, false); }
	public static Expression Int(int i) { return new ExpressionLiteral(Expression.INT, i); }
	public static Expression Double(double d) { return new ExpressionLiteral(Expression.DOUBLE, d); }
	public static Expression Not(Expression expr) { return new ExpressionUnaryOp(ExpressionUnaryOp.NOT, expr); }
	public static Expression Parenth(Expression expr) { return new ExpressionUnaryOp(ExpressionUnaryOp.PARENTH, expr); }
}

//------------------------------------------------------------------------------