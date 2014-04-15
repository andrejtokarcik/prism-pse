//==============================================================================
//	
//	Copyright (c) 2013-
//	Authors:
//	* Ernst Moritz Hahn <emhahn@cs.ox.ac.uk> (University of Oxford)
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

package param;

import java.util.ArrayList;
import java.util.List;

import parser.ast.Expression;
import parser.ast.ExpressionConstant;
import parser.visitor.ASTTraverse;
import prism.PrismException;
import prism.PrismLangException;

/**
 */
final class ExprFunctionFactory extends FunctionFactory {
	private ExprFunction zero;
	private ExprFunction one;
	private ExprFunction nan;
	private ExprFunction inf;
	private ExprFunction minf;
	private ExprFunction[] parameters;
	
	private static List<ExpressionConstant> containedParameters;

	/**
	 * Creates a new function factory.
	 * 
	 * @param parameterNames names of parameters
	 * @param lowerBounds lower bounds of parameters
	 * @param upperBounds upper bounds of parameters
	 */
	ExprFunctionFactory(String[] paramNames, BigRational[] lowerBounds, BigRational[] upperBounds)
	{
		super(paramNames, lowerBounds, upperBounds);

		Expression defaultParametersMultiplied = Expression.Double(1);
		one = new ExprFunction(this, Expression.Double(1), ExprFunction.NORMAL, defaultParametersMultiplied);
		zero = new ExprFunction(this, Expression.Double(0), ExprFunction.NORMAL, defaultParametersMultiplied);
		nan = new ExprFunction(this, Expression.Double(0), ExprFunction.NAN, defaultParametersMultiplied);
		inf = new ExprFunction(this, Expression.Double(0), ExprFunction.INF, defaultParametersMultiplied);
		minf = new ExprFunction(this, Expression.Double(0), ExprFunction.MINF, defaultParametersMultiplied);
		parameters = new ExprFunction[paramNames.length];
		for (int i = 0; i < paramNames.length; i++) {
			Expression paramExpr = Expression.ConstantDouble(paramNames[i]);
			parameters[i] = new ExprFunction(this, paramExpr, ExprFunction.NORMAL, defaultParametersMultiplied);
		}
	}	

	@Override
	public Function getOne()
	{
		return one;
	}
	
	@Override
	public Function getZero()
	{
		return zero;
	}

	@Override
	public Function getNaN()
	{
		return nan;
	}

	@Override
	public Function getInf()
	{
		return inf;
	}
	
	@Override
	public Function getMInf()
	{
		return minf;
	}

	@Override
	public Function fromBigRational(BigRational from)
	{
		throw new UnsupportedOperationException();
	}
	
	@Override
	Function getVar(int var)
	{
		return parameters[var];
	}
	
	/**
	 */
	public ExprFunction fromExpression(Expression expr)
	{
		try {
			computeContainedParameters(expr);
		} catch (PrismException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		Expression parametersMultiplied = ((ExprFunction) getOne()).getExpression();
		for (ExpressionConstant parameterExpr : containedParameters) {
			parametersMultiplied = Expression.Times(parametersMultiplied, parameterExpr);
		}
		return new ExprFunction(this, expr, ExprFunction.NORMAL, parametersMultiplied);
	}

	/**
	 */
	private static void computeContainedParameters(Expression expr) throws PrismException
	{
		containedParameters = new ArrayList<ExpressionConstant>();
		expr.accept(new ASTTraverse()
		{
			public Object visit(ExpressionConstant e) throws PrismLangException
			{
				containedParameters.add(e);
				return null;
			}
		});
	}
	
	/**
	 */
	public ExprFunction fromDouble(double d)
	{
		return fromExpression(Expression.Double(d));
	}
}
