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

import parser.ast.Expression;
import prism.PrismException;

/**
 */
final class ExprFunction extends Function {
	/** Expression the function is wrapping */
	private Expression expr;
	/** type of function (rational function, infinity, etc.) */
	int type;
	final static int NORMAL = 0;
	final static int INF = 1;
	final static int MINF = 2;
	final static int NAN = 3;
	/** */
	private Expression parametersMultipliedExpr;
	/** */
	private ExprFunction parametersMultiplied = null;
	private Double population = null;
	private Double lower = null;
	private Double upper = null;
	
	// constructors
	
	/**
	 * Creates a new Expression-based function.
	 */
	ExprFunction(ExprFunctionFactory functionContext, Expression expr, int type, Expression parametersMultipliedExpr)
	{
		super(functionContext);
		this.expr = expr;
		this.type = type;
		this.parametersMultipliedExpr = parametersMultipliedExpr;
	}

	@Override
	public String toString()
	{
		if (isNaN()) {
			return "NaN";
		} else if (isInf()) {
			return "Inf";
		} else if (isMInf()) {
			return "MInf";
		}
		return expr.toString();
	}
	
	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof ExprFunction)) {
			return false;
		}
		ExprFunction function = (ExprFunction) obj;
		if (isNaN()) {
			return function.isNaN();
		}
		if (isInf()) {
			return function.isInf();
		}
		if (isMInf()) {
			return function.isMInf();
		}
		return expr.equals(function.expr);
	}
	
	@Override
	public int hashCode() {
		return expr.hashCode();
	}
	
	/**
	 */
	Expression getExpression()
	{
		return expr;
	}

	/**
	 */
	public void computeParametersMultiplied(boolean force) throws PrismException
	{
		if (parametersMultiplied == null || force) {
			parametersMultiplied = ((ExprFunctionFactory) factory).fromExpression(parametersMultipliedExpr);
		}
	}
	
	ExprFunction getParametersMultiplied()
	{
		/*
		if (parametersMultiplied == null) {
			parametersMultiplied = ((ExprFunctionFactory) factory).fromExpression(parametersMultipliedExpr);
		}
		*/
		assert parametersMultiplied != null;
		return parametersMultiplied;
	}

	/**
	 */
	public void computePopulation(boolean force) throws PrismException
	{
		if (population == null || force) {
			// Could call evaluateAtLower() as well
			population = ((ExprFunctionFactory) factory).fromExpression(Expression.Divide(expr, parametersMultipliedExpr)).evaluateAtUpper();
		}
	}
	
	double getPopulation()
	{
		/*
		if (population == null) {
			// Could call evaluateAtLower() as well
			population = ((ExprFunctionFactory) factory).fromExpression(Expression.Divide(expr, parametersMultipliedExpr)).evaluateAtUpper();
		}
		*/
		assert population != null;
		return population;
	}

	@Override
	public Function add(Function other)
	{
		if (this.isNaN() || other.isNaN()) {
			return factory.getNaN();
		}
		if (this.isInf() || other.isInf()) {
			if (this.isMInf() || other.isMInf()) {
				return factory.getZero();
			}
			return factory.getInf();
		}
		if (this.isMInf() || other.isMInf()) {
			return factory.getMInf();
		}
		return ((ExprFunctionFactory) factory).fromExpression(Expression.Plus(expr, ((ExprFunction) other).expr));
	}

	@Override
	public Function negate()
	{
		if (this.isNaN()) {
			return factory.getNaN();			
		}
		if (this.isInf()) {
			return factory.getMInf();
		}
		if (this.isMInf()) {
			return factory.getMInf();			
		}
		return ((ExprFunctionFactory) factory).fromExpression(Expression.Minus(((ExprFunction) factory.getZero()).getExpression(), expr));
	}
	
	@Override
	public Function multiply(Function other)
	{
		if (this.isNaN() || other.isNaN()) {
			return factory.getNaN();
		}
		if (this.isZero() || other.isZero()) {
			return factory.getZero();
		}
		if (this.isInf() || other.isInf()) {
			if (this.isMInf() || other.isMInf()) {
				return factory.getMInf();
			} else {
				return factory.getInf();
			}
		}
		return ((ExprFunctionFactory) factory).fromExpression(Expression.Times(expr, ((ExprFunction) other).expr));
	}
	
	public Function multiply(double byNumber)
	{
		Function byFunction = ((ExprFunctionFactory) factory).fromDouble(byNumber);
		return multiply(byFunction);
	}
	
	@Override
	public Function divide(Function other)
	{
		if (this.isNaN() || other.isNaN()) {
			return factory.getNaN();
		}
		return ((ExprFunctionFactory) factory).fromExpression(Expression.Divide(expr, ((ExprFunction) other).expr));
	}
	
	public Function divide(double byNumber)
	{
		Function byFunction = ((ExprFunctionFactory) factory).fromDouble(byNumber);
		return divide(byFunction);
	}
	
	@Override
	public Function star()
	{
		throw new UnsupportedOperationException();
	}
	
	@Override
	public Function toConstraint() {
		return this;
	}

	@Override
	public BigRational evaluate(Point point, boolean cancel) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public BigRational evaluate(Point point) {
		return evaluate(point, true);
	}
	
	public double evaluateDouble(Point point) throws PrismException {
		parser.Values constantValues = new parser.Values();
		double[] doubleValues = point.doubleValues();
		for (int i = 0; i < factory.getNumVariables(); i++) {
			constantValues.addValue(factory.getParameterName(i), doubleValues[i]);
		}
		return expr.evaluateDouble(constantValues);
	}
	
	public double evaluateAtLower() throws PrismException {
		if (lower == null) {
			lower = evaluateDouble(new Point(factory.getLowerBounds()));
		}
		return lower;
	}
	
	public double evaluateAtUpper() throws PrismException {
		if (upper == null) {
			upper = evaluateDouble(new Point(factory.getUpperBounds()));
		}
		return upper;
	}
	
	@Override
	public boolean check(Point point, boolean strict)
	{
		BigRational value = evaluate(point, false);
		int compare = value.signum();
		return strict ? (compare > 0) : (compare >= 0);
	}

	@Override
	public BigRational asBigRational() {
		if (isNaN()) {
			return BigRational.NAN;
		} else if (isInf()) {
			return BigRational.INF;
		} else if (isMInf()) {
			return BigRational.MINF;
		}
		BigRational[] point = new BigRational[factory.getNumVariables()];
		for (int dim = 0; dim < factory.getNumVariables(); dim++) {
			point[dim] = new BigRational(0);
		}
		return evaluate(new Point(point));
	}

	@Override
	public boolean isNaN() {
		return type == NAN;
	}

	@Override
	public boolean isInf() {
		return type == INF;
	}

	@Override
	public boolean isMInf() {
		return type == MINF;
	}

	@Override
	public boolean isOne() {
		if (type != NORMAL) {
			return false;
		}
		return equals(factory.getOne());
	}

	@Override
	public boolean isZero() {
		if (type != NORMAL) {
			return false;
		}
		return equals(factory.getZero());
	}
}
