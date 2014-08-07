//==============================================================================
//	
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import param.SymbolicEngine;
import param.TransitionList;
import parser.State;
import parser.Values;
import parser.ast.Expression;
import parser.ast.ExpressionConstant;
import parser.ast.ModulesFile;
import parser.visitor.ASTTraverse;
import prism.Pair;
import prism.PrismException;
import explicit.ModelExplorer;

public final class PSEModelExplorer implements ModelExplorer<Expression>
{
	private SymbolicEngine engine;
	private State currentState;
	private TransitionList transitionList;

	private BoxRegionFactory regionFactory;

	private Map<Expression, RateParametersAndPopulation> rateDataCache = new HashMap<Expression, RateParametersAndPopulation>();
	private Map<State, TransitionList> transitionsCache = new HashMap<State, TransitionList>();

	protected static class RateParametersAndPopulation extends Pair<Expression, Double>
	{
		public RateParametersAndPopulation(Expression first, Double second) {
			super(first, second);
		}

		public Expression getParameters()
		{
			return first;
		}

		public double getPopulation()
		{
			return second;
		}
	}

	public PSEModelExplorer(ModulesFile modulesFile) throws PrismException
	{
		modulesFile = (ModulesFile) modulesFile.deepCopy().replaceConstants(modulesFile.getConstantValues()).simplify();
		engine = new SymbolicEngine(modulesFile);
	}

	public void setParameters(String[] paramNames, double[] lower, double[] upper)
	{
		Values lowerParams = new Values();
		Values upperParams = new Values();
		for (int i = 0; i < paramNames.length; i++) {
			lowerParams.addValue(paramNames[i], lower[i]);
			upperParams.addValue(paramNames[i], upper[i]);
		}
		regionFactory = new BoxRegionFactory(lowerParams, upperParams);
	}

	public BoxRegionFactory getRegionFactory()
	{
		return regionFactory;
	}

	protected RateParametersAndPopulation extractRateParametersAndPopulation(Expression rateExpression) throws PrismException
	{
		if (rateDataCache.containsKey(rateExpression)) {
			return rateDataCache.get(rateExpression);
		}

		final List<ExpressionConstant> containedParameters = new ArrayList<ExpressionConstant>();
		rateExpression.accept(new ASTTraverse()
		{
			// TODO: visit() for doubles to handle rate population directly.
			// Subsequently, the whole method could be made static and moved
			// into pse.RateUtils or something.
			public Object visit(ExpressionConstant e)
			{
				containedParameters.add(e);
				return null;
			}
		});

		Expression rateParameters = Expression.Double(1);
		for (ExpressionConstant parameterExpr : containedParameters) {
			rateParameters = Expression.Times(rateParameters, parameterExpr);
		}
		double ratePopulation = Expression.Divide(rateExpression, rateParameters).evaluateDouble(regionFactory.completeSpace().getUpperBounds());
		RateParametersAndPopulation result = new RateParametersAndPopulation(rateParameters, ratePopulation);
		rateDataCache.put(rateExpression, result);
		return result;
	}

	protected RateParametersAndPopulation[] extractRateParametersAndPopulation(Expression[] rateExpressions) throws PrismException
	{
		if (rateExpressions == null) {
			return null;
		}
		int n = rateExpressions.length;
		RateParametersAndPopulation[] result = new RateParametersAndPopulation[n];
		for (int i = 0; i < n; i++) {
			result[i] = extractRateParametersAndPopulation(rateExpressions[i]);
		}
		return result;
	}

	protected ModulesFile getModulesFile()
	{
		return engine.getModulesFile();
	}

	@Override
	public State getDefaultInitialState() throws PrismException
	{
		return getModulesFile().getDefaultInitialState();
	}

	@Override
	public void queryState(State state) throws PrismException
	{
		currentState = state;
		if (transitionsCache.containsKey(state)) {
			transitionList = transitionsCache.get(state);
		} else {
			transitionList = engine.calculateTransitions(state);
			transitionsCache.put(state, transitionList);
		}
	}

	@Override
	public void queryState(State state, double time) throws PrismException
	{
		queryState(state);
	}

	@Override
	public int getNumChoices() throws PrismException
	{
		return transitionList.getNumChoices();
	}

	@Override
	public int getNumTransitions() throws PrismException
	{
		return transitionList.getNumTransitions();
	}

	@Override
	public int getNumTransitions(int choiceNr) throws PrismException
	{
		return transitionList.getChoice(choiceNr).size();
	}

	protected int getTotalIndexOfTransition(int choiceNr, int offset)
	{
		return transitionList.getTotalIndexOfTransition(choiceNr, offset);
	}

	@Override
	public String getTransitionAction(int choiceNr, int offset) throws PrismException
	{
		int a = transitionList.getTransitionModuleOrActionIndex(getTotalIndexOfTransition(choiceNr, offset));
		return a < 0 ? null : getModulesFile().getSynch(a - 1);
	}

	@Override
	public String getTransitionAction(int succNr) throws PrismException
	{
		int a = transitionList.getTransitionModuleOrActionIndex(succNr);
		return a < 0 ? null : getModulesFile().getSynch(a - 1);
	}

	@Override
	public Expression getTransitionProbability(int choiceNr, int offset) throws PrismException
	{
		return getTransitionProbability(getTotalIndexOfTransition(choiceNr, offset));
	}

	@Override
	public Expression getTransitionProbability(int succNr) throws PrismException
	{
		return transitionList.getTransitionProbability(succNr);
	}

	@Override
	public State computeTransitionTarget(int choiceNr, int offset) throws PrismException
	{
		return computeTransitionTarget(getTotalIndexOfTransition(choiceNr, offset));
	}

	@Override
	public State computeTransitionTarget(int succNr) throws PrismException
	{
		return transitionList.computeTransitionTarget(succNr, currentState);
	}

	public int getReaction(int succNr)
	{
		return transitionList.getChoiceOfTransition(succNr).hashCode();
	}
	
}
