//==============================================================================
//	
//	Copyright (c) 2013-
//	Authors:
//	* Ernst Moritz Hahn <emhahn@cs.ox.ac.uk> (University of Oxford)
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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import param.ChoiceListFlexi;
import param.SymbolicEngine;
import param.TransitionList;
import parser.State;
import parser.Values;
import parser.ast.Expression;
import parser.ast.ExpressionConstant;
import parser.ast.ModulesFile;
import parser.visitor.ASTTraverse;
import prism.ModelType;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismLangException;
import explicit.IndexedSet;
import explicit.StateStorage;

public final class ModelBuilder extends PrismComponent
{
	/** {@code ModulesFile} to be transformed to a {@code ParamModel} */
	private ModulesFile modulesFile;
	/** parametric model constructed from {@code modulesFile} */
	private PSEModel model;
	/** bounds of parameters */
	private Values paramsLower;
	private Values paramsUpper;
	/** */
	private BoxRegionFactory regionFactory;
	/** */
	private Map<State, TransitionList> transitionsCache = new HashMap<State, TransitionList>();
	private Map<Expression, Double> ratePopulationsCache = new HashMap<Expression, Double>();
	private Map<Expression, Expression> rateParamsCache = new HashMap<Expression, Expression>();

	/**
	 * Constructor
	 */
	public ModelBuilder(PrismComponent parent) throws PrismException
	{
		super(parent);
	}
	
	// setters and getters

	/**
	 * Set modules file to be transformed to parametric Markov model.
	 * 
	 * @param modulesFile modules file to be transformed to parametric Markov model
	 */
	public void setModulesFile(ModulesFile modulesFile)
	{
		this.modulesFile = modulesFile;
	}

	/**
	 * Set parameter informations.
	 * Obviously, all of {@code paramNames}, {@code lower}, {@code} upper
	 * must have the same length, and {@code lower} bounds of parameters must
	 * not be higher than {@code upper} bounds.
	 * 
	 * @param paramNames names of parameters
	 * @param lower lower bounds of parameters
	 * @param upper upper bounds of parameters
	 */
	public void setParameters(String[] paramNames, double[] lower, double[] upper)
	{
		paramsLower = new Values();
		paramsUpper = new Values();
		for (int i = 0; i < paramNames.length; i++) {
			paramsLower.addValue(paramNames[i], lower[i]);
			paramsUpper.addValue(paramNames[i], upper[i]);
		}

		regionFactory = new BoxRegionFactory(paramsLower, paramsUpper);
	}

	/**
	 * Construct parametric Markov model.
	 * For this to work, module file, PRISM log, etc. must have been set
	 * beforehand.
	 * 
	 * @throws PrismException in case the model cannot be constructed
	 */
	public void build() throws PrismException
	{
		long time;

		mainLog.print("\nBuilding model...\n");

		// build model
		time = System.currentTimeMillis();
		modulesFile = (ModulesFile) modulesFile.deepCopy().replaceConstants(modulesFile.getConstantValues()).simplify();
		PSEModel modelExpl = constructModel(modulesFile);
		modelExpl.scaleParameterSpace(regionFactory.completeSpace());
		time = System.currentTimeMillis() - time;

		mainLog.println("\nTime for model construction: " + time / 1000.0 + " seconds.");
		model = modelExpl;
	}

	/**
	 * Returns the constructed parametric Markov model.
	 * 
	 * @return constructed parametric Markov model
	 */
	public explicit.Model getModel()
	{
		return model;
	}

	/**
	 */
	public BoxRegionFactory getRegionFactory()
	{
		return regionFactory;
	}

	/**
	 * Reserves memory needed for parametric model and reserves necessary space.
	 * Afterwards, transition probabilities etc. can be added. 
	 * 
	 * @param modulesFile modules file of which to explore states
	 * @param model model in which to reserve memory
	 * @param modelType type of the model to construct
	 * @param engine the engine used to compute state successors, etc.
	 * @param states list of states to be filled by this method
	 * @throws PrismException thrown if problems in underlying methods occur
	 */
	private void reserveMemoryAndExploreStates(ModulesFile modulesFile, PSEModel model, ModelType modelType, SymbolicEngine engine, StateStorage<State> states)
			throws PrismException
	{
		int numStates = 0;
		int numTotalSuccessors = 0;

		LinkedList<State> explore = new LinkedList<State>();

		State state = modulesFile.getDefaultInitialState();
		states.add(state);
		explore.add(state);
		numStates++;

		while (!explore.isEmpty()) {
			state = explore.removeFirst();
			TransitionList tranlist = engine.calculateTransitions(state);
			transitionsCache.put(state, tranlist);
			
			int numChoices = tranlist.getNumChoices();
			for (int choiceNr = 0; choiceNr < numChoices; choiceNr++) {
				int numSuccessors = tranlist.getChoice(choiceNr).size();
				numTotalSuccessors += numSuccessors;
				for (int succNr = 0; succNr < numSuccessors; succNr++) {
					State stateNew = tranlist.getChoice(choiceNr).computeTarget(succNr, state);
					if (states.add(stateNew)) {
						numStates++;
						explore.add(stateNew);
					}
				}
			}
			if (numChoices == 0) {
				numTotalSuccessors++;
			}
		}

		model.reserveMem(numStates, numTotalSuccessors);
	}

	/**
	 */
	private Expression getRateParams(Expression rate) throws PrismException
	{
		if (rateParamsCache.containsKey(rate)) {
			return rateParamsCache.get(rate);
		}

		final List<ExpressionConstant> containedParameters = new ArrayList<ExpressionConstant>();
		rate.accept(new ASTTraverse()
		{
			public Object visit(ExpressionConstant e) throws PrismLangException
			{
				containedParameters.add(e);
				return null;
			}
		});

		Expression parametersMultiplied = Expression.Double(1);
		for (ExpressionConstant parameterExpr : containedParameters) {
			parametersMultiplied = Expression.Times(parametersMultiplied, parameterExpr);
		}
		rateParamsCache.put(rate, parametersMultiplied);
		return parametersMultiplied;
	}
	
	private double getRatePopulation(Expression rate) throws PrismException
	{
		if (ratePopulationsCache.containsKey(rate)) {
			return ratePopulationsCache.get(rate);
		}
		Expression rateParams = getRateParams(rate);
		// It is less safe to evaluate the rate population with paramsLower
		// since a lower bound may be reasonably set to zero, which could in turn
		// result in lots of NaNs due to division by zero.
		double ratePopulation = Expression.Divide(rate, rateParams).evaluateDouble(paramsUpper);
		ratePopulationsCache.put(rate, ratePopulation);
		return ratePopulation;
	}

	/**
	 * Construct model once function factory etc. has been allocated.
	 * 
	 * @param modulesFile modules file of which to construct parametric model
	 * @return parametric model constructed
	 * @throws PrismException thrown if model cannot be constructed
	 */
	private PSEModel constructModel(ModulesFile modulesFile) throws PrismException
	{
		ModelType modelType;
		PSEModel model;

		if (modulesFile.getInitialStates() != null) {
			throw new PrismException("Cannot do explicit-state reachability if there are multiple initial states");
		}

		modelType = modulesFile.getModelType();
		if (modelType != ModelType.CTMC) {
			throw new PrismException("Unsupported model type: " + modelType);
		}

		mainLog.print("\nComputing reachable states...");
		mainLog.flush();
		long timer = System.currentTimeMillis();

		model = new PSEModel();
		model.setModelType(modelType);

		SymbolicEngine engine = new SymbolicEngine(modulesFile);

		if (modulesFile.getInitialStates() != null) {
			throw new PrismException("Explicit model construction does not support multiple initial states");
		}

		StateStorage<State> states = new IndexedSet<State>(true);
		reserveMemoryAndExploreStates(modulesFile, model, modelType, engine, states);
		int[] permut = states.buildSortingPermutation();
		List<State> statesList = states.toPermutedArrayList(permut);
		model.setStatesList(statesList);
		model.addInitialState(permut[0]);
		for (State state : statesList) {
			TransitionList tranlist = transitionsCache.get(state);
			int numChoices = tranlist.getNumChoices();
			Expression sumOut = Expression.Double(0);
			for (int choiceNr = 0; choiceNr < numChoices; choiceNr++) {
				ChoiceListFlexi choice = tranlist.getChoice(choiceNr);
				int a = tranlist.getTransitionModuleOrActionIndex(tranlist.getTotalIndexOfTransition(choiceNr, 0));
				String action = a < 0 ? null : modulesFile.getSynch(a - 1);

				// CTMCs should only have a single nondeterministic choice per state
				assert choice.size() == 1;

				State stateNew = choice.computeTarget(0, state);
				Expression rateExpr = choice.getProbability(0);
				model.addTransition(choice.hashCode(), permut[states.get(state)], permut[states.get(stateNew)], getRateParams(rateExpr), getRatePopulation(rateExpr), action);
				sumOut = Expression.Plus(sumOut, rateExpr);
			}
			model.setSumLeaving(sumOut.evaluateDouble(paramsUpper));
			model.finishState();
		}

		mainLog.println();

		mainLog.print("Reachable states exploration and model construction");
		mainLog.println(" done in " + ((System.currentTimeMillis() - timer) / 1000.0) + " secs.");

		return model;
	}
}
