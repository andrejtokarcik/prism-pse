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

import java.io.File;
import java.util.AbstractMap.SimpleImmutableEntry;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;

import parser.Values;
import prism.ModelType;
import prism.PrismException;
import prism.PrismLog;
import explicit.ModelExplicit;

/**
 * Represents a parametric Markov model.
 * This class is used to store all types of parametric model, both models
 * with and without nondeterminism, discrete- as well as continuous-time.
 * This turned out the be the most convenient way to implement model checking
 * for parametric models.
 */
final class ParamModel extends ModelExplicit
{
	/** total number of nondeterministic choices over all states */
	private int numTotalChoices;
	/** total number of probabilistic transitions over all states */
	private int numTotalTransitions;
	/** begin and end of state transitions */
	private int[] rows;
	/** begin and end of a distribution in a nondeterministic choice */
	private int[] choices;
	/** origins and targets of distribution branches */
	private int[] colsFrom;
	private int[] cols;
	/** hash codes of updates/reactions associated with distribution branches */
	private int[] reactions;
	/** rates of distribution branches */
	private Function[] rates;
	/** probabilities of distribution branches */
	private Function[] nonZeros;
	/** labels - per transition, <i>not</i> per action */
	private String[] labels;
	/** total sum of leaving rates for a given nondeterministic choice */
	private Function[] sumRates;
	/** model type */
	private ModelType modelType;
	/** function factory which manages functions used on transitions, etc. */
	private FunctionFactory functionFactory;
	/** */
	private Set<Integer> predecessorsViaReaction = new HashSet<Integer>();
	/** */
	// TODO convert into arrays like rows/choices/cols/... above?
	private Map<Integer, List<Integer>> inReactions;
	private Map<Integer, List<Entry<Integer, Integer>>> inoutReactions;
	private Map<Integer, List<Integer>> outReactions;
	
	/**
	 * Constructs a new parametric model.
	 */
	ParamModel()
	{
		numStates = 0;
		numTotalChoices = 0;
		numTotalTransitions = 0;
		initialStates = new LinkedList<Integer>();
		deadlocks = new TreeSet<Integer>();
	}

	/**
	 * Sets the type of the model.
	 * 
	 * @param modelType type the model shall have
	 */
	void setModelType(ModelType modelType)
	{
		this.modelType = modelType;
	}

	// Accessors (for Model)

	@Override
	public ModelType getModelType()
	{
		return modelType;
	}

	@Override
	public Values getConstantValues()
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public int getNumTransitions()
	{
		return numTotalTransitions;
	}

	@Override
	public Iterator<Integer> getSuccessorsIterator(int s)
	{
		throw new UnsupportedOperationException();
	}
	
	@Override
	public boolean isSuccessor(int s1, int s2)
	{
		for (int choice = stateBegin(s1); choice < stateEnd(s1); choice++) {
			for (int succ = choiceBegin(choice); succ < choiceEnd(choice); succ++) {
				if (succState(succ) == s2) {
					return true;
				}
			}
		}
		return false;
	}

	@Override
	public boolean allSuccessorsInSet(int s, BitSet set)
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean someSuccessorsInSet(int s, BitSet set)
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public void findDeadlocks(boolean fix) throws PrismException
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public void checkForDeadlocks() throws PrismException
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public void checkForDeadlocks(BitSet except) throws PrismException
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public void exportToPrismExplicit(String baseFilename) throws PrismException
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public void exportToPrismExplicitTra(String filename) throws PrismException
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public void exportToPrismExplicitTra(File file) throws PrismException
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public void exportToPrismExplicitTra(PrismLog log)
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public void exportToDotFile(String filename) throws PrismException
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public void exportToDotFile(String filename, BitSet mark) throws PrismException
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public void exportToPrismLanguage(String filename) throws PrismException
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public String infoString()
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public String infoStringTable()
	{
		String s = "";
		s += "States:      " + numStates + " (" + getNumInitialStates() + " initial)\n";
		s += "Transitions: " + getNumTransitions() + "\n";
		return s;
	}

	// Other
	
	public int getNumChoices(int state)
	{
		return stateEnd(state) - stateBegin(state);
	}

	public int getNumTotalChoices()
	{
		return numTotalChoices;
	}

	/**
	 * Allocates memory for subsequent construction of model. 
	 * 
	 * @param numStates number of states of the model
	 * @param numTotalChoices total number of nondeterministic choices of the model
	 * @param numTotalSuccessors total number of probabilistic transitions of the model
	 */
	void reserveMem(int numStates, int numTotalChoices, int numTotalSuccessors)
	{
		rows = new int[numStates + 1];
		choices = new int[numTotalChoices + 1];
		labels = new String[numTotalSuccessors];
		reactions = new int[numTotalSuccessors];
		rates = new Function[numTotalSuccessors];
		cols = new int[numTotalSuccessors];
		colsFrom = new int[numTotalSuccessors];
		nonZeros = new Function[numTotalSuccessors];
		sumRates = new Function[numTotalChoices];
	}

	/**
	 * Finish the current state.
	 * Starting with the 0th state, this function shall be called once all
	 * nondeterministic decisions of the current nth state have been added.
	 * Subsequent method calls of {@code finishChoice} and {@code addTransition}
	 * will then apply to the (n+1)th state. Notice that this method must be
	 * called for each state of the method, even the last one, once all its
	 * transitions have been added.
	 */
	void finishState()
	{
		rows[numStates + 1] = numTotalChoices;
		numStates++;
	}

	/**
	 * Finished the current nondeterministic choice.
	 * Subsequent calls of {@code addTransition} will in turn add probabilistic
	 * branches to the next nondeterministic choice. Notice that this method has
	 * to be called for all nondeterministic choices of a state, even the last
	 * one. Notice that DTMCs and CTMCs should only have a single
	 * nondeterministic choice per state.
	 */
	void finishChoice()
	{
		choices[numTotalChoices + 1] = numTotalTransitions;
		numTotalChoices++;
	}

	/**
	 * Adds a probabilistic branch to the current nondeterministic choice.
	 * Notice that by this function the probability to leave to a given state
	 * shall be specified, <i>not</i> the rate to this state. Instead, the
	 * sum of rates leaving a certain nondeterministic decision shall be
	 * specified using {@code setSumLeaving}.
	 * 
	 * @param reactionHash
	 * @param toState to which state the probabilistic choice leads
	 * @param rate
	 * @param probFn with which probability it leads to this state
	 * @param action action with which the choice is labelled
	 */
	void addTransition(int reactionHash, int fromState, int toState, Function rate, Function probFn, String action)
	{
		reactions[numTotalTransitions] = reactionHash;
		colsFrom[numTotalTransitions] = fromState;
		cols[numTotalTransitions] = toState;
		rates[numTotalTransitions] = rate;
		nonZeros[numTotalTransitions] = probFn;
		labels[numTotalTransitions] = action;
		numTotalTransitions++;
	}

	/**
	 * Sets the total sum of leaving rate of the current nondeterministic choice.
	 * For discrete-time models, this function shall always be called with
	 * {@code leaving = 1}.
	 * 
	 * @param leaving total sum of leaving rate of the current nondeterministic choice
	 */
	void setSumLeaving(Function leaving)
	{
		sumRates[numTotalChoices] = leaving;
	}

	/**
	 * Returns the number of the first nondeterministic choice of {@code state}.
	 * 
	 * @param state state to return number of first nondeterministic choice of
	 * @return number of first nondeterministic choice of {@code state}
	 */
	int stateBegin(int state)
	{
		return rows[state];
	}

	/**
	 * Returns the number of the last nondeterministic choice of {@code state} plus one.
	 * 
	 * @param state state to return number of last nondeterministic choice of
	 * @return number of first nondeterministic choice of {@code state} plus one
	 */
	int stateEnd(int state)
	{
		return rows[state + 1];
	}

	/**
	 * Returns the first probabilistic branch of the given nondeterministic decision.
	 * 
	 * @param choice choice of which to return the first probabilitic branch
	 * @return number of first probabilistic branch of {@choice}
	 */
	int choiceBegin(int choice)
	{
		return choices[choice];
	}

	/**
	 * Returns the last probabilistic branch of the given nondeterministic decision plus one.
	 * 
	 * @param choice choice of which to return the first probabilitic branch
	 * @return number of last probabilistic branch of {@choice} plus one
	 */
	int choiceEnd(int choice)
	{
		return choices[choice + 1];
	}

	/**
	 */
	int getReaction(int succNr)
	{
		return reactions[succNr];
	}

	/**
	 */
	Function succRate(int succNr)
	{
		return rates[succNr];
	}

	/**
	 * Returns the successor state of the given probabilistic branch.
	 * 
	 * @param succNr probabilistic branch to return successor state of
	 * @return state which probabilistic branch leads to
	 */
	int succState(int succNr)
	{
		return cols[succNr];
	}
	
	/**
	 */
	int currState(int succNr)
	{
		return colsFrom[succNr];
	}

	/**
	 * Returns the probability of the given probabilistic branch
	 * 
	 * @param succNr probabilistic branch to return probability of
	 * @return probability of given probabilistic branch
	 */
	Function succProb(int succNr)
	{
		return nonZeros[succNr];
	}

	/**
	 * Returns the label of the given probabilistic branch
	 * 
	 * @param succNr probabilistic branch to return label of
	 * @return label of given probabilistic branch
	 */
	String getLabel(int succNr)
	{
		return labels[succNr];
	}

	/**
	 * Returns the total sum of leaving rates of the given nondeterministic choice 
	 * 
	 * @param choice nondeterministic choice to return sum of leaving rates of
	 * @return sum of leaving rates of given nondeterministic choice
	 */
	Function sumLeaving(int choice)
	{
		return sumRates[choice];
	}

	/**
	 */
	double maxSumLeaving() throws PrismException
	{
		int i;
		double d, max = Double.NEGATIVE_INFINITY;
		for (i = 0; i < numTotalChoices; i++) {
			d = ((ExprFunction) sumLeaving(i)).evaluateAtUpper();
			if (d > max)
				max = d;
		}
		return max;
	}

	/**
	 * Instantiates the parametric model at a given point.
	 * All transition probabilities, etc. will be evaluated and set as
	 * probabilities of the concrete model at the given point. 
	 * 
	 * @param point point to instantiate model at
	 * @return nonparametric model instantiated at {@code point}
	 */
	ParamModel instantiate(Point point)
	{
		ParamModel result = new ParamModel();
		result.reserveMem(numStates, numTotalChoices, numTotalTransitions);
		result.initialStates = new LinkedList<Integer>(this.initialStates);
		for (int state = 0; state < numStates; state++) {
			for (int choice = stateBegin(state); choice < stateEnd(state); choice++) {
				for (int succ = choiceBegin(choice); succ < choiceEnd(choice); succ++) {
					result.addTransition(reactions[succ],
							colsFrom[succ],
							cols[succ],
							functionFactory.fromBigRational(rates[succ].evaluate(point)),
							functionFactory.fromBigRational(nonZeros[succ].evaluate(point)),
							labels[succ]);
				}
				result.setSumLeaving(functionFactory.fromBigRational(this.sumLeaving(choice).evaluate(point)));
				result.finishChoice();
			}
			result.finishState();
		}
		result.functionFactory = this.functionFactory;

		return result;
	}

	@Override
	public void buildFromPrismExplicit(String filename) throws PrismException
	{
		throw new UnsupportedOperationException();
	}

	/**
	 * Sets the function factory to be used in this parametric model.
	 * 
	 * @param functionFactory function factory to be used in this parametric model
	 */
	void setFunctionFactory(FunctionFactory functionFactory)
	{
		this.functionFactory = functionFactory;
	}

	/**
	 * Returns the function factory used in this parametric model.
	 * 
	 * @return function factory used in this parametric model
	 */
	FunctionFactory getFunctionFactory()
	{
		return functionFactory;
	}

	/**
	 */
	public boolean hasPredecessorViaReaction(int state, int reaction)
	{
		//if (!refreshCache) {
		if (predecessorsViaReaction.contains(state ^ reaction)) {
			return true;
		}
		//}
		
		for (int pred = 0; pred < getNumStates(); pred++) {
			for (int choice = stateBegin(pred); choice < stateEnd(pred); choice++) {
				for (int succ = choiceBegin(choice); succ < choiceEnd(choice); succ++) {
					if (getReaction(succ) != reaction)
						continue;
					if (succState(succ) != state)
						continue;
					predecessorsViaReaction.add(state ^ reaction);
					return true;
				}
			}
		}
		return false;
	}
	
	/**
	 * Computes the disjoint sets of in-flowing and/or out-flowing reactions.
	 */
	public void computeInOutReactions() throws PrismException
	{
		int pred, predChoice, predSucc, predReaction, state, choice, succ;
		boolean inout;
		
		assert inReactions == null;
		assert inoutReactions == null;
		assert outReactions == null;

		// Initialise the reaction sets
		inReactions = new HashMap<Integer, List<Integer>>(numStates);
		inoutReactions = new HashMap<Integer, List<Entry<Integer, Integer>>>(numStates);
		outReactions = new HashMap<Integer, List<Integer>>(numStates);
		for (state = 0; state < numStates; state++) {
			inReactions.put(state, new ArrayList<Integer>());
			inoutReactions.put(state, new ArrayList<Entry<Integer, Integer>>());
			outReactions.put(state, new ArrayList<Integer>());
		}
		
		// Populate the sets with transition indices
		for (pred = 0; pred < numStates; pred++) {
			for (predChoice = stateBegin(pred); predChoice < stateEnd(pred); predChoice++) {
				for (predSucc = choiceBegin(predChoice); predSucc < choiceEnd(predChoice); predSucc++) {
					inout = false;
					predReaction = getReaction(predSucc);
					state = succState(predSucc);

					for (choice = stateBegin(state); choice < stateEnd(state); choice++) {
						for (succ = choiceBegin(choice); succ < choiceEnd(choice); succ++) {
							if (getReaction(succ) == predReaction) {
								inout = true;
								((ExprFunction) succRate(predSucc)).computePopulation(false);
								((ExprFunction) succRate(predSucc)).computeParametersMultiplied(false);
								((ExprFunction) succRate(succ)).computePopulation(false);
								((ExprFunction) succRate(succ)).computeParametersMultiplied(false);
								inoutReactions.get(state).add(new SimpleImmutableEntry<Integer, Integer>(predSucc, succ));

								// TODO: Perhaps we can break the two innermost for-loops from here?
								// I.e., is `state` guaranteed not to have another succ with this reaction?
							}
						}
					}

					if (!inout) {
						inReactions.get(state).add(predSucc);
					}

					if (!hasPredecessorViaReaction(pred, predReaction)) {
						outReactions.get(pred).add(predSucc);
					}
				}
			}
		}		
	}
	
	public List<Integer> getInReactions(int state)
	{
		/*
		if (inReactions == null) {
			computeInOutReactions();
		}
		*/
		assert inReactions != null;
		return inReactions.get(state);
	}
	
	public List<Entry<Integer, Integer>> getInoutReactions(int state)
	{
		/*
		if (inoutReactions == null) {
			computeInOutReactions();
		}
		*/
		assert inoutReactions != null;
		return inoutReactions.get(state);
	}
	
	public List<Integer> getOutReactions(int state)
	{
		/*
		if (outReactions == null) {
			computeInOutReactions();
		}
		*/
		assert outReactions != null;
		return outReactions.get(state);
	}

	public void vmMult(double vectMin[], double resultMin[], double vectMax[], double resultMax[], double qmax) throws PrismException
	{
		int pred, state;
		ExprFunction predRate, rate;
		ExprFunction predRateParams;
		double predPopulation, population;
		double midSumNumeratorMin, midSumNumeratorMax;

		for (state = 0; state < numStates; state++) {
			// Initialise the result
			resultMin[state] = vectMin[state];
			resultMax[state] = vectMax[state];

			// Incoming reactions
			for (int trans : getInReactions(state)) {
				rate = (ExprFunction) succRate(trans);
				pred = currState(trans);
				resultMin[state] += rate.evaluateAtLower() * vectMin[pred] / qmax;
				resultMax[state] += rate.evaluateAtUpper() * vectMax[pred] / qmax;
			}

			// Outgoing reactions
			for (int trans : getOutReactions(state)) {
				rate = (ExprFunction) succRate(trans);
				resultMin[state] -= rate.evaluateAtUpper() * vectMin[state] / qmax;
				resultMax[state] -= rate.evaluateAtLower() * vectMax[state] / qmax;
			}

			// Both incoming and outgoing
			for (Entry<Integer, Integer> transs : getInoutReactions(state)) {
				pred = currState(transs.getKey());
				predRate = (ExprFunction) succRate(transs.getKey());
				predPopulation = predRate.getPopulation();
				predRateParams = predRate.getParametersMultiplied();

				rate = (ExprFunction) succRate(transs.getValue());
				population = rate.getPopulation();
				// The rate params assumed to be the same for both `pred` and `state`
				assert predRateParams.equals(rate.getParametersMultiplied());

				midSumNumeratorMin = vectMin[pred] * predPopulation - vectMin[state] * population;
				if (midSumNumeratorMin > 0) {
					resultMin[state] += predRateParams.evaluateAtLower() * midSumNumeratorMin / qmax;
				} else {
					resultMin[state] += predRateParams.evaluateAtUpper() * midSumNumeratorMin / qmax;
				}

				midSumNumeratorMax = vectMax[pred] * predPopulation - vectMax[state] * population;
				if (midSumNumeratorMax > 0) {
					resultMax[state] += predRateParams.evaluateAtUpper() * midSumNumeratorMax / qmax;
				} else {
					resultMax[state] += predRateParams.evaluateAtLower() * midSumNumeratorMax / qmax;
				}
			}
		}
	}
}
