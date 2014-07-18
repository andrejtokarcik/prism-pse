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

import java.io.File;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import parser.Values;
import parser.ast.Expression;
import prism.ModelType;
import prism.Pair;
import prism.PrismException;
import prism.PrismLog;
import explicit.ModelExplicit;

public final class PSEModel extends ModelExplicit
{
	/** total number of probabilistic transitions over all states */
	private int numTotalTransitions;
	/** begin and end of state transitions */
	private int[] rows;
	/** origins and targets of distribution branches */
	private int[] colsFrom;
	private int[] colsTo;
	/** */
	private Expression[] rateParams;
	private double[] rateParamsLowers;
	private double[] rateParamsUppers;
	/** */
	private double[] ratePopulations;
	/** */
	private boolean[] parametrisedTransitions;
	/** reactions - classes of transitions */
	private int[] reactions;
	/** labels - per transition, <i>not</i> per action */
	private String[] labels;
	/** total sum of leaving rates for a state */
	private double[] exitRates;
	/** model type */
	private ModelType modelType;
	/** */
	private Set<Integer> predecessorsViaReaction = new HashSet<Integer>();
	/** */
	private Map<Integer, List<Integer>> inReactions;
	private Map<Integer, List<Pair<Integer, Integer>>> inoutReactions;
	private Map<Integer, List<Integer>> outReactions;
	
	/**
	 * Constructs a new parametric model.
	 */
	PSEModel()
	{
		numStates = 0;
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
		for (int trans = stateBegin(s1); trans < stateEnd(s1); trans++) {
			if (toState(trans) == s2)
				return true;
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
	public void buildFromPrismExplicit(String filename) throws PrismException
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
	public void exportTransitionsToDotFile(int i, PrismLog out)
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

	/**
	 * Allocates memory for subsequent construction of model. 
	 * 
	 * @param numStates number of states of the model
	 * @param numTotalTransitions total number of probabilistic transitions of the model
	 */
	void reserveMem(int numStates, int numTotalTransitions)
	{
		rows = new int[numStates + 1];
		labels = new String[numTotalTransitions];
		reactions = new int[numTotalTransitions];
		rateParams = new Expression[numTotalTransitions];
		rateParamsLowers = new double[numTotalTransitions];
		rateParamsUppers = new double[numTotalTransitions];
		parametrisedTransitions = new boolean[numTotalTransitions];
		ratePopulations = new double[numTotalTransitions];
		colsTo = new int[numTotalTransitions];
		colsFrom = new int[numTotalTransitions];
		exitRates = new double[numStates];
	}

	/**
	 * Finish the current state.
	 * Starting with the 0th state, this function shall be called once all
	 * nondeterministic decisions of the current nth state have been added.
	 * Subsequent method calls of {@code addTransition}
	 * will then apply to the (n+1)th state. Notice that this method must be
	 * called for each state of the method, even the last one, once all its
	 * transitions have been added.
	 */
	void finishState()
	{
		rows[numStates + 1] = numTotalTransitions;
		numStates++;
	}

	/**
	 * Adds a probabilistic transition from the current state.
	 */
	void addTransition(int reaction, int fromState, int toState, Expression rateParamsExpr, double ratePopulation, String action)
	{
		reactions[numTotalTransitions] = reaction;
		colsFrom[numTotalTransitions] = fromState;
		colsTo[numTotalTransitions] = toState;
		rateParams[numTotalTransitions] = rateParamsExpr;
		ratePopulations[numTotalTransitions] = ratePopulation;
		labels[numTotalTransitions] = action;

		predecessorsViaReaction.add(toState ^ reaction);

		numTotalTransitions++;
	}

	/**
	 * Sets the total sum of leaving rates from the current state.
	 */
	void setSumLeaving(double leaving)
	{
		exitRates[numStates] = leaving;
	}

	/**
	 * Returns the number of the first transition of {@code state}.
	 */
	int stateBegin(int state)
	{
		return rows[state];
	}

	/**
	 * Returns the number of the last transition of {@code state} plus one.
	 */
	int stateEnd(int state)
	{
		return rows[state + 1];
	}

	/**
	 */
	boolean parametrised(int trans)
	{
		return parametrisedTransitions[trans];
	}

	/**
	 */
	int getReaction(int trans)
	{
		return reactions[trans];
	}

	/**
	 */
	int fromState(int trans)
	{
		return colsFrom[trans];
	}

	/**
	 * Returns the successor state of the given transition.
	 */
	int toState(int trans)
	{
		return colsTo[trans];
	}

	/**
	 * Returns the label of the given transition.
	 */
	String getLabel(int trans)
	{
		return labels[trans];
	}
	
	/**
	 */
	double getMaxExitRate()
	{
		BitSet allStates = new BitSet(numStates);
		allStates.set(0, numStates - 1);
		return getMaxExitRate(allStates);
	}

	double getMaxExitRate(BitSet subset)
	{
		double max = Double.NEGATIVE_INFINITY;
		for (int state = subset.nextSetBit(0); state >= 0; state = subset.nextSetBit(state + 1)) {
			if (exitRates[state] > max)
				max = exitRates[state];
		}
		return max;
	}

	double getDefaultUniformisationRate()
	{
		return 1.02 * getMaxExitRate();
	}

	double getDefaultUniformisationRate(BitSet nonAbs)
	{
		return 1.02 * getMaxExitRate(nonAbs);
	}

	/**
	 */
	public void computeInOutReactions() throws PrismException
	{
		if (inReactions != null && inoutReactions != null && outReactions != null)
			return;

		// Initialise the reaction sets
		inReactions = new HashMap<Integer, List<Integer>>(numStates);
		inoutReactions = new HashMap<Integer, List<Pair<Integer, Integer>>>(numStates);
		outReactions = new HashMap<Integer, List<Integer>>(numStates);
		for (int state = 0; state < numStates; state++) {
			inReactions.put(state, new LinkedList<Integer>());
			inoutReactions.put(state, new LinkedList<Pair<Integer, Integer>>());
			outReactions.put(state, new LinkedList<Integer>());
		}

		// Populate the sets with transition indices
		for (int pred = 0; pred < numStates; pred++) {
			for (int predTrans = stateBegin(pred); predTrans < stateEnd(pred); predTrans++) {
				if (!parametrised(predTrans))
					continue;

				boolean inout = false;
				int predReaction = getReaction(predTrans);
				int state = toState(predTrans);

				for (int trans = stateBegin(state); trans < stateEnd(state); trans++) {
					if (getReaction(trans) == predReaction) {
						inout = true;
						inoutReactions.get(state).add(new Pair<Integer, Integer>(predTrans, trans));
						break;
					}
				}

				if (!inout)
					inReactions.get(state).add(predTrans);

				if (!predecessorsViaReaction.contains(pred ^ predReaction))
					outReactions.get(pred).add(predTrans);
			}
		}
	}

	public void vmMult(double vectMin[], double resultMin[], double vectMax[], double resultMax[], double q)
			throws PrismException
	{
		for (int state = 0; state < numStates; state++) {
			// Initialise the result
			resultMin[state] = vectMin[state];
			resultMax[state] = vectMax[state];

			// Incoming reactions
			for (int trans : inReactions.get(state)) {
				int pred = fromState(trans);
				resultMin[state] += rateParamsLowers[trans] * ratePopulations[trans] * vectMin[pred] / q;
				resultMax[state] += rateParamsUppers[trans] * ratePopulations[trans] * vectMax[pred] / q;
			}

			// Outgoing reactions
			for (int trans : outReactions.get(state)) {
				resultMin[state] -= rateParamsUppers[trans] * ratePopulations[trans] * vectMin[state] / q;
				resultMax[state] -= rateParamsLowers[trans] * ratePopulations[trans] * vectMax[state] / q;
			}

			// Both incoming and outgoing
			for (Pair<Integer, Integer> transs : inoutReactions.get(state)) {
				int predTrans = transs.first;
				int trans = transs.second;

				int pred = fromState(predTrans);
				assert fromState(trans) == state;

				// The rate params of the two considered transitions must be identical
				assert rateParamsLowers[predTrans] == rateParamsLowers[trans] && rateParamsUppers[predTrans] == rateParamsUppers[trans];

				double midSumNumeratorMin = ratePopulations[predTrans] * vectMin[pred] - ratePopulations[trans] * vectMin[state];
				if (midSumNumeratorMin > 0.0) resultMin[state] += rateParamsLowers[trans] * midSumNumeratorMin / q;
				else resultMin[state] += rateParamsUppers[trans] * midSumNumeratorMin / q;

				double midSumNumeratorMax = ratePopulations[predTrans] * vectMax[pred] - ratePopulations[trans] * vectMax[state];
				if (midSumNumeratorMax > 0.0) resultMax[state] += rateParamsUppers[trans] * midSumNumeratorMax / q;
				else resultMax[state] += rateParamsLowers[trans] * midSumNumeratorMax / q;
			}
		}

		// Optimisation: Non-parametrised transitions
		for (int trans = 0; trans < numTotalTransitions; trans++) {
			if (parametrised(trans))
				continue;

			int pred = fromState(trans);
			int state = toState(trans);

			double rate = rateParamsLowers[trans] * ratePopulations[trans];

			resultMin[pred] -= rate * vectMin[pred] / q;
			resultMax[pred] -= rate * vectMax[pred] / q;

			resultMin[state] += rate * vectMin[pred] / q;
			resultMax[state] += rate * vectMax[pred] / q;
		}
	}

	private double mvMultMidSumEvalMin(int trans, double vectMinPred, double vectMinState, double q)
	{
		double midSumNumeratorMin = ratePopulations[trans] * vectMinPred - ratePopulations[trans] * vectMinState;
		if (midSumNumeratorMin > 0.0) return rateParamsLowers[trans] * midSumNumeratorMin / q;
		else return rateParamsUppers[trans] * midSumNumeratorMin / q;
	}

	private double mvMultMidSumEvalMax(int trans, double vectMaxPred, double vectMaxState, double q)
	{
		double midSumNumeratorMax = ratePopulations[trans] * vectMaxPred - ratePopulations[trans] * vectMaxState;
		if (midSumNumeratorMax > 0.0) return rateParamsUppers[trans] * midSumNumeratorMax / q;
		else return rateParamsLowers[trans] * midSumNumeratorMax / q;
	}

	/**
	 * NB: Semantics of mvMult() is *not* analogical to that of vmMult(), the difference
	 * is crucial:  result[k]_i in vmMult() is simply the probability of being in state k
	 * after i iterations starting from an initial state.  On the other hand, mvMult()'s
	 * result[k]_i denotes the probability that an absorbing state (i.e., a state
	 * in the complement of subset) is reached after i iterations starting from k.
	 */
	public void mvMult(double vectMin[], double resultMin[], double vectMax[], double resultMax[], BitSet subset, boolean complement, double q)
			throws PrismException
	{
		if (subset == null) {
			// Loop over all states
			subset = new BitSet(numStates);
			subset.set(0, numStates - 1);
		}

		if (complement) {
			subset.flip(0, numStates - 1);
		}

		for (int state = subset.nextSetBit(0); state >= 0; state = subset.nextSetBit(state + 1)) {
			// Initialise the result
			resultMin[state] = vectMin[state];
			resultMax[state] = vectMax[state];

			for (int trans : outReactions.get(state)) {
				int succ = toState(trans);
				resultMin[state] += mvMultMidSumEvalMin(trans, vectMin[succ], vectMin[state], q);
				resultMax[state] += mvMultMidSumEvalMax(trans, vectMax[succ], vectMax[state], q);
			}

			for (Pair<Integer, Integer> transs : inoutReactions.get(state)) {
				int trans = transs.first;
				int succTrans = transs.second;

				assert toState(trans) == state;
				int succ = toState(succTrans);

				if (!subset.get(fromState(trans))) {
					// Reduce to the case of an incoming reaction
					resultMin[state] += mvMultMidSumEvalMin(trans, vectMin[succ], vectMin[state], q);
					resultMax[state] += mvMultMidSumEvalMax(trans, vectMax[succ], vectMax[state], q);
					continue;
				}

				// The rate params of the two considered transitions must be identical
				assert rateParamsLowers[succTrans] == rateParamsLowers[trans] && rateParamsUppers[succTrans] == rateParamsUppers[trans];

				resultMin[state] += mvMultMidSumEvalMin(succTrans, vectMin[succ], vectMin[state], q);
				resultMax[state] += mvMultMidSumEvalMax(succTrans, vectMax[succ], vectMax[state], q);
			}
		}

		// Optimisation: Non-parametrised transitions
		for (int trans = 0; trans < numTotalTransitions; trans++) {
			if (parametrised(trans))
				continue;

			int state = fromState(trans);
			int succ = toState(trans);

			if (!subset.get(state))
				continue;

			double rate = rateParamsLowers[trans] * ratePopulations[trans];
			resultMin[state] += rate * (vectMin[succ] - vectMin[state]) / q;
			resultMax[state] += rate * (vectMax[succ] - vectMax[state]) / q;
		}
	}

	public void setRegion(BoxRegion region) throws PrismException
	{
		for (int trans = 0; trans < numTotalTransitions; trans++) {
			rateParamsLowers[trans] = rateParams[trans].evaluateDouble(region.getLowerBounds());
			rateParamsUppers[trans] = rateParams[trans].evaluateDouble(region.getUpperBounds());
			parametrisedTransitions[trans] = rateParamsLowers[trans] != rateParamsUppers[trans];

		}
	}
}
