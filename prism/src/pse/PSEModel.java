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

import java.io.File;
import java.util.ArrayList;
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
import prism.ModelType;
import prism.Pair;
import prism.PrismException;
import prism.PrismLog;
import explicit.ModelExplicit;

final class PSEModel extends ModelExplicit
{
	/** total number of probabilistic transitions over all states */
	private int numTotalTransitions;
	/** begin and end of state transitions */
	private int[] rows;
	/** origins and targets of distribution branches */
	private int[] colsFrom;
	private int[] colsTo;
	/** */
	private double[] basicRateParamsLowers;
	private double[] basicRateParamsUppers;
	private double[] rateParamsLowers;
	private double[] rateParamsUppers;
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
			if (succState(trans) == s2) {
				return true;
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
		basicRateParamsLowers = new double[numTotalTransitions];
		basicRateParamsUppers = new double[numTotalTransitions];
		rateParamsLowers = new double[numTotalTransitions];
		rateParamsUppers = new double[numTotalTransitions];
		ratePopulations = new double[numTotalTransitions];
		parametrisedTransitions = new boolean[numTotalTransitions];
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
	void addTransition(int reaction, int fromState, int toState, double rateParamsLower, double rateParamsUpper, double ratePopulation, String action)
	{
		reactions[numTotalTransitions] = reaction;
		colsFrom[numTotalTransitions] = fromState;
		colsTo[numTotalTransitions] = toState;
		basicRateParamsLowers[numTotalTransitions] = rateParamsLower;
		basicRateParamsUppers[numTotalTransitions] = rateParamsUpper;
		rateParamsLowers[numTotalTransitions] = rateParamsLower;
		rateParamsUppers[numTotalTransitions] = rateParamsUpper;
		ratePopulations[numTotalTransitions] = ratePopulation;
		labels[numTotalTransitions] = action;
		parametrisedTransitions[numTotalTransitions] = rateParamsLower != rateParamsUpper;

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
	int getReaction(int succNr)
	{
		return reactions[succNr];
	}

	/**
	 * Returns the successor state of the given probabilistic branch.
	 * 
	 * @param succNr probabilistic branch to return successor state of
	 * @return state which probabilistic branch leads to
	 */
	int succState(int succNr)
	{
		return colsTo[succNr];
	}
	
	/**
	 */
	int currState(int succNr)
	{
		return colsFrom[succNr];
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
		/*
		assert inReactions == null;
		assert inoutReactions == null;
		assert outReactions == null;
		*/
		if (inReactions != null && inoutReactions != null && outReactions != null)
			return;

		// Initialise the reaction sets
		inReactions = new HashMap<Integer, List<Integer>>(numStates);
		inoutReactions = new HashMap<Integer, List<Pair<Integer, Integer>>>(numStates);
		outReactions = new HashMap<Integer, List<Integer>>(numStates);
		for (int state = 0; state < numStates; state++) {
			inReactions.put(state, new ArrayList<Integer>());
			inoutReactions.put(state, new ArrayList<Pair<Integer, Integer>>());
			outReactions.put(state, new ArrayList<Integer>());
		}

		// Populate the sets with transition indices
		for (int pred = 0; pred < numStates; pred++) {
			for (int predTrans = stateBegin(pred); predTrans < stateEnd(pred); predTrans++) {
				if (!parametrisedTransitions[predTrans])
					continue;

				boolean inout = false;
				int predReaction = getReaction(predTrans);
				int state = succState(predTrans);

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
		int pred, state;
		double midSumNumeratorMin, midSumNumeratorMax;
		
		for (state = 0; state < numStates; state++) {
			// Initialise the result
			resultMin[state] = vectMin[state];
			resultMax[state] = vectMax[state];

			// Incoming reactions
			for (int succ : inReactions.get(state)) {
				pred = currState(succ);
				resultMin[state] += rateParamsLowers[succ] * ratePopulations[succ] * vectMin[pred] / q;
				resultMax[state] += rateParamsUppers[succ] * ratePopulations[succ] * vectMax[pred] / q;
			}

			// Outgoing reactions
			for (int succ : outReactions.get(state)) {
				resultMin[state] -= rateParamsUppers[succ] * ratePopulations[succ] * vectMin[state] / q;
				resultMax[state] -= rateParamsLowers[succ] * ratePopulations[succ] * vectMax[state] / q;
			}

			// Both incoming and outgoing
			for (Pair<Integer, Integer> succs : inoutReactions.get(state)) {
				int predSucc = succs.first;
				int succ = succs.second;

				pred = currState(predSucc);
				assert currState(succ) == state;

				// The rate params assumed to be the same for both `pred` and `state`
				assert rateParamsLowers[predSucc] == rateParamsLowers[succ] && rateParamsUppers[predSucc] == rateParamsUppers[succ];

				midSumNumeratorMin = ratePopulations[predSucc] * vectMin[pred] - ratePopulations[succ] * vectMin[state];
				if (midSumNumeratorMin > 0) resultMin[state] += rateParamsLowers[succ] * midSumNumeratorMin / q;
				else resultMin[state] += rateParamsUppers[succ] * midSumNumeratorMin / q;

				midSumNumeratorMax = ratePopulations[predSucc] * vectMax[pred] - ratePopulations[succ] * vectMax[state];
				if (midSumNumeratorMax > 0) resultMax[state] += rateParamsUppers[succ] * midSumNumeratorMax / q;
				else resultMax[state] += rateParamsLowers[succ] * midSumNumeratorMax / q;
			}
		}

		// Optimisation: Non-parametrised transitions
		for (int succ = 0; succ < numTotalTransitions; succ++) {
			if (parametrisedTransitions[succ])
				continue;

			pred = currState(succ);
			state = succState(succ);

			double rate = rateParamsLowers[succ] * ratePopulations[succ];

			resultMin[pred] -= rate * vectMin[pred] / q;
			resultMax[pred] -= rate * vectMax[pred] / q;

			resultMin[state] += rate * vectMin[pred] / q;
			resultMax[state] += rate * vectMax[pred] / q;
		}
	}

	private double mvMultMidSumEvalMin(int succ, double vectMinPred, double vectMinState, double q)
	{
		double midSumNumeratorMin = ratePopulations[succ] * vectMinPred - ratePopulations[succ] * vectMinState;
		if (midSumNumeratorMin > 0) return rateParamsLowers[succ] * midSumNumeratorMin / q;
		else return rateParamsUppers[succ] * midSumNumeratorMin / q;
	}

	private double mvMultMidSumEvalMax(int succ, double vectMaxPred, double vectMaxState, double q)
	{
		double midSumNumeratorMin = ratePopulations[succ] * vectMaxPred - ratePopulations[succ] * vectMaxState;
		if (midSumNumeratorMin > 0) return rateParamsUppers[succ] * midSumNumeratorMin / q;
		else return rateParamsLowers[succ] * midSumNumeratorMin / q;
	}
	
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

		int pred, state;

		for (state = subset.nextSetBit(0); state >= 0; state = subset.nextSetBit(state + 1)) {
			// Initialise the result
			resultMin[state] = vectMin[state];
			resultMax[state] = vectMax[state];

			// Incoming reactions (NB propagating backwards!)
			for (int succ : outReactions.get(state)) {
				// Note the exchange: succState(succ) is stored as pred
				pred = succState(succ);
				resultMin[state] += mvMultMidSumEvalMin(succ, vectMin[pred], vectMin[state], q);
				resultMax[state] += mvMultMidSumEvalMax(succ, vectMax[pred], vectMax[state], q);
			}
			
			// Outgoing reactions (taken backwards) are not considered
			// when computing backwards transient probabilities.

			// Both incoming and outgoing
			for (Pair<Integer, Integer> succs : inoutReactions.get(state)) {
				// Note the exchange: succs.second is stored as predSucc
				int predSucc = succs.second;
				int succ = succs.first;

				// Note the exchange: succState(predSucc) is stored as pred
				pred = succState(predSucc);
				assert succState(succ) == state;

				if (!subset.get(currState(succ))) {
					// Reduce to the case of an incoming reaction
					resultMin[state] += mvMultMidSumEvalMin(succ, vectMin[pred], vectMin[state], q);
					resultMax[state] += mvMultMidSumEvalMax(succ, vectMax[pred], vectMax[state], q);
					continue;
				}

				// The rate params assumed to be the same for both `pred` and `state`
				assert rateParamsLowers[predSucc] == rateParamsLowers[succ] && rateParamsUppers[predSucc] == rateParamsUppers[succ];

				resultMin[state] += mvMultMidSumEvalMin(predSucc, vectMin[pred], vectMin[state], q);
				resultMax[state] += mvMultMidSumEvalMax(predSucc, vectMax[pred], vectMax[state], q);
			}
		}

		// Optimisation: Non-parametrised transitions
		for (int succ = 0; succ < numTotalTransitions; succ++) {
			if (parametrisedTransitions[succ])
				continue;

			// Note the exchange: succState(succ) is stored as pred
			pred = succState(succ);
			state = currState(succ);

			if (!subset.get(state))
				continue;

			double rate = rateParamsLowers[succ] * ratePopulations[succ];
			resultMin[state] += rate * (vectMin[pred] - vectMin[state]) / q;
			resultMax[state] += rate * (vectMax[pred] - vectMax[state]) / q;
		}
	}

	/**
	 * Assumption: Transition rates involve at most one parameter
	 * (cf. the CAV 2013 article, p. 4).
	 * 
	 * For instance, the method below would *not* properly scale
	 * the _parameter_ space if a rate formula included k1 * k2 where
	 * k1, k2 are parameters. If the parameters were each scaled
	 * by a factor F, the rate would in turn need to be scaled by F^2.
	 * The method below, however, directly scales the induced _rate_
	 * space -- always by F.
	 * 
	 * It works correctly with at most one parameter because in that case
	 * scaling of the rate space corresponds to scaling of the parameter
	 * space.
	 */
	public void scaleParameterSpace(double scaleLower, double scaleUpper)
	{
		for (int succ = 0; succ < numTotalTransitions; succ++) {
			rateParamsLowers[succ] = basicRateParamsLowers[succ] + scaleLower * (basicRateParamsUppers[succ] - basicRateParamsLowers[succ]);
			rateParamsUppers[succ] = basicRateParamsLowers[succ] + scaleUpper * (basicRateParamsUppers[succ] - basicRateParamsLowers[succ]);
		}
	}
}
