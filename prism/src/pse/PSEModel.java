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
import java.util.*;

import parser.Values;
import parser.ast.Expression;
import prism.ModelType;
import prism.Pair;
import prism.PrismException;
import prism.PrismLog;
import explicit.ModelExplicit;

/**
 * Represents a parametrised CTMC model to be used for PSE-based
 * techniques of analysis.
 * 
 * @see PSEModelBuilder
 */
public final class PSEModel extends ModelExplicit
{
	/** complete parameter space */
	private BoxRegion completeSpace;
	/** total number of probabilistic transitions over all states */
	private int numTransitions;
	/** begin and end of state transitions */
	private int[] rows;
	/** origins of distribution branches */
	private int[] trStSrc;
	/** targets of distribution branches */
	private int[] trStTrg;
	/** all transitions' rate parameters, as expressions */
	private Expression[] rateParams;
	/** all transitions' rate parameters, evaluated with lower bounds of current region */
	private double[] trRateLower;
	/** all transitions' rate parameters, evaluated with upper bounds of current region */
	private double[] trRateUpper;
	/** species populations in all transitions' origin states */
	private double[] trRatePopul;
	/** indication on all transitions, whether their rate depends on parameters */
	private boolean[] parametrisedTransitions;
	/** all transitions' reactions, i.e. transition kinds */
	private int[] reactions;
	/** labels - per transition, <i>not</i> per action */
	private String[] labels;
	/** total sum of leaving rates for a state, as expressions */
	private Expression[] exitRatesExpr;
	/** total sum of leaving rates for a state, evaluated */
	private double[] exitRates;
	/** set of hash codes for deciding whether state has predecessors via reaction */
	private Set<Integer> predecessorsViaReaction;
	/** map from state to transitions coming into the state (and not outgoing) */
	private Map<Integer, List<Integer>> trsI;
	/** map from state to transitions both incoming in and outgoing from the state */
	private Map<Integer, List<Pair<Integer, Integer>>> trsIO;
	/** map from state to transitions going out from the state (and not incoming) */
	private Map<Integer, List<Integer>> trsO;
	/** map from state to transitions that are nopt parametrised */
	private Map<Integer, List<Integer>> trsNPBySrc;
	private Map<Integer, List<Integer>> trsNPByTrg;

	private PSEModelForVM modelVM;
	private PSEModelForMV modelMV;

	/**
	 * Constructs a new parametric model.
	 */
	PSEModel()
	{
		numStates = 0;
		numTransitions = 0;
		initialStates = new LinkedList<Integer>();
		deadlocks = new TreeSet<Integer>();
		predecessorsViaReaction = new HashSet<Integer>();
	}

	// Accessors (for Model)

	@Override
	public ModelType getModelType()
	{
		return ModelType.CTMC;
	}

	public BoxRegion getCompleteSpace()
	{
		return completeSpace;
	}

	@Override
	public Values getConstantValues()
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public int getNumTransitions()
	{
		return numTransitions;
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
	 * @param numTransitions total number of probabilistic transitions of the model
	 */
	void reserveMem(int numStates, int numTransitions)
	{
		rows = new int[numStates + 1];
		labels = new String[numTransitions];
		reactions = new int[numTransitions];
		rateParams = new Expression[numTransitions];
		trRateLower = new double[numTransitions];
		trRateUpper = new double[numTransitions];
		parametrisedTransitions = new boolean[numTransitions];
		trRatePopul = new double[numTransitions];
		trStTrg = new int[numTransitions];
		trStSrc = new int[numTransitions];
		exitRatesExpr = new Expression[numStates];
		exitRates = new double[numStates];
	}

	/**
	 * Finishes the current state.
	 * Starting with the 0th state, this function shall be called once all
	 * transitions outgoing from the current nth state have been added.
	 * Subsequent method calls of {@code addTransition}
	 * will then apply to the (n+1)th state. Notice that this method must be
	 * called for each state of the method, even the last one, once all its
	 * transitions have been added.
	 */
	void finishState()
	{
		rows[numStates + 1] = numTransitions;
		numStates++;
	}

	/**
	 * Adds a probabilistic transition from the current state.
	 * 
	 * @param reaction kind of the transition being added
	 * @param fromState from which state the transition goes
	 * @param toState to which state the transition leads
	 * @param rateParamsExpr the transition's rate parameters as expression
	 * @param ratePopulation the transition's origin state's species population
	 * @param action action with which the transition is labelled
	 */
	void addTransition(int reaction, int fromState, int toState, Expression rateParamsExpr, double ratePopulation, String action)
	{
		reactions[numTransitions] = reaction;
		trStSrc[numTransitions] = fromState;
		trStTrg[numTransitions] = toState;
		rateParams[numTransitions] = rateParamsExpr;
		trRatePopul[numTransitions] = ratePopulation;
		labels[numTransitions] = action;

		predecessorsViaReaction.add(toState ^ reaction);

		numTransitions++;
	}

	/**
	 * Sets the total sum of leaving rates from the current state.
	 * 
	 * @param leaving sum of leaving rates from the current state
	 */
	void setSumLeaving(Expression leaving)
	{
		exitRatesExpr[numStates] = leaving;
	}

	/**
	 * Returns the number of the first transition going from {@code state}.
	 * 
	 * 
	 * @param state state to return number of first transition of
	 * @return number of first transition going from {@code state}
	 */
	int stateBegin(int state)
	{
		return rows[state];
	}

	/**
	 * Returns the number of the last transition going from {@code state} plus one.
	 * 
	 * @param state state to return number of last transition of
	 * @return number of last transition going from {@code state} plus one
	 */
	int stateEnd(int state)
	{
		return rows[state + 1];
	}

	/**
	 * Returns whether the rate of the given transition depends on parameters.
	 * 
	 * @param trans transition about which to decide whether it's parametrised
	 * @return true iff the rate {@code trans} depends on parameters
	 */
	boolean isParametrised(int trans)
	{
		return parametrisedTransitions[trans];
	}

	/**
	 * Returns reaction to which the given transition belongs.
	 * 
	 * @param trans transition to return reaction for
	 * @return reaction of {@code trans}
	 */
	int getReaction(int trans)
	{
		return reactions[trans];
	}

	/**
	 * Returns the predecessor state of the given transition.
	 * 
	 * @param trans transition to return predecessor for
	 * @return predecessor state of {@code trans}
	 */
	int fromState(int trans)
	{
		return trStSrc[trans];
	}

	/**
	 * Returns the successor state of the given transition.
	 * 
	 * @param trans transition to return successor for
	 * @return successor state of {@code trans}
	 */
	int toState(int trans)
	{
		return trStTrg[trans];
	}

	/**
	 * Returns the label of the given transition.
	 * 
	 * @param trans transition to return label of
	 * @return label of {@code trans}
	 */
	String getLabel(int trans)
	{
		return labels[trans];
	}

	/**
	 * Computes the maximum exit rate over all states in the model,
	 * i.e. max_i { sum_j R(i,j) }.
	 * 
	 * @return maximum exit rate
	 */
	double getMaxExitRate()
	{
		return getMaxExitRate(null);
	}

	/**
	 * Computes the maximum exit rate over states in {@code subset},
	 * i.e. max_{i in subset} { sum_j R(i,j) }.
	 * 
	 * @param subset subset of states over which to compute maximum exit rate
	 * @return maximum exit rate over states in {@code subset}
	 */
	double getMaxExitRate(BitSet subset)
	{
		if (subset == null) {
			// Will loop over all states
			subset = new BitSet(numStates);
			subset.set(0, numStates);
		}
		double max = Double.NEGATIVE_INFINITY;
		for (int state = subset.nextSetBit(0); state >= 0; state = subset.nextSetBit(state + 1)) {
			if (exitRates[state] > max) {
				max = exitRates[state];
			}
		}
		return max;
	}

	/**
	 * Computes the default rate used to uniformise this parametrised CTMC.
	 */
	double getDefaultUniformisationRate()
	{
		return 1.02 * getMaxExitRate();
	}

	/**
	 * Computes the default rate used to uniformise this parametrised CTMC,
	 * assuming that all states *not* in {@code nonAbs} have been made absorbing.
	 */
	double getDefaultUniformisationRate(BitSet nonAbs)
	{
		return 1.02 * getMaxExitRate(nonAbs);
	}

	/**
	 * Analyses the model's transitions in order to divide them between exclusively
	 * incoming, exclusively outgoing or both incoming/outgoing from the perspective
	 * of particular states. The results are stored in {@code trsI},
	 * {@code trsO} and {@code trsIO}, respectively.
	 */
	public void computeInOutTransitions()
	{
		if (trsI != null)
			return;

		// Initialise the transition sets
		trsI = new HashMap<Integer, List<Integer>>(numStates);
		trsO = new HashMap<Integer, List<Integer>>(numStates);
		trsIO = new HashMap<Integer, List<Pair<Integer, Integer>>>(numStates);
		trsNPBySrc = new HashMap<Integer, List<Integer>>(numStates);
		trsNPByTrg = new HashMap<Integer, List<Integer>>(numStates);
		for (int state = 0; state < numStates; state++) {
			trsI.put(state, new LinkedList<Integer>());
			trsO.put(state, new LinkedList<Integer>());
			trsIO.put(state, new LinkedList<Pair<Integer, Integer>>());
			trsNPBySrc.put(state, new LinkedList<Integer>());
			trsNPByTrg.put(state, new LinkedList<Integer>());
		}

		// Populate the sets with transition indices
		for (int pred = 0; pred < numStates; pred++) {
			for (int predTrans = stateBegin(pred); predTrans < stateEnd(pred); predTrans++) {
				boolean inout = false;
				int predReaction = getReaction(predTrans);
				int state = toState(predTrans);
				if (!isParametrised(predTrans)) {
					trsNPBySrc.get(pred).add(predTrans);
					trsNPByTrg.get(state).add(predTrans);
					continue;
				}
				for (int trans = stateBegin(state); trans < stateEnd(state); trans++) {
					if (getReaction(trans) == predReaction) {
						inout = true;
						trsIO.get(state).add(new Pair<Integer, Integer>(predTrans, trans));
						break;
					}
				}
				if (!inout) {
					trsI.get(state).add(predTrans);
				}
				if (!predecessorsViaReaction.contains(pred ^ predReaction)) {
					trsO.get(pred).add(predTrans);
				}
			}
		}
	}

    final public void prepareForVM()
	{
	    final double qrec = 1.0 / getDefaultUniformisationRate();

	    VectorOfDouble matMinVal = new VectorOfDouble();
	    VectorOfInt matMinSrc = new VectorOfInt();
	    int[] matMinTrgBeg = new int [numStates + 1];
	    int matMinPos = 0;

	    VectorOfDouble matMaxVal = new VectorOfDouble();
	    VectorOfInt matMaxSrc = new VectorOfInt();
	    int[] matMaxTrgBeg = new int [numStates + 1];
	    int matMaxPos = 0;

	    VectorOfDouble matVal = new VectorOfDouble();
	    VectorOfInt matSrc = new VectorOfInt();
	    int[] matTrgBeg = new int [numStates + 1];
	    int matPos = 0;

	    double[] matMinDiagVal = new double[numStates];
	    double[] matMaxDiagVal = new double[numStates];
	    for (int i = 0; i < numStates; ++i)
	    {
		    matMinDiagVal[i] = 1;
		    matMaxDiagVal[i] = 1;
	    }


	    VectorOfDouble matIOLowerVal0 = new VectorOfDouble();
	    VectorOfDouble matIOLowerVal1 = new VectorOfDouble();
	    VectorOfDouble matIOUpperVal0 = new VectorOfDouble();
	    VectorOfDouble matIOUpperVal1 = new VectorOfDouble();
	    VectorOfInt matIOSrc = new VectorOfInt();
	    int[] matIOTrgBeg = new int [numStates + 1];
	    int matIOPos = 0;

	    for (int state = 0; state < numStates; ++state)
	    {
		    matMinTrgBeg[state] = matMinPos;
		    matMaxTrgBeg[state] = matMaxPos;
		    matTrgBeg[state] = matPos;
		    matIOTrgBeg[state] = matIOPos;

		    List<Integer> stTrsI = trsI.get(state);
		    List<Integer> stTrsO = trsO.get(state);
		    List<Pair<Integer, Integer>> stTrsIO = trsIO.get(state);
		    List<Integer> stTrsNP = trsNPByTrg.get(state);

		    for (Pair<Integer, Integer> p : stTrsIO)
		    {
			    final int t0 = p.first;
			    final int t1 = p.second;
			    final int v0 = trStSrc[t0];
			    final int v1 = trStTrg[t0]; // == trStSrc[t1]

			    final double valLower0 = trRateLower[t0] * trRatePopul[t0] * qrec;
			    final double valLower1 = trRateLower[t1] * trRatePopul[t1] * qrec;
			    final double valUpper0 = trRateUpper[t0] * trRatePopul[t0] * qrec;
			    final double valUpper1 = trRateUpper[t1] * trRatePopul[t1] * qrec;

			    // The rate params of t0 and t1 must be identical
			    // assert trRateLower[t0] == trRateLower[t1];
			    // assert trRateUpper[t0] == trRateUpper[t1];
			    //
			    // The lower rate == 0 iff upper rate == 0
			    // assert (trRateLower[t0] == 0 && trRateUpper[t0] == 0) ||
			    //        (trRateLower[t0] != 0 && trRateUpper[t0] != 0)
			    //

			    // if (valLower0 != 0) should be enough -- see above
			    if (!(valLower0 == 0 && valLower1 == 0 && valUpper0 == 0 && valUpper1 == 0))
			    {
				    matIOLowerVal0.pushBack(valLower0);
				    matIOLowerVal1.pushBack(valLower1);
				    matIOUpperVal0.pushBack(valUpper0);
				    matIOUpperVal1.pushBack(valUpper1);

				    matIOSrc.pushBack(v0);
				    ++matIOPos;
			    }
		    }

		    for (Integer t : stTrsI)
		    {
			    final double valMin = trRateLower[t] * trRatePopul[t] * qrec;
			    final double valMax = trRateUpper[t] * trRatePopul[t] * qrec;
			    if (valMin != 0)
			    {
				    matMinVal.pushBack(valMin);
				    matMinSrc.pushBack(trStSrc[t]);
				    ++matMinPos;
			    }
			    if (valMax != 0)
			    {
				    matMaxVal.pushBack(valMax);
				    matMaxSrc.pushBack(trStSrc[t]);
				    ++matMaxPos;
			    }
		    }
		    for (Integer t : stTrsO)
		    {
			    matMinDiagVal[trStSrc[t]] -= trRateUpper[t] * trRatePopul[t] * qrec;
			    matMaxDiagVal[trStSrc[t]] -= trRateLower[t] * trRatePopul[t] * qrec;
		    }

		    for (Integer t : stTrsNP)
		    {
			    final double val = trRateLower[t] * trRatePopul[t] * qrec;
			    matMinDiagVal[trStSrc[t]] -= val;
			    matMaxDiagVal[trStSrc[t]] -= val;
			    if (val != 0)
			    {
				    matVal.pushBack(val);
				    matSrc.pushBack(trStSrc[t]);
				    ++matPos;
			    }
		    }
	    }
	    matMinTrgBeg[numStates] = matMinPos;
	    matMaxTrgBeg[numStates] = matMaxPos;
	    matTrgBeg[numStates] = matPos;
	    matIOTrgBeg[numStates] = matIOPos;

	    modelVM = new PSEModelForVM
        ( numStates, numTransitions
        , matIOLowerVal0.data()
        , matIOLowerVal1.data()
        , matIOUpperVal0.data()
        , matIOUpperVal1.data()
        , matIOSrc.data()
        , matIOTrgBeg

        , matMinVal.data()
        , matMinSrc.data()
        , matMinTrgBeg

        , matMaxVal.data()
        , matMaxSrc.data()
        , matMaxTrgBeg

        , matMinDiagVal
        , matMaxDiagVal
        , matVal.data()
        , matSrc.data()
        , matTrgBeg
        );
	}

	/**
	 * @param subset only do multiplication for these rows (null means "all")
	 * @param complement if true, {@code subset} is taken to be its complement
	 */
	final public void prepareForMV(BitSet subset, boolean complement)
	{
		final double qrec = 1.0 / getDefaultUniformisationRate(subset);

		if (subset == null) {
			// Loop over all states
			subset = new BitSet(numStates);
			subset.set(0, numStates - 1);
		} else {
			subset = (BitSet)subset.clone();
		}

		if (complement) {
			subset.flip(0, numStates - 1);
		}

		VectorOfDouble matNPVal = new VectorOfDouble();
		VectorOfInt matNPCol = new VectorOfInt();
		int matNPRowCnt = subset.cardinality();
		int[] matNPRow = new int [matNPRowCnt];
		int[] matNPRowBeg = new int [matNPRowCnt + 1];
		int matNPPos = 0;

		int rowCntAll = 0;
		int matRowCnt = 0;
		int matValCnt = 0;
		ArrayList<TreeMap<Integer, Pair<Double,Double>>> matExplicit =
				new ArrayList<TreeMap<Integer, Pair<Double,Double>>>(numStates); // Array of rows.
		for (int i = 0; i < numStates; ++i) matExplicit.add(null);
		for (int state = subset.nextSetBit(0); state >= 0; state = subset.nextSetBit(state + 1))
		{
			TreeMap<Integer, Pair<Double,Double>> matRow = new TreeMap<Integer, Pair<Double,Double>>();
			matExplicit.set(state, matRow);
			matNPRow[rowCntAll] = state;
			matNPRowBeg[rowCntAll] = matNPPos;

			List<Integer> stTrsO = trsO.get(state);
			List<Pair<Integer, Integer>> stTrsIO = trsIO.get(state);
			List<Integer> stTrsNP = trsNPBySrc.get(state);

			for (int t : stTrsO)
			{
				final double valLower = trRateLower[t] * trRatePopul[t] * qrec;
				final double valUpper = trRateUpper[t] * trRatePopul[t] * qrec;
				final int col = trStTrg[t];
				if (!(valLower == 0 && valUpper == 0))
				{
					Pair<Double, Double> prev = matRow.get(col);
					double prevLowerVal = 0;
					if (prev != null) prevLowerVal = prev.first;
					double prevUpperVal = 0;
					if (prev != null) prevUpperVal = prev.second;
					matRow.put(col, new Pair<Double,Double>(prevLowerVal + valLower, prevUpperVal + valUpper));
				}
			}

			for (Pair<Integer, Integer> p : stTrsIO)
			{
				final int t0 = p.first;
				final int t1 = p.second;
				final int v0 = trStSrc[t0];
				final int v1 = trStTrg[t0]; // == trStSrc[t1] == state
				final int v2 = trStTrg[t1];

				double valLower = 0;
				double valUpper = 0;
				if (!subset.get(v0))
				{
					valLower = trRateLower[t0] * trRatePopul[t0] * qrec;
					valUpper = trRateUpper[t0] * trRatePopul[t0] * qrec;
				}
				else
				{
					valLower = trRateLower[t1] * trRatePopul[t1] * qrec;
					valUpper = trRateUpper[t1] * trRatePopul[t1] * qrec;
				}

				final int col = v2;
				if (!(valLower == 0 && valUpper == 0))
				{
					Pair<Double, Double> prev = matRow.get(col);
					double prevLowerVal = 0;
					if (prev != null) prevLowerVal = prev.first;
					double prevUpperVal = 0;
					if (prev != null) prevUpperVal = prev.second;
					matRow.put(col, new Pair<Double,Double>(prevLowerVal + valLower, prevUpperVal + valUpper));
				}
			}
			if (!matRow.isEmpty())
			{
				++matRowCnt;
				matValCnt += matRow.size();
			}

			int ps = matNPVal.size();
			for (int t : stTrsNP)
			{
				int v0 = trStSrc[t]; // == state
				int v1 = trStTrg[t];
				final double val = trRateLower[t] * trRatePopul[t] * qrec;

				if (val != 0)
				{
					matNPVal.pushBack(val);
					matNPCol.pushBack(v1);
					++matNPPos;
				}
			}
			++rowCntAll;
		}
		matNPRowBeg[rowCntAll] = matNPPos;

		double[] matLowerVal = new double[matValCnt];
		double[] matUpperVal = new double[matValCnt];
		int[] matCol = new int[matValCnt];
		int[] matRow = new int[matRowCnt];
		int[] matRowBeg = new int[matRowCnt + 1];
		int matPos = 0;

		matRowCnt = 0;
		for (int row = 0; row < numStates; ++row)
		{
			TreeMap<Integer, Pair<Double,Double>> matExplicitRow = matExplicit.get(row);
			if (matExplicitRow == null || matExplicitRow.isEmpty())
			{
				continue;
			}

			matRow[matRowCnt] = row;
			matRowBeg[matRowCnt] = matPos;

			for (Map.Entry<Integer, Pair<Double,Double>> e : matExplicitRow.entrySet())
			{
				matCol[matPos] = e.getKey();
				matLowerVal[matPos] = e.getValue().first;
				matUpperVal[matPos] = e.getValue().second;
				++matPos;
			}
			++matRowCnt;
		}
		matRowBeg[matRowCnt] = matPos;

		modelMV = new PSEModelForMV
			( numStates
			, matLowerVal
			, matUpperVal
			, matCol
			, matRow
			, matRowBeg
			, matRowCnt

			, matNPVal.data()
			, matNPCol.data()
			, matNPRow
			, matNPRowBeg
			, matNPRowCnt
			);
	}

	/**
	 * Does a vector-matrix multiplication for this parametrised CTMC's transition
	 * probability matrix (uniformised with rate {@code q}) and the vector's min/max
	 * components ({@code vectMin} and {@code vectMax}, respectively) passed in.
	 * The code follows closely the algorithm described in the article:
	 * <p>
	 * L. Brim‚ M. Češka‚ S. Dražan and D. Šafránek: Exploring Parameter Space
	 * of Stochastic Biochemical Systems Using Quantitative Model Checking
	 * In Computer Aided Verification (CAV'13): 107−123, 2013.
	 * 
	 * @param vectMin vector to multiply by when computing minimised result
	 * @param resultMin vector to store minimised result in
	 * @param vectMax vector to multiply by when computing maximised result
	 * @param resultMax vector to store maximised result in
	 * @see #mvMult(double[], double[], double[], double[])
	 */
	public void vmMult(double vectMin[], double resultMin[], double vectMax[], double resultMax[])
			throws PrismException
	{
		modelVM.vmMult(vectMin, resultMin, vectMax, resultMax);
	}

	/**
	 * Does a matrix-vector multiplication for this parametrised CTMC's transition
	 * probability matrix (uniformised with rate {@code q}) and the vector's min/max
	 * components ({@code vectMin} and {@code vectMax}, respectively) passed in.
	 * <p>
	 * NB: Semantics of {@code mvMult} is <i>not</i> analogical to that of {@code vmMult},
	 * the difference is crucial:  {@code result[k]_i} in {@code vmMult} is simply
	 * the probability of being in state {@code k} after {@code i} iterations starting
	 * from the initial state.  On the other hand, {@code mvMult}'s {@code result[k]_i}
	 * denotes the probability that an absorbing state (i.e., a state not in {@code subset})
	 * is reached after {@code i} iterations starting from {@code k}.
	 *
	 * @param vectMin vector to multiply by when computing minimised result
	 * @param resultMin vector to store minimised result in
	 * @param vectMax vector to multiply by when computing maximised result
	 * @param resultMax vector to store maximised result in
	 * @see #vmMult(double[], double[], double[], double[])
	 */
	public void mvMult(double vectMin[], double resultMin[], double vectMax[], double resultMax[])
			throws PrismException
	{
		modelMV.mvMult(vectMin, resultMin, vectMax, resultMax);
	}

	/**
	 * Updates the transition rates and other parametrised data
	 * of this parametrised CTMC according to the given parameter region.
	 * 
	 * @param region parameter region according to which configure the model's
	 * parameter space
	 * @throws PrismException thrown if rates cannot be evaluated with the new
	 * parameter region's bounds
	 */
	public void evaluateParameters(BoxRegion region) throws PrismException
	{
		for (int trans = 0; trans < numTransitions; trans++) {
			trRateLower[trans] = rateParams[trans].evaluateDouble(region.getLowerBounds());
			trRateUpper[trans] = rateParams[trans].evaluateDouble(region.getUpperBounds());
			parametrisedTransitions[trans] = trRateLower[trans] != trRateUpper[trans];
		}
		for (int state = 0; state < numStates; state++) {
			exitRates[state] = exitRatesExpr[state].evaluateDouble(region.getUpperBounds());
		}
		modelVM = null; // This marks the model as dirty (i.e. it needs to be rebuilt)
		modelMV = null;
	}

	/**
	 */
	public void setParameterSpace(BoxRegion region) throws PrismException
	{
		completeSpace = region;
		evaluateParameters(region);
	}
}
