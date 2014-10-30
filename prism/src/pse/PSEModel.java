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
import parser.ast.ModulesFile;
import prism.ModelType;
import prism.Pair;
import prism.PrismException;
import prism.PrismLog;
import explicit.CTMC;
import explicit.ModelExplicit;

/**
 * Represents a parametrised CTMC model to be used for PSE-based
 * techniques of analysis.
 * 
 * @see PSEModelBuilder
 */
public final class PSEModel extends ModelExplicit
{
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
	/** total sum of leaving rates for a state */
	private double[] exitRates;
	/** set of hash codes for deciding whether state has predecessors via reaction */
	private Set<Integer> predecessorsViaReaction;
	/** map from state to transitions coming into the state (and not outgoing) */
	private Map<Integer, List<Integer>> trsI;
	private int trsICnt;
	/** map from state to transitions both incoming in and outgoing from the state */
	private Map<Integer, List<Pair<Integer, Integer>>> trsIO;
	private int trsIOCnt;
	/** map from state to transitions going out from the state (and not incoming) */
	private Map<Integer, List<Integer>> trsO;
	private int trsOCnt;
	/** map from state to transitions that are nopt parametrised */
	private Map<Integer, List<Integer>> trsNP;

	private PSEModelForVM modelVM;

	/**
	 * Constructs a new parametric model.
	 */
	PSEModel()
	{
		numStates = 0;
		numTransitions = 0;
		trsICnt = 0;
		trsOCnt = 0;
		trsIOCnt = 0;
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
	void setSumLeaving(double leaving)
	{
		exitRates[numStates] = leaving;
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
			subset.set(0, numStates - 1);
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
		trsNP = new HashMap<Integer, List<Integer>>(numStates);
		for (int state = 0; state < numStates; state++) {
			trsI.put(state, new LinkedList<Integer>());
			trsO.put(state, new LinkedList<Integer>());
			trsIO.put(state, new LinkedList<Pair<Integer, Integer>>());
			trsNP.put(state, new LinkedList<Integer>());
		}

		// Populate the sets with transition indices
		for (int pred = 0; pred < numStates; pred++) {
			for (int predTrans = stateBegin(pred); predTrans < stateEnd(pred); predTrans++) {
				if (!isParametrised(predTrans)) {
					trsNP.get(pred).add(predTrans);
					continue;
				}
				boolean inout = false;
				int predReaction = getReaction(predTrans);
				int state = toState(predTrans);
				for (int trans = stateBegin(state); trans < stateEnd(state); trans++) {
					if (getReaction(trans) == predReaction) {
						inout = true;
						trsIO.get(state).add(new Pair<Integer, Integer>(predTrans, trans));
						++trsIOCnt;
						break;
					}
				}
				if (!inout) {
					trsI.get(state).add(predTrans);
					++trsICnt;
				}
				if (!predecessorsViaReaction.contains(pred ^ predReaction)) {
					trsO.get(pred).add(predTrans);
					++trsOCnt;
				}
			}
		}

		modelVM = buildModelForVM();
	}

  public PSEModelForVM buildModelForVM()
	{
		final double qrec = 1.0 / getDefaultUniformisationRate();
		final double[] trRatePopul_ = new double[trRatePopul.length];
		for (int i = 0; i < trRatePopul_.length; ++i) { trRatePopul_[i] = trRatePopul[i] * qrec; }

		int[] trsI_ = new int[trsICnt];
		int[] trsO_ = new int[trsOCnt];
		int[] trsIO_ = new int[trsIOCnt * 2];

		VectorOfDouble trsNPVal = new VectorOfDouble();
		VectorOfInt trsNPTrg = new VectorOfInt();
		int[] trsNPSrcBeg = new int[numStates + 1];
		int trsNPPos = 0;

		int trsIPos = 0;
		int trsOPos = 0;
		int trsIOPos = 0;
		for (int state = 0; state < numStates; ++state)
		{
			List<Integer> stTrsI = trsI.get(state);
			List<Integer> stTrsO = trsO.get(state);
			List<Pair<Integer, Integer>> stTrsIO = trsIO.get(state);
			List<Integer> stTrsNP = trsNP.get(state);

			for (Integer tr : stTrsI)
			{
				trsI_[trsIPos++] = tr;
			}
			for (Integer tr : stTrsO)
			{
				trsO_[trsOPos++] = tr;
			}
			for (Pair<Integer, Integer> p : stTrsIO)
			{
				trsIO_[trsIOPos++] = p.first;
				trsIO_[trsIOPos++] = p.second;
			}

			trsNPSrcBeg[state] = trsNPPos;
			for (Integer t : stTrsNP)
			{
				final double rate = trRateLower[t] * trRatePopul_[t];
				if (rate != 0)
				{
					trsNPVal.pushBack(rate);
					trsNPTrg.pushBack(trStTrg[t]);
					++trsNPPos;
				}
			}
		}
		trsNPSrcBeg[numStates] = trsNPPos;

		// Sort for better locality when accessing trRateLower/Upper/Popul.
		Arrays.sort(trsI_);
		Arrays.sort(trsO_);

		return new PSEModelForVM
			( numStates, numTransitions
			, trRateLower
			, trRateUpper
			, trRatePopul_
			, trStSrc
			, trStTrg
			, trsI_
			, trsO_
			, trsIO_
			, trsNPVal.data()
			, trsNPTrg.data()
			, trsNPSrcBeg
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
	 * @see #mvMult(double[], double[], double[], double[], BitSet, boolean, double)
	 */
	public void vmMult(double vectMin[], double resultMin[], double vectMax[], double resultMax[])
			throws PrismException
	{
		modelVM.vmMult(vectMin, resultMin, vectMax, resultMax);
	}

	private double mvMultMidSumEvalMin(int trans, double vectMinPred, double vectMinState, double q)
	{
		double midSumNumeratorMin = trRatePopul[trans] * vectMinPred - trRatePopul[trans] * vectMinState;
		if (midSumNumeratorMin > 0.0) {
			return trRateLower[trans] * midSumNumeratorMin / q;
		} else {
			return trRateUpper[trans] * midSumNumeratorMin / q;
		}
	}

	private double mvMultMidSumEvalMax(int trans, double vectMaxPred, double vectMaxState, double q)
	{
		double midSumNumeratorMax = trRatePopul[trans] * vectMaxPred - trRatePopul[trans] * vectMaxState;
		if (midSumNumeratorMax > 0.0) {
			return trRateUpper[trans] * midSumNumeratorMax / q;
		} else {
			return trRateLower[trans] * midSumNumeratorMax / q;
		}
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
	 * @param subset only do multiplication for these rows (null means "all")
	 * @param complement if true, {@code subset} is taken to be its complement
	 * @param q uniformisation rate
	 * @see #vmMult(double[], double[], double[], double[])
	 */
	public void mvMult(double vectMin[], double resultMin[], double vectMax[], double resultMax[], BitSet subset, boolean complement, double q)
			throws PrismException
	{
		if (subset == null) {
			// Will loop over all states
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

			for (int trans : trsO.get(state)) {
				int succ = toState(trans);
				resultMin[state] += mvMultMidSumEvalMin(trans, vectMin[succ], vectMin[state], q);
				resultMax[state] += mvMultMidSumEvalMax(trans, vectMax[succ], vectMax[state], q);
			}

			for (Pair<Integer, Integer> transs : trsIO.get(state)) {
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
				assert trRateLower[succTrans] == trRateLower[trans];
				assert trRateUpper[succTrans] == trRateUpper[trans];

				resultMin[state] += mvMultMidSumEvalMin(succTrans, vectMin[succ], vectMin[state], q);
				resultMax[state] += mvMultMidSumEvalMax(succTrans, vectMax[succ], vectMax[state], q);
			}
		}

		// Optimisation: Non-parametrised transitions
		for (int trans = 0; trans < numTransitions; trans++) {
			if (isParametrised(trans))
				continue;

			int state = fromState(trans);
			int succ = toState(trans);

			if (!subset.get(state))
				continue;

			double rate = trRateLower[trans] * trRatePopul[trans];
			resultMin[state] += rate * (vectMin[succ] - vectMin[state]) / q;
			resultMax[state] += rate * (vectMax[succ] - vectMax[state]) / q;
		}
	}

	/**
	 * Updates the transition rates of this parametrised CTMC according
	 * to the given parameter region.
	 * 
	 * @param region parameter region according to which configure the model's
	 * parameter space
	 * @throws PrismException thrown if rates cannot be evaluated with the new
	 * parameter region's bounds
	 */
	public void configureParameterSpace(BoxRegion region) throws PrismException
	{
		for (int trans = 0; trans < numTransitions; trans++) {
			trRateLower[trans] = rateParams[trans].evaluateDouble(region.getLowerBounds());
			trRateUpper[trans] = rateParams[trans].evaluateDouble(region.getUpperBounds());
			parametrisedTransitions[trans] = trRateLower[trans] != trRateUpper[trans];
		}
		if (modelVM != null)
		{
			modelVM = buildModelForVM();
		}
	}

	/**
	 * Returns a particular non-parametrised CTMC associated with
	 * the given point of the parameter space of this parametrised CTMC.
	 * 
	 * @param point point of parameter space determining the parameters' values
	 * @param modulesFile modules file
	 * @param constructModel object conducting construction of {@code explicit.CTMC}
	 * models
	 * @return non-parametrised CTMC obtained by substituting {@code point}
	 * for parameter ranges
	 * @throws PrismException thrown if an error occurred during construction
	 * of the non-parametrised CTMC
	 */
	public CTMC instantiate(Point point, ModulesFile modulesFile, explicit.ConstructModel constructModel)
			throws PrismException
	{
		modulesFile = (ModulesFile) modulesFile.deepCopy();
		// Add point dimensions to constants of the modules file
		modulesFile.getConstantValues().addValues(point.getDimensions());
		return (CTMC) constructModel.constructModel(modulesFile);
	}
}
