//==============================================================================
//	
//	Copyright (c) 2014-
//	Authors:
//	* Andrej Tokarcik <andrejtokarcik@gmail.com> (Masaryk University)
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
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

package pse;

import java.util.BitSet;
import java.util.LinkedList;
import java.util.Map.Entry;

import parser.Values;
import parser.ast.Expression;
import parser.ast.ExpressionFilter;
import parser.ast.ExpressionLabel;
import parser.ast.ExpressionProb;
import parser.ast.ExpressionTemporal;
import parser.ast.ExpressionUnaryOp;
import parser.ast.ModulesFile;
import parser.ast.PropertiesFile;
import parser.ast.RelOp;
import parser.ast.ExpressionFilter.FilterOperator;
import parser.type.TypeBool;
import parser.type.TypeDouble;
import prism.ModelType;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismLog;
import prism.PrismPrintStreamLog;
import prism.Result;
import explicit.FoxGlynn;
import explicit.StateModelChecker;
import explicit.StateValues;
import explicit.Utils;

public final class PSEModelChecker extends PrismComponent
{
	// Log for output (default to System.out)
	private PrismLog mainLog = new PrismPrintStreamLog(System.out);

	// Model file (for instantiating)
	private ModulesFile modulesFile = null;

	// Properties file (for labels, constants, etc.)
	private PropertiesFile propertiesFile = null;

	// Constants (extracted from model/properties)
	private Values constantValues;

	private StateModelChecker stateChecker;

	private BoxRegionFactory regionFactory;

	/**
	 * Constructor
	 */
	public PSEModelChecker(PrismComponent parent, BoxRegionFactory regionFactory) throws PrismException
	{
		super(parent);
		this.regionFactory = regionFactory;
		this.stateChecker = new StateModelChecker(this);
	}

	// Setters/getters

	/**
	 * Set the attached model file (for e.g. reward structures when model checking)
	 * and the attached properties file (for e.g. constants/labels when model checking)
	 */
	public void setModulesFileAndPropertiesFile(ModulesFile modulesFile, PropertiesFile propertiesFile)
	{
		this.modulesFile = modulesFile;
		this.propertiesFile = propertiesFile;
		// Get combined constant values from model/properties
		constantValues = new Values();
		constantValues.addValues(modulesFile.getConstantValues());
		if (propertiesFile != null)
			constantValues.addValues(propertiesFile.getConstantValues());

		stateChecker.setModulesFileAndPropertiesFile(modulesFile, propertiesFile);
	}

	public ModulesFile getModulesFile()
	{
		return modulesFile;
	}

	public Values getConstantValues()
	{
		return constantValues;
	}

	// Model checking functions

	/**
	 * Model check an expression, process and return the result.
	 * Information about states and model constants should be attached to the model.
	 */
	public Result check(PSEModel model, Expression expr, DecompositionProcedure decompositionProcedure) throws PrismException
	{
		long timer = 0;

		// Remove labels from property, using combined label list (on a copy of the expression)
		// This is done now so that we can handle labels nested below operators that are not
		// handled natively by the model checker yet (just evaluate()ed in a loop).
		expr = (Expression) expr.deepCopy().expandLabels(propertiesFile.getCombinedLabelList());

		// Also evaluate/replace any constants
		//expr = (Expression) expr.replaceConstants(constantValues);

		// Initialise decomposition procedure before model checking starts
		decompositionProcedure.initialise(this, model, expr);

		// Let decomposition procedure adjust the property, if needed
		expr = decompositionProcedure.getPropertyExpression();

		// Do model checking and store result vector
		timer = System.currentTimeMillis();
		BoxRegionValues regionValues = checkExpression(model, expr, decompositionProcedure);
		timer = System.currentTimeMillis() - timer;
		mainLog.println("\nTotal time for model checking: " + timer / 1000.0 + " seconds.");

		// Return result
		Result result = new Result();
		result.setResult(regionValues);
		decompositionProcedure.printSolution(mainLog);
		return result;
	}

	public BoxRegionValues checkExpression(PSEModel model, Expression expr, DecompositionProcedure decompositionProcedure) throws PrismException
	{
		if (expr instanceof ExpressionFilter) {
			return checkExpressionFilter(model, (ExpressionFilter) expr, decompositionProcedure);
		}
		if (expr instanceof ExpressionProb) {
			return checkExpressionProb(model, (ExpressionProb) expr, decompositionProcedure);
		}
		// Let explicit.StateModelChecker take care of other expressions
		StateValues vals = stateChecker.checkExpression(model, expr);
		return new BoxRegionValues(model, regionFactory.completeSpace(), vals, vals);
	}

	/**
	 * Model check a filter.
	 */
	protected BoxRegionValues checkExpressionFilter(PSEModel model, ExpressionFilter expr, DecompositionProcedure decompositionProcedure) throws PrismException
	{
		Object resObjMin, resObjMax;
		StateValues resRgnValsMin, resRgnValsMax;
		BitSet bsMin, bsMax;
		String resultExpl;
		BoxRegionValues resRgnVals = null;

		Expression filter = expr.getFilter();
		if (filter == null) {
			filter = Expression.True();
		}
		boolean filterTrue = Expression.isTrue(filter);
		String filterStatesString = filterTrue ? "all states" : "states satisfying filter";

		BoxRegionValues filterRgnVals = checkExpression(model, filter, decompositionProcedure);
		assert filterRgnVals.getNumRegions() == 1;
		BitSet bsFilterMin = filterRgnVals.getMin(regionFactory.completeSpace()).getBitSet();
		BitSet bsFilterMax = filterRgnVals.getMax(regionFactory.completeSpace()).getBitSet();
		boolean filterInit = (filter instanceof ExpressionLabel && ((ExpressionLabel) filter).getName().equals("init"));
		boolean filterInitSingle = filterInit & model.getNumInitialStates() == 1;
		if (bsFilterMin.isEmpty()) {
			throw new PrismException("Filter satisfies no states");
		}
		if (!filterInit) {
			mainLog.println("\nStates satisfying filter " + filter + ":");
			mainLog.println("min = " + bsFilterMin.cardinality());
			mainLog.println("max = " + bsFilterMax.cardinality());
		}

		BoxRegionValues subRgnVals = checkExpression(model, expr.getOperand(), decompositionProcedure);
		for (Entry<BoxRegion, BoxRegionValues.StateValuesPair> entry : subRgnVals) {
			BoxRegion region = entry.getKey();
			StateValues subRgnValsMin = entry.getValue().getMin();
			StateValues subRgnValsMax = entry.getValue().getMax();

			// Prepend all result strings with info about current region
			mainLog.println("\n== Region " + region + " ==");

			// Compute result according to filter type
			FilterOperator op = expr.getOperatorType();
			switch (op) {
			case PRINT:
			case PRINTALL:
				if (expr.getType() instanceof TypeBool) {
					mainLog.print("\nSatisfying states");
					mainLog.println(filterTrue ? ":" : " that are also in filter " + filter + ":");
					mainLog.println("\nmin\n");
					subRgnValsMin.printFiltered(mainLog, bsFilterMin);
					mainLog.println("\nmax\n");
					subRgnValsMax.printFiltered(mainLog, bsFilterMax);
				} else {
					if (op == FilterOperator.PRINT) {
						mainLog.println("\nResults (non-zero only) for filter " + filter + ":");
						mainLog.println("\nmin\n");
						subRgnValsMin.printFiltered(mainLog, bsFilterMin);
						mainLog.println("\nmax\n");
						subRgnValsMax.printFiltered(mainLog, bsFilterMax);
					} else {
						mainLog.println("\nResults (including zeros) for filter " + filter + ":");
						mainLog.println("\nmin\n");
						subRgnValsMin.printFiltered(mainLog, bsFilterMin, false, false, true, true);
						mainLog.println("\nmax\n");
						subRgnValsMax.printFiltered(mainLog, bsFilterMax, false, false, true, true);
					}
				}
				resRgnVals = subRgnVals;
				break;
			case MIN:
			case MAX:
			case ARGMIN:
			case ARGMAX:
				throw new UnsupportedOperationException();
			case COUNT:
				resObjMin = new Integer(subRgnValsMin.countOverBitSet(bsFilterMin));
				resObjMax = new Integer(subRgnValsMax.countOverBitSet(bsFilterMax));
				resRgnValsMin = new StateValues(expr.getType(), resObjMin, model);
				resRgnValsMax = new StateValues(expr.getType(), resObjMax, model);
				resRgnVals = new BoxRegionValues(model, region, resRgnValsMin, resRgnValsMax);

				resultExpl = filterTrue ? "Count of satisfying states" : "Count of satisfying states also in filter";
				mainLog.println("\n" + resultExpl + ":");
				mainLog.println("min = " + resObjMin);
				mainLog.println("max = " + resObjMax);
				break;
			case SUM:
			case AVG:
				throw new PrismException("Operation not implemented for parametric models");
			case FIRST:
				resObjMin = subRgnValsMin.firstFromBitSet(bsFilterMin);
				resObjMax = subRgnValsMax.firstFromBitSet(bsFilterMax);
				resRgnValsMin = new StateValues(expr.getType(), resObjMin, model);
				resRgnValsMax = new StateValues(expr.getType(), resObjMax, model);
				resRgnVals = new BoxRegionValues(model, region, resRgnValsMin, resRgnValsMax);

				resultExpl = "Value in ";
				if (filterInit) {
					resultExpl += filterInitSingle ? "the initial state" : "first initial state";
				} else {
					resultExpl += filterTrue ? "the first state" : "first state satisfying filter";
				}
				mainLog.println("\n" + resultExpl + ":");
				mainLog.println("min = " + resObjMin);
				mainLog.println("max = " + resObjMax);
				break;
			case RANGE:
				resObjMin = new prism.Interval(subRgnValsMin.minOverBitSet(bsFilterMin), subRgnValsMin.maxOverBitSet(bsFilterMin));
				resObjMax = new prism.Interval(subRgnValsMax.minOverBitSet(bsFilterMax), subRgnValsMax.maxOverBitSet(bsFilterMax));
				resRgnVals = subRgnVals;

				resultExpl = "Range of values over ";
				resultExpl += filterInit ? "initial states" : filterStatesString;
				mainLog.println("\n" + resultExpl + ":");
				mainLog.println("min = " + resObjMin);
				mainLog.println("max = " + resObjMax);
				break;
			case FORALL:
				bsMin = subRgnValsMin.getBitSet();
				bsMax = subRgnValsMax.getBitSet();

				mainLog.println("\nNumber of states satisfying " + expr.getOperand() + ":");
				mainLog.print("min = " + bsMin.cardinality());
				mainLog.println(bsMin.cardinality() == model.getNumStates() ? " (all in model)" : "");
				mainLog.print("max = " + bsMax.cardinality());
				mainLog.println(bsMax.cardinality() == model.getNumStates() ? " (all in model)" : "");

				resObjMin = subRgnValsMin.forallOverBitSet(bsFilterMin);
				resObjMax = subRgnValsMax.forallOverBitSet(bsFilterMax);
				resRgnValsMin = new StateValues(expr.getType(), resObjMin, model);
				resRgnValsMax = new StateValues(expr.getType(), resObjMax, model);
				resRgnVals = new BoxRegionValues(model, region, resRgnValsMin, resRgnValsMax);

				mainLog.print("\nProperty satisfied in {");
				mainLog.print("min = " + subRgnValsMin.countOverBitSet(bsFilterMin) + ", ");
				mainLog.print("max = " + subRgnValsMin.countOverBitSet(bsFilterMin) + "}");
				if (filterInit) {
					mainLog.println(" of " + model.getNumInitialStates() + " initial states.");
				} else {
					if (filterTrue) {
						mainLog.println(" of all " + model.getNumStates() + " states.");
					} else {
						mainLog.print(" of {");
						mainLog.print("min = " + bsFilterMin.cardinality() + ", ");
						mainLog.print("max = " + bsFilterMax.cardinality() + "}");
						mainLog.println(" filter states.");
					}
				}
				break;
			case EXISTS:
				bsMin = subRgnValsMin.getBitSet();
				bsMax = subRgnValsMax.getBitSet();
				resObjMin = subRgnValsMin.existsOverBitSet(bsFilterMin);
				resObjMax = subRgnValsMax.existsOverBitSet(bsFilterMax);
				resRgnValsMin = new StateValues(expr.getType(), resObjMin, model);
				resRgnValsMax = new StateValues(expr.getType(), resObjMax, model);
				resRgnVals = new BoxRegionValues(model, region, resRgnValsMin, resRgnValsMax);

				// Create explanation of result and print some details to log
				resultExpl = "min = Property satisfied in ";
				if (filterTrue) {
					resultExpl += ((Boolean) resObjMin) ? "at least one state" : "no states";
				} else {
					resultExpl += ((Boolean) resObjMin) ? "at least one filter state" : "no filter states";
				}
				resultExpl += "\nmax = Property satisfied in ";
				if (filterTrue) {
					resultExpl += ((Boolean) resObjMax) ? "at least one state" : "no states";
				} else {
					resultExpl += ((Boolean) resObjMax) ? "at least one filter state" : "no filter states";
				}
				mainLog.println("\n" + resultExpl);
				break;
			case STATE:
				// Check filter satisfied by exactly one state
				if (bsFilterMin.cardinality() != 1 || bsFilterMax.cardinality() != 1) {
					String s = "Filter should be satisfied in exactly 1 state";
					s += " (but \"" + filter + "\" is true in {";
					s += "min = " + bsFilterMin.cardinality() + ", ";
					s += "max = " + bsFilterMax.cardinality() + "} states)";
					throw new PrismException(s);
				}

				resObjMin = subRgnValsMin.firstFromBitSet(bsFilterMin);
				resObjMax = subRgnValsMax.firstFromBitSet(bsFilterMax);
				resRgnValsMin = new StateValues(expr.getType(), resObjMin, model);
				resRgnValsMax = new StateValues(expr.getType(), resObjMax, model);
				resRgnVals = new BoxRegionValues(model, region, resRgnValsMin, resRgnValsMax);

				// Create explanation of result and print some details to log
				resultExpl = "Value in ";
				if (filterInit) {
					resultExpl += "the initial state";
				} else {
					resultExpl += "the filter state";
				}
				mainLog.println("\n" + resultExpl + ":");
				mainLog.println("min = " + resObjMin);
				mainLog.println("max = " + resObjMax);
				break;
			default:
				throw new PrismException("Unrecognised filter type \"" + expr.getOperatorName() + "\"");
			}
		}

		return resRgnVals;
	}

	/**
	 * Model check a P operator expression and return the values for all states.
	 */
	protected BoxRegionValues checkExpressionProb(PSEModel model, ExpressionProb expr, DecompositionProcedure decompositionProcedure) throws PrismException
	{
		Expression pb; // Probability bound (expression)
		double p = 0; // Probability bound (actual value)
		RelOp relOp; // Relational operator
		ModelType modelType = model.getModelType();

		BoxRegionValues regionValues = null;

		// Get info from prob operator
		relOp = expr.getRelOp();
		pb = expr.getProb();
		if (pb != null) {
			p = pb.evaluateDouble(constantValues);
			if (p < 0 || p > 1)
				throw new PrismException("Invalid probability bound " + p + " in P operator");
		}

		// Compute probabilities
		switch (modelType) {
		case CTMC:
			regionValues = checkProbPathFormula(model, expr.getExpression(), decompositionProcedure);
			break;
		default:
			throw new PrismException("Cannot model check " + expr + " for a " + modelType);
		}

		// For =? properties, just return values
		if (pb == null) {
			return regionValues;
		}
		// Otherwise, compare against bound to get set of satisfying states
		else {
			for (Entry<BoxRegion, BoxRegionValues.StateValuesPair> entry : regionValues) {
				BitSet solnMin = entry.getValue().getMin().getBitSetFromInterval(relOp, p);
				BitSet solnMax = entry.getValue().getMax().getBitSetFromInterval(relOp, p);
				regionValues.put(entry.getKey(), solnMin, solnMax);
			}
			return regionValues;
		}
	}

	/**
	 * Compute probabilities for the contents of a P operator.
	 */
	protected BoxRegionValues checkProbPathFormula(PSEModel model, Expression expr, DecompositionProcedure decompositionProcedure) throws PrismException
	{
		// Test whether this is a simple path formula (i.e. PCTL)
		// and then pass control to appropriate method.
		if (expr.isSimplePathFormula()) {
			return checkProbPathFormulaSimple(model, expr, decompositionProcedure);
		} else {
			throw new PrismException("PSE supports unnested CSL formulae only");
		}
	}

	/**
	 * Compute probabilities for a simple, non-LTL path operator.
	 */
	protected BoxRegionValues checkProbPathFormulaSimple(PSEModel model, Expression expr, DecompositionProcedure decompositionProcedure) throws PrismException
	{
		BoxRegionValues regionValues = null;

		// Negation/parentheses
		if (expr instanceof ExpressionUnaryOp) {
			ExpressionUnaryOp exprUnary = (ExpressionUnaryOp) expr;
			// Parentheses
			if (exprUnary.getOperator() == ExpressionUnaryOp.PARENTH) {
				// Recurse
				regionValues = checkProbPathFormulaSimple(model, exprUnary.getOperand(), decompositionProcedure);
			}
			// Negation
			else if (exprUnary.getOperator() == ExpressionUnaryOp.NOT) {
				// Compute, then subtract from 1
				regionValues = checkProbPathFormulaSimple(model, exprUnary.getOperand(), decompositionProcedure);
				for (Entry<BoxRegion, BoxRegionValues.StateValuesPair> entry : regionValues) {
					entry.getValue().getMin().timesConstant(-1.0);
					entry.getValue().getMin().plusConstant(1.0);
					entry.getValue().getMax().timesConstant(-1.0);
					entry.getValue().getMax().plusConstant(1.0);
				}
			}
		}
		// Temporal operators
		else if (expr instanceof ExpressionTemporal) {
			ExpressionTemporal exprTemp = (ExpressionTemporal) expr;
			// Next
			if (exprTemp.getOperator() == ExpressionTemporal.P_X) {
				throw new PrismException("X operator requires the embedded-DTMC representation not supported by PSE");
			}
			// Until
			else if (exprTemp.getOperator() == ExpressionTemporal.P_U) {
				if (exprTemp.hasBounds()) {
					regionValues = checkProbBoundedUntil(model, exprTemp, decompositionProcedure);
				} else {
					throw new PrismException("PSE supports bounded until only");
				}
			}
			// Anything else - convert to until and recurse
			else {
				regionValues = checkProbPathFormulaSimple(model, exprTemp.convertToUntilForm(), decompositionProcedure);
			}
		}

		if (regionValues == null)
			throw new PrismException("Unrecognised path operator in P operator");

		return regionValues;
	}

	/**
	 * Model check a time-bounded until operator; return vector of probabilities for all states.
	 */
	protected BoxRegionValues checkProbBoundedUntil(PSEModel model, ExpressionTemporal expr, DecompositionProcedure decompositionProcedure)
			throws PrismException
	{
		double lTime, uTime; // time bounds
		Expression exprTmp;
		BitSet b1Min, b1Max, b2Min, b2Max, tmpMin, tmpMax;
		BoxRegionValues regionValues = null, tmpRegionValues = null;
		BoxRegionValues oldRegionValues = null, oldTmpRegionValues = null;
		StateValues probsMin, probsMax;

		// get info from bounded until

		// lower bound is 0 if not specified
		// (i.e. if until is of form U<=t)
		exprTmp = expr.getLowerBound();
		if (exprTmp != null) {
			lTime = exprTmp.evaluateDouble(constantValues);
			if (lTime < 0) {
				throw new PrismException("Invalid lower bound " + lTime + " in time-bounded until formula");
			}
		} else {
			lTime = 0;
		}
		// upper bound is -1 if not specified
		// (i.e. if until is of form U>=t)
		exprTmp = expr.getUpperBound();
		if (exprTmp != null) {
			uTime = exprTmp.evaluateDouble(constantValues);
			if (uTime < 0 || (uTime == 0 && expr.upperBoundIsStrict())) {
				String bound = (expr.upperBoundIsStrict() ? "<" : "<=") + uTime;
				throw new PrismException("Invalid upper bound " + bound + " in time-bounded until formula");
			}
			if (uTime < lTime) {
				throw new PrismException("Upper bound must exceed lower bound in time-bounded until formula");
			}
		} else {
			uTime = -1;
		}

		// model check operands first
		BoxRegionValues op1RgnVals = checkExpression(model, expr.getOperand1(), decompositionProcedure);
		BoxRegionValues op2RgnVals = checkExpression(model, expr.getOperand2(), decompositionProcedure);
		assert op1RgnVals.getNumRegions() == 1 && op2RgnVals.getNumRegions() == 1;
		b1Min = op1RgnVals.getMin(regionFactory.completeSpace()).getBitSet();
		b1Max = op1RgnVals.getMax(regionFactory.completeSpace()).getBitSet();
		b2Min = op2RgnVals.getMin(regionFactory.completeSpace()).getBitSet();
		b2Max = op2RgnVals.getMax(regionFactory.completeSpace()).getBitSet();

		// compute probabilities

		// a trivial case: "U<=0"
		if (lTime == 0 && uTime == 0) {
			// prob is 1 in b2 states, 0 otherwise
			probsMin = StateValues.createFromBitSetAsDoubles(b2Min, model);
			probsMax = StateValues.createFromBitSetAsDoubles(b2Max, model);
			regionValues = new BoxRegionValues(model, regionFactory.completeSpace(), probsMin, probsMax);
		} else {
			BoxRegionValues onesMultProbs = BoxRegionValues.createWithOnes(model, regionFactory.completeSpace());

			// >= lTime
			if (uTime == -1) {
				throw new PrismException("PSE supports bounded until only");
			}
			// <= uTime
			else if (lTime == 0) {
				// nb: uTime != 0 since would be caught above (trivial case)
				b1Min.andNot(b2Min);
				b1Max.andNot(b2Max);

				while (true) {
					try {
						regionValues = computeTransientBackwardsProbs(model, b2Min, b1Min, b2Max, b1Max, uTime, onesMultProbs, decompositionProcedure, oldRegionValues);
						break;
					} catch (DecompositionProcedure.DecompositionNeeded e) {
						e.getRegionsToDecompose().print(mainLog);
						for (BoxRegion region : e.getRegionsToDecompose())
							onesMultProbs.divideRegion(region);
						oldRegionValues = e.getExaminedRegionValues();
					}
				}
			}
			// [lTime,uTime] (including where lTime == uTime)
			else {
				tmpMin = (BitSet) b1Min.clone();
				tmpMin.andNot(b2Min);
				tmpMax = (BitSet) b1Max.clone();
				tmpMax.andNot(b2Max);

				while (true) {
					try {
						// Decomposition is always performed due to the failed condition in the second transient computation.
						// The first transient computation is executed with decomposing disabled.
						tmpRegionValues = computeTransientBackwardsProbs(model, b2Min, tmpMin, b2Max, tmpMax, uTime - lTime, onesMultProbs, SimpleDecompositionProcedure.NoDecomposing.getInstance(), oldTmpRegionValues);
						regionValues = computeTransientBackwardsProbs(model, b1Min, b1Min, b1Max, b1Max, lTime, tmpRegionValues, decompositionProcedure, oldRegionValues);
						break;
					} catch (DecompositionProcedure.DecompositionNeeded e) {
						e.getRegionsToDecompose().print(mainLog);
						for (BoxRegion region : e.getRegionsToDecompose())
							onesMultProbs.divideRegion(region);
						oldTmpRegionValues = tmpRegionValues;
						oldRegionValues = e.getExaminedRegionValues();
					}
				}
			}
		}

		mainLog.print("\nThe backwards transient probability computations produced ");
		mainLog.println(regionValues.getNumRegions() + " final regions.");
		return regionValues;
	}

	public BoxRegionValues computeTransientBackwardsProbs(PSEModel model,
			BitSet targetMin, BitSet nonAbsMin, BitSet targetMax, BitSet nonAbsMax,
			double t, BoxRegionValues multProbs, DecompositionProcedure decompositionProcedure)
					throws PrismException, DecompositionProcedure.DecompositionNeeded
	{
		return computeTransientBackwardsProbs(model, targetMin, nonAbsMin, targetMax, nonAbsMax, t, multProbs, decompositionProcedure, null);
	}

	/**
	 * NB: Decompositions of the parameter space must be performed explicitly,
	 * DecompositionNeeded is not handled within the method.
	 */
	public BoxRegionValues computeTransientBackwardsProbs(PSEModel model,
			BitSet targetMin, BitSet nonAbsMin, BitSet targetMax, BitSet nonAbsMax,
			double t, BoxRegionValues multProbs, DecompositionProcedure decompositionProcedure,
			BoxRegionValues previousResult)
					throws PrismException, DecompositionProcedure.DecompositionNeeded
	{
		BoxRegionValues regionValues = new BoxRegionValues(model);
		int i, n, iters, totalIters;
		double solnMin[], soln2Min[], sumMin[];
		double solnMax[], soln2Max[], sumMax[];
		double tmpsoln[];
		long timer;
		// Fox-Glynn stuff
		FoxGlynn fg;
		int left, right;
		double termCritParam, q, qt, acc, weights[], totalWeight;

		assert nonAbsMin.equals(nonAbsMax);
		BitSet nonAbs = nonAbsMin;
		if (previousResult == null)
			previousResult = new BoxRegionValues(model);

		// Store num states
		n = model.getNumStates();

		// Optimisations: If (nonAbs is empty or t = 0) and multProbs is null, this is easy.
		if (((nonAbsMin != null && nonAbsMin.isEmpty()) || t == 0) &&
				((nonAbsMax != null && nonAbsMax.isEmpty()) || t == 0) &&
				multProbs == null) {
			solnMin = Utils.bitsetToDoubleArray(targetMin, n);
			solnMax = Utils.bitsetToDoubleArray(targetMax, n);
			return new BoxRegionValues(model, regionFactory.completeSpace(), solnMin, solnMax);
		}

		// Start backwards transient computation
		timer = System.currentTimeMillis();
		mainLog.println("\nStarting backwards transient probability computation...");

		// Compute the in, out, inout sets of reactions
		model.computeInOutReactions();

		// Get uniformisation rate; do Fox-Glynn
		q = model.getDefaultUniformisationRate(nonAbs);
		qt = q * t;
		mainLog.println("\nUniformisation: q.t = " + q + " x " + t + " = " + qt);
		termCritParam = 1e-6;
		acc = termCritParam / 8.0;
		fg = new FoxGlynn(qt, 1e-300, 1e+300, acc);
		left = fg.getLeftTruncationPoint();
		right = fg.getRightTruncationPoint();
		if (right < 0)
			throw new PrismException("Overflow in Fox-Glynn computation (time bound too big?)");
		weights = fg.getWeights();
		totalWeight = fg.getTotalWeight();
		for (i = left; i <= right; i++)
			weights[i - left] /= totalWeight;
		mainLog.println("Fox-Glynn (" + acc + "): left = " + left + ", right = " + right);

		totalIters = 0;
		for (Entry<BoxRegion, BoxRegionValues.StateValuesPair> entry : multProbs) {
			BoxRegion region = entry.getKey();

			// If the previous region values contain probs for this region, i.e. the region
			// has not been decomposed, then just use the previous result directly.
			if (previousResult.hasRegion(region)) {
				regionValues.put(region, previousResult.getMin(region), previousResult.getMax(region));
				continue;
			}

			double[] multProbsMin = entry.getValue().getMin().getDoubleArray();
			double[] multProbsMax = entry.getValue().getMax().getDoubleArray();

			// Configure parameter space
			model.setRegion(region);
			mainLog.println("Computing probabilities for parameter region " + region);

			// Create solution vectors
			solnMin = new double[n];
			soln2Min = new double[n];
			sumMin = new double[n];
			solnMax = new double[n];
			soln2Max = new double[n];
			sumMax = new double[n];

			// Initialise solution vectors.
			// Vectors soln/soln2 are multProbs[i] for target states.
			// Vector sum is all zeros (done by array creation).
			for (i = 0; i < n; i++) {
				solnMin[i] = soln2Min[i] = targetMin.get(i) ? multProbsMin[i] : 0.0;
				solnMax[i] = soln2Max[i] = targetMax.get(i) ? multProbsMax[i] : 0.0;
			}

			// If necessary, do 0th element of summation (doesn't require any matrix powers)
			if (left == 0) {
				for (i = 0; i < n; i++) {
					sumMin[i] += weights[0] * solnMin[i];
					sumMax[i] += weights[0] * solnMax[i];
				}
			}

			// Start iterations
			iters = 1;
			while (iters <= right) {
				// Matrix-vector multiply				
				model.mvMult(solnMin, soln2Min, solnMax, soln2Max, nonAbs, false, q);

				// Swap vectors for next iter
				tmpsoln = solnMin;
				solnMin = soln2Min;
				soln2Min = tmpsoln;
				tmpsoln = solnMax;
				solnMax = soln2Max;
				soln2Max = tmpsoln;

				// Add to sum
				if (iters >= left) {
					for (i = 0; i < n; i++) {
						sumMin[i] += weights[iters - left] * solnMin[i];
						sumMax[i] += weights[iters - left] * solnMax[i];;
					}
					decompositionProcedure.examineSingleIteration(region, sumMin, sumMax);
				}

				iters++;
				totalIters++;
			}

			// Store result
			regionValues.put(region, sumMin, sumMax);
		}

		decompositionProcedure.examineWholeComputation(regionValues);

		// Finished bounded probabilistic reachability
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Backwards transient probability computation");
		mainLog.print(" took " + totalIters + " iters");
		mainLog.println(" and " + timer / 1000.0 + " seconds.");

		return regionValues;
	}

	// Transient analysis

	public BoxRegionValues doTransient(PSEModel model, double t, StateValues initDistMin, StateValues initDistMax, DecompositionProcedure decompositionProcedure) throws PrismException
	{
		StateValues initDistMinNew = null, initDistMaxNew = null;

		// Build initial distribution (if not specified)
		if (initDistMin == null) {
			initDistMinNew = new StateValues(TypeDouble.getInstance(), new Double(0.0), model);
			double initVal = 1.0 / model.getNumInitialStates();
			for (int in : model.getInitialStates()) {
				initDistMinNew.setDoubleValue(in, initVal);
			}
		} else {
			initDistMinNew = initDistMin;
		}
		if (initDistMax == null) {
			initDistMaxNew = new StateValues(TypeDouble.getInstance(), new Double(0.0), model);
			double initVal = 1.0 / model.getNumInitialStates();
			for (int in : model.getInitialStates()) {
				initDistMaxNew.setDoubleValue(in, initVal);
			}
		} else {
			initDistMaxNew = initDistMax;
		}

		// Compute transient probabilities
		return computeTransientProbs(model, t, initDistMinNew.getDoubleArray(), initDistMaxNew.getDoubleArray(), decompositionProcedure);
	}

	/**
	 * NB: Decompositions of the parameter space are performed implicitly.
	 */
	// TODO: Perform decompositions explicitly from doTransient analogically to backwards transient computation,
	// i.e. replace the double arrays initDist{Min,Max} with BoxRegionValues.
	public BoxRegionValues computeTransientProbs(PSEModel model, double t, double initDistMin[], double initDistMax[], DecompositionProcedure decompositionProcedure)
			throws PrismException
	{
		BoxRegionValues regionValues = new BoxRegionValues(model);
		int i, n, iters, totalIters;
		double solnMin[], soln2Min[], sumMin[];
		double solnMax[], soln2Max[], sumMax[];
		double tmpsoln[];
		long timer;
		// Fox-Glynn stuff
		FoxGlynn fg;
		int left, right;
		double termCritParam, q, qt, acc, weights[], totalWeight;

		// For decomposing the parameter space
		LinkedList<BoxRegion> regions = new LinkedList<BoxRegion>();
		regions.add(regionFactory.completeSpace());

		// Start bounded probabilistic reachability
		timer = System.currentTimeMillis();
		mainLog.println("\nStarting transient probability computation...");

		// Compute the in, out, inout sets of reactions
		model.computeInOutReactions();

		// Store num states
		n = model.getNumStates();

		// Get uniformisation rate; do Fox-Glynn
		q = model.getDefaultUniformisationRate();
		qt = q * t;
		mainLog.println("\nUniformisation: q.t = " + q + " x " + t + " = " + qt);
		termCritParam = 1e-6;
		acc = termCritParam / 8.0;
		fg = new FoxGlynn(qt, 1e-300, 1e+300, acc);
		left = fg.getLeftTruncationPoint();
		right = fg.getRightTruncationPoint();
		if (right < 0)
			throw new PrismException("Overflow in Fox-Glynn computation (time bound too big?)");
		weights = fg.getWeights();
		totalWeight = fg.getTotalWeight();
		for (i = left; i <= right; i++)
			weights[i - left] /= totalWeight;
		mainLog.println("Fox-Glynn (" + acc + "): left = " + left + ", right = " + right);

		totalIters = 0;
		while (regions.size() != 0) {
			// Create solution vectors
			solnMin = new double[n];
			soln2Min = new double[n];
			sumMin = new double[n];
			solnMax = new double[n];
			soln2Max = new double[n];
			sumMax = new double[n];

			// Initialise solution vectors.
			// Don't need to do soln2 since will be immediately overwritten.
			// Vector sum is all zeros (done by array creation).
			solnMin = initDistMin.clone();
			solnMax = initDistMax.clone();

			// If necessary, do 0th element of summation (doesn't require any matrix powers)
			if (left == 0) {
				for (i = 0; i < n; i++) {
					sumMin[i] += weights[0] * solnMin[i];
					sumMax[i] += weights[0] * solnMax[i];
				}
			}

			// Configure parameter space
			BoxRegion region = regions.remove();
			model.setRegion(region);
			mainLog.println("Computing probabilities for parameter region " + region);

			try {
				// Start iterations
				iters = 1;
				totalIters++;
				while (iters <= right) {
					// Vector-matrix multiply
					model.vmMult(solnMin, soln2Min, solnMax, soln2Max, q);

					// Swap vectors for next iter
					tmpsoln = solnMin;
					solnMin = soln2Min;
					soln2Min = tmpsoln;
					tmpsoln = solnMax;
					solnMax = soln2Max;
					soln2Max = tmpsoln;

					// Add to sum
					if (iters >= left) {
						for (i = 0; i < n; i++) {
							sumMin[i] += weights[iters - left] * solnMin[i];
							sumMax[i] += weights[iters - left] * solnMax[i];
						}
						decompositionProcedure.examineSingleIteration(region, sumMin, sumMax);
					}

					iters++;
					totalIters++;
				}

				// Store result
				regionValues.put(region, sumMin, sumMax);
			} catch (DecompositionProcedure.DecompositionNeeded e) {
				e.getRegionsToDecompose().print(mainLog);
				for (BoxRegion regionToDecompose : e.getRegionsToDecompose()) {
					regions.add(regionToDecompose.lowerHalf());
					regions.add(regionToDecompose.upperHalf());
				}
			}
		}

		// TODO call decompositionProcedure.examineWholeComputation()?

		// Finished bounded probabilistic reachability
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Transient probability computation");
		mainLog.print(" took " + totalIters + " iters");
		mainLog.print(" and " + timer / 1000.0 + " seconds");
		mainLog.println(" (producing " + regionValues.getNumRegions() + " final regions).");

		return regionValues;
	}
}
