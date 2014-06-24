//==============================================================================
//	
//	Copyright (c) 2013-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
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
import explicit.Model;
import explicit.StateModelChecker;
import explicit.StateValues;
import explicit.Utils;

final public class PSEModelChecker extends PrismComponent
{
	// Log for output (default to System.out)
	private PrismLog mainLog = new PrismPrintStreamLog(System.out);

	// Model file (for reward structures, etc.)
	private ModulesFile modulesFile = null;

	// Properties file (for labels, constants, etc.)
	private PropertiesFile propertiesFile = null;

	// Constants (extracted from model/properties)
	private Values constantValues;

	// The result of model checking will be stored here
	//private Result result;

	private StateModelChecker stateChecker;

	/**
	 * Constructor
	 */
	public PSEModelChecker(PrismComponent parent) throws PrismException
	{
		super(parent);
		stateChecker = new StateModelChecker(this);
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

	// Model checking functions

	/**
	 * Model check an expression, process and return the result.
	 * Information about states and model constants should be attached to the model.
	 */
	public Result check(Model model, Expression expr) throws PrismException
	{
		//PSEModel paramModel = (PSEModel) model;

		ExpressionFilter exprFilter = null;
		long timer = 0;

		// Remove labels from property, using combined label list (on a copy of the expression)
		// This is done now so that we can handle labels nested below operators that are not
		// handled natively by the model checker yet (just evaluate()ed in a loop).
		expr = (Expression) expr.deepCopy().expandLabels(propertiesFile.getCombinedLabelList());

		// Also evaluate/replace any constants
		//expr = (Expression) expr.replaceConstants(constantValues);

		// *** TODO replace the following blocks with ExpressionFilter.addDefaultFilterIfNeeded, cf. SVN 8394

		// The final result of model checking will be a single value. If the expression to be checked does not
		// already yield a single value (e.g. because a filter has not been explicitly included), we need to wrap
		// a new (invisible) filter around it. Note that some filters (e.g. print/argmin/argmax) also do not
		// return single values and have to be treated in this way.
		if (!expr.returnsSingleValue()) {
			// New filter depends on expression type and number of initial states.
			// Boolean expressions...
			if (expr.getType() instanceof TypeBool) {
				// Result is true iff true for all initial states
				exprFilter = new ExpressionFilter("forall", expr, new ExpressionLabel("init"));
			}
			// Non-Boolean (double or integer) expressions...
			else {
				// Result is for the initial state, if there is just one,
				// or the range over all initial states, if multiple
				if (model.getNumInitialStates() == 1) {
					exprFilter = new ExpressionFilter("state", expr, new ExpressionLabel("init"));
				} else {
					exprFilter = new ExpressionFilter("range", expr, new ExpressionLabel("init"));
				}
			}
		}
		// Even, when the expression does already return a single value, if the the outermost operator
		// of the expression is not a filter, we still need to wrap a new filter around it.
		// e.g. 2*filter(...) or 1-P=?[...{...}]
		// This because the final result of model checking is only stored when we process a filter.
		else if (!(expr instanceof ExpressionFilter)) {
			// We just pick the first value (they are all the same)
			exprFilter = new ExpressionFilter("first", expr, new ExpressionLabel("init"));
			// We stop any additional explanation being displayed to avoid confusion.
			exprFilter.setExplanationEnabled(false);
		}

		// For any case where a new filter was created above...
		if (exprFilter != null) {
			// Make it invisible (not that it will be displayed)
			exprFilter.setInvisible(true);
			// Compute type of new filter expression (will be same as child)
			exprFilter.typeCheck();
			// Store as expression to be model checked
			expr = exprFilter;
		}

		// *** end TODO above

		// Do model checking and store result vector
		timer = System.currentTimeMillis();
		BoxRegionValues regionValues = checkExpression(model, expr);
		timer = System.currentTimeMillis() - timer;
		mainLog.println("\nTime for model checking: " + timer / 1000.0 + " seconds.");

		// Store result
		Result result = new Result();
		result.setResult(regionValues);

		// Print result to log
		String resultString = "Result";
		if (!("Result".equals(expr.getResultName())))
			resultString += " (" + expr.getResultName().toLowerCase() + ")";
		resultString += ": " + result.getResultString();
		mainLog.print("\n" + resultString);

		return result;
	}

	public BoxRegionValues checkExpression(Model model, Expression expr) throws PrismException
	{
		if (expr instanceof ExpressionFilter) {
			return checkExpressionFilter(model, (ExpressionFilter) expr);
		}
		if (expr instanceof ExpressionProb) {
			return checkExpressionProb(model, (ExpressionProb) expr);
		}
		// Let explicit.StateModelChecker take care of other expressions
		StateValues vals = stateChecker.checkExpression(model, expr);
		return new BoxRegionValues(model, BoxRegion.completeSpace, vals, vals);
	}

	/**
	 * Model check a filter.
	 */
	protected BoxRegionValues checkExpressionFilter(Model model, ExpressionFilter expr) throws PrismException
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

		BitSet bsFilterMin = checkExpression(model, filter).getMin(BoxRegion.completeSpace).getBitSet();
		BitSet bsFilterMax = checkExpression(model, filter).getMax(BoxRegion.completeSpace).getBitSet();
		boolean filterInit = (filter instanceof ExpressionLabel && ((ExpressionLabel) filter).getName().equals("init"));
		boolean filterInitSingle = filterInit & model.getNumInitialStates() == 1;
		if (bsFilterMin.isEmpty()) {
			throw new PrismException("Filter satisfies no states");
		}
		if (!filterInit) {
			mainLog.println("\nStates satisfying filter " + filter + ":");
			mainLog.println("MIN= " + bsFilterMin.cardinality());
			mainLog.println("MAX = " + bsFilterMax.cardinality());
		}

		BoxRegionValues rgnVals = checkExpression(model, expr.getOperand());
		StateValues rgnValsMin = rgnVals.getMin(BoxRegion.completeSpace);
		StateValues rgnValsMax = rgnVals.getMax(BoxRegion.completeSpace);

		// Compute result according to filter type
		FilterOperator op = expr.getOperatorType();
		switch (op) {
		case PRINT:
		case PRINTALL:
			if (expr.getType() instanceof TypeBool) {
				mainLog.print("\nSatisfying states");
				mainLog.println(filterTrue ? ":" : " that are also in filter " + filter + ":");
				mainLog.println("\nMIN\n");
				rgnValsMin.printFiltered(mainLog, bsFilterMin);
				mainLog.println("\nMAX\n");
				rgnValsMax.printFiltered(mainLog, bsFilterMax);
			} else {
				if (op == FilterOperator.PRINT) {
					mainLog.println("\nResults (non-zero only) for filter " + filter + ":");
					mainLog.println("\nMIN\n");
					rgnValsMin.printFiltered(mainLog, bsFilterMin);
					mainLog.println("\nMAX\n");
					rgnValsMax.printFiltered(mainLog, bsFilterMax);
				} else {
					mainLog.println("\nResults (including zeros) for filter " + filter + ":");
					mainLog.println("\nMIN\n");
					rgnValsMin.printFiltered(mainLog, bsFilterMin, false, false, true, true);
					mainLog.println("\nMAX\n");
					rgnValsMax.printFiltered(mainLog, bsFilterMax, false, false, true, true);
				}
			}
			resRgnVals = rgnVals;
			break;
		case MIN:
		case MAX:
		case ARGMIN:
		case ARGMAX:
			throw new UnsupportedOperationException();
		case COUNT:
			resObjMin = new Integer(rgnValsMin.countOverBitSet(bsFilterMin));
			resObjMax = new Integer(rgnValsMax.countOverBitSet(bsFilterMax));
			resRgnValsMin = new StateValues(expr.getType(), resObjMin, model);
			resRgnValsMax = new StateValues(expr.getType(), resObjMax, model);
			resRgnVals = new BoxRegionValues(model, BoxRegion.completeSpace, resRgnValsMin, resRgnValsMax);

			resultExpl = filterTrue ? "Count of satisfying states" : "Count of satisfying states also in filter";
			mainLog.println("\n" + resultExpl + ":");
			mainLog.println("MIN = " + resObjMin);
			mainLog.println("MAX = " + resObjMax);
			break;
		case SUM:
		case AVG:
			throw new PrismException("Operation not implemented for parametric models");
		case FIRST:
			resObjMin = rgnValsMin.firstFromBitSet(bsFilterMin);
			resObjMax = rgnValsMax.firstFromBitSet(bsFilterMax);
			resRgnValsMin = new StateValues(expr.getType(), resObjMin, model);
			resRgnValsMax = new StateValues(expr.getType(), resObjMax, model);
			resRgnVals = new BoxRegionValues(model, BoxRegion.completeSpace, resRgnValsMin, resRgnValsMax);

			resultExpl = "Value in ";
			if (filterInit) {
				resultExpl += filterInitSingle ? "the initial state" : "first initial state";
			} else {
				resultExpl += filterTrue ? "the first state" : "first state satisfying filter";
			}
			mainLog.println("\n" + resultExpl + ":");
			mainLog.println("MIN = " + resObjMin);
			mainLog.println("MAX = " + resObjMax);
			break;
		case RANGE:
			resObjMin = new prism.Interval(rgnValsMin.minOverBitSet(bsFilterMin), rgnValsMin.maxOverBitSet(bsFilterMin));
			resObjMax = new prism.Interval(rgnValsMax.minOverBitSet(bsFilterMax), rgnValsMax.maxOverBitSet(bsFilterMax));
			resRgnVals = rgnVals;

			resultExpl = "Range of values over ";
			resultExpl += filterInit ? "initial states" : filterStatesString;
			mainLog.println("\n" + resultExpl + ":");
			mainLog.println("MIN = " + resObjMin);
			mainLog.println("MAX = " + resObjMax);
			break;
		case FORALL:
			bsMin = rgnValsMin.getBitSet();
			bsMax = rgnValsMax.getBitSet();

			mainLog.println("\nNumber of states satisfying " + expr.getOperand() + ":");
			mainLog.print("MIN = " + bsMin.cardinality());
			mainLog.println(bsMin.cardinality() == model.getNumStates() ? " (all in model)" : "");
			mainLog.print("MAX = " + bsMax.cardinality());
			mainLog.println(bsMax.cardinality() == model.getNumStates() ? " (all in model)" : "");

			resObjMin = new Boolean(rgnValsMin.forallOverBitSet(bsFilterMin));
			resObjMax = new Boolean(rgnValsMax.forallOverBitSet(bsFilterMax));
			resRgnValsMin = new StateValues(expr.getType(), resObjMin, model);
			resRgnValsMax = new StateValues(expr.getType(), resObjMax, model);
			resRgnVals = new BoxRegionValues(model, BoxRegion.completeSpace, resRgnValsMin, resRgnValsMax);

			mainLog.print("\nProperty satisfied in {");
			mainLog.print("MIN = " + rgnValsMin.countOverBitSet(bsFilterMin) + ", ");
			mainLog.print("MAX = " + rgnValsMin.countOverBitSet(bsFilterMin) + "}");
			if (filterInit) {
				mainLog.println(" of " + model.getNumInitialStates() + " initial states.");
			} else {
				if (filterTrue) {
					mainLog.println(" of all " + model.getNumStates() + " states.");
				} else {
					mainLog.print(" of {");
					mainLog.print("MIN = " + bsFilterMin.cardinality() + ", ");
					mainLog.print("MAX = " + bsFilterMax.cardinality() + "}");
					mainLog.println(" filter states.");
				}
			}
			break;
		case EXISTS:
			bsMin = rgnValsMin.getBitSet();
			bsMax = rgnValsMax.getBitSet();
			resObjMin = new Boolean(rgnValsMin.existsOverBitSet(bsFilterMin));
			resObjMax = new Boolean(rgnValsMax.existsOverBitSet(bsFilterMax));
			resRgnValsMin = new StateValues(expr.getType(), resObjMin, model);
			resRgnValsMax = new StateValues(expr.getType(), resObjMax, model);
			resRgnVals = new BoxRegionValues(model, BoxRegion.completeSpace, resRgnValsMin, resRgnValsMax);

			// Create explanation of result and print some details to log
			resultExpl = "MIN = Property satisfied in ";
			if (filterTrue) {
				resultExpl += ((boolean) resObjMin) ? "at least one state" : "no states";
			} else {
				resultExpl += ((boolean) resObjMin) ? "at least one filter state" : "no filter states";
			}
			resultExpl = "\nMAX = Property satisfied in ";
			if (filterTrue) {
				resultExpl += ((boolean) resObjMax) ? "at least one state" : "no states";
			} else {
				resultExpl += ((boolean) resObjMax) ? "at least one filter state" : "no filter states";
			}
			mainLog.println("\n" + resultExpl);
			break;
		case STATE:
			// Check filter satisfied by exactly one state
			if (bsFilterMin.cardinality() != 1 || bsFilterMax.cardinality() != 1) {
				String s = "Filter should be satisfied in exactly 1 state";
				s += " (but \"" + filter + "\" is true in {";
				s += "MIN = " + bsFilterMin.cardinality() + ", ";
				s += "MAX = " + bsFilterMax.cardinality() + "} states)";
				throw new PrismException(s);
			}

			resObjMin = rgnValsMin.firstFromBitSet(bsFilterMin);
			resObjMax = rgnValsMax.firstFromBitSet(bsFilterMax);
			resRgnValsMin = new StateValues(expr.getType(), resObjMin, model);
			resRgnValsMax = new StateValues(expr.getType(), resObjMax, model);
			resRgnVals = new BoxRegionValues(model, BoxRegion.completeSpace, resRgnValsMin, resRgnValsMax);

			// Create explanation of result and print some details to log
			resultExpl = "Value in ";
			if (filterInit) {
				resultExpl += "the initial state";
			} else {
				resultExpl += "the filter state";
			}
			mainLog.println("\n" + resultExpl + ":");
			mainLog.println("MIN = " + resObjMin);
			mainLog.println("MAX = " + resObjMax);
			break;
		default:
			throw new PrismException("Unrecognised filter type \"" + expr.getOperatorName() + "\"");
		}

		return resRgnVals;
	}

	/**
	 * Model check a P operator expression and return the values for all states.
	 */
	protected BoxRegionValues checkExpressionProb(Model model, ExpressionProb expr) throws PrismException
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
			regionValues = checkProbPathFormula(model, expr.getExpression());
			break;
		default:
			throw new PrismException("Cannot model check " + expr + " for a " + modelType);
		}

		// Print out probabilities
		/*
		if (settings.getBoolean(PrismSettings.PRISM_VERBOSE)) {
			mainLog.print("\nProbabilities (non-zero only) for all states:\n");
			res.print(mainLog);
		}
		*/

		// For =? properties, just return values
		if (pb == null) {
			return regionValues;
		}
		// Otherwise, compare against bound to get set of satisfying states
		else {
			for (Entry<BoxRegion, BoxRegionValues.StateValuesPair> entry : regionValues) {
				BitSet solMin = entry.getValue().getMin().getBitSetFromInterval(relOp, p);
				BitSet solMax = entry.getValue().getMax().getBitSetFromInterval(relOp, p);
				regionValues.put(entry.getKey(), solMin, solMax);
			}
			return regionValues;
		}
	}

	/**
	 * Compute probabilities for the contents of a P operator.
	 */
	protected BoxRegionValues checkProbPathFormula(Model model, Expression expr) throws PrismException
	{
		// Test whether this is a simple path formula (i.e. PCTL)
		// and then pass control to appropriate method.
		if (expr.isSimplePathFormula()) {
			return checkProbPathFormulaSimple(model, expr);
		} else {
			throw new PrismException("Operation not implemented for parametric models");
		}
	}

	/**
	 * Compute probabilities for a simple, non-LTL path operator.
	 */
	protected BoxRegionValues checkProbPathFormulaSimple(Model model, Expression expr) throws PrismException
	{
		BoxRegionValues regionValues = null;

		// Negation/parentheses
		if (expr instanceof ExpressionUnaryOp) {
			ExpressionUnaryOp exprUnary = (ExpressionUnaryOp) expr;
			// Parentheses
			if (exprUnary.getOperator() == ExpressionUnaryOp.PARENTH) {
				// Recurse
				regionValues = checkProbPathFormulaSimple(model, exprUnary.getOperand());
			}
			// Negation
			else if (exprUnary.getOperator() == ExpressionUnaryOp.NOT) {
				// TODO
				/*
				// Compute, then subtract from 1
				mcRes = checkProbPathFormulaSimple(model, exprUnary.getOperand());
				mcRes.timesConstant(-1.0);
				mcRes.plusConstant(1.0);
				*/
				throw new UnsupportedOperationException();
			}
		}
		// Temporal operators
		else if (expr instanceof ExpressionTemporal) {
			ExpressionTemporal exprTemp = (ExpressionTemporal) expr;
			// Next
			if (exprTemp.getOperator() == ExpressionTemporal.P_X) {
				// TODO
				//probs = checkProbNext(model, exprTemp);
				throw new UnsupportedOperationException();
			}
			// Until
			else if (exprTemp.getOperator() == ExpressionTemporal.P_U) {
				if (exprTemp.hasBounds()) {
					regionValues = checkProbBoundedUntil(model, exprTemp);
				} else {
					// TODO
					//probs = checkProbUntil(model, exprTemp);
					throw new UnsupportedOperationException();
				}
			}
			// Anything else - convert to until and recurse
			else {
				regionValues = checkProbPathFormulaSimple(model, exprTemp.convertToUntilForm());
			}
		}

		if (regionValues == null)
			throw new PrismException("Unrecognised path operator in P operator");

		return regionValues;
	}

	/**
	 * Model check a time-bounded until operator; return vector of probabilities for all states.
	 */
	protected BoxRegionValues checkProbBoundedUntil(Model model, ExpressionTemporal expr) throws PrismException
	{
		double lTime, uTime; // time bounds
		Expression exprTmp;
		BitSet b1Min, b1Max, b2Min, b2Max, tmpMin, tmpMax;
		BoxRegionValues regionValues, tmpRegionValues;
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
		b1Min = checkExpression(model, expr.getOperand1()).getMin(BoxRegion.completeSpace).getBitSet();
		b1Max = checkExpression(model, expr.getOperand1()).getMax(BoxRegion.completeSpace).getBitSet();
		b2Min = checkExpression(model, expr.getOperand2()).getMin(BoxRegion.completeSpace).getBitSet();
		b2Max = checkExpression(model, expr.getOperand2()).getMax(BoxRegion.completeSpace).getBitSet();

		// compute probabilities

		// a trivial case: "U<=0"
		if (lTime == 0 && uTime == 0) {
			// prob is 1 in b2 states, 0 otherwise
			probsMin = StateValues.createFromBitSetAsDoubles(b2Min, model);
			probsMax = StateValues.createFromBitSetAsDoubles(b2Max, model);
			regionValues = new BoxRegionValues(model, BoxRegion.completeSpace, probsMin, probsMax);
		} else {
			// break down into different cases to compute probabilities

			// >= lTime
			if (uTime == -1) {
				// check for special case of lTime == 0, this is actually an unbounded until
				if (lTime == 0) {
					// TODO
					// compute probs
					//regionValues = computeUntilProbs((PSEModel) model, b1, b2);
					throw new UnsupportedOperationException();
				} else {
					// TODO
					/*
					// compute unbounded until probs
					tmpRegionValues = computeUntilProbs((PSEModel) model, b1, b2);
					// compute bounded until probs
					regionValues = computeTransientBackwardsProbsRanged((PSEModel) model, b1, b1, lTime, tmpRegionValues);
					*/
					throw new UnsupportedOperationException();
				}
			}
			// <= uTime
			else if (lTime == 0) {
				// nb: uTime != 0 since would be caught above (trivial case)
				b1Min.andNot(b2Min);
				b1Max.andNot(b2Max);
				regionValues = computeTransientBackwardsProbsRanged((PSEModel) model, b2Min, b1Min, b2Max, b1Max, uTime, null, null);
				// set values to exactly 1 for target (b2) states
				// (these are computed inexactly during uniformisation)
				int n = model.getNumStates();
				for (int i = 0; i < n; i++) {
					if (b2Min.get(i))
						regionValues.getMin(BoxRegion.completeSpace).setDoubleValue(i, 1);
					if (b2Max.get(i))
						regionValues.getMax(BoxRegion.completeSpace).setDoubleValue(i, 1);
				}
			}
			// [lTime,uTime] (including where lTime == uTime)
			else {
				tmpMin = (BitSet) b1Min.clone();
				tmpMin.andNot(b2Min);
				tmpMax = (BitSet) b1Max.clone();
				tmpMax.andNot(b2Max);
				tmpRegionValues = computeTransientBackwardsProbsRanged((PSEModel) model, b2Min, tmpMin, b2Max, tmpMax, uTime - lTime, null, null);
				double[] multProbsMin = tmpRegionValues.getMin(BoxRegion.completeSpace).getDoubleArray();
				double[] multProbsMax = tmpRegionValues.getMax(BoxRegion.completeSpace).getDoubleArray();
				regionValues = computeTransientBackwardsProbsRanged((PSEModel) model, b1Min, b1Min, b1Max, b1Max, lTime, multProbsMin, multProbsMax);
			}
		}

		return regionValues;
	}

	public BoxRegionValues computeTransientBackwardsProbsRanged(PSEModel ctmcRanged,
			BitSet targetMin, BitSet nonAbsMin, BitSet targetMax, BitSet nonAbsMax,
			double t, double multProbsMin[], double multProbsMax[]) throws PrismException
	{
		BoxRegionValues regionValues = new BoxRegionValues(ctmcRanged);
		int i, n, iters;
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

		// Optimisations: If (nonAbs is empty or t = 0) and multProbs is null, this is easy.
		if ((((nonAbsMin != null && nonAbsMin.isEmpty()) || (t == 0)) && multProbsMin == null) &&
				(((nonAbsMax != null && nonAbsMax.isEmpty()) || (t == 0)) && multProbsMax == null)) {
			solnMin = Utils.bitsetToDoubleArray(targetMin, ctmcRanged.getNumStates());
			solnMax = Utils.bitsetToDoubleArray(targetMax, ctmcRanged.getNumStates());
			return new BoxRegionValues(ctmcRanged, BoxRegion.completeSpace, solnMin, solnMax);
		}

		// Start backwards transient computation
		timer = System.currentTimeMillis();
		mainLog.println("\nStarting backwards transient probability computation...");

		mainLog.println("\nComputing in, out, inout reactions...");
		ctmcRanged.computeInOutReactions();

		// Store num states
		n = ctmcRanged.getNumStates();

		// Get uniformisation rate; do Fox-Glynn
		q = ctmcRanged.getDefaultUniformisationRate(nonAbs);
		qt = q * t;
		mainLog.println("\nUniformisation: q.t = " + q + " x " + t + " = " + qt);
		termCritParam = 1e-6;
		acc = termCritParam / 8.0;
		fg = new FoxGlynn(qt, 1e-300, 1e+300, acc);
		left = fg.getLeftTruncationPoint();
		right = fg.getRightTruncationPoint();
		if (right < 0) {
			throw new PrismException("Overflow in Fox-Glynn computation (time bound too big?)");
		}
		weights = fg.getWeights();
		totalWeight = fg.getTotalWeight();
		for (i = left; i <= right; i++) {
			weights[i - left] /= totalWeight;
		}
		mainLog.println("Fox-Glynn (" + acc + "): left = " + left + ", right = " + right);

		// Create solution vector(s)
		solnMin = new double[n];
		soln2Min = new double[n];
		sumMin = new double[n];
		solnMax = new double[n];
		soln2Max = new double[n];
		sumMax = new double[n];

		// Initialise solution vectors.
		// Vectors soln/soln2 are 1 for target states, or multProbs[i] if supplied.
		// Vector sum is all zeros (done by array creation).
		// XXX: soln2 is immediately overwritten here as well, isn't it?
		if (multProbsMin != null) {
			for (i = 0; i < n; i++)
				solnMin[i] = soln2Min[i] = targetMin.get(i) ? multProbsMin[i] : 0.0;
		} else {
			for (i = 0; i < n; i++)
				solnMin[i] = soln2Min[i] = targetMin.get(i) ? 1.0 : 0.0;
		}
		if (multProbsMax != null) {
			for (i = 0; i < n; i++)
				solnMax[i] = soln2Max[i] = targetMax.get(i) ? multProbsMax[i] : 0.0;
		} else {
			for (i = 0; i < n; i++)
				solnMax[i] = soln2Max[i] = targetMax.get(i) ? 1.0 : 0.0;
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
			ctmcRanged.mvMult(solnMin, soln2Min, solnMax, soln2Max, nonAbs, false, q);

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
			}
			iters++;
		}

		// Store result
		regionValues.put(BoxRegion.completeSpace, sumMin, sumMax);

		// Finished bounded probabilistic reachability
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Backwards transient probability computation");
		mainLog.print(" took " + iters + " iters");
		mainLog.print(" and " + timer / 1000.0 + " seconds in total.\n");

		return regionValues;
	}

	/**
	 * Compute probabilities for an (unbounded) until operator.
	 */
	protected BoxRegionValues checkProbUntil(Model model, ExpressionTemporal expr) throws PrismException
	{
		BitSet b1, b2;
		BoxRegionValues probs = null;

		// TODO
		/*
		// model check operands first
		b1 = checkExpression(model, expr.getOperand1()).getBitSet();
		b2 = checkExpression(model, expr.getOperand2()).getBitSet();

		probs = computeUntilProbs((PSEModel) model, b1, b2);

		return probs;
		*/

		throw new UnsupportedOperationException();
	}

	public BoxRegionValues computeUntilProbs(PSEModel ctmcRanged, BitSet remain, BitSet target) throws PrismException
	{
		// TODO
		//return computeReachProbs(ctmcRanged, remain, target, null, null);
		throw new UnsupportedOperationException();
	}

	// Transient analysis

	public BoxRegionValues doTransientRanged(Model model, double t, double accuracy, explicit.StateValues initDistMin, explicit.StateValues initDistMax) throws PrismException
	{
		PSEModel ctmcRanged = (PSEModel) model;
		BoxRegionValues res = null;
		explicit.StateValues initDistMinNew = null, initDistMaxNew = null;

		// Build initial distribution (if not specified)
		if (initDistMin == null) {
			initDistMinNew = new explicit.StateValues(TypeDouble.getInstance(), new Double(0.0), ctmcRanged);
			double initVal = 1.0 / ctmcRanged.getNumInitialStates();
			for (int in : ctmcRanged.getInitialStates()) {
				initDistMinNew.setDoubleValue(in, initVal);
			}
		} else {
			initDistMinNew = initDistMin;
		}
		if (initDistMax == null) {
			initDistMaxNew = new explicit.StateValues(TypeDouble.getInstance(), new Double(0.0), ctmcRanged);
			double initVal = 1.0 / ctmcRanged.getNumInitialStates();
			for (int in : ctmcRanged.getInitialStates()) {
				initDistMaxNew.setDoubleValue(in, initVal);
			}
		} else {
			initDistMaxNew = initDistMax;
		}
		
		// Compute transient probabilities
		res = computeTransientProbsRanged(ctmcRanged, t, accuracy, initDistMinNew.getDoubleArray(), initDistMaxNew.getDoubleArray());

		return res;
	}

	public BoxRegionValues computeTransientProbsRanged(PSEModel ctmcRanged, double t, double accuracy, double initDistMin[], double initDistMax[]) throws PrismException
	{
		BoxRegionValues regionValues = new BoxRegionValues(ctmcRanged);
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
		regions.add(new BoxRegion(0.0, 1.0));
		int numDecompositions = 0;

		// Start bounded probabilistic reachability
		timer = System.currentTimeMillis();
		mainLog.println("\nStarting transient probability computation...");
		totalIters = 0;

		mainLog.println("\nComputing in, out, inout reactions...");
		ctmcRanged.computeInOutReactions();

		// Store num states
		n = ctmcRanged.getNumStates();

		// Get uniformisation rate; do Fox-Glynn
		q = ctmcRanged.getDefaultUniformisationRate();
		qt = q * t;
		mainLog.println("\nUniformisation: q.t = " + q + " x " + t + " = " + qt);
		termCritParam = 1e-6;
		acc = termCritParam / 8.0;
		fg = new FoxGlynn(qt, 1e-300, 1e+300, acc);
		left = fg.getLeftTruncationPoint();
		right = fg.getRightTruncationPoint();
		if (right < 0) {
			throw new PrismException("Overflow in Fox-Glynn computation (time bound too big?)");
		}
		weights = fg.getWeights();
		totalWeight = fg.getTotalWeight();
		for (i = left; i <= right; i++) {
			weights[i - left] /= totalWeight;
		}
		mainLog.println("Fox-Glynn (" + acc + "): left = " + left + ", right = " + right);

		while (regions.size() != 0) {
			// Create solution vector(s)
			solnMin = new double[n];
			soln2Min = new double[n];
			sumMin = new double[n];
			solnMax = new double[n];
			soln2Max = new double[n];
			sumMax = new double[n];

			// Initialise solution vectors
			// (don't need to do soln2 since will be immediately overwritten)
			solnMin = initDistMin.clone();
			solnMax = initDistMax.clone();
			for (i = 0; i < n; i++) {
				sumMin[i] = 0.0;
				sumMax[i] = 0.0;
			}

			// If necessary, do 0th element of summation (doesn't require any matrix powers)
			if (left == 0) {
				for (i = 0; i < n; i++) {
					sumMin[i] += weights[0] * solnMin[i];
					sumMax[i] += weights[0] * solnMax[i];
				}
			}

			// Shrink the parameter space
			BoxRegion region = regions.remove();
			ctmcRanged.scaleParameterSpace(region.getMinCoeff(), region.getMaxCoeff());
			numDecompositions++;

			try {
				// Start iterations
				iters = 1;
				totalIters++;
				while (iters <= right) {
					// Vector-matrix multiply
					ctmcRanged.vmMult(solnMin, soln2Min, solnMax, soln2Max, q);

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

							// Check whether the minimised/maximised probs are accurate enough
							if ((sumMax[i] - sumMin[i]) > accuracy) {
								throw new SignificantInaccuracy();
							}
						}
					}

					iters++;
					totalIters++;
				}

				// Store result
				regionValues.put(region, sumMin, sumMax);
			} catch (SignificantInaccuracy e) {
				// Decompose the current region!
				regions.add(region.lowerHalf());
				regions.add(region.upperHalf());
			}
		}

		// Finished bounded probabilistic reachability
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Transient probability computation");
		mainLog.print(" took " + totalIters + " iters, " + numDecompositions + " decompositions");
		mainLog.print(" (resulting in " + regionValues.getNumRegions() + " final regions)");
		mainLog.print(" and " + timer / 1000.0 + " seconds in total.\n");

		return regionValues;
	}
}
