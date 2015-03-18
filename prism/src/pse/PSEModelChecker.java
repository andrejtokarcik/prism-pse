//==============================================================================
//	
//	Copyright (c) 2014-
//	Authors:
//	* Andrej Tokarcik <andrejtokarcik@gmail.com> (Masaryk University)
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

import java.io.File;
import java.util.BitSet;
import java.util.Map.Entry;

import parser.Values;
import parser.ast.Expression;
import parser.ast.ExpressionFilter;
import parser.ast.ExpressionFilter.FilterOperator;
import parser.ast.ExpressionLabel;
import parser.ast.ExpressionProb;
import parser.ast.ExpressionReward;
import parser.ast.ExpressionTemporal;
import parser.ast.ExpressionUnaryOp;
import parser.ast.ModulesFile;
import parser.ast.PropertiesFile;
import parser.ast.RelOp;
import parser.ast.RewardStruct;
import parser.type.TypeBool;
import parser.type.TypeDouble;
import prism.PrismComponent;
import prism.PrismException;
import prism.Result;
import prism.PrismSettings;
import explicit.FoxGlynn;
import explicit.StateModelChecker;
import explicit.StateValues;
import explicit.Utils;
import explicit.rewards.ConstructRewards;
import explicit.rewards.MCRewards;

/**
 * Model checker for {@link PSEModel}.
 * Covering such functionality as forwards/backwards transient
 * computations and checking of unnested CSL formulae.
 * The results are examined for accuracy by invocations
 * of passed-in decomposition procedures.
 */
public final class PSEModelChecker extends PrismComponent
{
	/** modules file, for instantiating of the PSE model */
	private ModulesFile modulesFile = null;

	/** properties file (for labels, constants, etc.) */
	private PropertiesFile propertiesFile = null;

	/** constants (extracted from modules & properties) */
	private Values constantValues;

	/** state model checker, used as fall-back for checking
	 *  of those expressions that PSE isn't concerned with */
	private StateModelChecker stateChecker;

	/**
	 * Constructor.
	 */
	public PSEModelChecker(PrismComponent parent)
	{
		super(parent);
		this.stateChecker = new StateModelChecker(this);
	}

	// Setters/getters

	/**
	 * Sets the attached model file (for e.g. reward structures when model checking)
	 * and the attached properties file (for e.g. constants/labels when model checking)
	 */
	public void setModulesFileAndPropertiesFile(ModulesFile modulesFile, PropertiesFile propertiesFile)
	{
		this.modulesFile = modulesFile;
		this.propertiesFile = propertiesFile;
		// Get combined constant values from model/properties
		constantValues = new Values();
		constantValues.addValues(modulesFile.getConstantValues());
		if (propertiesFile != null) {
			constantValues.addValues(propertiesFile.getConstantValues());
		}

		stateChecker.setModulesFileAndPropertiesFile(modulesFile, propertiesFile);
	}

	public ModulesFile getModulesFile()
	{
		return modulesFile;
	}

	public PropertiesFile getPropertiesFile()
	{
		return propertiesFile;
	}

	public Values getConstantValues()
	{
		return constantValues;
	}

	// Model checking functions, including backwards transient analysis

	/**
	 * Prepares for and performs model checking of the given property expression.
	 * Prints and returns the result.
	 * 
	 * @param model model to check
	 * @param expr property to check
	 * @param decompositionProcedure decomposition procedure to verify accuracy
	 * of results
	 * @return very final result of model checking
	 * @throws PrismException
	 */
	public Result check(PSEModel model, Expression expr, DecompositionProcedure decompositionProcedure)
			throws PrismException
	{
		long timer = 0;

		// Remove labels from property, using combined label list (on a copy of the expression)
		// This is done now so that we can handle labels nested below operators that are not
		// handled natively by the model checker yet (just evaluate()ed in a loop).
		expr = (Expression) expr.deepCopy().expandLabels(propertiesFile.getCombinedLabelList());

		// Also evaluate/replace any constants
		//expr = (Expression) expr.replaceConstants(constantValues);

		// Initialise decomposition procedure before model checking starts
		decompositionProcedure.initialiseModelChecking(this, model, expr);

		// Let decomposition procedure adjust the property, if needed
		expr = decompositionProcedure.getPropertyExpression();

		// Do model checking and store result vector
		timer = System.currentTimeMillis();
		BoxRegionValues regionValues = checkExpression(model, expr, decompositionProcedure);
		timer = System.currentTimeMillis() - timer;
		mainLog.println("\nTotal time for model checking: " + timer / 1000.0 + " seconds.");

		// Print and return result
		decompositionProcedure.printSolution(mainLog, settings.getBoolean(PrismSettings.PRISM_VERBOSE));
		Result result = new Result();
		result.setResult(regionValues);
		return result;
	}

	/**
	 * Performs model checking of the given property expression.
	 * 
	 * @param model model to check
	 * @param expr property to check
	 * @param decompositionProcedure decomposition procedure to verify accuracy
	 * of results
	 * @return minimised & maximised values per state and region, as returned
	 * by model checking sub-procedures
	 * @throws PrismException
	 */
	public BoxRegionValues checkExpression(PSEModel model, Expression expr, DecompositionProcedure decompositionProcedure)
			throws PrismException
	{
		if (expr instanceof ExpressionFilter) {
			return checkExpressionFilter(model, (ExpressionFilter) expr, decompositionProcedure);
		}
		if (expr instanceof ExpressionProb) {
			return checkExpressionProb(model, (ExpressionProb) expr, decompositionProcedure);
		}
		if (expr instanceof ExpressionReward) {
			return checkExpressionReward(model, (ExpressionReward) expr, decompositionProcedure);
		}
		// Let explicit.StateModelChecker take care of other kinds of expressions
		StateValues vals = stateChecker.checkExpression(model, expr);
		return new BoxRegionValues(model, model.getCompleteSpace(), vals, vals.deepCopy());
	}

	/**
	 * Checks a filter.
	 */
	protected BoxRegionValues checkExpressionFilter(PSEModel model, ExpressionFilter expr, DecompositionProcedure decompositionProcedure)
			throws PrismException
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
		BitSet bsFilterMin = filterRgnVals.getMin(model.getCompleteSpace()).getBitSet();
		BitSet bsFilterMax = filterRgnVals.getMax(model.getCompleteSpace()).getBitSet();
		boolean filterInit = (filter instanceof ExpressionLabel && ((ExpressionLabel) filter).getName().equals("init"));
		boolean filterInitSingle = filterInit & model.getNumInitialStates() == 1;
		if (bsFilterMin.isEmpty()) {
			throw new PrismException("Filter satisfies no states");
		}
		if (!filterInit) {
			mainLog.println("\nNumber of states satisfying filter " + filter + ":");
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
					mainLog.println("\n=== Minimised ===\n");
					subRgnValsMin.printFiltered(mainLog, bsFilterMin);
					mainLog.println("\n=== Maximised ===\n");
					subRgnValsMax.printFiltered(mainLog, bsFilterMax);
				} else {
					if (op == FilterOperator.PRINT) {
						mainLog.println("\nResults (non-zero only) for filter " + filter + ":");
						mainLog.println("\n=== Minimised ===\n");
						subRgnValsMin.printFiltered(mainLog, bsFilterMin);
						mainLog.println("\n=== Maximised ===\n");
						subRgnValsMax.printFiltered(mainLog, bsFilterMax);
					} else {
						mainLog.println("\nResults (including zeros) for filter " + filter + ":");
						mainLog.println("\n=== Minimised ===\n");
						subRgnValsMin.printFiltered(mainLog, bsFilterMin, false, false, true, true);
						mainLog.println("\n=== Maximised ===\n");
						subRgnValsMax.printFiltered(mainLog, bsFilterMax, false, false, true, true);
					}
				}
				resRgnVals = subRgnVals;
				break;
			case MIN:
			case MAX:
			case ARGMIN:
			case ARGMAX:
				throw new PrismException("Operation not implemented for PSE models");
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
				throw new PrismException("Operation not implemented for PSE models");
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
	 * Checks a P operator expression and returns the values for all states.
	 */
	protected BoxRegionValues checkExpressionProb(PSEModel model, ExpressionProb expr, DecompositionProcedure decompositionProcedure) throws PrismException
	{
		// Get info from prob operator
		RelOp relOp = expr.getRelOp(); // Relational operator
		Expression pb = expr.getProb(); // Probability bound (expression)
		double p = 0;  // Probability bound (actual value)
		if (pb != null) {
			p = pb.evaluateDouble(constantValues);
			if (p < 0 || p > 1)
				throw new PrismException("Invalid probability bound " + p + " in P operator");
		}

		// Compute probabilities
		BoxRegionValues regionValues = checkProbPathFormula(model, expr.getExpression(), decompositionProcedure);

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
	 * Computes probabilities for the contents of a P operator.
	 */
	protected BoxRegionValues checkProbPathFormula(PSEModel model, Expression expr, DecompositionProcedure decompositionProcedure) throws PrismException
	{
		// Test whether this is a simple path formula (i.e. PCTL)
		// and then pass control to appropriate method.
		if (expr.isSimplePathFormula()) {
			return checkProbPathFormulaSimple(model, expr, decompositionProcedure);
		} else {
			throw new PrismException("PSE does not support nested CSL formulae");
		}
	}

	/**
	 * Computes probabilities for a simple, non-LTL path operator.
	 */
	protected BoxRegionValues checkProbPathFormulaSimple(PSEModel model, Expression expr, DecompositionProcedure decompositionProcedure)
			throws PrismException
	{
		return checkProbPathFormulaSimple(model, expr, decompositionProcedure, false);
	}

	/**
	 * Computes probabilities for a simple, non-LTL path operator,
	 * possibly with need to negate the results in a sub-procedure
	 * (as indicated by {@code negate}).
	 */
	protected BoxRegionValues checkProbPathFormulaSimple(PSEModel model, Expression expr, DecompositionProcedure decompositionProcedure, boolean negate)
			throws PrismException
	{
		BoxRegionValues regionValues = null;

		// Negation/parentheses
		if (expr instanceof ExpressionUnaryOp) {
			ExpressionUnaryOp exprUnary = (ExpressionUnaryOp) expr;
			// Parentheses
			if (exprUnary.getOperator() == ExpressionUnaryOp.PARENTH) {
				// Recurse
				regionValues = checkProbPathFormulaSimple(model, exprUnary.getOperand(), decompositionProcedure, negate);
			}
			// Negation
			// NB: Handling of negations depends on the formula's being unnested/simple,
			// viz. on the fact that the only top-level operators beyond the one temporal
			// operator can be negations and/or parentheses.
			else if (exprUnary.getOperator() == ExpressionUnaryOp.NOT) {
				regionValues = checkProbPathFormulaSimple(model, exprUnary.getOperand(), decompositionProcedure, !negate);
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
					regionValues = checkProbBoundedUntil(model, exprTemp, decompositionProcedure, negate);
				} else {
					throw new PrismException("PSE supports bounded U operator only");
				}
			}
			// Anything else - convert to until and recurse
			else {
				regionValues = checkProbPathFormulaSimple(model, exprTemp.convertToUntilForm(), decompositionProcedure, negate);
			}
		}

		if (regionValues == null)
			throw new PrismException("Unrecognised path operator in P operator");

		return regionValues;
	}

	protected BoxRegionValues checkProbBoundedUntil(PSEModel model, ExpressionTemporal expr, DecompositionProcedure decompositionProcedure)
			throws PrismException
	{
		return checkProbBoundedUntil(model, expr, decompositionProcedure, false);
	}

	/**
	 * Checks a time-bounded until operator, possibly with need
	 * to negate the results in a sub-procedure (as indicated
	 * by {@code negate}). 
	 */
	protected BoxRegionValues checkProbBoundedUntil(PSEModel model, ExpressionTemporal expr, DecompositionProcedure decompositionProcedure, boolean negate)
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
		b1Min = op1RgnVals.getMin(model.getCompleteSpace()).getBitSet();
		b1Max = op1RgnVals.getMax(model.getCompleteSpace()).getBitSet();
		b2Min = op2RgnVals.getMin(model.getCompleteSpace()).getBitSet();
		b2Max = op2RgnVals.getMax(model.getCompleteSpace()).getBitSet();

		// compute probabilities

		// a trivial case: "U<=0"
		if (lTime == 0 && uTime == 0) {
			// prob is 1 in b2 states, 0 otherwise
			probsMin = StateValues.createFromBitSetAsDoubles(b2Min, model);
			probsMax = StateValues.createFromBitSetAsDoubles(b2Max, model);
			regionValues = new BoxRegionValues(model, model.getCompleteSpace(), probsMin, probsMax);
		} else {
			BoxRegionValues onesMultProbs = BoxRegionValues.createWithAllOnes(model, model.getCompleteSpace());

			// >= lTime
			if (uTime == -1) {
				throw new PrismException("PSE requires until formula time-bounded from above");
			}
			// <= uTime
			else if (lTime == 0) {
				// nb: uTime != 0 since would be caught above (trivial case)
				b1Min.andNot(b2Min);
				b1Max.andNot(b2Max);

				while (true) {
					try {
						regionValues = computeTransientBackwardsProbs(
								model, b2Min, b1Min, b2Max, b1Max, uTime, onesMultProbs,
								decompositionProcedure, oldRegionValues, negate);
						break;
					} catch (DecompositionProcedure.DecompositionNeeded e) {
						e.printRegionsToDecompose(mainLog);
						for (BoxRegion region : e.getRegionsToDecompose()) {
							onesMultProbs.decomposeRegion(region);
						}
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
						// DecompositionNeeded is always thrown within the second transient computation,
						// since the first transient computation is executed with decomposing disabled.
						tmpRegionValues = computeTransientBackwardsProbs(
								model, b2Min, tmpMin, b2Max, tmpMax, uTime - lTime, onesMultProbs,
								SimpleDecompositionProcedure.NoDecomposing.getInstance(), oldTmpRegionValues, false);
						regionValues = computeTransientBackwardsProbs(
								model, b1Min, b1Min, b1Max, b1Max, lTime, tmpRegionValues,
								decompositionProcedure, oldRegionValues, negate);
						break;
					} catch (DecompositionProcedure.DecompositionNeeded e) {
						e.printRegionsToDecompose(mainLog);
						for (BoxRegion region : e.getRegionsToDecompose()) {
							onesMultProbs.decomposeRegion(region);
						}
						oldTmpRegionValues = tmpRegionValues;
						oldRegionValues = e.getExaminedRegionValues();
					}
				}
			}
		}

		mainLog.print("\nThe PSE backwards transient probability computations altogether produced ");
		mainLog.println(regionValues.getNumRegions() + " final regions.");
		return regionValues;
	}

	public BoxRegionValues computeTransientBackwardsProbs(PSEModel model,
			BitSet targetMin, BitSet nonAbsMin, BitSet targetMax, BitSet nonAbsMax,
			double t, BoxRegionValues multProbs, DecompositionProcedure decompositionProcedure)
					throws PrismException, DecompositionProcedure.DecompositionNeeded
	{
		return computeTransientBackwardsProbs(model, targetMin, nonAbsMin, targetMax, nonAbsMax, t, multProbs, decompositionProcedure, null, false);
	}

	/**
	 * Performs backwards transient probability computation as required for CSL
	 * model checking.
	 * Computes, for each state, the sum over {@code target} states of the probability
	 * of being in that state at time {@code t} multiplied by the corresponding
	 * probability in the vector {@code multProbs}, assuming that all states
	 * <i>not</i> in {@code nonAbs} are made absorbing. This holds for both
	 * {@code *Min} and {@code *Max} variant, as they are in fact computed
	 * in parallel independently of each other.
	 * <p>
	 * NB: Decompositions of the parameter space must be performed explicitly,
	 * {@code DecompositionNeeded} is not handled within the method.
	 * 
	 * @param model model to check
	 * @param targetMin target states when producing minimised probabilities
	 * @param nonAbsMin states <i>not</i> to be made absorbing when producing
	 * minimised probabilities (null means "all")
	 * @param targetMax target states when producing maximised probabilities
	 * @param nonAbsMax  states <i>not</i> to be made absorbing when producing
	 * maximised probabilities (null means "all")
	 * @param t time bound
	 * @param multProbs vector containing minimised & maximised multiplication
	 * factors per region
	 * @param decompositionProcedure decomposition procedure to verify accuracy
	 * of results
	 * @param previousResult previous result, from which finished regions are to
	 * be re-used
	 * @param negate if true, results are subtracted from 1.0 and swapped before
	 * returning
	 * @return minimised & maximised probabilities per state and region
	 * @throws PrismException
	 * @throws DecompositionProcedure.DecompositionNeeded thrown by decomposition
	 * procedure if it found the results to be too inaccurate
	 * @see PSEModel#mvMult
	 */
	public BoxRegionValues computeTransientBackwardsProbs(PSEModel model,
			BitSet targetMin, BitSet nonAbsMin, BitSet targetMax, BitSet nonAbsMax,
			double t, BoxRegionValues multProbs, DecompositionProcedure decompositionProcedure,
			BoxRegionValues previousResult, boolean negate)
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
		if (previousResult == null) {
			previousResult = new BoxRegionValues(model);
		}

		// Store num states
		n = model.getNumStates();

		// Optimisations: If (nonAbs is empty or t = 0) and multProbs is all ones, this is easy.
		if ((nonAbs != null && nonAbs.isEmpty()) || (t == 0)) {
			if (multProbs.isAllOnes()) {
				for (BoxRegion region : multProbs.keySet()) {
					if (previousResult.hasRegion(region)) {
						regionValues.put(region, previousResult.getMin(region), previousResult.getMax(region));
						continue;
					}
					solnMin = Utils.bitsetToDoubleArray(targetMin, n);
					solnMax = Utils.bitsetToDoubleArray(targetMax, n);
					regionValues.put(region, solnMin, solnMax);
				}
				return regionValues;
			}
		}

		// Start backwards transient computation
		timer = System.currentTimeMillis();
		mainLog.println("\nStarting PSE backwards transient probability computation...");

		// Compute the in, out, inout sets of reactions
		model.computeInOutTransitions();

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
		for (i = left; i <= right; i++) {
			weights[i - left] /= totalWeight;
		}
		mainLog.println("Fox-Glynn (" + acc + "): left = " + left + ", right = " + right);
		mainLog.println();

		// Get number of iterations for partial examination
		int numItersExaminePartial = settings.getInteger(PrismSettings.PRISM_PSE_EXAMINEPARTIAL);

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
			model.evaluateParameters(region);
			model.prepareForMV(nonAbs, false);
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
				model.mvMult(solnMin, soln2Min, solnMax, soln2Max);

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
					// After a number of iters (default 50), examine the partially computed result
					if (iters % numItersExaminePartial == 0) {
						decompositionProcedure.examinePartialComputation(regionValues, region, sumMin, sumMax);
					}
				}

				iters++;
				totalIters++;
			}

			// Examine this region's result after all the iters have been finished
			decompositionProcedure.examinePartialComputation(regionValues, region, sumMin, sumMax);

			// Store result
			regionValues.put(region, sumMin, sumMax);
		}

		// Negate if necessary
		if (negate) {
			// Subtract all min/max values from 1 and swap
			for (BoxRegionValues.StateValuesPair pair: regionValues.values()) {
				pair.getMin().timesConstant(-1.0);
				pair.getMin().plusConstant(1.0);
				pair.getMax().timesConstant(-1.0);
				pair.getMax().plusConstant(1.0);
				pair.swap();
			}
		}

		// Examine the whole computation after it's completely finished
		decompositionProcedure.examineWholeComputation(regionValues);

		// Finished bounded probabilistic reachability
		timer = System.currentTimeMillis() - timer;
		mainLog.print("\nPSE backwards transient probability computation");
		mainLog.print(" took " + totalIters + " iters");
		mainLog.println(" and " + timer / 1000.0 + " seconds.");

		return regionValues;
	}

	/**
	 * Model check an R operator expression and return the values for all states.
	 */
	protected BoxRegionValues checkExpressionReward(PSEModel model, ExpressionReward expr, DecompositionProcedure decompositionProcedure)
			throws PrismException
	{
		// Get info from R operator
		RewardStruct rewStruct = expr.getRewardStructByIndexObject(modulesFile, constantValues);
		RelOp relOp = expr.getRelOp(); // Relational operator
		Expression rb = expr.getReward(); // Reward bound (expression)
		double r = 0; // Reward bound (actual value)
		if (rb != null) {
			r = rb.evaluateDouble(constantValues);
			if (r < 0)
				throw new PrismException("Invalid reward bound " + r + " in R[] formula");
		}

		// Build rewards
		mainLog.println("Building reward structure...");
		ConstructRewards constructRewards = new ConstructRewards(mainLog);
		MCRewards rewards = constructRewards.buildMCRewardStructure(model, rewStruct, constantValues);

		// Compute rewards
		BoxRegionValues regionValues = checkRewardFormula(model, rewards, expr.getExpression(), decompositionProcedure);

		// For =? properties, just return values
		if (rb == null) {
			return regionValues;
		}
		// Otherwise, compare against bound to get set of satisfying states
		else {
			for (Entry<BoxRegion, BoxRegionValues.StateValuesPair> entry : regionValues) {
				BitSet solnMin = entry.getValue().getMin().getBitSetFromInterval(relOp, r);
				BitSet solnMax = entry.getValue().getMax().getBitSetFromInterval(relOp, r);
				regionValues.put(entry.getKey(), solnMin, solnMax);
			}
			return regionValues;
		}
	}

	/**
	 * Compute rewards for the contents of an R operator.
	 */
	protected BoxRegionValues checkRewardFormula(PSEModel model, MCRewards modelRewards, Expression expr, DecompositionProcedure decompositionProcedure)
			throws PrismException
	{
		assert expr instanceof ExpressionTemporal : "Unrecognized operator in R operator";

		ExpressionTemporal exprTemp = (ExpressionTemporal) expr;
		switch (exprTemp.getOperator()) {
		case ExpressionTemporal.R_C:
			return checkRewardCumulative(model, modelRewards, exprTemp, decompositionProcedure);
		default:
			throw new PrismException("PSE does not yet handle the " + exprTemp.getOperatorSymbol() + " reward operator");
		}
	}

	/**
	 * Compute rewards for a cumulative reward operator.
	 */
	protected BoxRegionValues checkRewardCumulative(PSEModel model, MCRewards modelRewards, ExpressionTemporal expr, DecompositionProcedure decompositionProcedure)
			throws PrismException
	{
		// Check that there is an upper time bound
		if (expr.getUpperBound() == null) {
			throw new PrismException("Cumulative reward operator without time bound (C) is only allowed for multi-objective queries");
		}

		// Get time bound
		assert model.getModelType().continuousTime();
		double time = expr.getUpperBound().evaluateDouble(constantValues);
		if (time < 0) {
			throw new PrismException("Invalid time bound " + time + " in cumulative reward formula");
		}

		BoxRegionValues oldRegionValues = null;
		BoxRegionValues onesMultProbs = BoxRegionValues.createWithAllOnes(model, model.getCompleteSpace());
		while (true) {
			try {
				return computeCumulativeRewards(
						model, modelRewards, time, onesMultProbs,
						decompositionProcedure, oldRegionValues);
			} catch (DecompositionProcedure.DecompositionNeeded e) {
				e.printRegionsToDecompose(mainLog);
				for (BoxRegion region : e.getRegionsToDecompose()) {
					onesMultProbs.decomposeRegion(region);
				}
				oldRegionValues = e.getExaminedRegionValues();
			}
		}
	}

	/**
	 */
	public BoxRegionValues computeCumulativeRewards(PSEModel model,
			MCRewards mcRewards, double t, BoxRegionValues multProbs,
			DecompositionProcedure decompositionProcedure, BoxRegionValues previousResult)
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

		if (previousResult == null) {
			previousResult = new BoxRegionValues(model);
		}

		// Store num states
		n = model.getNumStates();

		// Optimisations: If t = 0, this is easy.
		if (t == 0) {
			double[] zeros = new double[n];
			for (BoxRegion region : multProbs.keySet()) {
				if (previousResult.hasRegion(region)) {
					regionValues.put(region, previousResult.getMin(region), previousResult.getMax(region));
					continue;
				}
				regionValues.put(region, zeros, zeros);
			}
			return regionValues;
		}

		// Start backwards transient computation
		timer = System.currentTimeMillis();
		mainLog.println("\nStarting PSE backwards cumulative rewards computation...");

		// Compute the in, out, inout sets of reactions
		model.computeInOutTransitions();

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
		for (i = left; i <= right; i++) {
			weights[i - left] /= totalWeight;
		}
		
		// modify the poisson probabilities to what we need for this computation
		// first make the kth value equal to the sum of the values for 0...k
		for (i = left+1; i <= right; i++) {
			weights[i - left] += weights[i - 1 - left];
		}
		// then subtract from 1 and divide by uniformisation constant (q) to give mixed poisson probabilities
		for (i = left; i <= right; i++) {
			weights[i - left] = (1 - weights[i - left]) / q;
		}

		mainLog.println("Fox-Glynn (" + acc + "): left = " + left + ", right = " + right);
		mainLog.println();

		// Get number of iterations for partial examination
		int numItersExaminePartial = settings.getInteger(PrismSettings.PRISM_PSE_EXAMINEPARTIAL);

		totalIters = 0;
		for (Entry<BoxRegion, BoxRegionValues.StateValuesPair> entry : multProbs) {
			BoxRegion region = entry.getKey();

			// If the previous region values contain probs for this region, i.e. the region
			// has not been decomposed, then just use the previous result directly.
			if (previousResult.hasRegion(region)) {
				regionValues.put(region, previousResult.getMin(region), previousResult.getMax(region));
				continue;
			}

			/*
			double[] multProbsMin = entry.getValue().getMin().getDoubleArray();
			double[] multProbsMax = entry.getValue().getMax().getDoubleArray();
			*/

			// Configure parameter space
			model.evaluateParameters(region);
			model.prepareForMV(null, false);
			mainLog.println("Computing probabilities for parameter region " + region);

			// Create solution vectors
			solnMin = new double[n];
			soln2Min = new double[n];
			sumMin = new double[n];
			solnMax = new double[n];
			soln2Max = new double[n];
			sumMax = new double[n];

			// Initialise solution vectors.
			for (i = 0; i < n; i++) {
				solnMin[i] = mcRewards.getStateReward(i);
				solnMax[i] = mcRewards.getStateReward(i);
			}

			// do 0th element of summation (doesn't require any matrix powers)
			if (left == 0) {
				for (i = 0; i < n; i++) {
					sumMin[i] += weights[0] * solnMin[i];
					sumMax[i] += weights[0] * solnMax[i];
				}
			} else {
				for (i = 0; i < n; i++) {
					sumMin[i] += solnMin[i] / q;
					sumMax[i] += solnMax[i] / q;
				}
			}

			// Start iterations
			iters = 1;
			while (iters <= right) {
				// Matrix-vector multiply				
				model.mvMult(solnMin, soln2Min, solnMax, soln2Max);

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
				} else {
					for (i = 0; i < n; i++) {
						sumMin[i] += solnMin[i] / q;
						sumMax[i] += solnMax[i] / q;
					}
				}
				
				// After a number of iters (default 50), examine the partially computed result
				if (iters % numItersExaminePartial == 0) {
					decompositionProcedure.examinePartialComputation(regionValues, region, sumMin, sumMax);
				}

				iters++;
				totalIters++;
			}

			// Examine this region's result after all the iters have been finished
			decompositionProcedure.examinePartialComputation(regionValues, region, sumMin, sumMax);

			// Store result
			regionValues.put(region, sumMin, sumMax);
		}

		/*
		// Negate if necessary
		if (negate) {
			// Subtract all min/max values from 1 and swap
			for (BoxRegionValues.StateValuesPair pair: regionValues.values()) {
				pair.getMin().timesConstant(-1.0);
				pair.getMin().plusConstant(1.0);
				pair.getMax().timesConstant(-1.0);
				pair.getMax().plusConstant(1.0);
				pair.swap();
			}
		}
		*/

		// Examine the whole computation after it's completely finished
		decompositionProcedure.examineWholeComputation(regionValues);

		// Finished bounded probabilistic reachability
		timer = System.currentTimeMillis() - timer;
		mainLog.print("\nPSE backwards transient cumulative rewards computation");
		mainLog.print(" took " + totalIters + " iters");
		mainLog.println(" and " + timer / 1000.0 + " seconds.");

		return regionValues;
	}

	// Transient analysis (forwards)

	/**
	 * Computes transient probability distribution (forwards).
	 * Optionally, uses the passed in region values {@code initDist}
	 * as the initial probability distribution (time 0).
	 * If the initial distribution is null, starts from initial state
	 * (or uniform distribution over multiple initial states).
	 */
	public BoxRegionValues doTransient(PSEModel model, double t, BoxRegionValues initDist, DecompositionProcedure decompositionProcedure)
			throws PrismException
	{
		// Build initial distribution (if not specified)
		if (initDist == null) {
			initDist = buildInitialDistribution(model);
		}

		// Compute transient probabilities
		BoxRegionValues regionValues = null;
		BoxRegionValues previousResult = null;
		while (true) {
			try {
				regionValues = computeTransientProbs(model, t, initDist, decompositionProcedure, previousResult);
				break;
			} catch (DecompositionProcedure.DecompositionNeeded e) {
				e.printRegionsToDecompose(mainLog);
				for (BoxRegion region : e.getRegionsToDecompose()) {
					initDist.decomposeRegion(region);
				}
				previousResult = e.getExaminedRegionValues(); 
			}
		}

		return regionValues;
	}

	public BoxRegionValues computeTransientProbs(PSEModel model, double t, BoxRegionValues initDist, DecompositionProcedure decompositionProcedure)
			throws PrismException, DecompositionProcedure.DecompositionNeeded
	{
		return computeTransientProbs(model, t, initDist, decompositionProcedure, null);
	}

	/**
	 * Performs forwards transient probability computation.
	 * Computes the minimised & maximised probability of being in each state
	 * at time {@code t}, assuming the initial distribution {@code initDist}.
	 * <p>
	 * NB: Decompositions of the parameter space must be performed explicitly,
	 * {@code DecompositionNeeded} is not handled within the method.
	 * 
	 * @param model model to be analysed
	 * @param t time point
	 * @param initDist initial distribution
	 * @param decompositionProcedure decomposition procedure to verify accuracy
	 * of results
	 * @return minimised & maximised probabilities per state and region
	 * @throws PrismException
	 * @see PSEModel#vmMult
	 */
	public BoxRegionValues computeTransientProbs(PSEModel model, double t, BoxRegionValues initDist, DecompositionProcedure decompositionProcedure, BoxRegionValues previousResult)
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

		if (previousResult == null) {
			previousResult = new BoxRegionValues(model);
		}

		// Start bounded probabilistic reachability
		timer = System.currentTimeMillis();
		mainLog.println("\nStarting PSE transient probability computation...");

		// Compute the in, out, inout sets of reactions
		model.computeInOutTransitions();

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

		// Get number of iterations for partial examination
		int numItersExaminePartial = settings.getInteger(PrismSettings.PRISM_PSE_EXAMINEPARTIAL);

		totalIters = 0;
		for (Entry<BoxRegion, BoxRegionValues.StateValuesPair> entry : initDist) {
			BoxRegion region = entry.getKey();
			
			// If the previous region values contain probs for this region, i.e. the region
			// has not been decomposed, then just use the previous result directly.
			if (previousResult.hasRegion(region)) {
				regionValues.put(region, previousResult.getMin(region), previousResult.getMax(region));
				continue;
			}

			// Create and initialise solution vectors.
			// Don't need to do soln2 since will be immediately overwritten.
			// Vector sum is all zeros (done by array creation).
			solnMin = entry.getValue().getMin().getDoubleArray().clone();
			soln2Min = new double[n];
			sumMin = new double[n];
			solnMax = entry.getValue().getMax().getDoubleArray().clone();
			soln2Max = new double[n];
			sumMax = new double[n];

			// If necessary, do 0th element of summation (doesn't require any matrix powers)
			if (left == 0) {
				for (i = 0; i < n; i++) {
					sumMin[i] += weights[0] * solnMin[i];
					sumMax[i] += weights[0] * solnMax[i];
				}
			}

			// Configure parameter space
			model.evaluateParameters(region);
			model.prepareForVM();
			mainLog.println("Computing probabilities for parameter region " + region);

			// Start iterations
			iters = 1;
			totalIters++;
			while (iters <= right) {
				// Vector-matrix multiply
				model.vmMult(solnMin, soln2Min, solnMax, soln2Max);

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
					// After a number of iters (default 50), examine the partially computed result
					if (iters % numItersExaminePartial == 0) {
						decompositionProcedure.examinePartialComputation(regionValues, region, sumMin, sumMax);
					}
				}

				iters++;
				totalIters++;
			}

			// Examine this region's result after all the iters have been finished
			decompositionProcedure.examinePartialComputation(regionValues, region, sumMin, sumMax);

			// Store result
			regionValues.put(region, sumMin, sumMax);
		}

		// Examine the whole computation after it's completely finished
		decompositionProcedure.examineWholeComputation(regionValues);

		// Finished bounded probabilistic reachability
		timer = System.currentTimeMillis() - timer;
		mainLog.print("\nPSE transient probability computation");
		mainLog.print(" took " + totalIters + " iters");
		mainLog.print(" and " + timer / 1000.0 + " seconds");
		mainLog.println(" (producing " + regionValues.getNumRegions() + " final regions).");

		return regionValues;
	}
	
	// Utility methods

	/**
	 * Read a probability distribution, stored as a {@code BoxRegionValues} object,
	 * from a file, where the only initial parameter region is the complete parameter
	 * space.
	 * If {@code distFile} is null, so is the return value.
	 */
	public BoxRegionValues readDistributionFromFile(File distFile, PSEModel model) throws PrismException
	{
		if (distFile == null) {
			return null;
		}

		mainLog.println("\nImporting probability distribution from file \"" + distFile + "\"...");
		BoxRegionValues regionValues = new BoxRegionValues(model);
		regionValues.readFromFile(distFile);
		return regionValues;
	}

	/**
	 * Build a probability distribution, stored as a {@code BoxRegionValues} object,
	 * from the initial states info of the current model, with the only initial
	 * parameter region being the complete parameter space: state values have either
	 * probability 1 for the (single) initial state or equiprobable over multiple
	 * initial states.
	 */
	public BoxRegionValues buildInitialDistribution(PSEModel model) throws PrismException
	{
		StateValues uniformInit = new StateValues(TypeDouble.getInstance(), new Double(0.0), model);
		double initVal = 1.0 / model.getNumInitialStates();
		for (int in : model.getInitialStates()) {
			uniformInit.setDoubleValue(in, initVal);
		}
		return new BoxRegionValues(model, model.getCompleteSpace(), uniformInit, uniformInit.deepCopy());
	}
}
