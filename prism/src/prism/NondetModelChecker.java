//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford, formerly University of Birmingham)
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

package prism;

import jdd.*;
import odd.*;
import dv.*;
import mtbdd.*;
import sparse.*;
import hybrid.*;
import parser.*;
import parser.ast.*;

// Model checker for MDPs

public class NondetModelChecker implements ModelChecker
{
	// logs
	private PrismLog mainLog;
	private PrismLog techLog;

	// options

	// which engine to use
	private int engine;
	// termination criterion (iterative methods)
	private int termCrit;
	// parameter for termination criterion
	private double termCritParam;
	// max num iterations (iterative methods)
	private int maxIters;
	// flags
	private boolean verbose; // verbose output?
	private boolean precomp; // use 0,1 precomputation algorithms?
	private boolean fairness; // use fairness?
	// sparse bits info
	private int SBMaxMem;
	private int numSBLevels;

	// properties file
	private PropertiesFile propertiesFile;

	// constant values
	private Values constantValues;

	// Expression2MTBDD object for translating expressions
	private Expression2MTBDD expr2mtbdd;

	// class-wide storage for probability in the initial state
	private double numericalRes;

	// model info
	private NondetModel model;
	private JDDNode trans;
	private JDDNode trans01;
	private JDDNode start;
	private JDDNode reach;
	private ODDNode odd;
	private JDDNode nondetMask;
	private JDDVars allDDRowVars;
	private JDDVars allDDColVars;
	private JDDVars allDDNondetVars;
	private JDDVars[] varDDRowVars;

	// constructor - set some defaults

	public NondetModelChecker(PrismLog log1, PrismLog log2, Model m, PropertiesFile pf) throws PrismException
	{
		// initialise
		mainLog = log1;
		techLog = log2;
		if (!(m instanceof NondetModel)) {
			throw new PrismException("Error: wrong model type passed to NondetModelChecker.");
		}
		model = (NondetModel) m;
		propertiesFile = pf;
		trans = model.getTrans();
		trans01 = model.getTrans01();
		start = model.getStart();
		reach = model.getReach();
		odd = model.getODD();
		nondetMask = model.getNondetMask();
		allDDRowVars = model.getAllDDRowVars();
		allDDColVars = model.getAllDDColVars();
		allDDNondetVars = model.getAllDDNondetVars();
		varDDRowVars = model.getVarDDRowVars();

		// create list of all constant values needed
		constantValues = new Values();
		constantValues.addValues(model.getConstantValues());
		if (pf != null)
			constantValues.addValues(pf.getConstantValues());

		// create Expression2MTBDD object
		expr2mtbdd = new Expression2MTBDD(mainLog, techLog, model.getVarList(), varDDRowVars, constantValues);
		expr2mtbdd.setFilter(reach);

		// set up some default options
		// (although all should be overridden before model checking)
		engine = Prism.HYBRID;
		termCrit = Prism.RELATIVE;
		termCritParam = 1e-6;
		maxIters = 10000;
		verbose = false;
		precomp = true;
		fairness = true;
		SBMaxMem = 1024;
		numSBLevels = -1;
	}

	// set engine

	public void setEngine(int e)
	{
		engine = e;
	}

	// set options (generic)

	public void setOption(String option, boolean b)
	{
		if (option.equals("verbose")) {
			verbose = b;
		} else if (option.equals("precomp")) {
			precomp = b;
		} else if (option.equals("fairness")) {
			fairness = b;
		} else if (option.equals("compact")) {
			PrismSparse.setCompact(b);
			PrismHybrid.setCompact(b);
		} else {
			mainLog.println("Warning: option \"" + option + "\" not supported by NondetModelChecker.");
		}
	}

	public void setOption(String option, int i)
	{
		if (option.equals("termcrit")) {
			termCrit = i;
			PrismMTBDD.setTermCrit(i);
			PrismSparse.setTermCrit(i);
			PrismHybrid.setTermCrit(i);
		} else if (option.equals("maxiters")) {
			maxIters = i;
			PrismMTBDD.setMaxIters(i);
			PrismSparse.setMaxIters(i);
			PrismHybrid.setMaxIters(i);
		} else if (option.equals("sbmaxmem")) {
			SBMaxMem = i;
			PrismHybrid.setSBMaxMem(i);
		} else if (option.equals("numsblevels")) {
			numSBLevels = i;
			PrismHybrid.setNumSBLevels(i);
		} else {
			mainLog.println("Warning: option \"" + option + "\" not supported by NondetModelChecker.");
		}
	}

	public void setOption(String option, double d)
	{
		if (option.equals("termcritparam")) {
			termCritParam = d;
			PrismMTBDD.setTermCritParam(d);
			PrismSparse.setTermCritParam(d);
			PrismHybrid.setTermCritParam(d);
		} else {
			mainLog.println("Warning: option \"" + option + "\" not supported by NondetModelChecker.");
		}
	}

	public void setOption(String option, String s)
	{
		mainLog.println("Warning: option \"" + option + "\" not supported by NondetModelChecker.");
	}

	// compute number of places to round solution vector to
	// (given termination criteria, etc.)

	public int getPlacesToRoundBy()
	{
		return 1 - (int) (Math.log(termCritParam) / Math.log(10));
	}

	// -----------------------------------------------------------------------------------
	// Check a property, i.e. an expression
	// -----------------------------------------------------------------------------------

	public Object check(Expression expr) throws PrismException
	{
		long timer = 0;
		JDDNode dd;
		StateList states;
		boolean b;
		Object res;

		// start timer
		timer = System.currentTimeMillis();

		// do model checking and store result
		dd = checkExpression(expr);
		states = new StateListMTBDD(dd, model);

		// stop timer
		timer = System.currentTimeMillis() - timer;

		// print out model checking time
		mainLog.println("\nTime for model checking: " + timer / 1000.0 + " seconds.");

		// output and result of model checking depend on return type
		if (expr.getType() == Expression.BOOLEAN) {

			// print out number of satisfying states
			mainLog.print("\nNumber of satisfying states: ");
			mainLog.print(states.size());
			if (states.size() == model.getNumStates()) {
				mainLog.print(" (all)");
			} else if (states.includes(model.getStart())) {
				mainLog.print((model.getNumStartStates() == 1) ? " (including initial state)"
						: " (including all initial states)");
			} else {
				mainLog.print((model.getNumStartStates() == 1) ? " (initial state not satisfied)"
						: " (initial states not all satisfied)");
			}
			mainLog.print("\n");

			// if "verbose", print out satisfying states
			if (verbose) {
				mainLog.print("\nSatisfying states:");
				if (states.size() > 0) {
					mainLog.print("\n");
					states.print(mainLog);
				} else {
					mainLog.print(" (none)\n");
				}
			}

			// result is true if all states satisfy, false otherwise
			b = (states.size() == model.getNumStates());
			res = b ? new Boolean(true) : new Boolean(false);

			// print result
			mainLog.print("\nResult: " + b + " (property " + (b ? "" : "not ") + "satisfied in all states)\n");
		} else {
			// result is value stored earlier
			res = new Double(numericalRes);

			// print result
			mainLog.print("\nResult: " + res + "\n");
		}

		// finished with memory
		states.clear();

		// return result
		return res;
	}

	// Check expression (recursive)

	public JDDNode checkExpression(Expression expr) throws PrismException
	{
		JDDNode res;

		// P operator
		if (expr instanceof ExpressionProb) {
			res = checkExpressionProb((ExpressionProb) expr);
		}
		// R operator
		else if (expr instanceof ExpressionReward) {
			res = checkExpressionReward((ExpressionReward) expr);
		}
		// Label
		else if (expr instanceof ExpressionLabel) {
			res = checkExpressionLabel((ExpressionLabel) expr);
		}
		// Otherwise, must be an ordinary expression
		else {
			res = expr2mtbdd.translateExpression(expr);
		}

		// Filter out non-reachable states from solution
		JDD.Ref(reach);
		res = JDD.Apply(JDD.TIMES, res, reach);

		return res;
	}

	// -----------------------------------------------------------------------------------
	// Check method for each operator
	// -----------------------------------------------------------------------------------

	// P operator

	private JDDNode checkExpressionProb(ExpressionProb expr) throws PrismException
	{
		Expression pb; // probability bound (expression)
		double p = 0; // probability value (actual value)
		String relOp; // relational operator
		boolean min; // are we finding min (true) or max (false) probs
		PathExpression pe; // path expression

		JDDNode filter, sol, tmp;
		StateProbs probs = null;
		StateList states = null;
		double minRes = 0, maxRes = 0;

		// get info from prob operator
		relOp = expr.getRelOp();
		pb = expr.getProb();
		if (pb != null) {
			p = pb.evaluateDouble(constantValues, null);
			if (p < 0 || p > 1)
				throw new PrismException("Invalid probability bound " + p + " in P operator");
		}

		// check for trivial (i.e. stupid) cases
		if (pb != null) {
			if ((p == 0 && relOp.equals(">=")) || (p == 1 && relOp.equals("<="))) {
				mainLog.print("\nWarning: checking for probability " + relOp + " " + p
						+ " - formula trivially satisfies all states\n");
				JDD.Ref(reach);
				return reach;
			} else if ((p == 0 && relOp.equals("<")) || (p == 1 && relOp.equals(">"))) {
				mainLog.print("\nWarning: checking for probability " + relOp + " " + p
						+ " - formula trivially satisfies no states\n");
				return JDD.Constant(0);
			}
		}

		// determine whether min or max probabilities needed
		if (relOp.equals(">") || relOp.equals(">=") || relOp.equals("min=")) {
			// min
			min = true;
		} else if (relOp.equals("<") || relOp.equals("<=") || relOp.equals("max=")) {
			// max
			min = false;
		} else {
			throw new PrismException("Can't use \"P=?\" for MDPs; use \"Pmin=?\" or \"Pmax=?\"");
		}

		// translate filter (if present)
		filter = null;
		if (expr.getFilter() != null) {
			filter = checkExpression(expr.getFilter().getExpression());
		}

		// check if filter satisfies no states
		if (filter != null)
			if (filter.equals(JDD.ZERO)) {
				// for P=? properties, this is an error
				if (pb == null) {
					throw new PrismException("Filter {" + expr.getFilter().getExpression() + "} in P=?[] property satisfies no states");
				}
				// otherwise just print a warning
				else {
					mainLog.println("\nWarning: Filter {" + expr.getFilter().getExpression()
							+ "} satisfies no states and is being ignored");
					JDD.Deref(filter);
					filter = null;
				}
			}

		// compute probabilities
		pe = expr.getPathExpression();
		try {
			if (pe instanceof PathExpressionTemporal) {
				if (((PathExpressionTemporal) pe).hasBounds()) {
					switch (((PathExpressionTemporal) pe).getOperator()) {
					case PathExpressionTemporal.P_U:
						probs = checkProbBoundedUntil((PathExpressionTemporal) pe, min);
						break;
					case PathExpressionTemporal.P_F:
						probs = checkProbBoundedFuture((PathExpressionTemporal) pe, min);
						break;
					case PathExpressionTemporal.P_G:
						probs = checkProbBoundedGlobal((PathExpressionTemporal) pe, min);
						break;
					}
				} else {
					switch (((PathExpressionTemporal) pe).getOperator()) {
					case PathExpressionTemporal.P_X:
						probs = checkProbNext((PathExpressionTemporal) pe, min);
						break;
					case PathExpressionTemporal.P_U:
						probs = checkProbUntil((PathExpressionTemporal) pe, pb, p, min);
						break;
					case PathExpressionTemporal.P_F:
						probs = checkProbFuture((PathExpressionTemporal) pe, pb, p, min);
						break;
					case PathExpressionTemporal.P_G:
						probs = checkProbGlobal((PathExpressionTemporal) pe, pb, p, min);
						break;
					}
				}
			}
			if (probs == null)
				throw new PrismException("Unrecognised path operator in P operator");
		} catch (PrismException e) {
			if (filter != null)
				JDD.Deref(filter);
			throw e;
		}

		// round off probabilities
		// probs.roundOff(getPlacesToRoundBy());

		// print out probabilities
		if (verbose) {
			mainLog.print("\n" + (min ? "Minimum" : "Maximum") + " probabilities (non-zero only) for all states:\n");
			probs.print(mainLog);
		}

		// if a filter was provided, there is some additional output...
		if (filter != null) {

			// for non-"P=?"-type properties, print out some probs
			if (pb != null) {
				if (!expr.noFilterRequests()) {
					mainLog
							.println("\nWarning: \"{min}\", \"{max}\" only apply to \"P=?\" properties; they are ignored here.");
				}
				mainLog.print("\n" + (min ? "Minimum" : "Maximum")
						+ " probabilities (non-zero only) for states satisfying " + expr.getFilter().getExpression() + ":\n");
				probs.printFiltered(mainLog, filter);
			}

			// for "P=?"-type properties...
			if (pb == null) {
				// compute/print min info
				if (expr.filterMinRequested()) {
					minRes = probs.minOverBDD(filter);
					mainLog.print("\nMinimum " + (min ? "minimum" : "maximum") + " probability for states satisfying "
							+ expr.getFilter().getExpression() + ": " + minRes + "\n");
					tmp = probs.getBDDFromInterval(minRes - termCritParam, minRes + termCritParam);
					JDD.Ref(filter);
					tmp = JDD.And(tmp, filter);
					states = new StateListMTBDD(tmp, model);
					mainLog.print("There are " + states.size() + " states with this minimum probability (+/- "
							+ termCritParam + ")");
					if (!verbose && states.size() > 10) {
						mainLog
								.print(".\nThe first 10 states are displayed below. To view them all, use verbose mode.\n");
						states.print(mainLog, 10);
					} else {
						mainLog.print(":\n");
						states.print(mainLog);
					}
					JDD.Deref(tmp);
				}
				// compute/print min info
				if (expr.filterMaxRequested()) {
					maxRes = probs.maxOverBDD(filter);
					mainLog.print("\nMaximum " + (min ? "minimum" : "maximum") + " probability for states satisfying "
							+ expr.getFilter().getExpression() + ": " + maxRes + "\n");
					tmp = probs.getBDDFromInterval(maxRes - termCritParam, maxRes + termCritParam);
					JDD.Ref(filter);
					tmp = JDD.And(tmp, filter);
					states = new StateListMTBDD(tmp, model);
					mainLog.print("There are " + states.size() + " states with this maximum probability (+/- "
							+ termCritParam + ")");
					if (!verbose && states.size() > 10) {
						mainLog
								.print(".\nThe first 10 states are displayed below. To view them all, use verbose mode.\n");
						states.print(mainLog, 10);
					} else {
						mainLog.print(":\n");
						states.print(mainLog);
					}
					JDD.Deref(tmp);
				}
			}
		}

		// for P=? properties...
		// if there are multiple initial states or if there is a filter
		// satisfying
		// multiple states with no {min}/{max}/etc., print a warning...
		if (pb == null) {
			if (filter == null) {
				if (model.getNumStartStates() > 1) {
					mainLog
							.print("\nWarning: There are multiple initial states; the result of model checking is for the first one: ");
					model.getStartStates().print(mainLog, 1);
				}
			} else if (expr.noFilterRequests()) {
				StateListMTBDD filterStates = new StateListMTBDD(filter, model);
				if (filterStates.size() > 1) {
					mainLog.print("\nWarning: The filter {" + expr.getFilter().getExpression() + "} is satisfied by "
							+ filterStates.size() + " states.\n");
					mainLog.print("The result of model checking is for the first of these: ");
					filterStates.print(mainLog, 1);
				}
			}
		}

		// compute result of property
		// if there's a bound, get set of satisfying states
		if (pb != null) {
			sol = probs.getBDDFromInterval(relOp, p);
			// remove unreachable states from solution
			JDD.Ref(reach);
			sol = JDD.And(sol, reach);
		}
		// if there's no bound, result will be a probability
		else {
			// just store empty set for sol
			sol = JDD.Constant(0);
			// use filter if present
			if (filter != null) {
				if (expr.filterMinRequested())
					numericalRes = minRes;
				else if (expr.filterMaxRequested())
					numericalRes = maxRes;
				else
					numericalRes = probs.firstFromBDD(filter);
			}
			// otherwise use initial state
			else {
				numericalRes = probs.firstFromBDD(start);
			}
		}

		// free vector
		probs.clear();

		// derefs
		if (filter != null)
			JDD.Deref(filter);

		return sol;
	}

	// R operator

	private JDDNode checkExpressionReward(ExpressionReward expr) throws PrismException
	{
		Object rs; // reward struct index
		Expression rb; // reward bound (expression)
		double r = 0; // reward value (actual value)
		String relOp; // relational operator
		boolean min; // are we finding min (true) or max (false) rewards
		PathExpression pe; // path expression

		JDDNode stateRewards = null, transRewards = null, filter, sol, tmp;
		StateProbs rewards = null;
		StateList states = null;
		double minRes = 0, maxRes = 0;
		int i;

		// get info from reward operator
		rs = expr.getRewardStructIndex();
		relOp = expr.getRelOp();
		rb = expr.getReward();
		if (rb != null) {
			r = rb.evaluateDouble(constantValues, null);
			if (r < 0)
				throw new PrismException("Invalid reward bound " + r + " in R operator");
		}

		// get reward info
		if (model.getNumRewardStructs() == 0)
			throw new PrismException("Model has no rewards specified");
		if (rs == null) {
			stateRewards = model.getStateRewards(0);
			transRewards = model.getTransRewards(0);
		} else if (rs instanceof Expression) {
			i = ((Expression) rs).evaluateInt(constantValues, null);
			rs = new Integer(i); // for better error reporting below
			stateRewards = model.getStateRewards(i - 1);
			transRewards = model.getTransRewards(i - 1);
		} else if (rs instanceof String) {
			stateRewards = model.getStateRewards((String) rs);
			transRewards = model.getTransRewards((String) rs);
		}
		if (stateRewards == null || transRewards == null)
			throw new PrismException("Invalid reward structure index \"" + rs + "\"");

		// check for trivial (i.e. stupid) cases
		if (rb != null) {
			if (r == 0 && relOp.equals(">=")) {
				mainLog.print("\nWarning: checking for reward " + relOp + " " + r
						+ " - formula trivially satisfies all states\n");
				JDD.Ref(reach);
				return reach;
			} else if (r == 0 && relOp.equals("<")) {
				mainLog.print("\nWarning: checking for reward " + relOp + " " + r
						+ " - formula trivially satisfies no states\n");
				return JDD.Constant(0);
			}
		}

		// determine whether min or max rewards needed
		if (relOp.equals(">") || relOp.equals(">=") || relOp.equals("min=")) {
			// min
			min = true;
		} else if (relOp.equals("<") || relOp.equals("<=") || relOp.equals("max=")) {
			// max
			min = false;
		} else {
			throw new PrismException("Can't use \"R=?\" for MDPs; use \"Rmin=?\" or \"Rmax=?\"");
		}

		// translate filter (if present)
		filter = null;
		if (expr.getFilter() != null) {
			filter = checkExpression(expr.getFilter().getExpression());
		}

		// check if filter satisfies no states
		if (filter != null)
			if (filter.equals(JDD.ZERO)) {
				// for R=? properties, this is an error
				if (rb == null) {
					throw new PrismException("Filter {" + expr.getFilter().getExpression() + "} in R=?[] property satisfies no states");
				}
				// otherwise just print a warning
				else {
					mainLog.println("\nWarning: Filter {" + expr.getFilter().getExpression()
							+ "} satisfies no states and is being ignored");
					JDD.Deref(filter);
					filter = null;
				}
			}

		// compute rewards
		pe = expr.getPathExpression();
		try {
			if (pe instanceof PathExpressionTemporal) {
				switch (((PathExpressionTemporal) pe).getOperator()) {
				case PathExpressionTemporal.R_I:
					rewards = checkRewardInst((PathExpressionTemporal) pe, stateRewards, transRewards, min);
					break;
				case PathExpressionTemporal.R_F:
					rewards = checkRewardReach((PathExpressionTemporal) pe, stateRewards, transRewards, min);
					break;
				}
			}
			if (rewards == null)
				throw new PrismException("Unrecognised path operator in R operator");
		} catch (PrismException e) {
			if (filter != null)
				JDD.Deref(filter);
			throw e;
		}

		// round off rewards
		// rewards.roundOff(getPlacesToRoundBy());

		// print out rewards
		if (verbose) {
			mainLog.print("\n" + (min ? "Minimum" : "Maximum") + " rewards (non-zero only) for all states:\n");
			rewards.print(mainLog);
		}

		// if a filter was provided, there is some additional output...
		if (filter != null) {

			// for non-"R=?"-type properties, print out some rewards
			if (rb != null) {
				if (!expr.noFilterRequests()) {
					mainLog
							.println("\nWarning: \"{min}\", \"{max}\" only apply to \"R=?\" properties; they are ignored here.");
				}
				mainLog.print("\n" + (min ? "Minimum" : "Maximum") + " rewards (non-zero only) for states satisfying "
						+ expr.getFilter().getExpression() + ":\n");
				rewards.printFiltered(mainLog, filter);
			}

			// for "R=?"-type properties...
			if (rb == null) {
				// compute/print min info
				if (expr.filterMinRequested()) {
					minRes = rewards.minOverBDD(filter);
					mainLog.print("\nMinimum " + (min ? "minimum" : "maximum") + " reward for states satisfying "
							+ expr.getFilter().getExpression() + ": " + minRes + "\n");
					tmp = rewards.getBDDFromInterval(minRes - termCritParam, minRes + termCritParam);
					JDD.Ref(filter);
					tmp = JDD.And(tmp, filter);
					states = new StateListMTBDD(tmp, model);
					mainLog.print("There are " + states.size() + " states with this minimum reward (+/- "
							+ termCritParam + ")");
					if (!verbose && states.size() > 10) {
						mainLog
								.print(".\nThe first 10 states are displayed below. To view them all, use verbose mode.\n");
						states.print(mainLog, 10);
					} else {
						mainLog.print(":\n");
						states.print(mainLog);
					}
					JDD.Deref(tmp);
				}
				// compute/print min info
				if (expr.filterMaxRequested()) {
					maxRes = rewards.maxOverBDD(filter);
					mainLog.print("\nMaximum " + (min ? "minimum" : "maximum") + " reward for states satisfying "
							+ expr.getFilter().getExpression() + ": " + maxRes + "\n");
					tmp = rewards.getBDDFromInterval(maxRes - termCritParam, maxRes + termCritParam);
					JDD.Ref(filter);
					tmp = JDD.And(tmp, filter);
					states = new StateListMTBDD(tmp, model);
					mainLog.print("There are " + states.size() + " states with this maximum reward (+/- "
							+ termCritParam + ")");
					if (!verbose && states.size() > 10) {
						mainLog
								.print(".\nThe first 10 states are displayed below. To view them all, use verbose mode.\n");
						states.print(mainLog, 10);
					} else {
						mainLog.print(":\n");
						states.print(mainLog);
					}
					JDD.Deref(tmp);
				}
			}
		}

		// for R=? properties...
		// if there are multiple initial states or if there is a filter
		// satisfying
		// multiple states with no {min}/{max}/etc., print a warning...
		if (rb == null) {
			if (filter == null) {
				if (model.getNumStartStates() > 1) {
					mainLog
							.print("\nWarning: There are multiple initial states; the result of model checking is for the first one: ");
					model.getStartStates().print(mainLog, 1);
				}
			} else if (expr.noFilterRequests()) {
				StateListMTBDD filterStates = new StateListMTBDD(filter, model);
				if (filterStates.size() > 1) {
					mainLog.print("\nWarning: The filter {" + expr.getFilter().getExpression() + "} is satisfied by "
							+ filterStates.size() + " states.\n");
					mainLog.print("The result of model checking is for the first of these: ");
					filterStates.print(mainLog, 1);
				}
			}
		}

		// compute result of property
		// if there's a bound, get set of satisfying states
		if (rb != null) {
			sol = rewards.getBDDFromInterval(relOp, r);
			// remove unreachable states from solution
			JDD.Ref(reach);
			sol = JDD.And(sol, reach);
		}
		// if there's no bound, result will be a reward
		else {
			// just store empty set for sol
			sol = JDD.Constant(0);
			// use filter if present
			if (filter != null) {
				if (expr.filterMinRequested())
					numericalRes = minRes;
				else if (expr.filterMaxRequested())
					numericalRes = maxRes;
				else
					numericalRes = rewards.firstFromBDD(filter);
			}
			// otherwise use initial state
			else {
				numericalRes = rewards.firstFromBDD(start);
			}
		}

		// free vector
		rewards.clear();

		// derefs
		if (filter != null)
			JDD.Deref(filter);

		return sol;
	}

	// next

	private StateProbs checkProbNext(PathExpressionTemporal pe, boolean min) throws PrismException
	{
		Expression expr;
		JDDNode b;
		StateProbs probs = null;

		// get operand
		if (!(pe.getOperand2() instanceof PathExpressionExpr))
			throw new PrismException("Invalid path formula");
		expr = ((PathExpressionExpr) pe.getOperand2()).getExpression();

		// model check operand first
		b = checkExpression(expr);

		// print out some info about num states
		// mainLog.print("\nb = " + JDD.GetNumMintermsString(b,
		// allDDRowVars.n()) + " states\n");

		// compute probabilities
		probs = computeNextProbs(trans, b, min);

		// derefs
		JDD.Deref(b);

		return probs;
	}

	// bounded until

	private StateProbs checkProbBoundedUntil(PathExpressionTemporal pe, boolean min) throws PrismException
	{
		Expression expr1, expr2;
		int time;
		JDDNode b1, b2;
		StateProbs probs = null;

		// get operands
		if (!(pe.getOperand1() instanceof PathExpressionExpr) || !(pe.getOperand2() instanceof PathExpressionExpr))
			throw new PrismException("Invalid path formula");
		expr1 = ((PathExpressionExpr) pe.getOperand1()).getExpression();
		expr2 = ((PathExpressionExpr) pe.getOperand2()).getExpression();

		// get info from bounded until
		time = pe.getUpperBound().evaluateInt(constantValues, null);
		if (time < 0) {
			throw new PrismException("Invalid bound " + time + " in bounded until formula");
		}

		// model check operands first
		b1 = checkExpression(expr1);
		try {
			b2 = checkExpression(expr2);
		} catch (PrismException e) {
			JDD.Deref(b1);
			throw e;
		}

		// print out some info about num states
		// mainLog.print("\nb1 = " + JDD.GetNumMintermsString(b1,
		// allDDRowVars.n()));
		// mainLog.print(" states, b2 = " + JDD.GetNumMintermsString(b2,
		// allDDRowVars.n()) + " states\n");

		// compute probabilities

		// a trivial case: "U<=0"
		if (time == 0) {
			// prob is 1 in b2 states, 0 otherwise
			JDD.Ref(b2);
			probs = new StateProbsMTBDD(b2, model);
		} else {
			try {
				probs = computeBoundedUntilProbs(trans, trans01, b1, b2, time, min);
			} catch (PrismException e) {
				JDD.Deref(b1);
				JDD.Deref(b2);
				throw e;
			}
		}

		// derefs
		JDD.Deref(b1);
		JDD.Deref(b2);

		return probs;
	}

	// until (unbounded)

	private StateProbs checkProbUntil(PathExpressionTemporal pe, Expression pb, double p, boolean min)
			throws PrismException
	{
		Expression expr1, expr2;
		JDDNode b1, b2, splus, newb1, newb2;
		StateProbs probs = null;
		long l;

		// get operands
		if (!(pe.getOperand1() instanceof PathExpressionExpr) || !(pe.getOperand2() instanceof PathExpressionExpr))
			throw new PrismException("Invalid path formula");
		expr1 = ((PathExpressionExpr) pe.getOperand1()).getExpression();
		expr2 = ((PathExpressionExpr) pe.getOperand2()).getExpression();

		// model check operands first
		b1 = checkExpression(expr1);
		try {
			b2 = checkExpression(expr2);
		} catch (PrismException e) {
			JDD.Deref(b1);
			throw e;
		}

		// print out some info about num states
		// mainLog.print("\nb1 = " + JDD.GetNumMintermsString(b1,
		// allDDRowVars.n()));
		// mainLog.print(" states, b2 = " + JDD.GetNumMintermsString(b2,
		// allDDRowVars.n()) + " states\n");

		// compute probabilities

		// if doing min with fairness, we solve a different problem
		// (as in christel's habilitation)
		if (min && fairness) {

			// print out reminder that we have to do conversion for fairness
			mainLog.print("\nDoing conversion for fairness...\n");
			// start timer
			l = System.currentTimeMillis();
			// convert to new problem
			newb2 = PrismMTBDD.Prob0A(trans01, reach, allDDRowVars, allDDColVars, allDDNondetVars, b1, b2);
			JDD.Ref(newb2);
			splus = JDD.Not(newb2);
			JDD.Ref(splus);
			JDD.Ref(b2);
			newb1 = JDD.And(splus, JDD.Not(b2));
			JDD.Deref(b1);
			JDD.Deref(b2);
			JDD.Deref(splus);
			JDD.Ref(reach);
			b1 = JDD.And(reach, newb1);
			JDD.Ref(reach);
			b2 = JDD.And(reach, newb2);
			// stop timer
			l = System.currentTimeMillis() - l;
			// print out conversion info
			mainLog.print("\nTime for fairness conversion: " + l / 1000.0 + " seconds.\n");
			// mainLog.print("\nb1 = " + JDD.GetNumMintermsString(b1,
			// allDDRowVars.n()));
			// mainLog.print(" states, b2 = " + JDD.GetNumMintermsString(b2,
			// allDDRowVars.n()) + " states\n");
		}

		// if p is 0 or 1 and precomputation algorithms are enabled, compute
		// probabilities qualitatively
		if (pb != null && ((p == 0) || (p == 1)) && precomp) {
			mainLog.print("\nWarning: probability bound in formula is " + p
					+ " so exact probabilities may not be computed\n");
			// for fairness, we compute max here
			probs = computeUntilProbsQual(trans01, b1, b2, min && !fairness);
		}
		// otherwise actually compute probabilities
		else {
			// for fairness, we compute max here
			try {
				probs = computeUntilProbs(trans, trans01, b1, b2, min && !fairness);
			} catch (PrismException e) {
				JDD.Deref(b1);
				JDD.Deref(b2);
				throw e;
			}
		}

		// if we're doing min with fairness,
		// we need to subtract the probabilities from
		// one to get the actual answer
		if (min && fairness) {
			probs.subtractFromOne();
		}

		// derefs
		JDD.Deref(b1);
		JDD.Deref(b2);

		return probs;
	}

	// bounded future (eventually)
	// F<=k phi == true U<=k phi

	private StateProbs checkProbBoundedFuture(PathExpressionTemporal pe, boolean min) throws PrismException
	{
		PathExpressionTemporal pe2;
		pe2 = new PathExpressionTemporal(PathExpressionTemporal.P_U, new PathExpressionExpr(Expression.True()), pe
				.getOperand2(), pe.getLowerBound(), pe.getUpperBound());
		return checkProbBoundedUntil(pe2, min);
	}

	// future (eventually)
	// F phi == true U phi

	private StateProbs checkProbFuture(PathExpressionTemporal pe, Expression pb, double p, boolean min)
			throws PrismException
	{
		PathExpressionTemporal pe2;
		pe2 = new PathExpressionTemporal(PathExpressionTemporal.P_U, new PathExpressionExpr(Expression.True()), pe
				.getOperand2());
		return checkProbUntil(pe2, pb, p, min);
	}

	// bounded global (always)
	// F<=k phi == true U<=k phi
	// P(G<=k phi) == 1-P(true U<=k !phi)
	// Pmin(G<=k phi) == 1-Pmax(true U<=k !phi)

	private StateProbs checkProbBoundedGlobal(PathExpressionTemporal pe, boolean min) throws PrismException
	{
		PathExpressionTemporal pe2;
		StateProbs probs;
		pe2 = new PathExpressionTemporal(PathExpressionTemporal.P_U, new PathExpressionExpr(Expression.True()),
				new PathExpressionExpr(Expression.Not(((PathExpressionExpr)pe.getOperand2()).getExpression())), pe.getLowerBound(), pe.getUpperBound());
		probs = checkProbBoundedUntil(pe2, !min);
		probs.subtractFromOne();
		return probs;
	}

	// global (always)
	// G phi == !(true U !phi)
	// P(G phi) == 1-P(true U !phi)
	// Pmin(G phi) == 1-Pmax(true U !phi)

	private StateProbs checkProbGlobal(PathExpressionTemporal pe, Expression pb, double p, boolean min)
			throws PrismException
	{
		PathExpressionTemporal pe2;
		StateProbs probs;
		pe2 = new PathExpressionTemporal(PathExpressionTemporal.P_U, new PathExpressionExpr(Expression.True()),
				new PathExpressionExpr(Expression.Not(((PathExpressionExpr)pe.getOperand2()).getExpression())));
		probs = checkProbUntil(pe2, pb, p, !min);
		probs.subtractFromOne();
		return probs;
	}

	// inst reward

	private StateProbs checkRewardInst(PathExpressionTemporal pe, JDDNode stateRewards, JDDNode transRewards,
			boolean min) throws PrismException
	{
		int time; // time bound
		StateProbs rewards = null;

		// get info from bounded until
		time = pe.getUpperBound().evaluateInt(constantValues, null);
		if (time < 0) {
			throw new PrismException("Invalid bound " + time + " in instantaneous reward property");
		}

		// compute rewards
		rewards = computeInstRewards(trans, stateRewards, time, min);

		return rewards;
	}

	// reach reward

	private StateProbs checkRewardReach(PathExpressionTemporal pe, JDDNode stateRewards, JDDNode transRewards,
			boolean min) throws PrismException
	{
		Expression expr;
		JDDNode b;
		StateProbs rewards = null;

		// get operand
		if (!(pe.getOperand2() instanceof PathExpressionExpr))
			throw new PrismException("Invalid path formula");
		expr = ((PathExpressionExpr) pe.getOperand2()).getExpression();

		// model check operand first
		b = checkExpression(expr);

		// print out some info about num states
		// mainLog.print("\nb = " + JDD.GetNumMintermsString(b,
		// allDDRowVars.n()) + " states\n");

		// compute rewards
		try {
			rewards = computeReachRewards(trans, trans01, stateRewards, transRewards, b, min);
		} catch (PrismException e) {
			JDD.Deref(b);
			throw e;
		}

		// derefs
		JDD.Deref(b);

		return rewards;
	}

	// Check label

	private JDDNode checkExpressionLabel(ExpressionLabel expr) throws PrismException
	{
		LabelList ll;
		JDDNode dd;
		int i;

		// treat special cases
		if (expr.getName().equals("deadlock")) {
			dd = model.getFixedDeadlocks();
			JDD.Ref(dd);
			return dd;
		}
		else if (expr.getName().equals("init")) {
			dd = start;
			JDD.Ref(dd);
			return dd;
		}
		// get expression associated with label
		ll = propertiesFile.getCombinedLabelList();
		i = ll.getLabelIndex(expr.getName());
		if (i == -1)
			throw new PrismException("Unknown label \"" + expr.getName() + "\" in property");
		// check recursively
		return checkExpression(ll.getLabel(i));
	}

	// -----------------------------------------------------------------------------------
	// probability computation methods
	// -----------------------------------------------------------------------------------

	// compute probabilities for next

	private StateProbs computeNextProbs(JDDNode tr, JDDNode b, boolean min)
	{
		JDDNode tmp;
		StateProbs probs = null;

		// matrix multiply: trans * b
		JDD.Ref(b);
		tmp = JDD.PermuteVariables(b, allDDRowVars, allDDColVars);
		JDD.Ref(tr);
		tmp = JDD.MatrixMultiply(tr, tmp, allDDColVars, JDD.BOULDER);
		// (then min or max)
		if (min) {
			// min
			JDD.Ref(nondetMask);
			tmp = JDD.Apply(JDD.MAX, tmp, nondetMask);
			tmp = JDD.MinAbstract(tmp, allDDNondetVars);
		} else {
			// max
			tmp = JDD.MaxAbstract(tmp, allDDNondetVars);
		}
		probs = new StateProbsMTBDD(tmp, model);

		return probs;
	}

	// compute probabilities for bounded until

	private StateProbs computeBoundedUntilProbs(JDDNode tr, JDDNode tr01, JDDNode b1, JDDNode b2, int time, boolean min)
			throws PrismException
	{
		JDDNode yes, no, maybe;
		JDDNode probsMTBDD;
		DoubleVector probsDV;
		StateProbs probs = null;

		// compute yes/no/maybe states
		if (b2.equals(JDD.ZERO)) {
			yes = JDD.Constant(0);
			JDD.Ref(reach);
			no = reach;
			maybe = JDD.Constant(0);
		} else if (b1.equals(JDD.ZERO)) {
			JDD.Ref(b2);
			yes = b2;
			JDD.Ref(reach);
			JDD.Ref(b2);
			no = JDD.And(reach, JDD.Not(b2));
			maybe = JDD.Constant(0);
		} else {
			// yes
			JDD.Ref(b2);
			yes = b2;
			// no
			if (yes.equals(reach)) {
				no = JDD.Constant(0);
			} else if (precomp) {
				if (min) {
					// "min prob = 0" equates to "there exists a prob 0"
					no = PrismMTBDD.Prob0E(tr01, reach, nondetMask, allDDRowVars, allDDColVars, allDDNondetVars, b1,
							yes);
				} else {
					// "max prob = 0" equates to "all probs 0"
					no = PrismMTBDD.Prob0A(tr01, reach, allDDRowVars, allDDColVars, allDDNondetVars, b1, yes);
				}
			} else {
				JDD.Ref(reach);
				JDD.Ref(b1);
				JDD.Ref(b2);
				no = JDD.And(reach, JDD.Not(JDD.Or(b1, b2)));
			}
			// maybe
			JDD.Ref(reach);
			JDD.Ref(yes);
			JDD.Ref(no);
			maybe = JDD.And(reach, JDD.Not(JDD.Or(yes, no)));
		}

		// print out yes/no/maybe
		mainLog.print("\nyes = " + JDD.GetNumMintermsString(yes, allDDRowVars.n()));
		mainLog.print(", no = " + JDD.GetNumMintermsString(no, allDDRowVars.n()));
		mainLog.print(", maybe = " + JDD.GetNumMintermsString(maybe, allDDRowVars.n()) + "\n");

		// if maybe is empty, we have the probabilities already
		if (maybe.equals(JDD.ZERO)) {
			JDD.Ref(yes);
			probs = new StateProbsMTBDD(yes, model);
		}
		// otherwise explicitly compute the remaining probabilities
		else {
			// compute probabilities
			try {
				switch (engine) {
				case Prism.MTBDD:
					probsMTBDD = PrismMTBDD.NondetBoundedUntil(tr, odd, nondetMask, allDDRowVars, allDDColVars,
							allDDNondetVars, yes, maybe, time, min);
					probs = new StateProbsMTBDD(probsMTBDD, model);
					break;
				case Prism.SPARSE:
					probsDV = PrismSparse.NondetBoundedUntil(tr, odd, allDDRowVars, allDDColVars, allDDNondetVars, yes,
							maybe, time, min);
					probs = new StateProbsDV(probsDV, model);
					break;
				case Prism.HYBRID:
					probsDV = PrismHybrid.NondetBoundedUntil(tr, odd, allDDRowVars, allDDColVars, allDDNondetVars, yes,
							maybe, time, min);
					probs = new StateProbsDV(probsDV, model);
					break;
				default:
					throw new PrismException("Engine does not support this numerical method");
				}
			} catch (PrismException e) {
				JDD.Deref(yes);
				JDD.Deref(no);
				JDD.Deref(maybe);
				throw e;
			}
		}

		// derefs
		JDD.Deref(yes);
		JDD.Deref(no);
		JDD.Deref(maybe);

		return probs;
	}

	// compute probabilities for until (for qualitative properties)

	// note: this function doesn't need to know anything about fairness
	// it is just told whether to compute min or max probabilities

	private StateProbs computeUntilProbsQual(JDDNode tr01, JDDNode b1, JDDNode b2, boolean min)
	{
		JDDNode yes = null, no = null, maybe;
		StateProbs probs = null;

		// note: we know precomputation is enabled else this function wouldn't
		// have been called

		// compute yes/no/maybe states
		if (b2.equals(JDD.ZERO)) {
			yes = JDD.Constant(0);
			JDD.Ref(reach);
			no = reach;
			maybe = JDD.Constant(0);
		} else if (b1.equals(JDD.ZERO)) {
			JDD.Ref(b2);
			yes = b2;
			JDD.Ref(reach);
			JDD.Ref(b2);
			no = JDD.And(reach, JDD.Not(b2));
			maybe = JDD.Constant(0);
		} else {
			// min
			if (min) {
				// no: "min prob = 0" equates to "there exists an adversary prob
				// equals 0"
				no = PrismMTBDD.Prob0E(tr01, reach, nondetMask, allDDRowVars, allDDColVars, allDDNondetVars, b1, b2);
				// yes: "min prob = 1" equates to "for all adversaries prob
				// equals 1"
				if (no.equals(JDD.ZERO)) {
					// if there are no "no" states, all states must have min
					// prob 1
					JDD.Ref(reach);
					yes = reach;
				} else {
					// use precomp alg
					// note that prob1a (unlike prob1e) requires no/b2, not
					// b1/b2
					yes = PrismMTBDD.Prob1A(tr01, reach, nondetMask, allDDRowVars, allDDColVars, allDDNondetVars, no,
							b2);
				}
			}
			// max
			else {
				// yes: "max prob = 1" equates to "there exists an adversary
				// prob equals 1"
				yes = PrismMTBDD.Prob1E(tr01, allDDRowVars, allDDColVars, allDDNondetVars, b1, b2);
				// no: "max prob = 0" equates to "for all adversaries prob
				// equals 0"
				if (yes.equals(reach)) {
					// trivial case - because yes+no+maybe=all, so yes=all =>
					// no,maybe=empty
					no = JDD.Constant(0);
				} else {
					// use precomp alg
					no = PrismMTBDD.Prob0A(tr01, reach, allDDRowVars, allDDColVars, allDDNondetVars, b1, yes);
				}
			}
			// maybe
			JDD.Ref(reach);
			JDD.Ref(yes);
			JDD.Ref(no);
			maybe = JDD.And(reach, JDD.Not(JDD.Or(yes, no)));
		}

		// print out yes/no/maybe
		mainLog.print("\nyes = " + JDD.GetNumMintermsString(yes, allDDRowVars.n()));
		mainLog.print(", no = " + JDD.GetNumMintermsString(no, allDDRowVars.n()));
		mainLog.print(", maybe = " + JDD.GetNumMintermsString(maybe, allDDRowVars.n()) + "\n");

		// if maybe is empty, we have the answer already...
		if (maybe.equals(JDD.ZERO)) {
			JDD.Ref(yes);
			probs = new StateProbsMTBDD(yes, model);
		}
		// otherwise we set the probabilities for maybe states to be 0.5
		// (actual probabilities for these states are unknown but definitely >0
		// and <1)
		// (this is safe because the results of this function will only be used
		// to compare against 0/1 bounds)
		// (this is not entirely elegant but is simpler and less error prone
		// than
		// trying to work out whether to use 0/1 for all case of min/max,
		// fairness, future/global, etc.)
		else {
			JDD.Ref(yes);
			JDD.Ref(maybe);
			probs = new StateProbsMTBDD(JDD.Apply(JDD.PLUS, yes, JDD.Apply(JDD.TIMES, maybe, JDD.Constant(0.5))), model);
		}

		// derefs
		JDD.Deref(yes);
		JDD.Deref(no);
		JDD.Deref(maybe);

		return probs;
	}

	// compute probabilities for until (general case)

	// note: this function doesn't need to know anything about fairness
	// it is just told whether to compute min or max probabilities

	private StateProbs computeUntilProbs(JDDNode tr, JDDNode tr01, JDDNode b1, JDDNode b2, boolean min)
			throws PrismException
	{
		JDDNode yes, no, maybe;
		JDDNode probsMTBDD;
		DoubleVector probsDV;
		StateProbs probs = null;

		// compute yes/no/maybe states
		if (b2.equals(JDD.ZERO)) {
			yes = JDD.Constant(0);
			JDD.Ref(reach);
			no = reach;
			maybe = JDD.Constant(0);
		} else if (b1.equals(JDD.ZERO)) {
			JDD.Ref(b2);
			yes = b2;
			JDD.Ref(reach);
			JDD.Ref(b2);
			no = JDD.And(reach, JDD.Not(b2));
			maybe = JDD.Constant(0);
		} else {
			// if precomputation enabled
			if (precomp) {
				// min
				if (min) {
					// no: "min prob = 0" equates to "there exists an adversary
					// prob equals 0"
					no = PrismMTBDD
							.Prob0E(tr01, reach, nondetMask, allDDRowVars, allDDColVars, allDDNondetVars, b1, b2);
					// yes: "min prob = 1" equates to "for all adversaries prob
					// equals 1"
					if (no.equals(JDD.ZERO)) {
						// if there are no "no" states, all states must have min
						// prob 1
						JDD.Ref(reach);
						yes = reach;
					} else {
						// use precomp alg
						// note that prob1a (unlike prob1e) requires no/b2, not
						// b1/b2
						yes = PrismMTBDD.Prob1A(tr01, reach, nondetMask, allDDRowVars, allDDColVars, allDDNondetVars,
								no, b2);
					}
				}
				// max
				else {
					// yes: "max prob = 1" equates to "there exists an adversary
					// prob equals 1"
					yes = PrismMTBDD.Prob1E(tr01, allDDRowVars, allDDColVars, allDDNondetVars, b1, b2);
					// no: "max prob = 0" equates to "for all adversaries prob
					// equals 0"
					if (yes.equals(reach)) {
						// trivial case - because yes+no+maybe=all, so yes=all
						// => no,maybe=empty
						no = JDD.Constant(0);
					} else {
						// use precomp alg
						no = PrismMTBDD.Prob0A(tr01, reach, allDDRowVars, allDDColVars, allDDNondetVars, b1, yes);
					}
				}
			}
			// if precomputation not enabled
			else {
				// yes
				JDD.Ref(b2);
				yes = b2;
				// no
				JDD.Ref(reach);
				JDD.Ref(b1);
				JDD.Ref(b2);
				no = JDD.And(reach, JDD.Not(JDD.Or(b1, b2)));
			}
			// maybe
			JDD.Ref(reach);
			JDD.Ref(yes);
			JDD.Ref(no);
			maybe = JDD.And(reach, JDD.Not(JDD.Or(yes, no)));
		}

		// print out yes/no/maybe
		mainLog.print("\nyes = " + JDD.GetNumMintermsString(yes, allDDRowVars.n()));
		mainLog.print(", no = " + JDD.GetNumMintermsString(no, allDDRowVars.n()));
		mainLog.print(", maybe = " + JDD.GetNumMintermsString(maybe, allDDRowVars.n()) + "\n");

		// if maybe is empty, we have the answer already...
		if (maybe.equals(JDD.ZERO)) {
			JDD.Ref(yes);
			probs = new StateProbsMTBDD(yes, model);
		}
		// otherwise we compute the actual probabilities
		else {
			// compute probabilities
			mainLog.println("\nComputing remaining probabilities...");

			try {
				switch (engine) {
				case Prism.MTBDD:
					probsMTBDD = PrismMTBDD.NondetUntil(tr, odd, nondetMask, allDDRowVars, allDDColVars,
							allDDNondetVars, yes, maybe, min);
					probs = new StateProbsMTBDD(probsMTBDD, model);
					break;
				case Prism.SPARSE:
					probsDV = PrismSparse.NondetUntil(tr, odd, allDDRowVars, allDDColVars, allDDNondetVars, yes, maybe,
							min);
					probs = new StateProbsDV(probsDV, model);
					break;
				case Prism.HYBRID:
					probsDV = PrismHybrid.NondetUntil(tr, odd, allDDRowVars, allDDColVars, allDDNondetVars, yes, maybe,
							min);
					probs = new StateProbsDV(probsDV, model);
					break;
				default:
					throw new PrismException("Engine does not support this numerical method");
				}
			} catch (PrismException e) {
				JDD.Deref(yes);
				JDD.Deref(no);
				JDD.Deref(maybe);
				throw e;
			}

			// round off probabilities
			// probs.roundOff(getPlacesToRoundBy());
		}

		// derefs
		JDD.Deref(yes);
		JDD.Deref(no);
		JDD.Deref(maybe);

		return probs;
	}

	// compute rewards for inst reward

	private StateProbs computeInstRewards(JDDNode tr, JDDNode sr, int time, boolean min) throws PrismException
	{
		JDDNode rewardsMTBDD;
		DoubleVector rewardsDV;
		StateProbs rewards = null;

		// a trivial case: "=0"
		if (time == 0) {
			JDD.Ref(sr);
			rewards = new StateProbsMTBDD(sr, model);
		}
		// otherwise we compute the actual rewards
		else {
			// compute the rewards
			try {
				switch (engine) {
				case Prism.MTBDD:
					rewardsMTBDD = PrismMTBDD.NondetInstReward(tr, sr, odd, nondetMask, allDDRowVars, allDDColVars,
							allDDNondetVars, time, min, start);
					rewards = new StateProbsMTBDD(rewardsMTBDD, model);
					break;
				case Prism.SPARSE:
					rewardsDV = PrismSparse.NondetInstReward(tr, sr, odd, allDDRowVars, allDDColVars, allDDNondetVars,
							time, min, start);
					rewards = new StateProbsDV(rewardsDV, model);
					break;
				default:
					throw new PrismException("Engine does not support this numerical method");
				}
			} catch (PrismException e) {
				throw e;
			}
		}

		return rewards;
	}

	// compute rewards for reach reward

	private StateProbs computeReachRewards(JDDNode tr, JDDNode tr01, JDDNode sr, JDDNode trr, JDDNode b, boolean min)
			throws PrismException
	{
		JDDNode inf, maybe, prob1, no;
		JDDNode rewardsMTBDD;
		DoubleVector rewardsDV;
		StateProbs rewards = null;

		// compute states which can't reach goal with probability 1
		if (b.equals(JDD.ZERO)) {
			JDD.Ref(reach);
			inf = reach;
			maybe = JDD.Constant(0);
		} else if (b.equals(reach)) {
			inf = JDD.Constant(0);
			maybe = JDD.Constant(0);
		} else {
			if (!min) {
				// compute states for which some adversaries don't reach goal
				// with probability 1
				// note that prob1a (unlike prob1e) requires no/b2, not b1/b2
				// hence we have to call prob0e first
				no = PrismMTBDD
						.Prob0E(tr01, reach, nondetMask, allDDRowVars, allDDColVars, allDDNondetVars, JDD.ONE, b);
				prob1 = PrismMTBDD.Prob1A(tr01, reach, nondetMask, allDDRowVars, allDDColVars, allDDNondetVars, no, b);
				JDD.Deref(no);
				JDD.Ref(reach);
				inf = JDD.And(reach, JDD.Not(prob1));
			} else {
				// compute states for which all adversaries don't reach goal
				// with probability 1
				prob1 = PrismMTBDD.Prob1E(tr01, allDDRowVars, allDDColVars, allDDNondetVars, reach, b);
				JDD.Ref(reach);
				inf = JDD.And(reach, JDD.Not(prob1));
			}
			JDD.Ref(reach);
			JDD.Ref(inf);
			JDD.Ref(b);
			maybe = JDD.And(reach, JDD.Not(JDD.Or(inf, b)));
		}

		// need to deal with zero loops yet
		if (min) {
			mainLog.println("\nWarning: PRISM hasn't checked for zero-reward loops.");
			mainLog.println("Your minimum rewards may be too low...");
		}

		// print out yes/no/maybe
		mainLog.print("\ngoal = " + JDD.GetNumMintermsString(b, allDDRowVars.n()));
		mainLog.print(", inf = " + JDD.GetNumMintermsString(inf, allDDRowVars.n()));
		mainLog.print(", maybe = " + JDD.GetNumMintermsString(maybe, allDDRowVars.n()) + "\n");

		// if maybe is empty, we have the rewards already
		if (maybe.equals(JDD.ZERO)) {
			JDD.Ref(inf);
			rewards = new StateProbsMTBDD(JDD.ITE(inf, JDD.PlusInfinity(), JDD.Constant(0)), model);
		}
		// otherwise we compute the actual rewards
		else {
			// compute the rewards
			try {
				switch (engine) {
				case Prism.MTBDD:
					rewardsMTBDD = PrismMTBDD.NondetReachReward(tr, sr, trr, odd, nondetMask, allDDRowVars,
							allDDColVars, allDDNondetVars, b, inf, maybe, min);
					rewards = new StateProbsMTBDD(rewardsMTBDD, model);
					break;
				case Prism.SPARSE:
					rewardsDV = PrismSparse.NondetReachReward(tr, sr, trr, odd, allDDRowVars, allDDColVars,
							allDDNondetVars, b, inf, maybe, min);
					rewards = new StateProbsDV(rewardsDV, model);
					break;
				case Prism.HYBRID:
					throw new PrismException("This functionality is not yet supported for this engine");
					// rewardsDV = PrismHybrid.NondetReachReward(tr, sr, trr,
					// odd, allDDRowVars, allDDColVars, allDDNondetVars, b, inf,
					// maybe, min);
					// rewards = new StateProbsDV(rewardsDV, model);
					// break;
				default:
					throw new PrismException("Engine does not support this numerical method");
				}
			} catch (PrismException e) {
				JDD.Deref(inf);
				JDD.Deref(maybe);
				throw e;
			}
		}

		// derefs
		JDD.Deref(inf);
		JDD.Deref(maybe);

		return rewards;
	}
}

// ------------------------------------------------------------------------------
