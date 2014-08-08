//==============================================================================
//	
//	Copyright (c) 2013-
//	Authors:
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

import java.util.LinkedList;

import explicit.StateValues;
import parser.Values;
import parser.ast.Expression;
import parser.ast.ExpressionProb;
import parser.ast.ExpressionReward;
import parser.ast.ExpressionTemporal;
import parser.ast.LabelList;
import parser.ast.ModulesFile;
import parser.ast.PropertiesFile;
import parser.ast.RewardStruct;
import prism.PrismComponent;
import prism.PrismException;
import prism.Result;
import simulator.PrismModelExplorer;
import simulator.SimulatorEngine;

/**
 * CTMC model checker based on fast adaptive uniformisation.
 */
public final class PSEFastAdaptiveUniformisationModelChecker extends PrismComponent
{
	// Region factory
	private BoxRegionFactory regionFactory;
	// Model file
	private ModulesFile modulesFile;
	// Properties file
	private PropertiesFile propertiesFile;
	// Constants from model
	private Values constantValues;
	// Labels from the model
	private LabelList labelListModel;
	// Labels from the property file
	private LabelList labelListProp;
	
	/**
	 * Constructor.
	 */
	public PSEFastAdaptiveUniformisationModelChecker(PrismComponent parent, BoxRegionFactory regionFactory)
	{
		super(parent);
		this.regionFactory = regionFactory;
	}

	public void setModulesFileAndPropertiesFile(ModulesFile modulesFile, PropertiesFile propertiesFile)
	{
		this.modulesFile = modulesFile;
		this.propertiesFile = propertiesFile;
		// Get combined constant values from model/properties
		constantValues = new Values();
		constantValues.addValues(modulesFile.getConstantValues());
		if (propertiesFile != null)
			constantValues.addValues(propertiesFile.getConstantValues());
		this.labelListModel = modulesFile.getLabelList();
		this.labelListProp = propertiesFile.getLabelList();

		//stateChecker.setModulesFileAndPropertiesFile(modulesFile, propertiesFile);
	}

	/**
	 * Model check a property.
	 */
	public Result check(Expression expr) throws PrismException
	{
		Result res;
		String resultString;
		long timer;

		// Starting model checking
		timer = System.currentTimeMillis();

		// Do model checking
		res = checkExpression(expr);

		// Model checking complete
		timer = System.currentTimeMillis() - timer;
		mainLog.println("\nModel checking completed in " + (timer / 1000.0) + " secs.");

		// Print result to log
		resultString = "Result";
		if (!("Result".equals(expr.getResultName())))
			resultString += " (" + expr.getResultName().toLowerCase() + ")";
		resultString += ": " + res;
		mainLog.print("\n" + resultString + "\n");

		// Return result
		return res;
	}

	/**
	 * Model check an expression (used recursively).
	 */
	private Result checkExpression(Expression expr) throws PrismException
	{
		if (expr instanceof ExpressionProb) {
			return checkExpressionProb((ExpressionProb) expr);
		} else {
			throw new PrismException("PSE+FAU not yet supported for this operator");
		}
	}

	/**
	 * Model check a P operator.
	 */
	private Result checkExpressionProb(ExpressionProb expr) throws PrismException
	{
		return null;

		// TODO
		/*
		// Check whether P=? (only case allowed)
		if (expr.getProb() != null) {
			throw new PrismException("PSE+FAU model checking currently only supports P=? properties");
		}

		if (!(expr.getExpression() instanceof ExpressionTemporal)) {
			throw new PrismException("PSE+FAU model checking currently only supports simple path operators");
		}
		ExpressionTemporal exprTemp = (ExpressionTemporal) expr.getExpression();
		if (!exprTemp.isSimplePathFormula()) {
			throw new PrismException("PSE+FAU window model checking currently only supports simple until operators");
		}

		double timeLower = 0.0;
		if (exprTemp.getLowerBound() != null) {
			timeLower = exprTemp.getLowerBound().evaluateDouble(constantValues);
		}
		if (exprTemp.getUpperBound() == null) {
			throw new PrismException("PSE+FAU window model checking currently requires an upper time bound");
		}
		double timeUpper = exprTemp.getUpperBound().evaluateDouble(constantValues);

		if (!exprTemp.hasBounds()) {
			throw new PrismException("PSE+FAU window model checking currently only supports timed properties");
		}

		mainLog.println("Starting transient probability computation using PSE+FAU...");
		PSEModelExplorer modelExplorer = new PSEModelExplorer(modulesFile);

		PSEFastAdaptiveUniformisation fau = new PSEFastAdaptiveUniformisation(this, modelExplorer);

		Expression op1 = exprTemp.getOperand1();
		if (op1 == null) {
			op1 = Expression.True();
		}
		Expression op2 = exprTemp.getOperand2();
		op1 = (Expression) op1.expandPropRefsAndLabels(propertiesFile, labelListModel);
		op1 = (Expression) op1.expandPropRefsAndLabels(propertiesFile, labelListProp);
		op2 = (Expression) op2.expandPropRefsAndLabels(propertiesFile, labelListModel);
		op2 = (Expression) op2.expandPropRefsAndLabels(propertiesFile, labelListProp);
		int operator = exprTemp.getOperator();

		Expression sink = null;
		Expression target = null;
		switch (operator) {
		case ExpressionTemporal.P_U:
		case ExpressionTemporal.P_F:
			sink = Expression.Not(op1);
			break;
		case ExpressionTemporal.P_G:
			sink = Expression.False();
			break;
		case ExpressionTemporal.P_W:
		case ExpressionTemporal.P_R:
		default:
			throw new PrismException("Operator currently not supported for PSE+FAU");
		}
		
		fau.setSink(sink);
		fau.computeTransientProbsAdaptive(timeLower);
		switch (operator) {
		case ExpressionTemporal.P_U:
		case ExpressionTemporal.P_F:
			sink = Expression.Or(Expression.Not(op1), op2);
			target = op2;
			break;
		case ExpressionTemporal.P_G:
			sink = Expression.Not(op2);
			target = op2;
			break;
		case ExpressionTemporal.P_W:
		case ExpressionTemporal.P_R:
		default:
			throw new PrismException("Operator currently not supported for PSE+FAU");
		}
		Values varValues = new Values();
		varValues.addValue("deadlock", "true");
		sink.replaceVars(varValues);
		fau.setAnalysisType(FastAdaptiveUniformisation.AnalysisType.REACH);
		fau.setSink(sink);
		fau.setTarget(target);
		fau.computeTransientProbsAdaptive(timeUpper - timeLower);
		mainLog.println("\nTotal probability lost is : " + fau.getTotalDiscreteLoss());
		mainLog.println("Maximal number of states stored during analysis : " + fau.getMaxNumStates());

		return new Result(new Double(fau.getValue()));
		*/
	}
	
	// Transient analysis
	
	public BoxRegionValues doTransient(PSEModelExplorer modelExplorer, double t, StateValues initDistMin, StateValues initDistMax, DecompositionProcedure decompositionProcedure)
			throws PrismException
	{
		// TODO: initial distributions

		// For decomposing of the parameter space
		LinkedList<BoxRegion> regions = new LinkedList<BoxRegion>();
		regions.add(regionFactory.completeSpace());

		// Start bounded probabilistic reachability
		long timer = System.currentTimeMillis();
		mainLog.println("\nStarting PSE+FAU transient probability computation...");

		BoxRegionValues regionValues = new BoxRegionValues();
		PSEFastAdaptiveUniformisation fau = new PSEFastAdaptiveUniformisation(this, modelExplorer);
		while (regions.size() != 0) {
			BoxRegion region = regions.remove();
			fau.configureParameterSpace(region);
			mainLog.println("Computing probabilities for parameter region " + region);

			StateValues probsMin, probsMax;
			probsMin = fau.doTransient(t, initDistMin, PSEFastAdaptiveUniformisation.VectorType.MIN);
			probsMax = fau.doTransient(t, initDistMax, PSEFastAdaptiveUniformisation.VectorType.MAX);
			try {
				decompositionProcedure.examinePartialComputation(regionValues, region, probsMin.getDoubleArray(), probsMax.getDoubleArray());
				regionValues.put(region, probsMin, probsMax);
			} catch (DecompositionProcedure.DecompositionNeeded e) {
				e.printRegionsToDecompose(mainLog);
				for (BoxRegion regionToDecompose : e.getRegionsToDecompose()) {
					regions.addAll(regionToDecompose.decompose());
				}
			}
		}

		// Finished bounded probabilistic reachability
		timer = System.currentTimeMillis() - timer;
		mainLog.print("\nPSE+FAU transient probability computation");
		mainLog.print(" took " + timer / 1000.0 + " seconds");
		mainLog.println(" (producing " + regionValues.getNumRegions() + " final regions).");

		return regionValues;
	}
}
