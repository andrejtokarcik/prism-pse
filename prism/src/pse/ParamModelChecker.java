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

import parser.type.TypeDouble;
import prism.PrismComponent;
import prism.PrismException;
import prism.PrismLog;
import prism.PrismPrintStreamLog;
import explicit.FoxGlynn;
import explicit.Model;
import explicit.ModelCheckerResult;
import explicit.ranged.ModelCheckerResultRanged;

/**
 * Model checker for parametric Markov models.
 * 
 * @author Ernst Moritz Hahn <emhahn@cs.ox.ac.uk> (University of Oxford)
 */
final public class ParamModelChecker extends PrismComponent
{
	// Log for output (default to System.out)
	private PrismLog mainLog = new PrismPrintStreamLog(System.out);
	
	/**
	 * Constructor
	 */
	public ParamModelChecker(PrismComponent parent) throws PrismException
	{
		super(parent);
	}
	
	// Setters/getters

	/**
	 * Set the attached model file (for e.g. reward structures when model checking)
	 * and the attached properties file (for e.g. constants/labels when model checking)
	 */
	/*
	public void setModulesFileAndPropertiesFile(ModulesFile modulesFile, PropertiesFile propertiesFile)
	{
		this.modulesFile = modulesFile;
		this.propertiesFile = propertiesFile;
		// Get combined constant values from model/properties
		constantValues = new Values();
		constantValues.addValues(modulesFile.getConstantValues());
		if (propertiesFile != null)
			constantValues.addValues(propertiesFile.getConstantValues());
	}
	*/

	// Model checking functions

	public ModelCheckerResultRanged doTransientRanged(Model model, double t, explicit.StateValues initDistMin, explicit.StateValues initDistMax) throws PrismException
	{
		ParamModel ctmcRanged = (ParamModel) model;
		ModelCheckerResultRanged res = null;
		explicit.StateValues initDistMinNew = null, initDistMaxNew = null; //probs = null;

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
		res = computeTransientProbsRanged(ctmcRanged, t, initDistMinNew.getDoubleArray(), initDistMaxNew.getDoubleArray());
		//probs = StateValues.createFromDoubleArray(res.soln, ctmcRanged);

		return res;
	}

	public ModelCheckerResultRanged computeTransientProbsRanged(ParamModel ctmcRanged, double t, double initDistMin[], double initDistMax[]) throws PrismException
	{
		ModelCheckerResult min = null, max = null;
		int i, n, iters;
		double solnMin[], soln2Min[], sumMin[];
		double solnMax[], soln2Max[], sumMax[];
		double tmpsoln[];
		long timer;
		// Fox-Glynn stuff
		FoxGlynn fg;
		int left, right;
		double termCritParam, q, qt, acc, weights[], totalWeight;

		// Start bounded probabilistic reachability
		timer = System.currentTimeMillis();
		mainLog.println("\nStarting transient probability computation...");

		mainLog.println("\nComputing in, out, inout reactions sets...");
		ctmcRanged.computeInOutReactions();

		// Store num states
		n = ctmcRanged.getNumStates();

		// Get uniformisation rate; do Fox-Glynn
		// XXX the choice of `q` may also cause the rounding differences in results?
		q = 1.02 * ctmcRanged.maxSumLeaving();
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

		// Build (implicit) uniformised DTMC
		//dtmc = ctmcRanged.buildImplicitUniformisedDTMC(q);

		// Create solution vector(s)
		// For soln, we just use init (since we are free to modify this vector)
		solnMin = initDistMin;
		soln2Min = new double[n];
		sumMin = new double[n];
		solnMax = initDistMax;
		soln2Max = new double[n];
		sumMax = new double[n];

		// Initialise solution vectors
		// (don't need to do soln2 since will be immediately overwritten)
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

		// Start iterations
		iters = 1;
		while (iters <= right) {
			// Matrix-vector multiply
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
				}
			}

			iters++;
		}

		// Finished bounded probabilistic reachability
		timer = System.currentTimeMillis() - timer;
		mainLog.print("Transient probability computation");
		mainLog.println(" took " + iters + " iters and " + timer / 1000.0 + " seconds.");

		// Return results
		min = new ModelCheckerResult();
		min.soln = sumMin;
		min.lastSoln = soln2Min;
		min.numIters = iters;
		min.timeTaken = timer / 1000.0;
		min.timePre = 0.0;

		max = new ModelCheckerResult();
		max.soln = sumMax;
		max.lastSoln = soln2Max;
		max.numIters = iters;
		max.timeTaken = timer / 1000.0;
		max.timePre = 0.0;

		return new ModelCheckerResultRanged(min, max);
	}
}
