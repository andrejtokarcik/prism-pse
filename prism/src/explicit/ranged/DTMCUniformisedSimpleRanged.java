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

package explicit.ranged;

import java.util.Map.Entry;

import explicit.Distribution;

public class DTMCUniformisedSimpleRanged extends DTMCExplicitRanged
{
	protected CTMCSimpleRanged ctmc;
	protected double q;

	public DTMCUniformisedSimpleRanged(CTMCSimpleRanged ctmc, double q)
	{
		this.ctmc = ctmc;
		this.numStates = ctmc.getNumStates();
		this.q = q;
	}

	public DTMCUniformisedSimpleRanged(CTMCSimpleRanged ctmc)
	{
		this(ctmc, ctmc.getDefaultUniformisationRate());
	}

	public void vmMultMax(double vect[], double result[])
	{
		int i, j;
		double prob;
		Distribution distrPredMax, distrSuccMin;
		
		// Initialise result
		for (i = 0; i < numStates; i++) {
			result[i] = vect[i];
		}

		for (i = 0; i < numStates; i++) {
			distrPredMax = ((CTMCSimpleRanged) ctmc).getTransitionsPredMax(i);
			for (Entry<Integer, Double> e : distrPredMax) {
				j = (Integer) e.getKey();
				prob = (Double) e.getValue();
				result[i] += (prob / q) * vect[j];
			}
			
			distrSuccMin = ((CTMCSimpleRanged) ctmc).getTransitionsSuccMin(i);		
			for (Entry<Integer, Double> e : distrSuccMin) {
				j = (Integer) e.getKey();
				prob = (Double) e.getValue();
				result[i] -= (prob / q) * vect[i];
			}
		}
	}
}
