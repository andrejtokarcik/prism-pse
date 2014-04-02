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

import explicit.Distribution;

public class CTMCSimpleRanged extends DTMCSimpleRanged implements CTMCRanged
{
	// Constructors

	/**
	 * Constructor: empty CTMC.
	 */
	public CTMCSimpleRanged()
	{
		super();
	}

	/**
	 * Constructor: new CTMC with fixed number of states.
	 */
	public CTMCSimpleRanged(int numStates)
	{
		super(numStates);
	}

	// Accessors (for CTMC)
	
	@Override
	public double getMaxExitRate()
	{
		int i;
		double d, max = Double.NEGATIVE_INFINITY;
		for (i = 0; i < numStates; i++) {
			d = transSuccMax.get(i).sum();
			if (d > max)
				max = d;
		}
		return max;
	}

	@Override
	public double getDefaultUniformisationRate()
	{
		return 1.02 * getMaxExitRate(); 
	}

	@Override
	public void uniformise(double q)
	{
		Distribution distr;
		int i;
		for (i = 0; i < numStates; i++) {
			// XXX correct?
			distr = transSuccMin.get(i);
			distr.set(i, q - distr.sumAllBut(i));
			distr = transSuccMax.get(i);
			distr.set(i, q - distr.sumAllBut(i));
			distr = transPredMin.get(i);
			distr.set(i, q - distr.sumAllBut(i));
			distr = transPredMax.get(i);
			distr.set(i, q - distr.sumAllBut(i));
		}
	}

	public DTMCRanged buildImplicitUniformisedDTMC(double q)
	{
		return new DTMCUniformisedSimpleRanged(this, q);
	}
}
