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

import java.util.Map.Entry;

import prism.PrismException;

abstract class MinSynthesis extends AbstractMinMaxSynthesis
{
	public MinSynthesis(double probTolerance, int initState)
	{
		super(probTolerance, initState);
		captionForOptimising = "minimising";
	}

	protected abstract double getMinimalUpperBound(BoxRegionValues regionValues) throws PrismException;

	@Override
	protected void determineOptimisingRegions(BoxRegionValues regionValues) throws PrismException
	{
		// Determine the minimal upper bound
		double minimalUpperBound = getMinimalUpperBound(regionValues);
		demarcationProbBounds.add(minimalUpperBound);
		
		// Determine the (non-)minimising regions
		regionsOptimising.clear();
		for (Entry<BoxRegion, BoxRegionValues.StateValuesPair> entry : regionValues) {
			if (regionsNonoptimising.contains(entry.getKey()))
				continue;

			if ((Double) entry.getValue().getMin().getValue(initState) > minimalUpperBound)
				regionsNonoptimising.add(entry.getKey());
			else
				regionsOptimising.add(entry.getKey());
		}
	}

	@Override
	public String toString()
	{
		return "min synthesis";
	}
}
