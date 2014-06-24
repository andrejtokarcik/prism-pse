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

import prism.Pair;

public class BoxRegion extends Pair.ComparablePair<Double, Double>
{
	public static BoxRegion completeSpace = new BoxRegion(0, 1);

	public BoxRegion(double min, double max)
	{
		super(min, max);
		assert 0 <= min && min <= 1;
		assert 0 <= max && max <= 1;
		assert min <= max;
	}

	public double getMinCoeff()
	{
		return first;
	}

	public double getMaxCoeff()
	{
		return second;
	}
	
	public BoxRegion lowerHalf()
	{
		return new BoxRegion(first, first + 0.5 * (second - first));
	}
	
	public BoxRegion upperHalf()
	{
		return new BoxRegion(first + 0.5 * (second - first), second);
	}
	
	@Override
	public String toString() {
		// TODO: Multiply by actual param bounds (e.g. supplied via BoxRegionFactory as in param?)
		return "Region [" + first + ", " + second + "]";
	}
}
