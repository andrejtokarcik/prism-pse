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

import parser.Values;

final class BoxRegion implements Comparable<BoxRegion>
{
	private Values boundsLower;
	private Values boundsUpper;

	public BoxRegion(Values boundsLower, Values boundsUpper)
	{
		assert boundsLower.compareTo(boundsUpper) <= 0;

		this.boundsLower = boundsLower;
		this.boundsUpper = boundsUpper;
	}

	public Values getLowerBounds()
	{
		return boundsLower;
	}

	public Values getUpperBounds()
	{
		return boundsUpper;
	}

	public BoxRegion lowerHalf()
	{
		Values newUpper = new Values();
		for (int i = 0; i < boundsLower.getNumValues(); i++) {
			double lowerValue = (Double) boundsLower.getValue(i);
			double upperValue = (Double) boundsUpper.getValue(i);
			newUpper.addValue(boundsLower.getName(i),
					lowerValue + 0.5 * (upperValue - lowerValue));
		}
		return new BoxRegion(boundsLower, newUpper);
	}

	public BoxRegion upperHalf()
	{
		Values newLower = new Values();
		for (int i = 0; i < boundsLower.getNumValues(); i++) {
			double lowerParamValue = (Double) boundsLower.getValue(i);
			double upperParamValue = (Double) boundsUpper.getValue(i);
			newLower.addValue(boundsLower.getName(i),
					lowerParamValue + 0.5 * (upperParamValue - lowerParamValue));
		}
		return new BoxRegion(newLower, boundsUpper);
	}

	public int compareTo(BoxRegion r)
	{
		return boundsLower.compareTo(r.boundsLower);
	}

	@Override
	public String toString() {
		StringBuilder builder = new StringBuilder("Region ");
		for (int i = 0; i < boundsLower.getNumValues(); i++) {
			if (i != 0) builder.append(",");
			builder.append(boundsLower.getName(i));
			builder.append("=");
			builder.append((Double) boundsLower.getValue(i));
			builder.append(":");
			builder.append((Double) boundsUpper.getValue(i));
		}
		return builder.toString();
	}
}
