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
	private Values boundsMid;

	public BoxRegion(Values boundsLower, Values boundsUpper)
	{
		assert boundsLower.compareTo(boundsUpper) <= 0;

		this.boundsLower = boundsLower;
		this.boundsUpper = boundsUpper;
		computeMidBounds();
	}

	private void computeMidBounds()
	{
		boundsMid = new Values();
		for (int i = 0; i < boundsLower.getNumValues(); i++) {
			double lowerValue = (Double) boundsLower.getValue(i);
			double upperValue = (Double) boundsUpper.getValue(i);
			boundsMid.addValue(boundsLower.getName(i), lowerValue + 0.5 * (upperValue - lowerValue));
		}
	}

	public Values getLowerBounds()
	{
		return boundsLower;
	}

	public Values getUpperBounds()
	{
		return boundsUpper;
	}

	public BoxRegion getLowerHalf()
	{
		return new BoxRegion(boundsLower, boundsMid);
	}

	public BoxRegion getUpperHalf()
	{
		return new BoxRegion(boundsMid, boundsUpper);
	}

	public int compareTo(BoxRegion r)
	{
		return boundsLower.compareTo(r.boundsLower);
	}

	@Override
	public String toString() {
		StringBuilder builder = new StringBuilder();
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
