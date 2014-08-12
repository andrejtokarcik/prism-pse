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

import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import parser.Values;
import explicit.Utils;

final class BoxRegion implements Comparable<BoxRegion>
{
	private Values lowerBounds;
	private Values upperBounds;

	private double volume = 0.0;

	public BoxRegion(Values boundsLower, Values boundsUpper)
	{
		assert boundsLower.compareTo(boundsUpper) <= 0;
		this.lowerBounds = boundsLower;
		this.upperBounds = boundsUpper;
	}


	public Values getLowerBounds()
	{
		return lowerBounds;
	}

	public Values getUpperBounds()
	{
		return upperBounds;
	}

	private Values computeMidBounds()
	{
		Values midBounds = new Values();
		for (int i = 0; i < lowerBounds.getNumValues(); i++) {
			double lowerValue = (Double) lowerBounds.getValue(i);
			double upperValue = (Double) upperBounds.getValue(i);
			midBounds.addValue(lowerBounds.getName(i), lowerValue + 0.5 * (upperValue - lowerValue));
		}
		return midBounds;
	}

	public Set<BoxRegion> decompose()
	{
		Set<BoxRegion> subregions = new HashSet<BoxRegion>();
		Values midBounds = computeMidBounds();
		Set<Integer> allIndices = new HashSet<Integer>();
		for (int i = 0; i < midBounds.getNumValues(); i++) {
			allIndices.add(i);
		}

		for (Set<Integer> indices : Utils.powerSet(allIndices)) {
			Values newLowerBounds = new Values();
			Values newUpperBounds = new Values();
			for (int i = 0; i < midBounds.getNumValues(); i++) {
				String name = midBounds.getName(i);
				double midValue = (Double) midBounds.getValue(i);
				if (indices.contains(i)) {
					newLowerBounds.addValue(name, (Double) lowerBounds.getValue(i));
					newUpperBounds.addValue(name, midValue);
				} else {
					newLowerBounds.addValue(name, midValue);
					newUpperBounds.addValue(name, (Double) upperBounds.getValue(i));
				}
			}
			subregions.add(new BoxRegion(newLowerBounds, newUpperBounds));
		}
		return subregions;
	}

	public double volume()
	{
		if (volume > 0.0)
			return volume;

		volume = 1.0;
		for (int i = 0; i < lowerBounds.getNumValues(); i++) {
			double lowerValue = (Double) lowerBounds.getValue(i);
			double upperValue = (Double) upperBounds.getValue(i);
			if (lowerValue != upperValue) {
				volume *= upperValue - lowerValue;
			}
		}
		return volume;
	}

	public Set<Point> generateSamplePoints()
	{
		return generateSamplePoints(2);
	}

	public Set<Point> generateSamplePoints(int numSamples)
	{
		Set<Point> samples = new HashSet<Point>();
		Random r = new Random();
		while (samples.size() != numSamples) {
			Values dimensions = new Values();
			for (int i = 0; i < lowerBounds.getNumValues(); i++) {
				double lowerValue = (Double) lowerBounds.getValue(i);
				double upperValue = (Double) upperBounds.getValue(i);
				double randomValue = lowerValue + r.nextDouble() * (upperValue - lowerValue);
				dimensions.addValue(lowerBounds.getName(i), randomValue);
			}
			samples.add(new Point(dimensions));
		}
		return samples;
	}

	@Override
	public int compareTo(BoxRegion r)
	{
		int lowerRes = lowerBounds.compareTo(r.lowerBounds);
		int upperRes = upperBounds.compareTo(r.upperBounds);
		if (lowerRes == upperRes)
			return lowerRes;

		int min = Math.min(lowerRes, upperRes);
		int max = Math.max(lowerRes, upperRes);
		if (min == 0)
			return max;
		if (max == 0)
			return min;
		return lowerRes;
	}

	@Override
	public int hashCode()
	{
		final int prime = 31;
		int result = 1;
		result = prime * result
				+ ((lowerBounds == null) ? 0 : lowerBounds.hashCode());
		result = prime * result
				+ ((upperBounds == null) ? 0 : upperBounds.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj)
	{
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		BoxRegion other = (BoxRegion) obj;
		if (lowerBounds == null) {
			if (other.lowerBounds != null)
				return false;
		} else if (!lowerBounds.equals(other.lowerBounds))
			return false;
		if (upperBounds == null) {
			if (other.upperBounds != null)
				return false;
		} else if (!upperBounds.equals(other.upperBounds))
			return false;
		return true;
	}

	@Override
	public String toString()
	{
		StringBuilder builder = new StringBuilder();
		for (int i = 0; i < lowerBounds.getNumValues(); i++) {
			if (i != 0) builder.append(",");
			builder.append(lowerBounds.getName(i));
			builder.append("=");
			builder.append((Double) lowerBounds.getValue(i));
			builder.append(":");
			builder.append((Double) upperBounds.getValue(i));
		}
		return builder.toString();
	}
}
