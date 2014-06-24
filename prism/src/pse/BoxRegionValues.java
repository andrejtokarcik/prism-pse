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

import java.util.TreeMap;
import java.util.Iterator;
import java.util.Map.Entry;

import explicit.Model;
import explicit.StateValues;
import prism.Pair;

public class BoxRegionValues implements Iterable<Entry<BoxRegion, BoxRegionValues.StateValuesPair>>
{
	private Model model;
	private TreeMap<BoxRegion, BoxRegionValues.StateValuesPair> valuesPairs;

	public static class StateValuesPair extends Pair<StateValues, StateValues>
	{
		public StateValues min;
		public StateValues max;

		public StateValuesPair(StateValues min, StateValues max) {
			super(min, max);
			this.min = min;
			this.max = max;
		}

		public StateValues getMin() {
			return min;
		}

		public StateValues getMax() {
			return max;
		}
	}

	public BoxRegionValues(Model model)
	{
		this.model = model;
		valuesPairs = new TreeMap<BoxRegion, BoxRegionValues.StateValuesPair>();
	}

	public BoxRegionValues(Model model, BoxRegion region, StateValues minValues, StateValues maxValues)
	{
		this(model);
		add(region, minValues, maxValues);
	}

	public BoxRegionValues(Model model, BoxRegion region, double[] minValues, double[] maxValues)
	{
		this(model);
		add(region, minValues, maxValues);
	}

	public void add(BoxRegion region, StateValues minValues, StateValues maxValues)
	{
		valuesPairs.put(region, new StateValuesPair(minValues, maxValues));
	}

	public void add(BoxRegion region, double[] min, double[] max)
	{
		StateValues minValues = StateValues.createFromDoubleArray(min, model);
		StateValues maxValues = StateValues.createFromDoubleArray(max, model);
		add(region, minValues, maxValues);
	}

	public StateValues getMin(BoxRegion region)
	{
		return valuesPairs.get(region).getMin();
	}

	public StateValues getMax(BoxRegion region)
	{
		return valuesPairs.get(region).getMax();
	}

	public int getNumRegions()
	{
		return valuesPairs.keySet().size();
	}

	@Override
	public Iterator<Entry<BoxRegion, BoxRegionValues.StateValuesPair>> iterator()
	{
		return valuesPairs.entrySet().iterator();
	}
}
