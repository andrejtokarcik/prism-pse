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

import java.util.BitSet;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.TreeMap;

import explicit.Model;
import explicit.StateValues;
import parser.type.TypeDouble;
import prism.Pair;
import prism.PrismException;
import prism.PrismLog;

@SuppressWarnings("serial")
public class BoxRegionValues extends TreeMap<BoxRegion, BoxRegionValues.StateValuesPair> implements Iterable<Entry<BoxRegion, BoxRegionValues.StateValuesPair>>
{
	private Model model;

	public static class StateValuesPair extends Pair<StateValues, StateValues>
	{
		public StateValuesPair(StateValues min, StateValues max)
		{
			super(min, max);
		}

		public StateValues getMin()
		{
			return first;
		}

		public StateValues getMax()
		{
			return second;
		}

		public void swap() {
			StateValues tmpFirst = first; 
			first = second;
			second = tmpFirst;
		}
	}

	public BoxRegionValues(Model model)
	{
		this.model = model;
	}

	public BoxRegionValues(Model model, BoxRegion region, StateValues minValues, StateValues maxValues)
	{
		this(model);
		put(region, minValues, maxValues);
	}

	public BoxRegionValues(Model model, BoxRegion region, double[] minValues, double[] maxValues)
	{
		this(model);
		put(region, minValues, maxValues);
	}

	public static BoxRegionValues createWithOnes(Model model, BoxRegion region) throws PrismException
	{
		StateValues ones = new StateValues(TypeDouble.getInstance(), new Double(1.0), model);
		return new BoxRegionValues(model, region, ones, ones);
	}
	
	public StateValuesPair put(BoxRegion region, StateValues minValues, StateValues maxValues)
	{
		return put(region, new StateValuesPair(minValues, maxValues));
	}

	public StateValuesPair put(BoxRegion region, double[] min, double[] max)
	{
		StateValues minValues = StateValues.createFromDoubleArray(min, model);
		StateValues maxValues = StateValues.createFromDoubleArray(max, model);
		return put(region, minValues, maxValues);
	}

	public StateValuesPair put(BoxRegion region, BitSet min, BitSet max)
	{
		StateValues minValues = StateValues.createFromBitSet(min, model);
		StateValues maxValues = StateValues.createFromBitSet(max, model);
		return put(region, minValues, maxValues);
	}

	public StateValuesPair remove(BoxRegion region)
	{
		return remove(region);
	}

	public void decomposeRegion(BoxRegion region)
	{
		StateValuesPair oldValuesPair = remove(region);
		for (BoxRegion subregion : region.decompose()) {
			put(subregion, oldValuesPair);
		}
	}

	public boolean hasRegion(BoxRegion region)
	{
		return containsKey(region);
	}

	public StateValues getMin(BoxRegion region)
	{
		return get(region).getMin();
	}

	public StateValues getMax(BoxRegion region)
	{
		return get(region).getMax();
	}

	public int getNumRegions()
	{
		return keySet().size();
	}

	@Override
	public Iterator<Entry<BoxRegion, BoxRegionValues.StateValuesPair>> iterator()
	{
		return entrySet().iterator();
	}

	public void print(PrismLog log)
	{
		for (Entry<BoxRegion, BoxRegionValues.StateValuesPair> entry : entrySet()) {
			log.println("\n== Region " + entry.getKey() + " ==");
			log.println("\n=== Minimised state values ===\n");
			entry.getValue().getMin().print(log);
			log.println("\n=== Maximised state values ===\n");
			entry.getValue().getMax().print(log);
		}
	}

	// TODO to facilitate garbage-collecting
	//public void clear()
}
