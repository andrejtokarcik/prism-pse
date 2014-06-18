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

public class BoxRegionValues implements Iterable<Entry<BoxRegion, Pair<StateValues, StateValues>>>
{
	private Model model;
	private TreeMap<BoxRegion, Pair<StateValues, StateValues>> valuesPairs;

	public BoxRegionValues(Model model)
	{
		this.model = model;
		valuesPairs = new TreeMap<BoxRegion, Pair<StateValues, StateValues>>();
	}

	public void add(BoxRegion region, StateValues minValues, StateValues maxValues)
	{
		valuesPairs.put(region, new Pair<StateValues, StateValues>(minValues, maxValues));
	}
	
	public void add(BoxRegion region, double[] min, double[] max)
	{
		StateValues minValues = StateValues.createFromDoubleArray(min, model);
		StateValues maxValues = StateValues.createFromDoubleArray(max, model);
		add(region, minValues, maxValues);
	}

	public StateValues getMin(BoxRegion region)
	{
		return valuesPairs.get(region).first;
	}

	public StateValues getMax(BoxRegion region)
	{
		return valuesPairs.get(region).second;
	}

	public int getNumRegions()
	{
		return valuesPairs.keySet().size();
	}

	@Override
	public Iterator<Entry<BoxRegion, Pair<StateValues, StateValues>>> iterator()
	{
		return valuesPairs.entrySet().iterator();
	}
}
