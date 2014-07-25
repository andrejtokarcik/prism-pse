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

import java.util.Iterator;
import java.util.LinkedList;
import java.util.Collection;
import java.util.Map;
import java.util.HashMap;

import prism.PrismLog;

@SuppressWarnings("serial")
final class BoxRegionsToDecompose extends HashMap<BoxRegion, Collection<String>> implements Iterable<BoxRegion> 
{
	public BoxRegionsToDecompose()
	{
		super();
	}

	public BoxRegionsToDecompose(BoxRegion region, String explanation)
	{
		super();
		initialiseExplanations(region);
		get(region).add(explanation);
	}

	private void initialiseExplanations(BoxRegion region)
	{
		put(region, new LinkedList<String>());
	}

	public Collection<String> addRegion(BoxRegion region, String explanation)
	{
		if (!containsKey(region))
			initialiseExplanations(region);
		
		Collection<String> explanations = get(region);
		explanations.add("has " + explanation);
		return explanations;
	}

	public void print(PrismLog log)
	{
		assert size() > 0;
		log.println("There are " + size() + " regions to be decomposed:");
		for (Map.Entry<BoxRegion, Collection<String>> entry : entrySet()) {
			log.println(" * " + entry.getKey() + " " + entry.getValue());
		}
	}

	@Override
	public Iterator<BoxRegion> iterator() {
		return keySet().iterator();
	}
}
