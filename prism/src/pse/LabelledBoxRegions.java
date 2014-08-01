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

import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;

import prism.PrismLog;

@SuppressWarnings("serial")
final class LabelledBoxRegions extends HashMap<BoxRegion, Collection<String>> 
{
	public LabelledBoxRegions()
	{
		super();
	}

	public LabelledBoxRegions(BoxRegion region, String label)
	{
		super();
		add(region, label);
	}

	public Collection<String> add(BoxRegion region)
	{
		return add(region, null);
	}

	public Collection<String> add(BoxRegion region, String label)
	{
		if (!containsKey(region)) {
			put(region, new LinkedList<String>());
		}
		Collection<String> labels = get(region);
		if (label != null) {
			labels.add(label);
		}
		return labels;
	}

	public boolean contains(BoxRegion region)
	{
		return containsKey(region);
	}

	public void print(PrismLog log)
	{
		if (isEmpty()) {
			log.println(" * [none]");
		} else {
			for (Map.Entry<BoxRegion, Collection<String>> entry : entrySet()) {
				log.print(" * " + entry.getKey());
				if (!entry.getValue().isEmpty()) {
					log.print(" " + entry.getValue());
				}
				log.println();
			}
		}
	}
}
