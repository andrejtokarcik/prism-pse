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

package explicit.ranged;

import java.util.BitSet;
import java.util.Iterator;

import explicit.ModelExplicit;
import prism.ModelType;
import prism.PrismException;
import prism.PrismLog;

public abstract class DTMCExplicitRanged extends ModelExplicit implements DTMCRanged
{
	@Override
	public ModelType getModelType()
	{
		return ModelType.DTMC;
	}

	@Override
	public String infoString()
	{
		throw new RuntimeException("Not implemented");
	}

	@Override
	public String infoStringTable()
	{
		throw new RuntimeException("Not implemented");
	}

	@Override
	public void exportToPrismLanguage(String filename) throws PrismException
	{
		throw new RuntimeException("Not implemented");
	}

	@Override
	public void exportToDotFile(String filename, BitSet mark) throws PrismException
	{
		throw new RuntimeException("Not implemented");
	}

	@Override
	public void exportToPrismExplicitTra(PrismLog out)
	{
		throw new RuntimeException("Not implemented");
	}

	@Override
	public void buildFromPrismExplicit(String filename) throws PrismException
	{
		throw new RuntimeException("Not implemented");
	}

	@Override
	public int getNumTransitions()
	{
		throw new RuntimeException("Not implemented");
	}

	@Override
	public Iterator<Integer> getSuccessorsIterator(final int s)
	{
		throw new RuntimeException("Not implemented");
	}
	
	@Override
	public boolean isSuccessor(int s1, int s2)
	{
		throw new RuntimeException("Not implemented");
	}

	@Override
	public boolean allSuccessorsInSet(int s, BitSet set)
	{
		throw new RuntimeException("Not implemented");
	}

	@Override
	public boolean someSuccessorsInSet(int s, BitSet set)
	{
		throw new RuntimeException("Not implemented");
	}

	@Override
	public void findDeadlocks(boolean fix) throws PrismException
	{
		throw new RuntimeException("Not implemented");
	}

	@Override
	public void checkForDeadlocks(BitSet except) throws PrismException
	{
		throw new RuntimeException("Not implemented");
	}
}
