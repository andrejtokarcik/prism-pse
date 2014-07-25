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

public final class SimpleDecompositionProcedure extends DecompositionProcedure {
	private double accuracy;
	private int numStates;

	public static final class NoDecomposing extends DecompositionProcedure {
		private static final NoDecomposing INSTANCE = new NoDecomposing();

		private NoDecomposing() {}

		public static NoDecomposing getInstance()
		{
			return INSTANCE;
		}

		@Override
		public void examineSingleIteration(BoxRegion region, double probsMin[], double probsMax[]) throws DecompositionNeeded {}
	}
	
	public SimpleDecompositionProcedure(double accuracy, int numStates)
	{
		this.accuracy = accuracy;
		this.numStates = numStates;
	}

	@Override
	public void examineSingleIteration(BoxRegion region, double probsMin[], double probsMax[]) throws DecompositionNeeded
	{
		for (int state = 0; state < numStates; state++) {
			if (probsMax[state] - probsMin[state] > accuracy)
				throw new DecompositionNeeded(region, "significant inaccuracy");
		}
	}
}
