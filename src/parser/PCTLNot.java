//==============================================================================
//	
//	Copyright (c) 2002-2004, Dave Parker, Andrew Hinton
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

package parser;

import apmc.*;
import simulator.*;

public class PCTLNot extends PCTLFormulaUnary
{
	// constructor
	
	public PCTLNot(PCTLFormula f)
	{
		super(f);
	}

	// convert to apmc data structures
	
	public int toApmc(Apmc apmc) throws ApmcException
	{
		return apmc.newUnaryOperand( apmc.NOT, operand.toApmc(apmc));
	}

	/**
	 *	Convert and build simulator data structures
	 */
	public int toSimulator(SimulatorEngine sim) throws SimulatorException
	{
		return SimulatorEngine.createNot(operand.toSimulator(sim));
	}

	// convert to string
	
	public String toString()
	{
		return "!" + operand;
	}
}

//------------------------------------------------------------------------------
