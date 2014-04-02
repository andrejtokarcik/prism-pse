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

import java.util.ArrayList;
import java.util.List;

import explicit.ModelSimple;
import explicit.Distribution;

public class DTMCSimpleRanged extends DTMCExplicitRanged implements ModelSimple
{
	// Transition matrices (distribution lists)
	protected List<Distribution> transSuccMin;
	protected List<Distribution> transSuccMax;
	protected List<Distribution> transPredMin;
	protected List<Distribution> transPredMax;

	// Other statistics
	protected int numTransitionsMin;
	protected int numTransitionsMax;

	// Constructors

	/**
	 * Constructor: empty DTMC.
	 */
	public DTMCSimpleRanged()
	{
		initialise(0);
	}

	/**
	 * Constructor: new DTMC with fixed number of states.
	 */
	public DTMCSimpleRanged(int numStates)
	{
		initialise(numStates);
	}

	// Mutators (for ModelSimple)

	@Override
	public void initialise(int numStates)
	{
		super.initialise(numStates);
		transSuccMin = new ArrayList<Distribution>(numStates);
		transSuccMax = new ArrayList<Distribution>(numStates);
		transPredMin = new ArrayList<Distribution>(numStates);
		transPredMax = new ArrayList<Distribution>(numStates);
		for (int i = 0; i < numStates; i++) {
			transSuccMin.add(new Distribution());
			transSuccMax.add(new Distribution());
			transPredMin.add(new Distribution());
			transPredMax.add(new Distribution());
		}
	}
	
	@Override
	public void clearState(int i)
	{
		// Do nothing if state does not exist
		if (i >= numStates || i < 0)
			return;
		// Clear data structures and update stats
		numTransitionsMin -= transSuccMin.get(i).size();
		numTransitionsMax -= transSuccMax.get(i).size();
		transSuccMin.get(i).clear();
		transSuccMax.get(i).clear();
		transPredMin.get(i).clear();
		transPredMax.get(i).clear();
	}

	@Override
	public int addState()
	{
		addStates(1);
		return numStates - 1;
	}

	@Override
	public void addStates(int numToAdd)
	{
		for (int i = 0; i < numToAdd; i++) {
			transSuccMin.add(new Distribution());
			transSuccMax.add(new Distribution());
			transPredMin.add(new Distribution());
			transPredMin.add(new Distribution());
			numStates++;
		}
	}

	// Mutators (other)

	/**
	 * Set the probability for a transition. 
	 */
	public void setProbabilityMin(int i, int j, double probNew)
	{
		double probCur = transSuccMin.get(i).get(j);
		if (probCur > 0.0)
			numTransitionsMin--;
		if (probNew > 0.0)
			numTransitionsMin++;
		transSuccMin.get(i).set(j, probNew);
		transPredMin.get(j).set(i, probNew);
	}

	public void setProbabilityMax(int i, int j, double probNew)
	{
		double probCur = transSuccMax.get(i).get(j);
		if (probCur > 0.0)
			numTransitionsMax--;
		if (probNew > 0.0)
			numTransitionsMax++;
		transSuccMax.get(i).set(j, probNew);
		transPredMax.get(j).set(i, probNew);
	}

	/**
	 * Add to the probability for a transition. 
	 */
	public void addToProbabilityMin(int i, int j, double prob)
	{
		if (!transSuccMin.get(i).add(j, prob) && !transPredMin.get(j).add(i, prob)) {
			if (prob > 0.0)
				numTransitionsMin++;
		}
	}
	public void addToProbabilityMax(int i, int j, double prob)
	{
		if (!transSuccMax.get(i).add(j, prob) && !transPredMax.get(j).add(i, prob)) {
			if (prob > 0.0)
				numTransitionsMax++;
		}
	}


	// Accessors (for DTMCRanged)

	@Override
	public void vmMultMin(double vect[], double result[])
	{
		// TODO?
		throw new RuntimeException("Not implemented");
	}

	@Override
	public void vmMultMax(double vect[], double result[])
	{
		// TODO?
		throw new RuntimeException("Not implemented");
	}

	// Accessors (other)

	/**
	 * Get the transitions (a distribution) for state s.
	 */
	public Distribution getTransitionsSuccMin(int s)
	{
		return transSuccMin.get(s);
	}
	public Distribution getTransitionsSuccMax(int s)
	{
		return transSuccMax.get(s);
	}
	public Distribution getTransitionsPredMin(int s)
	{
		return transPredMin.get(s);
	}
	public Distribution getTransitionsPredMax(int s)
	{
		return transPredMax.get(s);
	}
}
