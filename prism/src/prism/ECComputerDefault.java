//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Carlos S. Bederian (Universidad Nacional de Cordoba)
//	* Dave Parker <d.a.parker@cs.bham.ac.uk> (University of Birmingham/Oxford)
//	* Mark Kattenbelt <mark.kattenbelt@comlab.ox.ac.uk> (University of Oxford)
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

package prism;

import java.util.List;
import java.util.Stack;
import java.util.Vector;

import jdd.JDD;
import jdd.JDDNode;
import jdd.JDDVars;

/**
 * Symbolic maximal end component computer for a nondeterministic model such as MDP.
 */
public class ECComputerDefault extends ECComputer
{
	/**
	 * Build (M)EC computer for a given model.
	 */
	public ECComputerDefault(PrismComponent parent, JDDNode reach, JDDNode trans, JDDNode trans01, JDDVars allDDRowVars, JDDVars allDDColVars,
			JDDVars allDDNondetVars) throws PrismException
	{
		super(parent, reach, trans, trans01, allDDRowVars, allDDColVars, allDDNondetVars);
	}

	// Methods for ECComputer interface

	@Override
	public void computeMECStates() throws PrismException
	{
		mecs = findEndComponents(reach, null);
	}

	@Override
	public void computeMECStates(JDDNode states) throws PrismException
	{
		mecs = findEndComponents(states, null);
	}

	@Override
	public void computeMECStates(JDDNode states, JDDNode filter) throws PrismException
	{
		mecs = findEndComponents(states, filter);
	}

	// Computation

	/**
	 * Find all accepting maximal end components (MECs) contained within {@code states},
	 * where acceptance is defined as those which intersect with {@code filter}.
	 * (If {@code filter} is null, the acceptance condition is trivially satisfied.)
	 * @param states BDD of the set of containing states
	 * @param filter BDD for the set of accepting states
	 * @return a vector of (referenced) BDDs representing the ECs
	 */
	private Vector<JDDNode> findEndComponents(JDDNode states, JDDNode filter) throws PrismException
	{
		Stack<JDDNode> candidates = new Stack<JDDNode>();
		JDD.Ref(states);
		candidates.push(states);
		Vector<JDDNode> ecs = new Vector<JDDNode>();
		SCCComputer sccComputer;

		while (!candidates.isEmpty()) {
			JDDNode candidate = candidates.pop();
			// Compute the stable set
			JDD.Ref(candidate);
			JDDNode stableSet = findMaximalStableSet(candidate);
			// Drop empty sets
			if (stableSet.equals(JDD.ZERO)) {
				JDD.Deref(stableSet);
				JDD.Deref(candidate);
				continue;
			}

			if (stableSet.equals(candidate) && JDD.GetNumMinterms(stableSet, allDDRowVars.n()) == 1) {
				ecs.add(candidate);
				JDD.Deref(stableSet);
				continue;
			}

			// Filter bad transitions
			JDD.Ref(stableSet);
			JDDNode stableSetTrans = maxStableSetTrans(stableSet);

			// now find the maximal SCCs in (stableSet, stableSetTrans)
			List<JDDNode> sccs;
			sccComputer = SCCComputer.createSCCComputer(this, stableSet, stableSetTrans, allDDRowVars, allDDColVars);
			if (filter != null)
				sccComputer.computeSCCs(filter);
			else
				sccComputer.computeSCCs();
			JDD.Deref(stableSet);
			JDD.Deref(stableSetTrans);
			sccs = sccComputer.getSCCs();
			JDD.Deref(sccComputer.getNotInSCCs());
			if (sccs.size() > 0) {
				if (sccs.size() > 1 || !sccs.get(0).equals(candidate)) {
					candidates.addAll(sccs);
					JDD.Deref(candidate);
				} else {
					ecs.add(candidate);
					JDD.Deref(sccs.get(0));
				}
			} else
				JDD.Deref(candidate);
		}
		return ecs;
	}

	/**
	 * Returns a stable set of states contained in candidateStates
	 * 
	 * @param candidateStates set of candidate states S x H_i (dereferenced after calling this function)
	 * @return a referenced BDD with the maximal stable set in c
	 */
	private JDDNode findMaximalStableSet(JDDNode candidateStates)
	{
		JDDNode old = JDD.Constant(0);
		JDDNode current = candidateStates;

		while (!current.equals(old)) {
			JDD.Deref(old);
			JDD.Ref(current);
			old = current;

			JDD.Ref(current);
			JDD.Ref(trans);
			// Select transitions starting in current
			JDDNode currTrans = JDD.Apply(JDD.TIMES, trans, current);
			// Select transitions starting in current and ending in current
			JDDNode tmp = JDD.PermuteVariables(current, allDDRowVars, allDDColVars);
			tmp = JDD.Apply(JDD.TIMES, currTrans, tmp);
			// Sum all successor probabilities for each (state, action) tuple
			tmp = JDD.SumAbstract(tmp, allDDColVars);
			// If the sum for a (state,action) tuple is 1,
			// there is an action that remains in the stable set with prob 1
			tmp = JDD.GreaterThan(tmp, 1 - sumRoundOff);
			// Without fairness, we just need one action per state
			current = JDD.ThereExists(tmp, allDDNondetVars);
		}
		JDD.Deref(old);
		return current;

	}

	/**
	 * Returns the transition relation of a stable set
	 * 
	 * @param b BDD of a stable set (dereferenced after calling this function)
	 * @return referenced BDD of the transition relation restricted to the stable set
	 */
	private JDDNode maxStableSetTrans(JDDNode b)
	{

		JDD.Ref(b);
		JDD.Ref(trans);
		// Select transitions starting in b
		JDDNode currTrans = JDD.Apply(JDD.TIMES, trans, b);
		JDDNode mask = JDD.PermuteVariables(b, allDDRowVars, allDDColVars);
		// Select transitions starting in current and ending in current
		mask = JDD.Apply(JDD.TIMES, currTrans, mask);
		// Sum all successor probabilities for each (state, action) tuple
		mask = JDD.SumAbstract(mask, allDDColVars);
		// If the sum for a (state,action) tuple is 1,
		// there is an action that remains in the stable set with prob 1
		mask = JDD.GreaterThan(mask, 1 - sumRoundOff);
		// select the transitions starting in these tuples
		JDD.Ref(trans01);
		JDDNode stableTrans01 = JDD.And(trans01, mask);
		// Abstract over actions
		return JDD.ThereExists(stableTrans01, allDDNondetVars);
	}
}
