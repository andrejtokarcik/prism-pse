//==============================================================================
//	
//	Copyright (c) 2005, Andrew Hinton
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

package settings;

/**
 *
 * @author  Andrew Hinton
 */
public abstract class FontColorConstraint implements SettingConstraint
{
	
	/** Creates a new instance of FontColorConstraint */
	public FontColorConstraint()
	{
	}
	
	public void checkValue(Object value) throws SettingException
	{
		if(value instanceof FontColorPair)
		{
			checkValueFontColor((FontColorPair)value);
		}
		else
		{
			throw new SettingException("Invalid type for property, should be a Font Colour pair.");
		}
	}
	
	public abstract void checkValueFontColor(FontColorPair col) throws SettingException;
	
}
