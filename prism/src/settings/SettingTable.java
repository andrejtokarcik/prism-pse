//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Andrew Hinton <ug60axh@cs.bham.uc.uk> (University of Birmingham)
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

package settings;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
import javax.swing.table.*;
import javax.swing.border.*;
import javax.swing.text.*;

/**
 *
 * @author  ug60axh
 */
public class SettingTable extends JPanel implements ListSelectionListener, TableModelListener, ItemListener, SettingDisplay
{
	private Component parent;
	private SettingTableModel theModel;
	
	private int lineWidth;
	
	private boolean shouldRemove;
	/** Creates new form PropertyTable */
	public SettingTable(Component parent)
	{
		super();
		this.parent = parent;
		theModel = new SettingTableModel();
		initComponents();
		
		theModel.setJTable(theTable);
		theModel.addTableModelListener(this);
		lineWidth = theTable.getRowHeight();
		theTable.setModel(theModel);
		theTable.setRowSelectionAllowed(false);
		theTable.setColumnSelectionAllowed(false);
		//theTable.setCellSelectionEnabled(true);
		theTable.getSelectionModel().addListSelectionListener(this);
		theTable.setAutoResizeMode(JTable.AUTO_RESIZE_LAST_COLUMN);
		theCombo.setModel(theModel.getComboModel());
		theCombo.addItemListener(this);
		theTable.setSurrendersFocusOnKeystroke(true);
		
		
		theTable.getColumnModel().getColumn(0).setMinWidth(30);
		
		theTable.setRequestFocusEnabled(false);
		
		TableColumn column = theTable.getColumnModel().getColumn(1);
		column.setCellRenderer(new SettingCellRenderer());
		column.setCellEditor(new SettingCellEditor());
		
		TableResizer tr = new TableResizer(theTable);
		
		theTable.addMouseListener(tr);
		theTable.addMouseMotionListener(tr);
		
		doChoiceBox();
		
		commentLabel.setFont(new Font("serif", Font.BOLD, 12));
		
		shouldRemove = true;
		
	}
	
	public void setOwners(ArrayList owners)
	{
		theModel.setOwners(owners);
	}
	
	public void refreshGroupNames()
	{
		theModel.refreshGroupNames();
		
		//doChoiceBox();
		//repaint();
	}
	
	public void setNameColumnWidth(int width)
	{
		//theTable.getColumnModel().getColumn(0).setMinWidth(width);
		//theTable.getColumnModel().getColumn(0).setMaxWidth(width);
		//theTable.getColumnModel().getColumn(0).setPreferredWidth(width);
		//theTable.getColumnModel().getColumn(0).setMaxWidth(width);
		//theTable.repaint();
	}
	
	private void doChoiceBox()
	{
		//System.out.println("calling doChoiceBox()");
		if(theModel.getNumGroups() == 0)
		{
			//System.out.println("0 groups");
			topPanel.removeAll();
			JLabel lab = new JLabel("");
			topPanel.setLayout(new BorderLayout());
			topPanel.add(lab, BorderLayout.CENTER);
		}
		else if(theModel.getNumGroups() == 1)
		{
			//System.out.println("1 group");
			topPanel.removeAll();
			JLabel lab = new JLabel(theCombo.getModel().getElementAt(0).toString());
			//System.out.println(theCombo.getModel().getElementAt(0).toString());
			topPanel.setLayout(new BorderLayout());
			topPanel.add(lab, BorderLayout.CENTER);
		}
		else
		{
			//System.out.println("2 group");
			topPanel.removeAll();
			topPanel.setLayout(new BorderLayout());
			topPanel.add(theCombo, BorderLayout.CENTER);
		}
		this.revalidate();
	}
	
	public void stopEditing()
	{
		if(theTable.getCellEditor() != null) theTable.removeEditor();
		//        if(ce != null)ce.stopEditing();
	}
	
	//    public void setCurrEditor(SettingEditor ce)
	//    {
	//		TableColumn column = theTable.getColumnModel().getColumn(1);
	//		column.setCellEditor(ce);
	//		//if(this.ce != ce) System.out.println("THE CURREDITOR HAS CHANGED");
	//
	//        this.ce = ce;
	//    }
	//
	//    private SettingEditor ce;
	
	/** This method is called from within the constructor to
	 * initialize the form.
	 * WARNING: Do NOT modify this code. The content of this method is
	 * always regenerated by the Form Editor.
	 */
    private void initComponents()//GEN-BEGIN:initComponents
    {
        javax.swing.JPanel jPanel1;
        javax.swing.JPanel jPanel3;
        javax.swing.JSplitPane jSplitPane1;

        jPanel1 = new javax.swing.JPanel();
        jSplitPane1 = new javax.swing.JSplitPane();
        jScrollPane1 = new javax.swing.JScrollPane();
        jScrollPane1.getViewport().setBackground(Color.white);
        theTable = new JTable()
        {
            public void editingStopped(ChangeEvent e)
            {
                // Take in the new value
                TableCellEditor editor = getCellEditor();
                if (editor != null)
                {
                    Object value = editor.getCellEditorValue();
                    setValueAt(value, editingRow, editingColumn);
                    if(shouldRemove)
                    {
                        removeEditor();
                        getSelectionModel().setSelectionInterval(editingRow, editingRow);
                        getColumnModel().getSelectionModel().setSelectionInterval(editingColumn, editingColumn);
                    }
                    else
                    {
                        getSelectionModel().setSelectionInterval(editingRow, editingRow);
                        getColumnModel().getSelectionModel().setSelectionInterval(editingColumn, editingColumn);

                    }
                    shouldRemove = false;
                }
            }

            //This method is a fix from http://www.codeguru.com/java/articles/180.shtml by Zafir Anjum, cheers!
            //this is required because there is a bug in JTable where the
            //just saying tableScroll.setColumnHeader(null); does not work as it should
            //Unfortunately, it overrides a deprecated API, so let's hope they
            //sort it out by Java 5.0, nice one Sun...
            protected void configureEnclosingScrollPane()
            {
                Container p = getParent();
                if (p instanceof JViewport)
                {
                    Container gp = p.getParent();
                    if (gp instanceof JScrollPane)
                    {
                        JScrollPane scrollPane = (JScrollPane)gp;
                        // Make certain we are the viewPort's view and not, for
                        // example, the rowHeaderView of the scrollPane -
                        // an implementor of fixed columns might do this.
                        JViewport viewport = scrollPane.getViewport();
                        if (viewport == null || viewport.getView() != this)
                        {
                            return;
                        }
                        //                scrollPane.setColumnHeaderView(getTableHeader());
                        scrollPane.getViewport().setBackingStoreEnabled(true);
                        scrollPane.setBorder(UIManager.getBorder("Table.scrollPaneBorder"));
                    }
                }
            }

        };
        theTable.setModel(theModel);
        theTable.setSelectionMode(DefaultListSelectionModel.SINGLE_SELECTION);
        theTable.setRowSelectionAllowed(false);
        theTable.setColumnSelectionAllowed(false);
        theTable.setCellSelectionEnabled(true);
        jPanel3 = new javax.swing.JPanel();
        commentText = new javax.swing.JTextArea();
        commentLabel = new javax.swing.JLabel();
        topPanel = new javax.swing.JPanel();
        theCombo = new javax.swing.JComboBox();

        setLayout(new java.awt.BorderLayout());

        jPanel1.setLayout(new java.awt.BorderLayout());

        jSplitPane1.setBackground(new java.awt.Color(255, 255, 255));
        jSplitPane1.setBorder(null);
        jSplitPane1.setDividerSize(3);
        jSplitPane1.setOrientation(javax.swing.JSplitPane.VERTICAL_SPLIT);
        jSplitPane1.setResizeWeight(1.0);
        jSplitPane1.setOneTouchExpandable(true);
        jScrollPane1.setBackground(new java.awt.Color(255, 255, 255));
        jScrollPane1.setBorder(new javax.swing.border.LineBorder(java.awt.SystemColor.textInactiveText));
        jScrollPane1.setHorizontalScrollBarPolicy(javax.swing.JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        theTable.setModel(new javax.swing.table.DefaultTableModel(
            new Object [][]
            {
                {null, null, null, null},
                {null, null, null, null},
                {null, null, null, null},
                {null, null, null, null}
            },
            new String []
            {
                "Title 1", "Title 2", "Title 3", "Title 4"
            }
        ));
        theTable.setDoubleBuffered(true);
        theTable.setGridColor(new java.awt.Color(198, 197, 197));
        jScrollPane1.setViewportView(theTable);

        jSplitPane1.setLeftComponent(jScrollPane1);

        jPanel3.setLayout(new java.awt.BorderLayout());

        jPanel3.setBorder(new javax.swing.border.LineBorder(java.awt.SystemColor.inactiveCaption));
        jPanel3.setMinimumSize(new java.awt.Dimension(10, 75));
        jPanel3.setPreferredSize(new java.awt.Dimension(100, 75));
        commentText.setBackground(javax.swing.UIManager.getDefaults().getColor("Panel.background"));
        commentText.setColumns(1);
        commentText.setEditable(false);
        commentText.setLineWrap(true);
        commentText.setWrapStyleWord(true);
        commentText.setBorder(null);
        commentText.setDoubleBuffered(true);
        commentText.setFocusable(false);
        commentText.setMinimumSize(new java.awt.Dimension(100, 75));
        commentText.setPreferredSize(new java.awt.Dimension(100, 75));
        jPanel3.add(commentText, java.awt.BorderLayout.CENTER);

        jPanel3.add(commentLabel, java.awt.BorderLayout.NORTH);

        jSplitPane1.setRightComponent(jPanel3);

        jPanel1.add(jSplitPane1, java.awt.BorderLayout.CENTER);

        add(jPanel1, java.awt.BorderLayout.CENTER);

        topPanel.setLayout(new java.awt.BorderLayout());

        topPanel.add(theCombo, java.awt.BorderLayout.NORTH);

        add(topPanel, java.awt.BorderLayout.NORTH);

    }//GEN-END:initComponents
	
	public void valueChanged(ListSelectionEvent e)
	{
		//System.out.println("list VALUE CHANGED");
		Setting selected = theModel.getSelectedProperty(theTable.getSelectedRow());
		
		if(selected != null)
		{
			commentLabel.setText(selected.getName());
			commentText.setText(selected.getComment());
		}
		else
		{
			commentLabel.setText("");
			commentText.setText("");
		}
		
		//                for(int i = 0; i < theModel.getRowCount(); i++)
		//                {
		//                    ////System.out.println("Row "+i);
		//                    String value = theModel.getValueAt(i, 1).toString();
		//                    int lines = 1;
		//
		//                    if(theModel.getValueAt(i, 1) instanceof FontColorProperty)
		//                    {
		//                int height = ((FontColorProperty)theModel.getValueAt(i,1)).getFontColorPair().f.getSize();
		//                height = Math.max(height, (lineWidth-2));
		//                theTable.setRowHeight(i, (height*lines)+4);
		//            }
		//            else if(theModel.getValueAt(i, 1) instanceof SingleProperty)
		//            {
		//                //lines = getNumLines(value);
		//                //int heightWanted =
		//                //theTable.setRowHeight(i, (lineWidth*lines)+2);
		//            }
		//        }
	}
	
	public void tableChanged(TableModelEvent e)
	{
		Setting selected = theModel.getSelectedProperty(theTable.getSelectedRow());
		
		if(selected != null)
		{
			commentLabel.setText(selected.getName());
			commentText.setText(selected.getComment());
		}
		else
		{
			commentLabel.setText("");
			commentText.setText("");
		}
		//System.out.println("TABLE CHANGED");
		//       CellEditor ce = theTable.getCellEditor();
		//        if(ce != null) ce.cancelCellEditing();
		//        theCombo.setModel(theModel.getComboModel());
		//        for(int i = 0; i < theModel.getRowCount(); i++)
		//        {
		//            ////System.out.println("Row "+i);
		//            String value = theModel.getValueAt(i, 1).toString();
		//            int lines = 1;
		//
		//            if(theModel.getValueAt(i, 1) instanceof FontColorProperty)
		//            {
		//                int height = ((FontColorProperty)theModel.getValueAt(i,1)).getFontColorPair().f.getSize();
		//
		//                height = Math.max(height, (lineWidth-2));
		//                theTable.setRowHeight(i, (height*lines)+4);
		//            }
		//            else if(theModel.getValueAt(i, 1) instanceof SingleProperty)
		//            {
		//                //lines = getNumLines(value);
		//                //theTable.setRowHeight(i, (lineWidth*lines)+2);
		//                //int heightWanted = (int)area.getPreferredSize().getHeight();
		//                  ///  if(heightWanted != theTable.getRowHeight(row));
		//                   // theTable.setRowHeight(row, heightWanted);
		//            }
		//            else if(theModel.getValueAt(i, 1) instanceof MultipleProperty)
		//            {
		//                lines = getNumLines(value);
		//                theTable.setRowHeight(i, (lineWidth*lines)+2);
		//            }
		//
		//
		//
		//        }
		//        doChoiceBox();
		//theTable.s
		//System.out.println("getting the combomodel");
		theCombo.setModel(theModel.getComboModel());
		doChoiceBox();
	}
	
	//    public static int getNumLines(String str)
	//    {
	//        int count = 1;
	//        for(int i = 0; i < str.length(); i++)
	//        {
	//            char curr = str.charAt(i);
	//            if(curr=='\n') count++;
	//        }
	//        ////System.out.println("count = "+count);
	//        return count;
	//    }
	
	public void itemStateChanged(ItemEvent e)
	{
		theModel.setCurrentGroup(theCombo.getSelectedIndex());
	}
	
	public void redisplaySetting(Setting setting)
	{
		Setting selected = theModel.getSelectedProperty(theTable.getSelectedRow());
		
		if(selected != null)
		{
			commentLabel.setText(selected.getName());
			commentText.setText(selected.getComment());
		}
		else
		{
			commentLabel.setText("");
			commentText.setText("");
		}
		theModel.fireTableDataChanged();
		theTable.repaint();
	}
	
    // Variables declaration - do not modify//GEN-BEGIN:variables
    javax.swing.JLabel commentLabel;
    javax.swing.JTextArea commentText;
    private javax.swing.JScrollPane jScrollPane1;
    javax.swing.JComboBox theCombo;
    javax.swing.JTable theTable;
    javax.swing.JPanel topPanel;
    // End of variables declaration//GEN-END:variables
	
	
	
	class SettingCellRenderer implements TableCellRenderer
	{
		
		/**
		 *  In this case value will be instanceof Setting
		 */
		public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column)
		{
			
			if(value instanceof Setting)
			{
				Setting setting = (Setting)value;
				
				//the getSettingRenderer() method returns a static SettingRenderer object for each Setting class, this is called with the actual
				//renderable data.
				return setting.getSettingRenderer().getTableCellRendererComponent(table, setting, setting.getValue(), isSelected, hasFocus, setting.isEnabled(), row, column);
				
			}
			else if(value instanceof ArrayList)
			{
				ArrayList settings = (ArrayList)value;
				ArrayList values = new ArrayList();
				boolean enabled = true;
				Setting first = null;
				for(int i = 0; i < settings.size();i++)
				{
					if(settings.get(i) instanceof Setting)
					{
						Setting setting = (Setting)settings.get(i);
						if(i == 0)first = setting;
						if(!setting.isEnabled()) enabled = false;
						values.add(setting.getValue());
					}
				}
				if(first != null)
					return first.getSettingRenderer().getTableCellRendererComponent(table, first, values, isSelected, hasFocus, enabled, row, column);
				else
					return new JLabel("ERRORRRRR!!!!");
			}
			else return new JLabel("ERRRORRRRRR!!!!");
		}
		
	}
	
	class SettingCellEditor extends AbstractCellEditor implements TableCellEditor
	{
		private SettingEditor currentEditor;
		
		public Object getCellEditorValue()
		{
			try
			{
				return currentEditor.getEditorValue();
			}
			catch(SettingException e)
			{
				return e; //actually return the exception for display
			}
		}
		
		public boolean stopCellEditing() { 
			fireEditingStopped(); 
			return shouldRemove;
		}
		
		public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int row, int column)
		{
			
			
			if(value instanceof Setting)
			{
				Setting setting = (Setting)value;
				currentEditor = setting.getSettingEditor();
				
				//the getSettingRenderer() method returns a static SettingRenderer object for each Setting class, this is called with the actual
				//renderable data.
				return currentEditor.getTableCellEditorComponent(table, setting, setting.getValue(), isSelected, row, column);
				
			}
			else if(value instanceof ArrayList)
			{
				ArrayList settings = (ArrayList)value;
				ArrayList values = new ArrayList();
				boolean enabled = true;
				Setting first = null;
				for(int i = 0; i < settings.size();i++)
				{
					if(settings.get(i) instanceof Setting)
					{
						Setting setting = (Setting)settings.get(i);
						if(i == 0)first = setting;
						if(!setting.isEnabled()) enabled = false;
						values.add(setting.getValue());
					}
				}
				if(first != null)
				{
					currentEditor = first.getSettingEditor();
					return currentEditor.getTableCellEditorComponent(table, first, values, isSelected, row, column);
				}
				else
					return new JLabel("ERRORRRRR!!!!");
			}
			else
			{
				return new JLabel("NEVER!!!!");
			}
		}
		
	}
	
	class SettingTableModel extends AbstractTableModel
	{
		//All of the data
		private ArrayList owners;
		
		//Current sorted data
		private ArrayList groupNames;
		private ArrayList groupStarts;
		private ArrayList groupSizes;
		private int currentGroup;
		
		private DefaultComboBoxModel comboModel;
		
		private JTable theTable;
		
		public SettingTableModel()
		{
			super();
			this.theTable = null;
			groupNames = new ArrayList();
			groupStarts = new ArrayList();
			groupSizes = new ArrayList();
			owners = new ArrayList();
			comboModel = new DefaultComboBoxModel();
		}
		
		public void setJTable(JTable tab)
		{
			this.theTable = tab;
		}
		
		public void setOwners(ArrayList owners)
		{
			this.owners = owners;
			Collections.sort(owners);
			
			Iterator it = owners.iterator();
			SettingOwner last = null;
			int currGroupCount = 0;
			String tempName = "";
			groupNames = new ArrayList();
			groupStarts = new ArrayList();
			groupSizes = new ArrayList();
			int index = 0;
			
			String ownerList = "";
			while(it.hasNext())
			{
				
				
				SettingOwner po = (SettingOwner)it.next();
				//            for(int i = 0; i < po.getNumProperties(); i++)
				//            {
				//                po.getProperty(i).addObserver(this);
				//                po.getProperty(i).setOwningModel(this);
				//            }
				if(last == null)
				{
					//this is the first group
					currGroupCount++;
					if(!po.getSettingOwnerName().equals(""))ownerList += "\'"+po.getSettingOwnerName()+"\'";
					tempName = po.getSettingOwnerClassName();
					groupStarts.add(new Integer(0));
				}
				else if(po.getSettingOwnerID() == last.getSettingOwnerID())
				{
					//this is for the second or after in the sequence
					currGroupCount++;
					//tempName = ""+currGroupCount+" "+po.getClassDescriptor()+"s";
					if(!po.getSettingOwnerClassName().endsWith("s"))
						tempName = po.getSettingOwnerClassName()+"s";
					if(!po.getSettingOwnerName().equals(""))ownerList += ", \'"+po.getSettingOwnerName()+"\'";
				}
				else
				{
					//this starts a new group
					tempName+=" "+ownerList+"";
					ownerList = "";
					groupNames.add(tempName);
					groupSizes.add(new Integer(currGroupCount));
					currGroupCount = 0;
					currGroupCount++;
					ownerList += "\'"+po.getSettingOwnerName()+"\'";
					if(!po.getSettingOwnerName().equals(""))tempName = po.getSettingOwnerClassName()+" \'"+po.getSettingOwnerName()+"\'";
					groupStarts.add(new Integer(index));
				}
				last = po;
				index++;
			}
			if(owners.size() != 0)
			{
				tempName += " "+ownerList+"";
				groupNames.add(tempName);
				groupSizes.add(new Integer(currGroupCount));
			}
			if(currentGroup > owners.size()-1) currentGroup = 0;
			comboModel = new DefaultComboBoxModel(groupNames.toArray());
			fireTableDataChanged();
		}
		
		public void refreshGroupNames()
		{
			//System.out.println("refreshing group names");
			Iterator it = owners.iterator();
			SettingOwner last = null;
			int currGroupCount = 0;
			String tempName = "";
			groupNames = new ArrayList();
			int index = 0;
			
			String ownerList = "";
			while(it.hasNext())
			{
				
				
				SettingOwner po = (SettingOwner)it.next();
				if(last == null)
				{
					//this is the first group
					currGroupCount++;
					if(!po.getSettingOwnerName().equals(""))ownerList += "\'"+po.getSettingOwnerName()+"\'";
					tempName = po.getSettingOwnerClassName();
					//groupStarts.add(new Integer(0));
				}
				else if(po.getSettingOwnerID() == last.getSettingOwnerID())
				{
					//this is for the second or after in the sequence
					currGroupCount++;
					//tempName = ""+currGroupCount+" "+po.getClassDescriptor()+"s";
					if(!po.getSettingOwnerClassName().endsWith("s"))
						tempName = po.getSettingOwnerClassName()+"s";
					if(!po.getSettingOwnerName().equals(""))ownerList += ", \'"+po.getSettingOwnerName()+"\'";
				}
				else
				{
					//this starts a new group
					tempName+=" "+ownerList+"";
					ownerList = "";
					groupNames.add(tempName);
					//System.out.println("adding: "+tempName);
					//groupSizes.add(new Integer(currGroupCount));
					currGroupCount = 0;
					currGroupCount++;
					ownerList += "\'"+po.getSettingOwnerName()+"\'";
					if(!po.getSettingOwnerName().equals(""))tempName = po.getSettingOwnerClassName()+" \'"+po.getSettingOwnerName()+"\'";
					//groupStarts.add(new Integer(index));
				}
				last = po;
				index++;
			}
			if(owners.size() != 0)
			{
				tempName += " "+ownerList+"";
				groupNames.add(tempName);
				//System.out.println("adding "+tempName);
				//groupSizes.add(new Integer(currGroupCount));
			}
			//if(currentGroup > owners.size()-1) currentGroup = 0;
			comboModel = new DefaultComboBoxModel(groupNames.toArray());
			
			
			fireTableDataChanged();
		}
		
		public String getGroupName(int i)
		{
			return (String)groupNames.get(i);
		}
		
		public int getNumGroupNames()
		{
			return groupNames.size();
		}
		
		
		public int getRowCount()
		{
			if(groupNames.size() == 0) return 0;
			SettingOwner firstInGroup = (SettingOwner)owners.get(((Integer)groupStarts.get(currentGroup)).intValue());
			return firstInGroup.getNumSettings();
		}
		
		public int getColumnCount()
		{
			return 2;
		}
		
		public String getColumnName(int column)
		{
			if(column == 0) return "Property";
			else return "Value";
		}
		
		public Object getValueAt(int row, int column)
		{
			if(column == 0)
			{
				SettingOwner firstInGroup = (SettingOwner)owners.get(((Integer)groupStarts.get(currentGroup)).intValue());
				//System.out.println("firstInGroup = "+firstInGroup);
				return firstInGroup.getSetting(row).getName();
			}
			else
			{
				
				//Simple if the selected owner group has only 1 member
				if(getCurrentGroupSize() == 1)
				{
					SettingOwner firstInGroup = getOwner(getCurrentGroupStart());
					return firstInGroup.getSetting(row);
				}
				else
				{
					ArrayList currProps = new ArrayList();
					for(int i = getCurrentGroupStart(); i < getCurrentGroupStart()+getCurrentGroupSize(); i++)
					{
						SettingOwner prop = getOwner(i);
						currProps.add(prop.getSetting(row));
					}
					
					return currProps;
				}
				
			}
		}
		
		public boolean isCellEditable(int row, int column)
		{
			if(column == 0)
			{
				return false;
			}
			else
			{
				return getSelectedProperty(row).isEnabled();
			}
		}
		
		public Setting getSelectedProperty(int listIndex)
		{
			if(listIndex < 0) return null;
			SettingOwner firstInGroup = getOwner(getCurrentGroupStart());
			return firstInGroup.getSetting(listIndex);
			
		}
		
		public void setValueAt(Object obj, int row, int column)
		{
			try
			{
				//Simple if the selected owner group has only 1 member
				if(column == 1)
				{
					if(getCurrentGroupSize() == 1)
					{
						//System.out.println("setting called obj class is "+obj.getClass().toString());
						//System.out.println("value = "+obj.toString());
						SettingOwner firstInGroup = getOwner(getCurrentGroupStart());
						if(!obj.equals(SettingEditor.NOT_CHANGED_VALUE))
							firstInGroup.getSetting(row).editValue(obj);
					}
					else
					{
						
						for(int i = getCurrentGroupStart(); i < getCurrentGroupStart()+getCurrentGroupSize(); i++)
						{
							SettingOwner prop = getOwner(i);
							if(!obj.equals(SettingEditor.NOT_CHANGED_VALUE))
								prop.getSetting(row).editValue(obj);
						}
					}
				}
				shouldRemove = true;
			}
			catch(SettingException e)
			{
				String message;
				if(obj instanceof SettingException) //if the message has been passed by the editor
				{
					message = ((SettingException)obj).getMessage();
				}
				else
				{
					message = e.getMessage();
				}
				JOptionPane.showMessageDialog(parent,message, "Error", JOptionPane.ERROR_MESSAGE);
				shouldRemove = false;
			}
		}
		
		/** Getter for property currentGroup.
		 * @return Value of property currentGroup.
		 *
		 */
		public int getCurrentGroup()
		{
			return currentGroup;
		}
		
		protected int getCurrentGroupSize()
		{
			return ((Integer)groupSizes.get(currentGroup)).intValue();
		}
		
		protected int getCurrentGroupStart()
		{
			return ((Integer)groupStarts.get(currentGroup)).intValue();
		}
		
		protected String getCurrentGroupName()
		{
			return (String)groupNames.get(currentGroup);
		}
		
		protected SettingOwner getOwner(int i)
		{
			return (SettingOwner)owners.get(i);
		}
		
		public int getNumGroups()
		{
			return groupNames.size();
		}
		
		
		/** Setter for property currentGroup.
		 * @param currentGroup New value of property currentGroup.
		 *
		 */
		public void setCurrentGroup(int currentGroup)
		{
			this.currentGroup = currentGroup;
			fireTableDataChanged();
		}
		
		
		
		
		/** Getter for property comboModel.
		 * @return Value of property comboModel.
		 *
		 */
		public javax.swing.DefaultComboBoxModel getComboModel()
		{
			return comboModel;
		}
		
	}
	
	public static void printArray(ArrayList a)
	{
		System.out.print("(");
		for(int i = 0; i < a.size(); i++)
			System.out.print(a.get(i)+" ");
		//System.out.println(")");
	}
	
	
	
	
	
	
}
