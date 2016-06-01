package Manananggal;

import java.util.TreeSet;

public class TextFilter
{
	private TreeSet<String> m_vcStrings;
	private TreeSet<String> m_vcSelectedStrings;
	
	TextFilter(TreeSet<String> vcStrings)
	{			
		m_vcStrings			= vcStrings;
		m_vcSelectedStrings = new TreeSet<String>();

		for(String strString : vcStrings)
		{
			m_vcStrings.add(strString);
			m_vcSelectedStrings.add(strString);
		}
	}
	
	public void SelectString(String strString)
	{
		m_vcSelectedStrings.add(strString);
	}
	
	public void UnselectString(String strString)
	{
		m_vcSelectedStrings.remove(strString);
	}
	
	public void SelectAll()
	{
		for(String strString : m_vcStrings)
			m_vcSelectedStrings.add(strString);
	}
	
	public void UnselectAll()
	{
		m_vcSelectedStrings.clear();
	}
	
	public boolean IsStringSelected(String strString)
	{
		if(m_vcSelectedStrings.contains(strString))
			return true;
		
		return false;
	}

	public TreeSet<String> GetStrings()
	{
		return m_vcStrings;
	}
}
