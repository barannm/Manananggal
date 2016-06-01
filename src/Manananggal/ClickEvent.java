package Manananggal;

import java.util.TreeSet;

import BioKit.Exon;

public class ClickEvent
{
	String m_strClickedA;	// name of first clicked element
	String m_strClickedB;	// name of second clicked element
	
	TreeSet<Exon> m_vcSelectedExons;
	TreeSet<String> m_vcTargetIsoforms;
	
	int m_nExonGroupStartA;
	int m_nExonGroupStartB;
	int m_nExonGroupEndA;
	int m_nExonGroupEndB;
	
	ClickEvent()
	{
		m_strClickedA = null;
		m_strClickedB = null;

		m_nExonGroupStartA	= -1;
		m_nExonGroupStartB	= -1;
		m_nExonGroupEndA	= -1;
		m_nExonGroupEndB	= -1;
		
		m_vcSelectedExons = new TreeSet<Exon>();
		m_vcTargetIsoforms = new TreeSet<String>();
	}
	
	void clear()
	{
		m_strClickedA = null;
		m_strClickedB = null;

		m_nExonGroupStartA	= -1;
		m_nExonGroupStartB	= -1;
		m_nExonGroupEndA	= -1;
		m_nExonGroupEndB	= -1;
		
		m_vcSelectedExons.clear();
		m_vcTargetIsoforms.clear();
	}
}
