/*
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * @author Matthias Barann
 */

package Manananggal;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;

public class AlternativeSplicingExon implements Comparable<AlternativeSplicingExon>
{
	protected int	   	m_nType;		// 1 = exon skipping event, 2 = alternative transcript start, 3 = alternative transcript end
	protected String 	m_strID;
	protected int	   	m_nStart;
	protected int	   	m_nEnd;
	protected String 	m_strExonGroup;
	protected String	m_strConditionA;
	protected String 	m_strConditionB;
	
	protected double 	m_fAltExonCovPerBaseA;	// expression value of the alternatively spliced exon
	protected double	 m_fAltExonCovPerBaseB;	// expression value of the alternatively spliced exon
	
	protected double 	m_fFractionChangeAbsolute;
	protected double 	m_fFractionChangeRelative;
	protected double	m_fPValue;
	
	protected double[] 	m_pFractionTestedExon;
	protected double[] 	m_pFractionOtherExons;
	
	AlternativeSplicingExon()
	{
		m_nType			= -1;
		m_nStart 		= -1;
		m_nEnd 			= -1;
		m_strExonGroup	= "?";
		m_strConditionA	= "?";
		m_strConditionB = "?";
		
		m_fFractionChangeAbsolute = -1.0;
		m_fFractionChangeRelative = -1.0;
		m_fPValue				  = -1.0;
		
		m_pFractionTestedExon = null;
		m_pFractionOtherExons = null;
		
		m_strID = "?";
	}
	
	AlternativeSplicingExon(int nType, String strExonGroup, int nStart, int nEnd, String strConditionA, String strConditionB, double fAltExonCovPerBaseA, double fAltExonCovPerBaseB, double fFractionChangeAbsolute, double fFractionChangeRelative, double fPValue, double[] pFractionTestedExon, double[] pFractionOtherExons)
	{
		m_nType			= nType;
		m_nStart		= nStart;
		m_nEnd			= nEnd;
		m_strExonGroup	= strExonGroup;
		m_strConditionA	= strConditionA;
		m_strConditionB = strConditionB;
		
		m_fFractionChangeAbsolute = fFractionChangeAbsolute;
		m_fFractionChangeRelative = fFractionChangeRelative;
		m_fPValue				  = fPValue;
		
		m_pFractionTestedExon = pFractionTestedExon;
		m_pFractionOtherExons = pFractionOtherExons;
		
		m_fAltExonCovPerBaseA = fAltExonCovPerBaseA;
		m_fAltExonCovPerBaseB = fAltExonCovPerBaseB;
		
		m_strID = m_nType + "_" + m_strExonGroup + "_" + m_strConditionA + "_" + m_strConditionB;
	}
	
	public int GetStart()
	{
		return m_nStart;
	}
	
	public int GetEnd()
	{
		return m_nEnd;
	}
	
	public String GetID()
	{
		return m_strID;
	}
	
	public String GetCondition(boolean bFirstCondition)
	{
		if(bFirstCondition)
			return m_strConditionA;
		
		return m_strConditionB;
	}

	@Override
	public int compareTo(AlternativeSplicingExon other)
	{
		return(m_strID.compareTo(other.m_strID));
	}
	
	public boolean equals(AlternativeSplicingExon other)
	{
		return(m_strID.equals(other.m_strID));
	}
	
	public double GetAbsoluteChange()
	{
		return m_fFractionChangeAbsolute;
	}
	
	public double GetRelativeChange()
	{
		return m_fFractionChangeRelative;
	}
	
	public double GetPValue()
	{
		return m_fPValue;
	}
	
	public int GetType()
	{
		return m_nType;
	}
	
	public String GetTypeAsString()
	{
		String strType = "?";
		switch(m_nType)
		{
			case SplicingWebApp.AS_TYPE_EXON_SKIPPING: 					strType = "exn_skipping"; 					break;
			case SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN: 			strType = "alt_start_unique_jun"; 			break;
			case SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN: 			strType = "alt_end_unique_jun"; 			break;
			case SplicingWebApp.AS_TYPE_ALT_START_SHARED_JUN: 			strType = "alt_start_shared_jun"; 			break;
			case SplicingWebApp.AS_TYPE_ALT_END_SHARED_JUN: 			strType = "alt_end_shared_jun"; 			break;
			case SplicingWebApp.AS_TYPE_RETAINED_INTRON:				strType = "retained_intron";				break;
			case SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN_DOUBLE: 	strType = "alt_start_unique_jun_double"; 	break;
			case SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN_DOUBLE: 		strType = "alt_end_unique_jun_double"; 		break;
			case SplicingWebApp.AS_TYPE_ALT_START_SHARED_JUN_DOUBLE: 	strType = "alt_start_shared_jun_double"; 	break;
			case SplicingWebApp.AS_TYPE_ALT_END_SHARED_JUN_DOUBLE: 		strType = "alt_end_shared_jun_double"; 		break;
			default:													strType = "" + m_nType; 					break;
		}
		
		return strType;
	}
	
	public void SetType(int nType)
	{
		m_nType = nType;
	}

	public void WriteToFile(FileOutputStream pOut) throws IOException
	{
		SplicingWebApp.WriteStringToFileOutputStream(m_strID, pOut);
		SplicingWebApp.WriteStringToFileOutputStream(m_strExonGroup, pOut);
		SplicingWebApp.WriteStringToFileOutputStream(m_strConditionA, pOut);
		SplicingWebApp.WriteStringToFileOutputStream(m_strConditionB, pOut);

		ByteBuffer bb = ByteBuffer.allocate(Integer.BYTES * 5 + Double.BYTES * 5);
		bb.putInt(m_nType);
		bb.putInt(m_nStart);
		bb.putInt(m_nEnd);
		bb.putDouble(m_fAltExonCovPerBaseA);
		bb.putDouble(m_fAltExonCovPerBaseB);
		bb.putDouble(m_fFractionChangeAbsolute);
		bb.putDouble(m_fFractionChangeRelative);
		bb.putDouble(m_fPValue);
		bb.putInt(m_pFractionTestedExon.length);
		bb.putInt(m_pFractionOtherExons.length);
		pOut.write(bb.array());
		
		bb = ByteBuffer.allocate(Double.BYTES * m_pFractionTestedExon.length);
		for(double fVal : m_pFractionTestedExon)
			bb.putDouble(fVal);
		pOut.write(bb.array());
		
		bb = ByteBuffer.allocate(Double.BYTES * m_pFractionOtherExons.length);
		for(double fVal : m_pFractionOtherExons)
			bb.putDouble(fVal);
		pOut.write(bb.array());
	}
	
	public void ReadFromFile(FileInputStream pIn) throws IOException
	{
		m_strID 		= SplicingWebApp.ReadStringFromFileInputStream(pIn);
		m_strExonGroup 	= SplicingWebApp.ReadStringFromFileInputStream(pIn);
		m_strConditionA	= SplicingWebApp.ReadStringFromFileInputStream(pIn);
		m_strConditionB	= SplicingWebApp.ReadStringFromFileInputStream(pIn);
		
		byte pBytes[] = new byte[Integer.BYTES * 5 + Double.BYTES * 5];
		if(pIn.read(pBytes) == -1) return;
		ByteBuffer bb = ByteBuffer.wrap(pBytes);
		m_nType 					= bb.getInt();
		m_nStart 					= bb.getInt();
		m_nEnd 						= bb.getInt();
		m_fAltExonCovPerBaseA 		= bb.getDouble();
		m_fAltExonCovPerBaseB 		= bb.getDouble();
		m_fFractionChangeAbsolute 	= bb.getDouble();
		m_fFractionChangeRelative 	= bb.getDouble();
		m_fPValue 					= bb.getDouble();
		int nTestedExonLength		= bb.getInt();
		int nOtherExonCount			= bb.getInt();
		
		m_pFractionTestedExon = new double[nTestedExonLength];
		pBytes = new byte[Double.BYTES * nTestedExonLength];
		if(pIn.read(pBytes) == -1) return;
		bb = ByteBuffer.wrap(pBytes);
		
		for(int i=0; i<nTestedExonLength; i++)
			m_pFractionTestedExon[i] = bb.getDouble();
		
		m_pFractionOtherExons = new double[nOtherExonCount];
		pBytes = new byte[Double.BYTES * nOtherExonCount];
		if(pIn.read(pBytes) == -1) return;
		bb = ByteBuffer.wrap(pBytes);
		
		for(int i=0; i<nOtherExonCount; i++)
			m_pFractionOtherExons[i] = bb.getDouble();
	}
};
