package Manananggal;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;

import org.apache.commons.math3.stat.StatUtils;

import BioKit.Exon;
import BioKit.Gene;

public class PSI_score_result implements Comparable<PSI_score_result>
{
	String 					m_strID;
	TreeSet<CountElement> 	m_vcPathA;	// stores all junctions of path A
	TreeSet<CountElement> 	m_vcPathB;	// stores all junctions of path B
	double 					m_fPValue;
	double					m_fMaxEffect;
	boolean					m_bNovelJunction;
	BioKit.Exon 					m_pUniqueExonsPathA[];
	Exon 					m_pUniqueExonsPathB[];
	Exon 					m_pSharedExons[];
	TreeSet<String> 		m_vcValidIsoforms;
	Integer 				m_nRating;
	TreeMap<String, double[]> m_mapScoresToConditions;
	
	public PSI_score_result(TreeSet<CountElement> vcPathA, TreeSet<CountElement> vcPathB, double fPValue, Gene gene, Vector<String> vcConditionsInOrder, Vector<double[]> vcScoresPerCondition)
	{
		m_vcPathA = vcPathA;
		m_vcPathB = vcPathB;
		m_fPValue = fPValue;
		
		m_bNovelJunction 	= true;
		m_strID				= "";
		m_nRating = 0;
		
		m_pUniqueExonsPathA		= null;
		m_pUniqueExonsPathB		= null;
		m_pSharedExons			= null;
		
		m_mapScoresToConditions = new TreeMap<String, double[]>();

		double m_fMinScore = Double.MAX_VALUE;
		double m_fMaxScore = 0.0;
		for(int i=0; i<vcConditionsInOrder.size(); i++)
		{
			m_mapScoresToConditions.put(vcConditionsInOrder.get(i), vcScoresPerCondition.get(i));
			
			m_fMinScore = Math.min(m_fMinScore, StatUtils.mean(vcScoresPerCondition.get(i)));
			m_fMaxScore = Math.max(m_fMaxScore, StatUtils.mean(vcScoresPerCondition.get(i)));
		}
		
		m_fMaxEffect = m_fMaxScore - m_fMinScore;
		
		m_vcValidIsoforms = new TreeSet<String>();
		
		ComputeAlternativelySplicedExons(gene);
	}
	
	public PSI_score_result()
	{
		m_strID 				= "?";
		m_vcPathA 				= new TreeSet<CountElement>();
		m_vcPathB 				= new TreeSet<CountElement>();
		m_fPValue 				= -1.0;
		m_pUniqueExonsPathA 	= null;
		m_pUniqueExonsPathB 	= null;
		m_pSharedExons 			= null;
		m_vcValidIsoforms	 	= new TreeSet<String>();
		m_nRating				= 0;
		m_bNovelJunction		= true;
		
		m_mapScoresToConditions = new TreeMap<String, double[]>();
	}

	public void Print()
	{
		System.out.println(m_strID);
		
		System.out.println("path A: "  + m_vcPathA);
		System.out.println("path B: "  + m_vcPathB);
		System.out.println("p-value: " + m_fPValue);
		
		System.out.println("valid isoforms: " + m_vcValidIsoforms);
		
		System.out.println("unique exons path A: ");
		System.out.println(Arrays.toString(m_pUniqueExonsPathA));
		
		System.out.println("unique exons path B: ");
		System.out.println(Arrays.toString(m_pUniqueExonsPathB));
		
		System.out.println("shared exons: ");
		System.out.println(Arrays.toString(m_pSharedExons));
	}
	
	public void ComputeAlternativelySplicedExons(Gene gene)
	{		
		int nPathStart	= Integer.MAX_VALUE;
		int nPathEnd	= Integer.MIN_VALUE;
		for(CountElement e : m_vcPathA)
		{
			nPathStart = Math.min(nPathStart, e.m_nStart);
			nPathEnd   = Math.max(nPathEnd,   e.m_nEnd);
		}
		for(CountElement e : m_vcPathB)
		{
			nPathStart = Math.min(nPathStart, e.m_nStart);
			nPathEnd   = Math.max(nPathEnd,   e.m_nEnd);
		}

		//############################################################################
		//     generate a unique list of exons that lie within the junction paths
		//############################################################################
		TreeSet<Exon> vcUniqueExonsPathA = new TreeSet<Exon>();
		TreeSet<Exon> vcUniqueExonsPathB = new TreeSet<Exon>();
		
		// get a list of all exons for pathA
		for(CountElement exon : m_vcPathA)
		{
			// only check exons
			if(!exon.m_bIsExon) continue;
			
			for(Exon ex : gene.getSortedExons())
			{
				if(ex.getCodingStart() == exon.m_nStart && ex.getCodingStop() == exon.m_nEnd)
					vcUniqueExonsPathA.add(ex);
			}
		}
		
		// now the same for the second path
		for(CountElement exon : m_vcPathB)
		{
			// only check exons
			if(!exon.m_bIsExon) continue;
			
			for(Exon ex : gene.getSortedExons())
			{
				if(ex.getCodingStart() == exon.m_nStart && ex.getCodingStop() == exon.m_nEnd)
					vcUniqueExonsPathB.add(ex);
			}
		}
		
		// get list of shared exons
		TreeSet<Exon> vcSharedExons = new TreeSet<Exon>();
		for(Exon ex : vcUniqueExonsPathA)
		{
			if(vcUniqueExonsPathB.contains(ex))
				vcSharedExons.add(ex);
		}
		
		// remove shared exons from exon vectors
		for(Exon ex : vcSharedExons)
		{
			vcUniqueExonsPathA.remove(ex);
			vcUniqueExonsPathB.remove(ex);
		}
		
		m_pSharedExons = new Exon[vcSharedExons.size()];
		int nIdx = 0;
		for(Exon ex : vcSharedExons)
		{
			m_pSharedExons[nIdx] = ex;
			nIdx++;
		}
		
		m_pUniqueExonsPathA = new Exon[vcUniqueExonsPathA.size()];
		nIdx = 0;
		for(Exon ex : vcUniqueExonsPathA)
		{
			m_pUniqueExonsPathA[nIdx] = ex;
			nIdx++;
		}
		
		m_pUniqueExonsPathB = new Exon[vcUniqueExonsPathB.size()];
		nIdx = 0;
		for(Exon ex : vcUniqueExonsPathB)
		{
			m_pUniqueExonsPathB[nIdx] = ex;
			nIdx++;
		}

		//###########################################
		//      get a list of isoforms that fit
		//###########################################
		Iterator<String> it = gene.getGeneProductNames();
		boolean bIsoformFitA = false;
		boolean bIsoformFitB = false;
		while(it.hasNext())
		{
			String strIsoform = it.next();
			
			// make sure the current isoform fits either path
			TreeSet<String> vcJunctions = gene.getSpliceJunctionInformationForGeneProduct(strIsoform);
			
			boolean bValidPathA = true;
			boolean bValidPathB = true;

			// isoform must contain all junctions from path A or B
			for(CountElement jun : m_vcPathA)
			{
				// skip exons
				if(jun.m_bIsExon) continue;
				
				String strJun = jun.m_nStart + "-" + jun.m_nEnd;
				if(!vcJunctions.contains(strJun))
				{
					bValidPathA = false;
					break;
				}
			}
			
			for(CountElement jun : m_vcPathB)
			{
				// skip exons
				if(jun.m_bIsExon) continue;
				
				String strJun = jun.m_nStart + "-" + jun.m_nEnd;
				if(!vcJunctions.contains(strJun))
				{
					bValidPathB = false;
					break;
				}
			}
			
			// skip invalid isoforms
			if(!(bValidPathA || bValidPathB))
				continue;
			
			if(bValidPathA)
				bIsoformFitA = true;
			
			if(bValidPathB)
				bIsoformFitB = true;

			m_vcValidIsoforms.add(strIsoform.split("\\.")[0]);
		}
		
		if(bIsoformFitA && bIsoformFitB)
			m_bNovelJunction = false;

		m_strID = "ID";
		for(CountElement e : m_vcPathA)
		{
			// skip exons
			if(e.m_bIsExon) continue;
			
			m_strID += "_" + e.m_nStart + "_" + e.m_nEnd;
		}
		
		for(CountElement e : m_vcPathB)
		{
			// skip exons
			if(e.m_bIsExon) continue;
						
			m_strID += "_" + e.m_nStart + "_" + e.m_nEnd;
		}
	}

	@Override
	public int compareTo(PSI_score_result other)
	{
		return m_strID.compareTo(other.m_strID);
	}

	@Override
	public boolean equals(Object o)
	{
		
		if(o instanceof PSI_score_result)
		{
			PSI_score_result other = (PSI_score_result) o;
			return m_strID.equals(other.m_strID);
		}
		else
			return super.equals(o);
	}

	public void WriteToFile(RandomAccessFile pOut) throws IOException
	{
		pOut.writeUTF(m_strID);
		
		pOut.writeBoolean(m_bNovelJunction);
		pOut.writeDouble(m_fMaxEffect);
		
		int nJunctionsA = 0;
		for(CountElement e : m_vcPathA)
		{
			// skip exons
			if(e.m_bIsExon) continue;
			
			nJunctionsA++;
		}
		
		int nJunctionsB = 0;
		for(CountElement e : m_vcPathB)
		{
			// skip exons
			if(e.m_bIsExon) continue;
			
			nJunctionsB++;
		}
		
		pOut.writeInt(nJunctionsA);
		for(CountElement e : m_vcPathA)
		{
			// skip exons
			if(e.m_bIsExon) continue;
			
			pOut.writeInt(e.m_nStart);
			pOut.writeInt(e.m_nEnd);
		}
	
		pOut.writeInt(nJunctionsB);
		for(CountElement e : m_vcPathB)
		{
			// skip exons
			if(e.m_bIsExon) continue;
						
			pOut.writeInt(e.m_nStart);
			pOut.writeInt(e.m_nEnd);
		}
		
		pOut.writeDouble(m_fPValue);
		
		if(m_pUniqueExonsPathA == null)
		{
			pOut.writeInt(0);
		}
		else
		{
			pOut.writeInt(m_pUniqueExonsPathA.length);
			for(Exon ex : m_pUniqueExonsPathA)
			{
				pOut.writeInt(ex.getExonID());
				pOut.writeInt(ex.getCodingStart());
				pOut.writeInt(ex.getCodingStop());
			}
		}
		
		if(m_pUniqueExonsPathB == null)
		{
			pOut.writeInt(0);
		}
		else
		{
			pOut.writeInt(m_pUniqueExonsPathB.length);
			for(Exon ex : m_pUniqueExonsPathB)
			{
				pOut.writeInt(ex.getExonID());
				pOut.writeInt(ex.getCodingStart());
				pOut.writeInt(ex.getCodingStop());
			}
		}
		
		if(m_pSharedExons == null)
		{
			pOut.writeInt(0);
		}
		else
		{
			pOut.writeInt(m_pSharedExons.length);
			for(Exon ex : m_pSharedExons)
			{
				pOut.writeInt(ex.getExonID());
				pOut.writeInt(ex.getCodingStart());
				pOut.writeInt(ex.getCodingStop());
			}
		}
		
		pOut.writeInt(m_vcValidIsoforms.size());
		for(String strIsoform : m_vcValidIsoforms)
			pOut.writeUTF(strIsoform);
		
		pOut.writeInt(m_nRating);
		
		pOut.writeInt(m_mapScoresToConditions.size());
		for(Map.Entry<String, double[]> e : m_mapScoresToConditions.entrySet())
		{
			pOut.writeUTF(e.getKey());
			double[] pValues = e.getValue();
			pOut.writeInt(pValues.length);
			
			for(double fVal : pValues)
				pOut.writeDouble(fVal);
		}
	}
	
	public void ReadFromFile(RandomAccessFile pIn) throws IOException
	{
		m_strID = pIn.readUTF();
		
		m_bNovelJunction = pIn.readBoolean();
		m_fMaxEffect	 = pIn.readDouble();
		
		int nJunctionsPathA = pIn.readInt();
		for(int i=0; i<nJunctionsPathA; i++)
		{
			CountElement e	= new CountElement();
			e.m_nStart		= pIn.readInt();
			e.m_nEnd		= pIn.readInt();
			m_vcPathA.add(e);
		}
		
		int nJunctionsPathB = pIn.readInt();
		for(int i=0; i<nJunctionsPathB; i++)
		{
			CountElement e	= new CountElement();
			e.m_nStart		= pIn.readInt();
			e.m_nEnd		= pIn.readInt();
			m_vcPathB.add(e);
		}

		m_fPValue = pIn.readDouble();
		
		int nExons = pIn.readInt();
		if(nExons == 0)
		{
			m_pUniqueExonsPathA = null;
		}
		else
		{
			m_pUniqueExonsPathA = new Exon[nExons];
			for(int i=0; i<nExons; i++)
			{
				int nID		= pIn.readInt();
				int nStart	= pIn.readInt();
				int nEnd	= pIn.readInt();
				
				Exon ex = new Exon(nStart, nEnd);
				ex.setID(nID);
				m_pUniqueExonsPathA[i] = ex;
			}
		}
		
		nExons = pIn.readInt();
		if(nExons == 0)
		{
			m_pUniqueExonsPathB = null;
		}
		else
		{
			m_pUniqueExonsPathB = new Exon[nExons];
			for(int i=0; i<nExons; i++)
			{
				int nID		= pIn.readInt();
				int nStart	= pIn.readInt();
				int nEnd	= pIn.readInt();
				
				Exon ex = new Exon(nStart, nEnd);
				ex.setID(nID);
				m_pUniqueExonsPathB[i] = ex;
			}
		}
		
		nExons = pIn.readInt();
		if(nExons == 0)
		{
			m_pSharedExons = null;
		}
		else
		{
			m_pSharedExons = new Exon[nExons];
			for(int i=0; i<nExons; i++)
			{
				int nID		= pIn.readInt();
				int nStart	= pIn.readInt();
				int nEnd	= pIn.readInt();
				
				Exon ex = new Exon(nStart, nEnd);
				ex.setID(nID);
				m_pSharedExons[i] = ex;
			}
		}
		
		int nIsoforms = pIn.readInt();
		for(int i=0; i<nIsoforms; i++)
		{
			m_vcValidIsoforms.add(pIn.readUTF());
		}
		
		m_nRating = pIn.readInt();
		
		int nConditions = pIn.readInt();
		for(int i=0; i<nConditions; i++)
		{
			String strCondition = pIn.readUTF();
			int nValues = pIn.readInt();
			
			double pValues[] = null;
			
			if(nValues > 0)
			{
				pValues = new double[nValues];
				for(int j=0; j<nValues; j++)
					pValues[j] = pIn.readDouble();
			}
			
			m_mapScoresToConditions.put(strCondition, pValues);
		}
	}

}