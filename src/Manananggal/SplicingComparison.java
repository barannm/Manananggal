package Manananggal;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.*;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Scanner;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;

import org.apache.commons.math3.stat.descriptive.moment.GeometricMean;
import org.apache.commons.math3.stat.descriptive.rank.Median;

import BioKit.Gene;

public class SplicingComparison {

	Connection m_SQLConnection;
	
	SplicingComparison()
	{
		m_SQLConnection = null;
	}
	
	// open database or create new database if necessary
	boolean InitSQLConnection(String strDatabase)
	{
		// create SQL data base		
		String currentDir = System.getProperty("user.dir");
        System.out.println("Current dir using System:" +currentDir);
		
		try
		{
			Class.forName("org.sqlite.JDBC");
		}
		catch(ClassNotFoundException e)
		{
			System.out.println(e.getMessage());
			System.out.println("ClassNotFoundException");
			return false;
		}
		
		try
		{
			m_SQLConnection = DriverManager.getConnection("jdbc:sqlite:" + strDatabase);
		}
		catch(SQLException e)
		{
			System.out.println(e.getMessage());
			System.out.println("SQLException");
			return false;
		}
		
		System.out.println("Opened database successfully");	
		return true;
	}

	boolean IsReady()
	{
		if(m_SQLConnection == null)
			return false;
		
		return true;
	}
	
	void CloseSQLConnection()
	{
		try
		{
			m_SQLConnection.close();
		}
		catch (SQLException e)
		{
			e.printStackTrace();
		}
	}

	@SuppressWarnings("unused")
	private class MyASEResult
	{
		String m_strGeneID;
		String m_strRef;
		char m_chStrand;
		int m_nFixedPos;
		int m_nMinorStart;
		int m_nMinorEnd;
		int m_nMajorStart;
		int m_nMajorEnd;
		TreeMap<String, int[]> m_vcMajorCountsPerCondition;	// key = condition, values = counts per sample
		TreeMap<String, int[]> m_vcMinorCountsPerCondition;	// key = condition, values = counts per sample
		TreeMap<String, double[]>  m_vcQuotients;			// key = condition, values = quotient (major/minor) per sample
		double m_fPValue;
		
		MyASEResult()
		{
			m_chStrand		= '?';
			m_strGeneID		= "?";
			m_strRef 		= "?";
			m_nFixedPos 	= -1;
			m_nMinorStart 	= -1;
			m_nMinorEnd 	= -1;
			m_nMajorStart	= -1;
			m_nMajorEnd 	= -1;
			m_fPValue 		= 0.0f;
			
			m_vcMajorCountsPerCondition = new TreeMap<String, int[]>();
			m_vcMinorCountsPerCondition = new TreeMap<String, int[]>();
			m_vcQuotients               = new TreeMap<String, double[]>();
		}
		
		public String toString(String strConditionA, String strConditionB)
		{		
			return(m_strGeneID + "\t" + m_strRef + "\t" + m_nFixedPos + "\t" + m_chStrand + "\t" +
					m_nMajorStart + "\t" + m_nMajorEnd  + "\t" + m_nMinorStart + "\t" + m_nMinorEnd + "\t" + m_fPValue + "\t" +
					Arrays.toString(m_vcMajorCountsPerCondition.get(strConditionA)) + "\t" + Arrays.toString(m_vcMinorCountsPerCondition.get(strConditionA)) + "\t" +
					Arrays.toString(m_vcMajorCountsPerCondition.get(strConditionB)) + "\t" + Arrays.toString(m_vcMinorCountsPerCondition.get(strConditionB)) + "\t" +
					Arrays.toString(m_vcQuotients.get(strConditionA)) + "\t" + Arrays.toString(m_vcQuotients.get(strConditionB)) + "\n");
		}
		
		public void ParseLine(String strLine, String strConditionA, String strConditionB)
		{
			String split[] = strLine.split("\t");
			if(split.length != 15)
			{
				System.out.println("invalid string for parsing");
				return;
			}
						
			m_strGeneID		= split[0];
			m_strRef 		= split[1];
			m_nFixedPos 	= Integer.parseInt(split[2]);
			m_chStrand		= split[3].charAt(0);			
			m_nMinorStart 	= Integer.parseInt(split[4]);
			m_nMinorEnd 	= Integer.parseInt(split[5]);
			m_nMajorStart	= Integer.parseInt(split[6]);
			m_nMajorEnd 	= Integer.parseInt(split[7]);
			m_fPValue 		= Double.parseDouble(split[8]);
			
			split[9].replace('[', ' ');
			split[9].replace(']', ' ');
			split[9].replace(',', ' ');
			
			split[10].replace('[', ' ');
			split[10].replace(']', ' ');
			split[10].replace(',', ' ');
			
			split[11].replace('[', ' ');
			split[11].replace(']', ' ');
			split[11].replace(',', ' ');
			
			split[12].replace('[', ' ');
			split[12].replace(']', ' ');
			split[12].replace(',', ' ');
			
			String splitMajorCountsA[] = split[ 9].split("\\s+");
			String splitMajorCountsB[] = split[10].split("\\s+");
			String splitMinorCountsA[] = split[11].split("\\s+");
			String splitMinorCountsB[] = split[12].split("\\s+");
			
			int nSamplesA = splitMajorCountsA.length;
			int nSamplesB = splitMajorCountsB.length;
			
			m_vcMajorCountsPerCondition = new TreeMap<String, int[]>();
			m_vcMinorCountsPerCondition = new TreeMap<String, int[]>();
			
			int pValuesMajor[] = new int[nSamplesA];
			int pValuesMinor[] = new int[nSamplesA];
			for(int i=0; i<nSamplesA; i++)
			{
				pValuesMajor[i] = Integer.parseInt(splitMajorCountsA[i]);
				pValuesMinor[i] = Integer.parseInt(splitMinorCountsA[i]);
			}
			m_vcMajorCountsPerCondition.put(strConditionA, pValuesMajor);
			m_vcMinorCountsPerCondition.put(strConditionA, pValuesMinor);
			
			pValuesMajor = new int[nSamplesB];
			pValuesMinor = new int[nSamplesB];
			for(int i=0; i<nSamplesB; i++)
			{
				pValuesMajor[i] = Integer.parseInt(splitMajorCountsB[i]);
				pValuesMinor[i] = Integer.parseInt(splitMinorCountsB[i]);
			}
			m_vcMajorCountsPerCondition.put(strConditionB, pValuesMajor);
			m_vcMinorCountsPerCondition.put(strConditionB, pValuesMinor);
			
			m_vcQuotients               = new TreeMap<String, double[]>();
			
			for(String strCondition : m_vcMajorCountsPerCondition.keySet())
			{
				int[] pMajorCounts = m_vcMajorCountsPerCondition.get(strCondition);
				int[] pMinorCounts = m_vcMajorCountsPerCondition.get(strCondition);
				
				double[] pQuotients = new double[pMajorCounts.length];
				
				for(int i=0; i<pMajorCounts.length; i++)
				{
					pQuotients[i] = (double)pMajorCounts[i] / (double)pMinorCounts[i];
				}
			}
		}
	};
	
	private class CountElement implements Comparable<CountElement>
	{
		String 	m_strRef;
		int 	m_nStart;
		int 	m_nEnd;
		String 	m_strGeneID;
		boolean m_bKnown;
		char 	m_chStrand;
		
		Vector<Integer> m_vcCounts;
		
		@SuppressWarnings("unused")
		CountElement()
		{
			m_strRef 	= "?";
			m_nStart 	= -1;
			m_nEnd		= -1;
			m_strGeneID	= "?";
			m_bKnown	= false;
			m_chStrand	= '?';
			
			m_vcCounts = new Vector<Integer>();
		}
		
		@Override
		public boolean equals(Object o)
		{
			if(o instanceof CountElement)
			{
				CountElement other = (CountElement) o;
				if(!m_strRef.equals(other.m_strRef)) return false;
				if(m_nStart != other.m_nStart || m_nEnd != other.m_nEnd) return false;
				if(m_chStrand != other.m_chStrand) return false;
				
				return true;
			}
			else
				return super.equals(o);
		}
		
		// use this function when reading from merged tables
		@SuppressWarnings("unused")
		void ParseLine(String strLine, Vector<Integer> vcSamplesToKeep)
		{
			String vcSplit[] = strLine.split("\t");
			if(vcSplit.length < 7)
			{
				System.out.println("Could not parse splice junction: invalid formated line: " + strLine);
				return;
			}
			
			m_strRef	= vcSplit[0];
			m_nStart	= Integer.parseInt(vcSplit[1]);
			m_nEnd		= Integer.parseInt(vcSplit[2]);
			m_chStrand	= vcSplit[3].charAt(0);
			if(vcSplit[4].equals("known"))
				m_bKnown = true;
			else
				m_bKnown = false;
			m_strGeneID	= vcSplit[5];
			
			for(int i=6; i<vcSplit.length; i++)
			{
				if(vcSamplesToKeep.contains(i))
					m_vcCounts.addElement(Integer.parseInt(vcSplit[i]));
			}
		}
		
		// use this function when reading from single counts files
		@SuppressWarnings("unused")
		void ParseLine(String strLine)
		{
			String vcSplit[] = strLine.split("\t");
			if(vcSplit.length < 10)
			{
				System.out.println("Could not parse splice junction: invalid formated line: " + strLine);
				return;
			}
			
			m_strGeneID	= vcSplit[1];
			
			if(vcSplit[2].equals("true"))
				m_chStrand = '+';
			else
				m_chStrand = '-';
			
			m_strRef	= vcSplit[4];
			m_nStart	= Integer.parseInt(vcSplit[5]);
			m_nEnd		= Integer.parseInt(vcSplit[6]);			

			m_chStrand	= vcSplit[3].charAt(0);
			if(vcSplit[0].equals("Novel_Junction"))
				m_bKnown = false;
			else
				m_bKnown = true;
			
			m_vcCounts.addElement(Integer.parseInt(vcSplit[8]));
		}
		
		@SuppressWarnings("unused")
		void ParseMergedData(String strInfo, Vector<Integer> vcValues)
		{
			String vcSplit[] = strInfo.split("\t");
			if(vcSplit.length < 6)
			{
				System.out.println("Could not parse splice junction: invalid formated line: " + strInfo);
				return;
			}
			
			m_strRef	= vcSplit[0];
			m_nStart	= Integer.parseInt(vcSplit[1]);
			m_nEnd		= Integer.parseInt(vcSplit[2]);
			m_chStrand	= vcSplit[3].charAt(0);
			if(vcSplit[4].equals("known"))
				m_bKnown = true;
			else
				m_bKnown = false;
			m_strGeneID	= vcSplit[5];
			
			m_vcCounts.addAll(vcValues);
		}
		
		int GetStart() {return m_nStart;}
		int GetEnd() {return m_nEnd;}
		@SuppressWarnings("unused")
		int GetLength() {return m_nEnd - m_nStart +1;}
		String GetRef() {return m_strRef;}
		@SuppressWarnings("unused")
		String GetGeneID() {return m_strGeneID;}

		public int compareTo(CountElement other)
		{
			int nDiff = this.m_nStart - other.GetStart();
			if(nDiff != 0)
				return nDiff;
			
			return this.m_nEnd - other.GetEnd();
		}
		
		public String toString()
		{
			return(m_strRef + "\t" + m_nStart + "\t" + m_nEnd + "\t" + m_bKnown + "\t" + m_vcCounts);
		}
	}
	
	private class NetworkNode
	{
		int m_nFirstPos;				// fixed position of the splice junction
		Vector<Integer> m_vcSecond;		// vector with indices to variable splice junction positions
		int m_nMajorIdx;				// index of major splice junction
		
		@SuppressWarnings("unused")
		NetworkNode()
		{
			m_nFirstPos	= -1; 
			m_nMajorIdx = -1;
			m_vcSecond 	= new Vector<Integer>();
		}
		
		NetworkNode(int nFirstPos, int nSecondIdx)
		{
			m_nFirstPos = nFirstPos;
			m_nMajorIdx = -1;
			m_vcSecond 	= new Vector<Integer>();
			m_vcSecond.add(nSecondIdx);
		}
		
		void AddConnection(int nIdx)
		{
			// only add new connections
			if(!m_vcSecond.contains(nIdx))
				m_vcSecond.add(nIdx);
			else
			{
				System.out.println("Warning: tried to add junction multiple times -> verify that your input matrix does not contain duplicate positions");
			}
		}
		
		@SuppressWarnings("unused")
		int GetPosition()
		{
			return m_nFirstPos;
		}
		
		Vector<Integer> GetConnections()
		{
			return m_vcSecond;
		}
		
		@SuppressWarnings("unused")
		int GetMajorIdx()
		{
			return m_nMajorIdx;
		}
	}
	
	@SuppressWarnings("unused")
	private class MyJunction implements Comparable<MyJunction>
	{
		int m_nExIDA;
		int m_nExIDB;
		
		@Override
		public boolean equals(Object o)
		{
			if(o instanceof MyJunction)
			{
				MyJunction other = (MyJunction) o;
				return(this.m_nExIDA == other.m_nExIDA && this.m_nExIDB == other.m_nExIDB);
			}
			return false;
		}

		@Override
		public int compareTo(MyJunction other)
		{
			int nDiff = this.m_nExIDA - other.m_nExIDA;
			if(nDiff != 0)
				return nDiff;
			
			return this.m_nExIDB - other.m_nExIDB;
		}
		
		@Override
		public String toString()
		{
			return("exon A:\t" + m_nExIDA + "\t" + "exon B:\t" + m_nExIDB);
		}
	}
	
	public void MergeTables(String strInputFolder, String strOutputFolder) throws IOException
	{
		if(m_SQLConnection == null)
		{
			System.out.println("ERROR: SQL connection invalid");
			return;
		}
		
		// add table for junction counts
		try
		{
			// disable auto commit for batch processing (=faster)
			m_SQLConnection.setAutoCommit(false);
			
			// add junction count table
			MergeJunctionTable(strInputFolder);
			
			// add exon count table
			MergeExonTable(strInputFolder);

			// re-enable auto commit
			m_SQLConnection.setAutoCommit(true);

			// add indices
			Statement stmt = m_SQLConnection.createStatement();
			String sql = "CREATE INDEX junction_gene_idx ON junctions (GENE);";
			stmt.executeUpdate(sql);
			
			sql = "CREATE INDEX exon_gene_idx ON exons (GENE);";
			stmt.executeUpdate(sql);
			stmt.close();
		}
		catch (SQLException e)
		{
			e.printStackTrace();
			System.out.println("SQL command failed");
		}
	}
	
	public void MergeJunctionTable(String strInputFolder) throws IOException
	{
		if(m_SQLConnection == null)
		{
			System.out.println("ERROR: SQL connection invalid");
			return;
		}
		
		System.out.println("Adding junction data to sqlite data base from: " + strInputFolder);
		File pIn = new File(strInputFolder);
		File pInputFiles[] = pIn.listFiles();
		
		TreeSet<String> vcSampleFiles = new TreeSet<String>();
		TreeSet<String> vcSampleNames = new TreeSet<String>();
		
		// get junction file names for all samples
		for(int i=0; i<pInputFiles.length; i++)
		{
			if(pInputFiles[i].getName().endsWith("junction_cnts.tsv"))
			{
				// get sample name
				String strSample = pInputFiles[i].getName().split(".spliced")[0];
				vcSampleNames.add(strSample);
				
				vcSampleFiles.add(pInputFiles[i].getAbsolutePath());
			}
		}
		
		int nCols = 7;
		String strColNames = "ID,REF,START,END,STRAND,TYPE,GENE";
		
		try
		{
			Statement stmt = m_SQLConnection.createStatement();
		    String sql = "CREATE TABLE junctions " +
		    			 "(ID		INT		PRIMARY KEY	NOT NULL," +
		    			 " REF		TEXT	NOT NULL," +
		    			 " START	INT 	NOT NULL," + 
		                 " END		INT		NOT NULL," +
		                 " STRAND	CHAR(1) NOT NULL," + 
		                 " TYPE		TEXT	NOT NULL," +
		                 " GENE		TEXT	NOT NULL";
		    
		    for(String strSample : vcSampleNames)
		    {
		    	sql += ", '" + strSample + "' INT NOT NULL";
		    	strColNames += ",'" + strSample + "'";
		    	nCols += 1;
		    }
		    sql += ")";
			stmt.executeUpdate(sql);
			stmt.close();
		}
		catch(SQLException e)
		{
			e.printStackTrace();
			System.out.println("SQLException");
		}
		
		// generate merged table
		TreeMap<String, Vector<Integer>> mapData = new TreeMap<String, Vector<Integer>>();
		
		Iterator<String> it = vcSampleFiles.iterator();
		int nSamples = vcSampleFiles.size();
		int nSampleIdx = 0;
		while(it.hasNext())
		{
			String strFile = it.next();
			File pFile = new File(strFile);
			try
			{
				Scanner pScanner = new Scanner(pFile);
				
				while(pScanner.hasNextLine())
				{
					String strLine = pScanner.nextLine();
					String vcSplit[] = strLine.split("\t");
					
					if(vcSplit.length != 10)
					{
						System.out.println("Invalid number of fields (!=10): " + strLine);
						continue;
					}
					
					// prepare junction information
					String strID = vcSplit[4] + "\t" + vcSplit[5] + "\t" + vcSplit[6];
					if(vcSplit[2].equals("true"))
						strID += "\t" + "+";
					else
						strID += "\t" + "-";
					
					if(strLine.contains("Novel_Junction"))
						strID += "\t" + "novel";
					else
						strID += "\t" + "known";
					
					strID += "\t" + vcSplit[1];
					
					// prepare junction counts
					if(mapData.containsKey(strID))
					{
						// update counts
						Vector<Integer> vcValues = mapData.get(strID);
						vcValues.set(nSampleIdx, Integer.parseInt(vcSplit[8]));
						mapData.put(strID, vcValues);
						
					}
					else
					{
						// add new junction with 0 counts
						Vector<Integer> vcValues = new Vector<Integer>();
						for(int i=0; i<nSamples; i++)
						{
							if(i == nSampleIdx)
							{
								vcValues.add(Integer.parseInt(vcSplit[8]));
							}
							else
								vcValues.add(0);
						}
						mapData.put(strID, vcValues);
					}
				}
				
				pScanner.close();
			}
			catch(FileNotFoundException e)
			{
				System.out.println(e.getMessage());
			}
			
			nSampleIdx += 1;
		}
		
		try
		{
			String sql = "INSERT INTO junctions (" +  strColNames + ") VALUES (";

			for(int i=0; i<nCols; i++)
			{
				if(i == 0)
					sql += "?";
				else
					sql += ",?";
			}
			sql += ")";
	
			PreparedStatement stmt = m_SQLConnection.prepareStatement(sql);
			
			// add count data
			int nID = 0;
			for(Map.Entry<String, Vector<Integer>> e : mapData.entrySet())
			{
				String strID = e.getKey();
				String split[] = strID.split("\t");
								
				stmt.setInt(1, nID);
				stmt.setString(2, split[0]);
				stmt.setInt(3, Integer.parseInt(split[1]));
				stmt.setInt(4, Integer.parseInt(split[2]));
				stmt.setString(5, split[3]);
				stmt.setString(6, split[4]);
				stmt.setString(7, split[5]);

				int nMaxCount = Collections.max(e.getValue());
				if(nMaxCount > 0)
				{
					int nIdx = 8;
					for(int nValue : e.getValue())
					{
						stmt.setInt(nIdx, nValue);
						nIdx += 1;
					}
					
					stmt.addBatch();
				}
				
				if(nID % 1000 == 0)
				{
					stmt.executeBatch();
					stmt.clearParameters();
					m_SQLConnection.commit();
				}
				
				nID += 1;
			}
			
			stmt.executeBatch();
			m_SQLConnection.commit();
			stmt.close();
		}
		catch(Exception e)
		{
			e.printStackTrace();
			System.out.println("SQL command failed");
		}
	}
	
	public void MergeExonTable(String strInputFolder)
	{
		if(m_SQLConnection == null)
		{
			System.out.println("ERROR: SQL connection invalid");
			return;
		}
		
		System.out.println("Merging exon count files contained in: " + strInputFolder);
		File pIn = new File(strInputFolder);
		File pInputFiles[] = pIn.listFiles();
		
		TreeSet<String> vcSampleFiles = new TreeSet<String>();
		TreeSet<String> vcSampleNames = new TreeSet<String>();
		
		// get junction file names for all samples
		for(int i=0; i<pInputFiles.length; i++)
		{
			if(pInputFiles[i].getName().endsWith("exon_cnts.tsv"))
			{
				// get sample name
				String strSample = pInputFiles[i].getName().split(".spliced")[0];
				vcSampleNames.add(strSample);
				
				vcSampleFiles.add(pInputFiles[i].getAbsolutePath());
			}
		}
		
		// generate merged table
		TreeMap<String, Vector<Integer>> mapData = new TreeMap<String, Vector<Integer>>();
		
		Iterator<String> it = vcSampleFiles.iterator();
		int nSamples = vcSampleFiles.size();
		int nSampleIdx = 0;
		while(it.hasNext())
		{
			String strFile = it.next();
			File pFile = new File(strFile);
			try
			{
				Scanner pScanner = new Scanner(pFile);
				
				while(pScanner.hasNextLine())
				{
					String strLine = pScanner.nextLine();
					String vcSplit[] = strLine.split("\t");
					
					if(vcSplit.length != 10)
					{
						System.out.println("Invalid number of fields (!=10): " + strLine);
						continue;
					}
					
					// prepare junction information
					String strID = vcSplit[4] + "\t" + vcSplit[5] + "\t" + vcSplit[6];
					if(vcSplit[2].equals("true"))
						strID += "\t" + "+";
					else
						strID += "\t" + "-";
					
					strID += "\t" + "known";
					
					strID += "\t" + vcSplit[1];
					
					// prepare junction counts
					if(mapData.containsKey(strID))
					{
						// update counts
						Vector<Integer> vcValues = mapData.get(strID);
						vcValues.set(nSampleIdx, Integer.parseInt(vcSplit[8]));
						mapData.put(strID, vcValues);
						
					}
					else
					{
						// add new junction with 0 counts
						Vector<Integer> vcValues = new Vector<Integer>();
						for(int i=0; i<nSamples; i++)
						{
							if(i == nSampleIdx)
							{
								vcValues.add(Integer.parseInt(vcSplit[8]));
							}
							else
								vcValues.add(0);
						}
						mapData.put(strID, vcValues);
					}
				}
				
				pScanner.close();
			}
			catch(FileNotFoundException e)
			{
				System.out.println(e.getMessage());
			}
			
			nSampleIdx += 1;
		}
		
		try
		{
			int nCols = 7;
			String strColNames = "ID,REF,START,END,STRAND,TYPE,GENE";
			
			// create the table
			Statement statement = m_SQLConnection.createStatement();
		    String sql = "CREATE TABLE exons " +
		    			 "(ID		INT		PRIMARY KEY	NOT NULL," +
		    			 " REF		TEXT	NOT NULL," +
		    			 " START	INT 	NOT NULL," + 
		                 " END		INT		NOT NULL," +
		                 " STRAND	CHAR(1) NOT NULL," + 
		                 " TYPE		TEXT	NOT NULL," +
		                 " GENE		TEXT	NOT NULL";
		    
		    for(String strSample : vcSampleNames)
		    {
		    	sql += ", '" + strSample + "' INT NOT NULL";
		    	strColNames += ",'" + strSample + "'";
		    	nCols += 1;
		    }
		    sql += ")";
		    statement.executeUpdate(sql);
		    statement.close();
		    
			// prepare SQL statement
			sql = "INSERT INTO exons (" +  strColNames + ") VALUES (";

			for(int i=0; i<nCols; i++)
			{
				if(i == 0)
					sql += "?";
				else
					sql += ",?";
			}
			sql += ")";
	
			PreparedStatement stmt = m_SQLConnection.prepareStatement(sql);

			// add count data
			int nID = 0;
			for(Map.Entry<String, Vector<Integer>> e : mapData.entrySet())
			{
				String strID = e.getKey();
				String split[] = strID.split("\t");
				
				stmt.setInt(1, nID);
				stmt.setString(2, split[0]);				
				stmt.setInt(3, Integer.parseInt(split[1]));
				stmt.setInt(4, Integer.parseInt(split[2]));
				stmt.setString(5, split[3]);
				stmt.setString(6, split[4]);
				stmt.setString(7, split[5]);
				
				// only add data if at least one sample had a read
				int nMaxCount = Collections.max(e.getValue());
				if(nMaxCount > 0)
				{
					int nIdx = 8;
					for(int nValue : e.getValue())
					{
						stmt.setInt(nIdx, nValue);
						nIdx += 1;
					}

					stmt.addBatch();
				}
				
				if(nID % 1000 == 0)
				{
					stmt.executeBatch();
					stmt.clearParameters();
					m_SQLConnection.commit();
				}
				
				nID += 1;
			}
			
			stmt.executeBatch();
			m_SQLConnection.commit();
			stmt.close();
		}
		catch(Exception e)
		{
			e.printStackTrace();
			System.out.println("SQL command failed");
		}
	}

	private boolean LoadSampleDesc(String strFile, String strDatabase)
	{
		if(m_SQLConnection == null)
		{
			System.out.println("ERROR: SQL connection invalid");
			return false;
		}
		
		//############################################
		//    create table with sample information
		//############################################
		try
		{
			Statement stmt = m_SQLConnection.createStatement();
		    String sql = "CREATE TABLE sample_info " +
		    			 "(ID			INT		NOT NULL," +
		    			 " SAMPLE_ID	TEXT	PRIMARY KEY	NOT NULL," +
		    			 " TISSUE		TEXT 	NOT NULL," + 
		                 " CONDITION1	TEXT	NOT NULL," +
		                 " CONDITION2	TEXT	NOT NULL," + 
		                 " CONDITION3	TEXT	NOT NULL," +
		                 " PAIRED_END	INT		NOT NULL," +
		                 " STRANDED		INT		NOT NULL," +
		                 " SIZE_FACTOR	REAL	NOT NULL," +
		    			 " BAM_FILE		TEXT," 			   +
		    			 " BIGWIG_FILE	TEXT";
		    sql += ")";
			stmt.executeUpdate(sql);
			stmt.close();
		}
		catch(SQLException e)
		{
			e.printStackTrace();
			System.out.println("SQLException");
		}
		
		File pIn = new File(strFile);
		try
		{
			Scanner pScanner = new Scanner(pIn);
			
			TreeMap<String, Integer> mapIndexToColnames = new TreeMap<String, Integer>();	// keys = column name, values = column index
			
			boolean bProcessedHeader = false;
			int nID = 0;
			while(pScanner.hasNextLine())
			{
				String strLine = pScanner.nextLine();
				
				String split[] = strLine.split("\t");

				// skip empty lines and lines
				if(split.length == 0)
					continue;
				
				if(!bProcessedHeader)
				{
					if(strLine.startsWith("#"))
						strLine = strLine.substring(1);
					
					String pHeader[] = strLine.split("\t");
					
					for(int i=0; i<pHeader.length; i++)
					{
						pHeader[i] = pHeader[i].toLowerCase();
						mapIndexToColnames.put(pHeader[i], i);
					}
					bProcessedHeader = true;
				}
				else
				{
					String strSample	 = null;
					String strCondition1 = null;
					String strCondition2 = "???";
					String strCondition3 = "???";
					String strTissue	 = "???";
					int nIsPairedEnd	 = -1;
					int nIsStranded	 	 = -1;
					String strBamFile	 = "???";
					String strBigWigFile = "???";
					
					if((mapIndexToColnames.containsKey("sample") && mapIndexToColnames.get("sample") < split.length) || (mapIndexToColnames.containsKey("sample_id") && mapIndexToColnames.get("sample_id") < split.length))
					{
						if(mapIndexToColnames.containsKey("sample"))
							strSample = split[mapIndexToColnames.get("sample")];
						else
							strSample = split[mapIndexToColnames.get("sample_id")];
					}
					else
						strSample = split[0];
					
					if(mapIndexToColnames.containsKey("bam") && mapIndexToColnames.get("bam") < split.length)
						strBamFile = split[mapIndexToColnames.get("bam")];
					
					if(mapIndexToColnames.containsKey("bigwig") && mapIndexToColnames.get("bigwig") < split.length)
						strBigWigFile = split[mapIndexToColnames.get("bigwig")];
					
					if(mapIndexToColnames.containsKey("condition1") && mapIndexToColnames.get("condition1") < split.length)
						strCondition1 = split[mapIndexToColnames.get("condition1")];
					
					if(mapIndexToColnames.containsKey("condition2") && mapIndexToColnames.get("condition2") < split.length)
						strCondition2 = split[mapIndexToColnames.get("condition2")];
					
					if(mapIndexToColnames.containsKey("condition3") && mapIndexToColnames.get("condition3") < split.length)
						strCondition3 = split[mapIndexToColnames.get("condition3")];
					
					if((mapIndexToColnames.containsKey("paired_end") && mapIndexToColnames.get("paired_end") < split.length) || (mapIndexToColnames.containsKey("paired") && mapIndexToColnames.get("paired") < split.length))
					{
						String strVal = "";
						if(mapIndexToColnames.containsKey("paired_end"))
							strVal = split[mapIndexToColnames.get("paired_end")];
						else
							strVal = split[mapIndexToColnames.get("paired")];
						
						if(strVal.equals("true"))
							nIsPairedEnd = 1;
						else if(strVal.equals("false"))
							nIsPairedEnd = 0;
						else
						{
							int nVal = Integer.parseInt(strVal);
							if(nVal == 0 || nVal == 1)
								nIsPairedEnd = nVal;
							else
							{
								System.out.println("invalid value for \"paired_end\"");
							}
						}
					}
					
					if(mapIndexToColnames.containsKey("stranded") && mapIndexToColnames.get("stranded") < split.length)
					{
						String strVal = split[mapIndexToColnames.get("stranded")];
						if(strVal.equals("true"))
							nIsPairedEnd = 1;
						else if(strVal.equals("false"))
							nIsPairedEnd = 0;
						else
						{
							int nVal = Integer.parseInt(strVal);
							if(nVal == 0 || nVal == 1)
								nIsPairedEnd = nVal;
							else
							{
								System.out.println("invalid value for \"stranded\"");
							}
						}
					}
					
					if(mapIndexToColnames.containsKey("tissue") && mapIndexToColnames.get("tissue") < split.length)
						strTissue = split[mapIndexToColnames.get("tissue")];
					
					if(strSample == null || strCondition1 == null)
					{
						System.out.println("Invalid sample description: " + strLine);
						System.out.println("A valid sample description must contain at least one condition.");
						pScanner.close();
						return false;
					}
					else
					{
						String sql = "INSERT INTO sample_info (ID,SAMPLE_ID,TISSUE,CONDITION1,CONDITION2,CONDITION3,PAIRED_END,STRANDED,SIZE_FACTOR,BAM_FILE,BIGWIG_FILE) VALUES(";
						
						sql += "'" + nID 			+ "'";
						sql += ",'" + strSample 	+ "'";
						sql += ",'" + strTissue 	+ "'";
						sql += ",'" + strCondition1 + "'";
						sql += ",'" + strCondition2 + "'";
						sql += ",'" + strCondition3 + "'";
						sql += ",'" + nIsPairedEnd 	+ "'";
						sql += ",'" + nIsStranded 	+ "'";
						sql += ",'-1'";
						sql += ",'" + strBamFile 	+"'";
						sql += ",'" + strBigWigFile +"'";
						sql += ");";
						
						try
						{
							Statement stmt = m_SQLConnection.createStatement();
							stmt.executeUpdate(sql);
							stmt.close();
						}
						catch(SQLException e)
						{
							System.out.println(e.getMessage());
							System.out.println(e.getStackTrace());
						}
					}
					nID += 1;
				}				
			}
			pScanner.close();
		}
		catch(FileNotFoundException e)
		{
			System.out.println(e.getMessage());
		}

		return true;
	}
	
	@SuppressWarnings("unused")
	private void GenerateSplicingNetwork(Vector<CountElement> vcJunctions, int nMinRowSum, int nMinCnt, HashMap<Integer, NetworkNode> mapNetwork)
	{
		//###############################################################
		//   generate a map of junctions, one end of the junction will
		//   be used as fixed position and the indices to all possible
		//   end positions will be stored. The major index corresponds
		//   to the highest covered splice junction across all samples
		//###############################################################
		int nJunctionIdx = 0;
		for(CountElement s : vcJunctions)
		{
			// network keys are the junction positions
			int nStart 	= s.GetStart();
			int nEnd 	= s.GetEnd();
			
			// check if junction has been added to the network already
			if(mapNetwork.containsKey(nStart))
			{
				NetworkNode node = mapNetwork.get(nStart);
				node.AddConnection(nJunctionIdx);
				mapNetwork.put(nStart, node);
			}
			else if(mapNetwork.containsKey(nEnd))
			{
				NetworkNode node = mapNetwork.get(nEnd);
				node.AddConnection(nJunctionIdx);
				mapNetwork.put(nEnd, node);
			}
			else
			{
				NetworkNode node1 = new NetworkNode(nStart, nJunctionIdx);
				mapNetwork.put(nStart, node1);
				NetworkNode node2 = new NetworkNode(nEnd, nJunctionIdx);
				mapNetwork.put(nEnd, node2);
			}
			
			nJunctionIdx += 1;
		}
		
		//#################################
		//     identify major junction
		//#################################
		for(Map.Entry<Integer, NetworkNode> e : mapNetwork.entrySet())
		{
			// get row sum and maximum count value
			NetworkNode node = e.getValue();
			Vector<Integer> vcConnections = node.GetConnections();
			int nMajorIdx = -1;
			double fMajorMean = -1.0f;
			for(int idx : vcConnections)
			{
				CountElement jun = vcJunctions.get(idx);
				if(jun.m_vcCounts.isEmpty())
				{
					System.out.println("skipping junction: " + jun.GetRef() + ":" + jun.GetStart() + "-" + jun.GetEnd() + " (no count data available)");
					continue;
				}
				
				int nMax = Collections.max(jun.m_vcCounts);
				
				// skip junctions which maximum read count is below the threshold
				if(nMax < nMinCnt)
					continue;
				
				// skip junctions which read count sum is below the threshold
				int nSum = 0;
				for(int val : jun.m_vcCounts)
					nSum += val;
				
				if(nSum < nMinRowSum)
					continue;
				
				// identify major isoform
				double fMean = (double) nSum / (double) jun.m_vcCounts.size();
				if(fMean > fMajorMean)
				{
					fMajorMean = fMean;
					nMajorIdx = idx;
				}
			}
			
			if(nMajorIdx == -1)
			{
				System.out.println(vcConnections.toString());
			}
			
			node.m_nMajorIdx = nMajorIdx;
		}
	}
	
	@SuppressWarnings("unused")
	private void OutputSplicingNetwork(String strOutputFolder, String strKey, HashMap<Integer, NetworkNode> mapNetwork, Vector<CountElement> vcJunctions)
	{
		try
		{
			PrintWriter p = new PrintWriter(strOutputFolder + "/network_" + strKey + ".txt");
			
			for(Map.Entry<Integer, NetworkNode> e : mapNetwork.entrySet())
			{
				NetworkNode node = e.getValue();
				
				if(node.m_vcSecond.size() > 1)
				{
					p.write(node.m_nFirstPos + "\nmajor\t " + vcJunctions.get(node.m_nMajorIdx).toString() + "\n");
					
					for(int nIdx : node.m_vcSecond)
					{
						p.write("\t" + vcJunctions.get(nIdx).toString() + "\n");
					}
				}
			}
			
			p.close();
		}
		catch(FileNotFoundException e)
		{
			System.out.println(e.getMessage());
		}
	}
	
	//#########################################################################################
	// size factors are calculated similar to DESeq size factors:
	//
	// Source: Simon Anders (http://seqanswers.com/forums/showpost.php?p=16468&postcount=13)
	//
	// - Construct a "reference sample" by taking, for each gene, the geometric mean of the
	//   counts in all samples.
	//
	// - To get the sequencing depth of a sample relative to the reference, calculate for each
	//   gene the quotient of the counts in your sample divided by the counts of the reference
	//   sample. Now you have, for each gene, an estimate of the depth ratio.
	//
	// - Simply take the median of all the quotients to get the relative depth of the library.
	//#########################################################################################
	private boolean CalculateSizeFactors(String strDatabase)
	{	
		try
		{
			Statement stmt = m_SQLConnection.createStatement();
			
			Vector<String> vcSamples = new Vector<String>();
			int nSamples = 0;
			
			String strSQL = "SELECT SAMPLE_ID FROM sample_info";
			ResultSet pRes = stmt.executeQuery(strSQL);
			while(pRes.next())
			{
				vcSamples.add(pRes.getString("SAMPLE_ID"));
				nSamples += 1;
			}
			
			System.out.println("number of samples: " + nSamples);

			Vector<Vector<Double>> vcQuotient = new Vector<Vector<Double>>();
			for(int i=0; i<nSamples; i++)
				vcQuotient.add(new Vector<Double>());
			
			GeometricMean gm = new GeometricMean();
			
			strSQL = "SELECT * FROM exons ORDER BY GENE";			
			pRes = stmt.executeQuery(strSQL);
			
			int nGenes = 0;
			while(pRes.next())
			{		
				double vcGeneReadCounts[] = new double[nSamples];
				
				int nIdx = 0;
				for(String strSample : vcSamples)
				{
					vcGeneReadCounts[nIdx] = pRes.getInt(strSample);
					nIdx += 1;
				}
				
				double fMean = gm.evaluate(vcGeneReadCounts);
				
				Vector<Double> vcCurrentQuotients = new Vector<Double>();
				
				// calculate per sample difference to geometric mean
				for(int i=0; i<nSamples; i++)
				{
					if(fMean != 0)
						vcCurrentQuotients.add(vcGeneReadCounts[i] / fMean);
				}
				
				if(vcCurrentQuotients.size() == nSamples)
				{
					for(int i=0; i<nSamples; i++)
					vcQuotient.get(i).add(vcCurrentQuotients.get(i));
					nGenes +=1;
				}
			}
			
			for(int i=0; i<nSamples; i++)
			{
				Median med = new Median();
				double[] pVals = new double[nGenes];
				for(int j=0; j<nGenes; j++)
					pVals[j] = vcQuotient.get(i).get(j);
				
				double fSizeFactor = med.evaluate(pVals);
				
				int nIdx = 0;
				String strSampleName = vcSamples.get(i);
				
				if(strSampleName == null)
					System.out.println("ERROR: invalid sample " + i + " " +  nIdx);
	
				strSQL = "UPDATE sample_info SET SIZE_FACTOR = '" + fSizeFactor + "' WHERE SAMPLE_ID = '" + strSampleName + "'";
				stmt.executeUpdate(strSQL);
				
				System.out.println(strSampleName + ": " + fSizeFactor);
			}
		}
		catch(SQLException e)
		{
			e.printStackTrace();
		}
		
		return true;
	}

	public void InitWithData(String strInputFolder, String strSampleDescFile, String strDataBaseName, String strOutputFolder) throws IOException
	{			
		String strDatabase =  strOutputFolder + "/" + strDataBaseName + "_counts.sqlite";
		File pDatabase     = new File(strDatabase);		

		//  check if database exists, otherwise create it
		if(!pDatabase.exists())
		{
			InitSQLConnection(strDatabase);
			// add sample information to data base
			System.out.println("loading sample description");
			if(!LoadSampleDesc(strSampleDescFile, strDatabase)) return;
			
			// create junction and exon read counts tables
			System.out.println("preparing count tables");
			MergeTables(strInputFolder, strOutputFolder);
			
			// calculate size factors and save them in the database
			System.out.println("calculating size factors");
			if(!CalculateSizeFactors(strDatabase)) return;
		}
		else
		{
			InitSQLConnection(strDatabase);
		}
	}

	/*
	 * Returns a tree map in the format: <strJunction, <strSample, nCount>>
	 */
	public TreeMap<String, TreeMap<String, Integer>> GetJunctionCountsForGene(Gene gene)
	{					
		// for each condition (1st key = string) store the junction read counts (values) per junction (2nd key)
		TreeMap<String, TreeMap<String, Integer>> vcResult = new TreeMap<String, TreeMap<String, Integer>>();

		//#############################################
		//    get junction counts from the database
		//##############################################
		try
		{
			//############################ 
			//    get all sample names
			//############################
			TreeSet<String> vcSamples = GetSamples();
			
			//##########################
			//      get count data
			//##########################
			
			String strGeneID = gene.getGeneID();
			if(strGeneID.contains("."))
			{
				strGeneID = strGeneID.split("\\.")[0];
			}
			
			String strSQL = "SELECT * FROM junctions WHERE GENE='" + strGeneID + "';";
			Statement stmt = m_SQLConnection.createStatement();
			ResultSet pRes = stmt.executeQuery(strSQL);

			while(pRes.next())
			{
//				System.out.println(pRes.getString(2) + ":" + pRes.getInt(3) + "-" + pRes.getInt(4) + " " + pRes.getString(7));
				
				TreeMap<String, Integer> mapCounts = new TreeMap<String, Integer>();
				for(String strSample : vcSamples)
				{
					mapCounts.put(strSample, pRes.getInt(strSample));
				}
				
				vcResult.put(pRes.getInt(3) + "_" + pRes.getInt(4), mapCounts);
			}
			
			stmt.close();
		}
		catch(SQLException e)
		{
			System.out.println(e.getMessage());
			e.printStackTrace();
		}

		return vcResult;
	}

	public TreeMap<String, TreeMap<String, Integer>> GetExonCountsForGene(Gene gene)
	{				
		// for each sample (2nd key = string) store the read counts (values) per exon (1st key)
		TreeMap<String, TreeMap<String, Integer>> vcResult = new TreeMap<String, TreeMap<String, Integer>>();

		//#############################################
		//    get junction counts from the database
		//##############################################
		try
		{
			//############################ 
			//    get all sample names
			//############################
			TreeSet<String> vcSamples = GetSamples();
			
			//##########################
			//      get count data
			//########################## 
			String strSQL = "SELECT * FROM exons WHERE GENE='" + gene.getGeneID().split("\\.")[0] + "';";
			Statement stmt = m_SQLConnection.createStatement();
			ResultSet pRes = stmt.executeQuery(strSQL);

			while(pRes.next())
			{	
				TreeMap<String, Integer> mapCounts = new TreeMap<String, Integer>();
				for(String strSample : vcSamples)
				{
					mapCounts.put(strSample, pRes.getInt(strSample));
				}
				
				vcResult.put(pRes.getInt(3) + "_" + pRes.getInt(4), mapCounts);
			}
			
			stmt.close();
		}
		catch(SQLException e)
		{
			System.out.println(e.getMessage());
			e.printStackTrace();
		}
	
		return vcResult;
	}
	
	public TreeMap<String, String> GetBamFilesForSamples()
	{
		if(m_SQLConnection == null)
			return null;
		
		TreeMap<String, String> mapBamFilesForSamples = new TreeMap<String, String>();
		
		//#############################################
		//    get junction counts from the database
		//##############################################
		try
		{
			Statement stmt = m_SQLConnection.createStatement();
			String strSQL = "SELECT SAMPLE_ID,BAM_FILE FROM sample_info;";
			
			ResultSet pRes = stmt.executeQuery(strSQL);
			while(pRes.next())
			{
				mapBamFilesForSamples.put(pRes.getString(1), pRes.getString(2));
			}
			
			stmt.close();
		}
		catch(SQLException e)
		{
			System.out.println(e.getMessage());
			e.printStackTrace();
		}
		
		return mapBamFilesForSamples;
	}

	public TreeMap<String, String> GetBigWigFilesForSamples()
	{
		if(m_SQLConnection == null)
			return null;
		
		TreeMap<String, String> mapBigwigFilesForSamples = new TreeMap<String, String>();
		
		//#############################################
		//    get junction counts from the database
		//##############################################
		try
		{
			Statement stmt = m_SQLConnection.createStatement();
			String strSQL = "SELECT SAMPLE_ID,BIGWIG_FILE FROM sample_info;";
			
			ResultSet pRes = stmt.executeQuery(strSQL);
			while(pRes.next())
			{
				mapBigwigFilesForSamples.put(pRes.getString(1), pRes.getString(2));
			}
			
			stmt.close();
		}
		catch(SQLException e)
		{
			System.out.println(e.getMessage());
			e.printStackTrace();
		}
		
		return mapBigwigFilesForSamples;
	}

	public TreeSet<String> GetSamples()
	{
		if(m_SQLConnection == null)
			return null;
		
		TreeSet<String> vcSamples = new TreeSet<String>();
		
		//#############################################
		//    get junction counts from the database
		//##############################################
		try
		{
			Statement stmt = m_SQLConnection.createStatement();
			String strSQL = "SELECT SAMPLE_ID FROM sample_info;";
			
			ResultSet pRes = stmt.executeQuery(strSQL);
			while(pRes.next())
			{
				vcSamples.add(pRes.getString(1));
			}
			
			stmt.close();
		}
		catch(SQLException e)
		{
			System.out.println(e.getMessage());
			e.printStackTrace();
		}
		
		return vcSamples;
	}

	public TreeSet<String> GetConditions(String strConditionType)
	{
		if(m_SQLConnection == null)
			return null;
		
		TreeSet<String> vcConditions = new TreeSet<String>();
		
		//#############################################
		//    get junction counts from the database
		//##############################################
		try
		{
			Statement stmt = m_SQLConnection.createStatement();
			String strSQL = "SELECT "+ strConditionType +" FROM sample_info;";
					
			ResultSet pRes = stmt.executeQuery(strSQL);
			while(pRes.next())
			{
				vcConditions.add(pRes.getString(1));
			}
			
			stmt.close();
		}
		catch(SQLException e)
		{
			System.out.println(e.getMessage());
			e.printStackTrace();
		}
		
		return vcConditions;
	}
	
 	public TreeMap<String, Double> GetSizeFactors()
	{
		if(m_SQLConnection == null)
			return null;
		
		TreeMap<String, Double> vcSizeFactorsForSamples = new TreeMap<String, Double>();
		
		//#############################################
		//    get junction counts from the database
		//##############################################
		try
		{
			Statement stmt = m_SQLConnection.createStatement();
			String strSQL = "SELECT SAMPLE_ID,SIZE_FACTOR FROM sample_info;";
			
			ResultSet pRes = stmt.executeQuery(strSQL);
			while(pRes.next())
			{
				vcSizeFactorsForSamples.put(pRes.getString(1), pRes.getDouble(2));
			}
			
			stmt.close();
		}
		catch(SQLException e)
		{
			System.out.println(e.getMessage());
			e.printStackTrace();
		}
		
		return vcSizeFactorsForSamples;
	}

	public TreeMap<String, TreeSet<String>> GetSamplesPerCondition(String strQuery)
	{
		if(m_SQLConnection == null)
			return null;
		
		// key = condition, value = list of samples
		TreeMap<String, TreeSet<String>> vcSamplesAndConditions = new TreeMap<String, TreeSet<String>>();	
		try
		{
			Statement stmt = m_SQLConnection.createStatement();
			String strSQL = "SELECT SAMPLE_ID," + strQuery + " FROM sample_info;";
			
			ResultSet pRes = stmt.executeQuery(strSQL);
			while(pRes.next())
			{
				String strSample	= pRes.getString(1);
				String strCondition = pRes.getString(2);
				
				if(strCondition.isEmpty())
				{
					System.out.println("sample " + strSample + " has no information for: " + strQuery);
					continue;
				}
				
				if(vcSamplesAndConditions.containsKey(strCondition))
				{
					TreeSet<String> vcSamples = vcSamplesAndConditions.get(strCondition);
					vcSamples.add(strSample);
					vcSamplesAndConditions.put(strCondition, vcSamples);
				}
				else
				{
					TreeSet<String> vcSamples = new TreeSet<String>();
					vcSamples.add(strSample);
					vcSamplesAndConditions.put(strCondition, vcSamples);
				}
			}
			stmt.close();
		}
		catch(SQLException e)
		{
			System.out.println(e.getMessage());
			e.printStackTrace();
		}
		
		return vcSamplesAndConditions;
	}

	/*
	 * Return 1 for bigwig files, 2 for bam files. Does not report samples that miss bigwig and bam files.
	 */
	public TreeMap<String, Integer> GetCoverageModePerSample()
	{
		TreeSet<String> vcSamples 					= GetSamples();
		TreeMap<String, String> mapBigWigFiles 		= GetBigWigFilesForSamples();
		TreeMap<String, String> mapBamFiles 		= GetBigWigFilesForSamples();
		
		TreeMap<String, Integer> mapModePerSample = new TreeMap<String, Integer>();	// 1 = use bigwig file, 2 = use bam file
		
		for(String strSample : vcSamples)
		{
			// check if bigwig file exists, otherwise search for the bam file
			if(mapBigWigFiles.containsKey(strSample))
			{
				File pFile = new File(mapBigWigFiles.get(strSample));
				if(!pFile.exists())
				{
					System.out.println("Warning: specified bigwig file does not exist: " + strSample + " -> " + mapBigWigFiles.get(strSample));
					continue;
				}

				mapModePerSample.put(strSample, 1);
			}
			else if(mapBamFiles.containsKey(strSample))
			{
				File pFile = new File(mapBamFiles.get(strSample));
				if(!pFile.exists())
				{
					System.out.println("Warning: specified bigwig file does not exist: " + strSample + " -> " + mapBamFiles.get(strSample));
					continue;
				}
				mapModePerSample.put(strSample, 2);
			}
			else
			{
				System.out.println("Warning: no valid bigwig or bam file detected for sample: " + strSample);
			}			
		}
		
		return mapModePerSample;
	}
	
	void AddIsoformExpressionToDataBase(TreeMap<String, Double> mapRelativeExpressionToIsoform, String strSample) throws SQLException
	{
		if(m_SQLConnection == null)
			return;
		
		Statement stmt = m_SQLConnection.createStatement();
		
		String strSQL = "CREATE TABLE IF NOT EXISTS isoform_expression " +
   			 			"(ISOFORM_ID	TEXT	NOT NULL," +
   			 			" SAMPLE_ID		TEXT	NOT NULL," +
   			 			" VALUE			REAL 	NOT NULL," +
   			 			" PRIMARY KEY (ISOFORM_ID, SAMPLE_ID)";
		strSQL += ")";
		stmt.executeUpdate(strSQL);

		for(String strIsoform : mapRelativeExpressionToIsoform.keySet())
		{
			strSQL = "INSERT INTO isoform_expression VALUES ('";
			strSQL += strIsoform +"'";
			strSQL += ", '" + strSample + "'";
			strSQL += ", '" + mapRelativeExpressionToIsoform.get(strIsoform);
			strSQL += "');";
			
			try
			{
				stmt.executeUpdate(strSQL);
			}
			catch(SQLException e)
			{
				System.out.println("entry already existing: " + strSQL);
				System.out.println(e.getMessage());
			}
		}
		
		stmt.close();
	}
	
	TreeMap<String, Double> GetIsoformExpression(String strIsoform) throws SQLException
	{
		TreeMap<String, Double> mapExpressionToSample = new TreeMap<String, Double>();
		
		if(m_SQLConnection == null)
			return null;
		
		Statement stmt = m_SQLConnection.createStatement();
		
		String strSQL = "SELECT * FROM sqlite_master WHERE type='table' AND name='isoform_expression' LIMIT 1"; 
		ResultSet pRes = stmt.executeQuery(strSQL);
		
		if(!pRes.next())
		{
			// table has not been created yet
			return null;
		}
		
		strSQL = "SELECT * FROM isoform_expression WHERE ISOFORM_ID=='" + strIsoform + "';";
		pRes = stmt.executeQuery(strSQL);
		
		boolean bEmpty = true;
		while(pRes.next())
		{
			mapExpressionToSample.put(pRes.getString(2), pRes.getDouble(3));
			bEmpty = false;
		}
		
		stmt.close();
		
		if(bEmpty)
			return null;
		else
			return mapExpressionToSample;
	}
}