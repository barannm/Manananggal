package Manananggal;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.IntBuffer;
import java.util.Scanner;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.zkoss.zul.Messagebox;

import BioKit.Exon;
import BioKit.Gene;

public class ProjectModel
{
	private String m_strProjectFilePath;
	private String m_strProjectName;
	private String m_strJunctionCountPath;
	private TreeMap<String, TreeSet<String>> m_mapConditonsToConditionTypes;
	private TreeSet<Sample>	m_vcSamples;
	private boolean bIsReady;
	private TreeMap<String, RandomJunctionReader> m_mapJunctionReader; // one junction reader per condition
	
	private class Sample implements Comparable<Sample>
	{
		String m_strName;
		String m_strBigWigFile;
		String m_strJunctionFile;
		double m_fSizeFactor;
		
		// map all conditions to the corresponding condition type for this sample
		TreeMap<String, String> m_mapConditionTypeToConditions;
		
		Sample()
		{
			m_strName			= "???";
			m_strBigWigFile  	= "???";
			m_strJunctionFile 	= "???";
			m_fSizeFactor		= 0.0;
			m_mapConditionTypeToConditions = new TreeMap<String, String>();
		}
		
		@Override
		public int compareTo(Sample other)
		{
			return m_strName.compareTo(other.m_strName);
		}
	}
	
	ProjectModel()
	{
		m_strProjectFilePath			= "?";
		m_strProjectName				= "?";
		m_strJunctionCountPath			= "?";
		m_mapConditonsToConditionTypes	= new TreeMap<String, TreeSet<String>>();
		m_vcSamples 					= new TreeSet<Sample>();
		bIsReady 						= false;
		m_mapJunctionReader				= new TreeMap<String, RandomJunctionReader>();
	}

	public void clear()
	{
		m_strProjectFilePath			= "?";
		m_strProjectName				= "?";
		m_strJunctionCountPath			= "?";
		bIsReady 						= false;
		m_mapJunctionReader				= new TreeMap<String, RandomJunctionReader>();
		
		m_mapConditonsToConditionTypes.clear();
		m_vcSamples.clear();
	}
	
	public boolean Init(String strProjectFile, int nMergingJunctions_minSamples, int nMergingJunctions_minSampleCoverage, double fMergingJunctions_MinAverageCount, int nMergingJunctions_minCoverageSingleSample, boolean bLoadJunctionCounts) throws IOException
	{	
		m_mapConditonsToConditionTypes.clear();
		m_vcSamples.clear();
		
		// get project name
		File pProject 			= new File(strProjectFile);
		m_strProjectName 		= pProject.getName();
		m_strProjectFilePath 	= pProject.getAbsolutePath().substring(0, pProject.getAbsolutePath().lastIndexOf(File.separator));
		
		if(!pProject.exists())
		{
			System.out.println("file does not exist: " + strProjectFile);
			return false;
		}
		
		if(!pProject.canRead())
		{
			System.out.println("cannot read from file: " + strProjectFile);
			return false;
		}
		
		System.out.println("Initializing project: " + m_strProjectName);
		System.out.println("Project path: " + m_strProjectFilePath);
		
		Scanner pScanner = new Scanner(new File(strProjectFile));
		
		// check if a header is present
		if(!pScanner.hasNextLine())
		{
			pScanner.close();
			return false;
		}
		
		// process header
		String strHeader = pScanner.nextLine();
		
		// map columns
		TreeMap<String, Integer> mapColumns = new TreeMap<String, Integer>();
		String[] pSplit = strHeader.split("\\s+");
		for(int i=0; i<pSplit.length; i++)
		{
			String strColumn = pSplit[i];
			mapColumns.put(strColumn, i);
		}
		
		// header must contain certain key-words
		if(!mapColumns.containsKey("sample") || !mapColumns.containsKey("bigwig_files") || !mapColumns.containsKey("junction_count_files"))
		{
			System.out.println("missing required column. Project file must contain the following columns:");
			System.out.println("sample");
			System.out.println("bigwig_files");
			System.out.println("junction_count_files");
			pScanner.close();
			return false;
		}
		
		// now add the samples
		while(pScanner.hasNextLine())
		{			
			String strLine = pScanner.nextLine();
			pSplit = strLine.split("\t");
			
			Sample sample = new Sample();

			for(String strColumn : mapColumns.keySet())
			{
				switch(strColumn)
				{
					case "sample": 					sample.m_strName 			= pSplit[mapColumns.get(strColumn)].trim(); break;
					case "bigwig_files": 			sample.m_strBigWigFile 		= pSplit[mapColumns.get(strColumn)].trim(); break;
					case "junction_count_files": 	sample.m_strJunctionFile 	= pSplit[mapColumns.get(strColumn)].trim(); break;
					case "size_factors":
					{
						String strNumber = pSplit[mapColumns.get(strColumn)].trim();
						if(!strNumber.isEmpty())
							sample.m_fSizeFactor = Double.parseDouble(strNumber);
						break;
					}
					default:
					{
						sample.m_mapConditionTypeToConditions.put(strColumn, pSplit[mapColumns.get(strColumn)]);

						if(m_mapConditonsToConditionTypes.containsKey(strColumn))
						{
							m_mapConditonsToConditionTypes.get(strColumn).add(pSplit[mapColumns.get(strColumn)].trim());
						}
						else
						{
							TreeSet<String> vcTmp = new TreeSet<String>();
							vcTmp.add(pSplit[mapColumns.get(strColumn)].trim());
							m_mapConditonsToConditionTypes.put(strColumn, vcTmp);
						}
						break;
					}
				}
			}
			
			if(m_strJunctionCountPath.equals("?"))
			{
				File pFile = new File(sample.m_strJunctionFile);
				String strFullPath = pFile.getAbsolutePath();
				m_strJunctionCountPath = strFullPath.substring(0, strFullPath.lastIndexOf(File.separatorChar));
				System.out.println("junction count path: "+ m_strJunctionCountPath);
			}
			
			m_vcSamples.add(sample);
		}

		if(bLoadJunctionCounts)
		{
			// use first condition type by default
			String strConditionType = m_mapConditonsToConditionTypes.firstKey();
			
			for(String strCondition : this.m_mapConditonsToConditionTypes.get(strConditionType))
			{
				// check if the merged and indexed junction count table exists				
				String strMergedJunctionCountFile = m_strJunctionCountPath + "/" + strCondition + ".merged_junction_counts.dat";

				File pFile = new File(strMergedJunctionCountFile);
				if(!pFile.exists())
				{
					System.out.println(strMergedJunctionCountFile + " does not exist. Creating merged table.");
					
					// use default settings if not otherwise specified
					if(nMergingJunctions_minSamples < 0)
						nMergingJunctions_minSamples = 10;
					
					if(nMergingJunctions_minSampleCoverage < 0)
					nMergingJunctions_minSampleCoverage = 1;
					
					if(fMergingJunctions_MinAverageCount < 0)
						fMergingJunctions_MinAverageCount = 1;
					
					if(nMergingJunctions_minCoverageSingleSample < 0)
						nMergingJunctions_minCoverageSingleSample = 10;
					
					// create the merged file and index it
					MergeJunctionTable(strCondition, strConditionType, strMergedJunctionCountFile, nMergingJunctions_minSamples, nMergingJunctions_minSampleCoverage, fMergingJunctions_MinAverageCount, nMergingJunctions_minCoverageSingleSample);
				}
				
				// if the file still does not exist, we got a problem
				if(!pFile.exists())
				{
					Messagebox.show("Failed to get junction count table");
					System.out.println("Failed to get junction count table");
					pScanner.close();
					return false;
				}
				
				// read index file
				RandomJunctionReader reader = new RandomJunctionReader();
				reader.Init(strMergedJunctionCountFile);
				m_mapJunctionReader.put(strCondition, reader);
			}
		}
		
		bIsReady = true;
		
		pScanner.close();
		return true;
	}

	public boolean IsReady()
	{
		return bIsReady;
	}
	
	public String GetProjectName()
	{
		return m_strProjectName;
	}
	
	public void MergeJunctionTable(String strCondition, String strConditionType, String strMergedFile, int nMergingJunctions_minSamples, int nMergingJunctions_minSampleCoverage, double fMergingJunctions_MinAverageCount, int nMergingJunctions_minCoverageSingleSample) throws IOException
	{		
		TreeMap<CountElement, int[]> mapCounts = new TreeMap<CountElement, int[]>();

		DataOutputStream pOutData = null;
//		RandomAccessFile pOutData;
		try
		{
//			pOutData = new RandomAccessFile(strMergedFile, "rw");
			pOutData = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(strMergedFile)));
		}
		catch (FileNotFoundException e1)
		{
			System.out.println("failed to open merged file: " + strMergedFile);
			e1.printStackTrace();
			return;
		}
		
		RandomAccessFile pOutIdx;
		try
		{
			pOutIdx = new RandomAccessFile(strMergedFile + ".idx", "rw");
		}
		catch (FileNotFoundException e1)
		{
			System.out.println("failed to open index file: " + strMergedFile + ".idx");
			e1.printStackTrace();
			pOutData.close();
			return;
		}
		
		System.out.println("merging junction counts for " + m_vcSamples.size() + " samples.");
		
		// get list of valid samples
		TreeSet<String> vcValidSamples = new TreeSet<String>();	
		TreeSet<String> vcSamples = GetSamplesPerCondition(strConditionType).get(strCondition);
		
		int nIdx = 0;
		for(Sample sample : m_vcSamples)
		{
			// skip samples of other indications
			if(!vcSamples.contains(sample.m_strName))
				continue;
			
			Scanner pScanner = null;
			try
			{
				pScanner = new Scanner(new File(sample.m_strJunctionFile));
			}
			catch (FileNotFoundException e1)
			{
				System.out.println("failed to open input file: " + sample.m_strJunctionFile);
				e1.printStackTrace();
				pOutData.close();
				pOutIdx.close();
				
				File tmp = new File(strMergedFile);
				tmp.delete();
				
				tmp = new File(strMergedFile + ".idx");
				tmp.delete();
				return;
			}
			
			System.out.println("processing file: " + sample.m_strJunctionFile);
			
			// track which samples were used in the output file
			vcValidSamples.add(sample.m_strName);
			
			while(pScanner.hasNextLine())
			{
				String strLine = pScanner.nextLine();
				
				// skip gene and exon entries
				if(strLine.startsWith("Gene") || strLine.startsWith("Exon"))
					continue;
				
				CountElement e = new CountElement();
				e.ParseLine(strLine);
				
				if(mapCounts.containsKey(e))
				{
					int pValues[] = mapCounts.get(e);
					pValues[nIdx] = e.m_vcCounts.firstElement();
				}
				else
				{
					int pValues[] = new int[vcSamples.size()];
					pValues[nIdx] = e.m_vcCounts.firstElement();
					mapCounts.put(e, pValues);
				}
			}
			
			nIdx++;
			pScanner.close();
		}
		
		System.out.println("done reading input files");
		
		// write number of samples and sample name to output file
		pOutData.writeInt(vcValidSamples.size());
		for(String strSample : vcValidSamples)
		{
			pOutData.writeUTF(strSample);
		}
		
		int nItemsWritten = 0;
		String strCurrentReference = "?";
		
		int nFirstItemStart = -1;
		int nLastItemEnd	= 0;
		long nOffset		= 0;
		
		for(CountElement e : mapCounts.keySet())
		{
			if(nFirstItemStart == -1)
			{
				nOffset 				= pOutData.size();
//				nOffset = pOutData.getFilePointer();
				strCurrentReference 	= e.m_strRef;
				nFirstItemStart 		= e.m_nStart;
			}
			
			//#######################################
			//     get array with coverage counts
			//#######################################
			int nSamplesWithEnoughCoverage = 0;
			int pValues[] = mapCounts.get(e);
			
			int nMax = 0;
			int nSum = 0;
			nIdx = 0;
			
//			long nStartTime = System.currentTimeMillis();
			
			for(int i=0; i<pValues.length; i++)
			{
				int nValue = pValues[i];
				
				if(nValue > nMergingJunctions_minSampleCoverage)
					nSamplesWithEnoughCoverage++;
				
				if(nValue > nMax)
					nMax = nValue;
				nSum += nValue;
			}
	
			//#############################################################
			//    test whether overall junction coverage is sufficient
			//#############################################################
			double fMean = nSum / pValues.length;
			
			// valid junctions require at least one sample with 10 reads
			if(nMax < nMergingJunctions_minCoverageSingleSample)
			{
				// or an average of at least 1 across all samples
				if(fMean < fMergingJunctions_MinAverageCount)
				{
					// or 10 samples with at least 1 read
					if(pValues.length >= nMergingJunctions_minSamples)
					{
						if(nSamplesWithEnoughCoverage < nMergingJunctions_minSamples)
							continue;
						// else it is okay
					}
					else
					{
						// skip this junction
						continue;
					}
				}
				// else it is okay
			}

			//#################################################################
			//    write index every 500 items or upon reference name change
			//#################################################################
			if((nItemsWritten % 500 == 0 || !strCurrentReference.equals(e.m_strRef) || e == mapCounts.lastEntry()) && nItemsWritten != 0)
			{				
				byte pBytes[] = strCurrentReference.getBytes("UTF-8");
				pOutIdx.writeInt(pBytes.length);
				pOutIdx.write(pBytes);
				pOutIdx.writeInt(nFirstItemStart);
				pOutIdx.writeInt(nLastItemEnd);
				pOutIdx.writeLong(nOffset);
	
				nOffset = pOutData.size();
//				nOffset = pOutData.getFilePointer();
				strCurrentReference = e.m_strRef;
				nFirstItemStart = e.m_nStart;
				nItemsWritten = 0;
			}
			nLastItemEnd	= e.m_nStart;

			//###################
			//    write data
			//###################
			
			// junction information
			pOutData.writeUTF(e.m_strRef);		// x byte
			pOutData.writeInt(e.m_nStart);		// 4 bytes
			pOutData.writeInt(e.m_nEnd);		// 4 bytes
			pOutData.writeBoolean(e.m_bKnown);	// 1 byte
			pOutData.writeUTF(e.GetGeneID());	// x byte
			pOutData.writeChar(e.m_chStrand);	// 2 bytes

			ByteBuffer	pBuffer 	= ByteBuffer.allocate(4*pValues.length);
			IntBuffer 	pIntBuffer	= pBuffer.asIntBuffer();
			pIntBuffer.put(pValues);	
			pOutData.write(pBuffer.array());
			
			nItemsWritten++;
		}
		
		pOutData.close();
		pOutIdx.close();
	}

	public TreeSet<String> GetSamples()
	{
		TreeSet<String> vcRes = new TreeSet<String>();
		
		for(Sample sample : m_vcSamples)
			vcRes.add(sample.m_strName);
		
		return vcRes;
	}
	
	public TreeMap<String, TreeSet<String>> GetSamplesPerCondition(String strConditionType)
	{
		TreeMap<String, TreeSet<String>> mapSamplesToConditions = new TreeMap<String, TreeSet<String>>();

		for(Sample sample : m_vcSamples)
		{
			String strCondition = sample.m_mapConditionTypeToConditions.get(strConditionType);
			
			if(mapSamplesToConditions.containsKey(strCondition))
			{
				TreeSet<String> vcTmp = mapSamplesToConditions.get(strCondition);
				vcTmp.add(sample.m_strName);
			}
			else
			{
				TreeSet<String> vcTmp = new TreeSet<String>();
				vcTmp.add(sample.m_strName);
				mapSamplesToConditions.put(strCondition, vcTmp);
			}
		}
		
		return mapSamplesToConditions;
	}

	public TreeMap<String, TreeSet<String>> GetConditionsToConditionTypes()
	{
		return m_mapConditonsToConditionTypes;
	}
	
	public TreeSet<String> GetConditions(String strConditionType)
	{
		if(!m_mapConditonsToConditionTypes.containsKey(strConditionType))
			return null;
		
		return m_mapConditonsToConditionTypes.get(strConditionType);
	}
	
	public TreeMap<String, Double> GetSizeFactors()
	{
		TreeMap<String, Double> mapSizeFactorsToSamples = new TreeMap<String, Double>();
		for(Sample sample : m_vcSamples)
			mapSizeFactorsToSamples.put(sample.m_strName, sample.m_fSizeFactor);
		
		return mapSizeFactorsToSamples;
	}
	
	public void WriteSizeFactors() throws FileNotFoundException
	{
		File pIn = new File(m_strProjectFilePath + "/" + m_strProjectName);
		Scanner pScanner = new Scanner(pIn);
		PrintWriter pOut = new PrintWriter(new File(m_strProjectFilePath + "/" + m_strProjectName + ".tmp"));
		
		// check if a header is present
		if(!pScanner.hasNextLine())
		{
			pScanner.close();
			return;
		}
		
		// process header
		String strHeader = pScanner.nextLine();
		
		// map columns
		TreeMap<String, Integer> mapColumns = new TreeMap<String, Integer>();
		String[] pSplit = strHeader.split("\t");
		for(int i=0; i<pSplit.length; i++)
		{
			String strColumn = pSplit[i];
			mapColumns.put(strColumn, i);
		}

		// check whether there is a sizeFactor column
		int nIdx = -1;
		if(mapColumns.containsKey("size_factors"))
		{
			nIdx = mapColumns.get("size_factors");
			pOut.println(StringUtils.join(pSplit, '\t'));
		}
		else
		{
			pOut.println(StringUtils.join(pSplit, '\t') + "\t" + "size_factors");
		}
		
		while(pScanner.hasNextLine())
		{
			String strLine = pScanner.nextLine();
			
			pSplit = strLine.split("\t");
			
			String strSample = pSplit[mapColumns.get("sample")];
			
			double fSizeFactor = 0.0;
			for(Sample s : m_vcSamples)
			{
				if(s.m_strName.equals(strSample))
				{
					fSizeFactor = s.m_fSizeFactor;
					break;
				}
			}
			
			if(nIdx != -1)
			{				
				pSplit[nIdx] = Double.toString(fSizeFactor);
				pOut.println(StringUtils.join(pSplit, '\t'));
			}
			else
			{
				pOut.println(StringUtils.join(pSplit, '\t') + "\t" + fSizeFactor);
			}
		}
		
		pScanner.close();
		pOut.close();
		
		File pFileOld = new File(m_strProjectFilePath + "/" + m_strProjectName + ".tmp");
		File pFileNew = new File(m_strProjectFilePath + "/" + m_strProjectName);
		
		if(pFileNew.exists())
			pFileNew.delete();

		if(!pFileOld.renameTo(pFileNew))
		{
			System.out.println("failed to rename size factor file");
		}
	}
	
	public void SetSizeFactor(String strSample, double fSizeFactor)
	{
		for(Sample s : m_vcSamples)
		{
			if(s.m_strName.equals(strSample))
				s.m_fSizeFactor = fSizeFactor;
		}
	}
	
	public TreeMap<String, String> GetBigWigFilesForSamples()
	{
		TreeMap<String, String> mapBigWigFilesToSamples = new TreeMap<String, String>();
		for(Sample sample : m_vcSamples)
			mapBigWigFilesToSamples.put(sample.m_strName, sample.m_strBigWigFile);
		
		return mapBigWigFilesToSamples;
	}

	public TreeMap<String, TreeMap<String, Integer>> GetJunctionCountsForGene(Gene gene) throws IOException
	{
		TreeSet<String> vcConditions = m_mapConditonsToConditionTypes.get(m_mapConditonsToConditionTypes.firstKey());
		
		TreeMap<String, TreeMap<String, Integer>> mapRes = new TreeMap<String, TreeMap<String, Integer>>(); 
		
		for(String strCondition : vcConditions)
		{
			TreeMap<String, TreeMap<String, Integer>> mapTmp = m_mapJunctionReader.get(strCondition).GetJunctionsForGene(gene);
			
			for(String strJunction : mapTmp.keySet())
			{
				if(mapRes.containsKey(strJunction))
				{
					mapRes.get(strJunction).putAll(mapTmp.get(strJunction));
				}
				else
				{
					mapRes.put(strJunction, mapTmp.get(strJunction));
				}
			}
		}
		
		return(mapRes);
	}
	
	public TreeMap<String, TreeMap<String, Integer>> GetJunctionsForRange(String strRef, int nStart, int nEnd) throws IOException
	{
		TreeSet<String> vcConditions = m_mapConditonsToConditionTypes.get(m_mapConditonsToConditionTypes.firstKey());
		
		TreeMap<String, TreeMap<String, Integer>> mapRes = new TreeMap<String, TreeMap<String, Integer>>(); 
		
		for(String strCondition : vcConditions)
		{
			TreeMap<String, TreeMap<String, Integer>> mapTmp = m_mapJunctionReader.get(strCondition).GetJunctionsForRange(strRef, nStart, nEnd, 0);
			
			for(String strJunction : mapTmp.keySet())
			{
				if(mapRes.containsKey(strJunction))
				{
					mapRes.get(strJunction).putAll(mapTmp.get(strJunction));
				}
				else
				{
					mapRes.put(strJunction, mapTmp.get(strJunction));
				}
			}
		}
		
		return(mapRes);
	}
	
	public TreeMap<CountElement, TreeMap<String, Integer>> GetJunctionsForRangeAsCountElements(String strRef, int nStart, int nEnd) throws IOException
	{
		TreeSet<String> vcConditions = m_mapConditonsToConditionTypes.get(m_mapConditonsToConditionTypes.firstKey());
		
		TreeMap<CountElement, TreeMap<String, Integer>> mapRes = new TreeMap<CountElement, TreeMap<String, Integer>>();
		
		for(String strCondition : vcConditions)
		{
			TreeMap<CountElement, TreeMap<String, Integer>> mapTmp = m_mapJunctionReader.get(strCondition).GetJunctionsForRangeAsCountElements(strRef, nStart, nEnd, 0);
			
			for(CountElement junction : mapTmp.keySet())
			{
				if(mapRes.containsKey(junction))
				{
					mapRes.get(junction).putAll(mapTmp.get(junction));
				}
				else
				{
					mapRes.put(junction, mapTmp.get(junction));
				}
			}
		}

		return(mapRes);
	}
	
	// returns a junction tree with all junctions and exons connecting the given range
	public JunctionTree GetJunctionsForRangeAsTree(Exon[] pExons, String strRef, int nStart, int nEnd, int nMinJunCoverage, String strConditionType) throws IOException
	{
		TreeMap<String, TreeSet<String>> mapSamplesToConditions = GetSamplesPerCondition(strConditionType);
		
		// first of all, get all junctions in the range		
		TreeMap<CountElement, TreeMap<String, Integer>> mapJunctionCounts = new TreeMap<CountElement, TreeMap<String, Integer>>();
		
		for(String strCondition : mapSamplesToConditions.keySet())
		{
			TreeMap<CountElement, TreeMap<String, Integer>> mapTmp = m_mapJunctionReader.get(strCondition).GetJunctionsForRangeAsCountElements(strRef, nStart, nEnd, 0);
			
			for(CountElement junction : mapTmp.keySet())
			{
				if(mapJunctionCounts.containsKey(junction))
				{
					mapJunctionCounts.get(junction).putAll(mapTmp.get(junction));
				}
				else
				{
					mapJunctionCounts.put(junction, mapTmp.get(junction));
				}
			}
		}
		
		//##################### prepare coverage array ################
		//TODO valid samples
		TreeMap<String, double[]> mapCoverageToConditions = new TreeMap<String, double[]>();
		for(String strCondition : mapSamplesToConditions.keySet())
		{
			int nSamples = mapSamplesToConditions.get(strCondition).size();
			mapCoverageToConditions.put(strCondition, new double[nSamples]);
		}

		// find invalid junctions
		TreeSet<CountElement> vcInvalidJunctions = new TreeSet<CountElement>();
		for(CountElement jun : mapJunctionCounts.keySet())
		{
			int nInvalid = 0;
			for(String strCondition : mapSamplesToConditions.keySet())
			{							
				int nSampleIdx = 0;
				for(String strSample : mapSamplesToConditions.get(strCondition))
				{
					mapCoverageToConditions.get(strCondition)[nSampleIdx] = mapJunctionCounts.get(jun).get(strSample);			
					nSampleIdx++;
				}
				double fMean = StatUtils.mean(mapCoverageToConditions.get(strCondition));
				
				if(fMean < nMinJunCoverage)
					nInvalid++;
			}
			
			if(nInvalid == mapSamplesToConditions.size())
			{
//				System.out.println("invalid junction: " + jun + " " + mapJunctionCounts.get(jun));
				vcInvalidJunctions.add(jun);
			}
		}
		
		for(CountElement jun : vcInvalidJunctions)
		{
			mapJunctionCounts.remove(jun);
		}
		
		// now combine the exons and junction counts to one tree
		JunctionTree junctionTree = new JunctionTree(mapJunctionCounts, pExons);
		
		return junctionTree;
	}

	public TreeMap<String, TreeMap<String, Double>> GetIsoformExpressions(String strGene, String strOutputPath) throws IOException
	{
		// <isoform, <sample, value>>
		TreeMap<String, TreeMap<String, Double>> mapExpressions = new TreeMap<String, TreeMap<String, Double>>();
		
		// check if mmseq results exist
		String strMMseqFile = strOutputPath + "/" + m_strProjectName + ".mmseq";
		File pFileData = new File(strMMseqFile);
		File pFileIdx  = new File((strMMseqFile + ".idx"));
		
		// check whether the input files exist
		if(!pFileIdx.exists())
		{
			System.out.println("Failed to open mmseq index file -> " + strMMseqFile +".idx");
			return null;
		}
		
		if(!pFileData.exists())
		{
			System.out.println("Failed to open mmseq data file -> " + strMMseqFile);
			return null;
		}
		
		RandomAccessFile pInData = new RandomAccessFile(pFileData, "r");
		Scanner pInIdx = new Scanner(pFileIdx);
		
		// read header
		TreeMap<String, Integer> mapSampleToColumn = new TreeMap<String, Integer>();
		String strHeader = pInData.readUTF();
		
		String pSplit[] = strHeader.split("\t");
		for(int i=2; i<pSplit.length; i++)
		{
			mapSampleToColumn.put(pSplit[i], i);
		}
		
		// get index
		System.out.println("getting index");
		long nFileOffset = -1;
		while(pInIdx.hasNextLine())
		{
			String strLine = pInIdx.nextLine();
			pSplit = strLine.split("\t");
			
			if(pSplit[0].equals(strGene))
			{
				System.out.println(strLine);
				nFileOffset = Long.parseLong(pSplit[1]);
				break;
			}
		}
		
		if(nFileOffset == -1)
		{
			pInData.close();
			pInIdx.close();
			return null;
		}
		
		// get data
		System.out.println("getting data");
		pInData.seek(nFileOffset);
		while(pInData.getFilePointer() < pInData.length())
		{
			String strLine = pInData.readUTF();
			pSplit = strLine.split("\t");
			
			if(!pSplit[0].equals(strGene))
				break;
			
			TreeMap<String, Double> mapExpressionToSamples = new TreeMap<String, Double>();
			
			for(String strSample : mapSampleToColumn.keySet())
			{
				mapExpressionToSamples.put(strSample, Double.parseDouble(pSplit[mapSampleToColumn.get(strSample)]));
			}
			mapExpressions.put(pSplit[1], mapExpressionToSamples);
		}
		
		pInIdx.close();
		pInData.close();
		return mapExpressions;
	}

	public void AddIsoformExpressionToDataBase(TreeMap<String, TreeMap<String, Double>> mapIsoformExpression, String strGeneID, String strOutputPath) throws IOException
	{
		RandomAccessFile pOutData = null;
		PrintWriter pOutIdx = null;
		
		// check if mmseq data for the gene exists already
		File pInIdx = new File(strOutputPath + "/" + m_strProjectName + ".mmseq.idx");
		if(pInIdx.exists())
		{
			Scanner pScanner = new Scanner(pInIdx);
			while(pScanner.hasNextLine())
			{
				String strLine = pScanner.nextLine();
				String pSplit[] = strLine.split("\t");
				if(pSplit[0].equals(strGeneID))
				{
					pScanner.close();
					return;
				}
			}
			pScanner.close();
		}
		
		File pData = new File(strOutputPath + "/" + m_strProjectName + ".mmseq");
		
		boolean bRequiresHeader = false;
		if(!pData.exists())
			bRequiresHeader = true;
		
		try
		{
			pOutData = new RandomAccessFile(pData, "rw");
			pOutIdx  = new PrintWriter(new FileOutputStream(new File(strOutputPath + "/" + m_strProjectName + ".mmseq.idx"), true));
		}
		catch (FileNotFoundException e)
		{
			System.out.println(e.getMessage());
			e.printStackTrace();
			return;
		}
		
		// go to end of file
		pOutData.seek(pOutData.length());
		long nPos = pOutData.getFilePointer();

		for(String strIsoform : mapIsoformExpression.keySet())
		{
			// write header
			String strLine = "#gene_id\ttrnascript_id";
			if(bRequiresHeader)
			{
				for(String strSample : mapIsoformExpression.get(strIsoform).keySet())
				{
					strLine += "\t" + strSample;
				}				
				bRequiresHeader = false;
				pOutData.writeUTF(strLine);
				nPos = pOutData.getFilePointer();
			}
			
			strLine = strGeneID + "\t" + strIsoform;
			for(String strSample : mapIsoformExpression.get(strIsoform).keySet())
			{
				strLine += "\t" + mapIsoformExpression.get(strIsoform).get(strSample);
			}
			pOutData.writeUTF(strLine);
		}

		pOutIdx.println(strGeneID + "\t" + nPos);
		
		pOutData.close();
		pOutIdx.close();
	}

	public TreeSet<AlternativeSplicingHit> LoadHitListOld(String strPath) throws FileNotFoundException
	{
		// open input file
		String strFile = strPath += "/" + m_strProjectName + "_hits.txt";
		
		File pFile = new File(strFile);
		if(!pFile.exists())
		{
//			Messagebox.show("File not found: "+ strFile);
			System.out.println("File not found: "+ strFile);
			return null;
		}
		
		// parse hits
		TreeSet<AlternativeSplicingHit> vcHits = new TreeSet<AlternativeSplicingHit>();
		Scanner pScanner = new Scanner(new File(strFile));
		while(pScanner.hasNextLine())
		{
			String strLine = pScanner.nextLine();
			
			AlternativeSplicingHit hit = new AlternativeSplicingHit();
			hit.ParseString(strLine);
			
			vcHits.add(hit);
		}
		pScanner.close();
		
		return vcHits;
	}
	
	public TreeSet<AlternativeSplicingHit> LoadHitList(String strPath) throws IOException
	{
		// open input file
		String strFile = strPath += "/" + m_strProjectName + "_hits.dat";
		
		File pFile = new File(strFile);
		if(!pFile.exists())
			return null;
		
		TreeSet<AlternativeSplicingHit> vcHits = new TreeSet<AlternativeSplicingHit>();
		
		RandomAccessFile pIn = new RandomAccessFile(strFile, "r");
		while(pIn.getFilePointer() < pIn.length())
		{
			AlternativeSplicingHit hit = new AlternativeSplicingHit();
			if(hit.ReadFromFile(pIn))
			{
				vcHits.add(hit);
			}
			else
			{
				System.out.println("failed to read hit from hit list -> potential file error");
			}
		}
		pIn.close();
		
		return vcHits;
	}
	
	public void AddHitToHitList(int nRating, GeneIdentifier gid, int nMinJunctionReads, int nMinCovPerBase, double fMinCoveredBases, double fVariableExonThreshold, TreeSet<String> vcSelectedIsoforms, String strFileGTF, String strComment, String strPath) throws FileNotFoundException
	{
		String strFile = strPath + "/" + m_strProjectName + "_hits.txt";
		File pFileIn = new File(strFile);
		
		// open output file
		String strFileTmp = strPath + "/" + m_strProjectName + "_hits.tmp";
		PrintWriter pTmp = new PrintWriter(new FileOutputStream(strFileTmp, true));
		
		// add all previous hits to the output file, except for hits matching the current gene id
		if(pFileIn.exists())
		{
			Scanner pIn = new Scanner(pFileIn);
			
			while(pIn.hasNextLine())
			{
				String strLine = pIn.nextLine();
				
				String pSplit[] = strLine.split("\t");
				String strID = pSplit[1];
				
				if(strID.equals(gid.m_strEnsemblGeneID))
					continue;
				
				pTmp.println(strLine);
			}
			
			pIn.close();
		}
		
		// add new result
		pFileIn.delete();
		String strOut = nRating + "\t" + gid.m_strEnsemblGeneID + "\t" + strFileGTF + "\t" + nMinJunctionReads + "\t" + nMinCovPerBase + "\t" + fMinCoveredBases + "\t" + fVariableExonThreshold + "\t" + strComment;
		for(String strIsoform : vcSelectedIsoforms)
			strOut += "\t" + strIsoform;
		
		pTmp.println(strOut);
		pTmp.close();
		
		File pOut = new File(strFileTmp);
		pOut.renameTo(pFileIn);
	}
	
	public void SaveHitListToFile(String strPath, TreeSet<AlternativeSplicingHit> vcHits) throws IOException
	{
		if(vcHits == null)
			return;
		
		String strFile = strPath + "/" + m_strProjectName + "_hits.dat";
		File pFile = new File(strFile);
		
		if(pFile.exists())
			pFile.delete();
		
		RandomAccessFile pOut = new RandomAccessFile(strFile, "rw");
//		PrintWriter pOut = new PrintWriter(new FileOutputStream(strFile, true));
		
		for(AlternativeSplicingHit hit : vcHits)
		{
			hit.WriteToFile(pOut);
		}

		pOut.close();
	}

	public void SaveSplicingHitsToFile(String strPath, TreeSet<SimpleSpliceScore> vcHits) throws IOException
	{
		String strFileData = strPath + "/" + m_strProjectName + "_splicing_hits.dat";
		String strFileIdx  = strPath + "/" + m_strProjectName + "_splicing_hits.idx";
		
		//###################################################################
		//    check if the results for this gene have been stored already
		//###################################################################
		
		String strGeneID = vcHits.first().m_strGeneID.split("\\.")[0];
		
		long nFileOffset = -1;
		File pFileIdx = new File(strFileIdx);		
		if(pFileIdx.exists())
		{
			Scanner pScanner = new Scanner(pFileIdx);
			while(pScanner.hasNext())
			{
				String strLine = pScanner.nextLine();
				
				String pSplit[] = strLine.split("\t");

				if(pSplit[0].equals(strGeneID))
				{
					nFileOffset = Long.parseLong(pSplit[1]);
					break;
				}
			}
			pScanner.close();
		}
		
		if(nFileOffset == -1)
		{
			for(SimpleSpliceScore hit : vcHits)
			{
				//################################
				//        add new result
				//################################
	
				// open data file for appending
				RandomAccessFile pOutData = new RandomAccessFile(strFileData, "rw");
				
				// override existing entry
				if(nFileOffset != -1)
				{
					pOutData.seek(nFileOffset);
				}
				else // otherwise add to end of file
				{
					pOutData.seek(pOutData.length());
					
					PrintWriter pOutIdx = new PrintWriter(new FileOutputStream(strFileIdx, true));
					pOutIdx.println(hit.m_strGeneID.split("\\.")[0] + "\t" + pOutData.length());
					pOutIdx.close();
				}
				
				// write gene information
				hit.WriteToFile(pOutData);
				
				pOutData.close();
			}
		}
	}
	
	public TreeSet<SimpleSpliceScore> GetSplicingHitForGene(String strPath, Gene gene) throws IOException
	{
		TreeSet<SimpleSpliceScore> vcResults = new TreeSet<SimpleSpliceScore>();

		String strFileData = strPath + "/" + m_strProjectName + "_splicing_hits.dat";
		String strFileIdx  = strPath + "/" + m_strProjectName + "_splicing_hits.idx";
		
		String strGeneID = gene.getGeneID().split("\\.")[0];
		
		// get file offset
		long nFileOffset = -1;
		File pFileIdx = new File(strFileIdx);		
		if(pFileIdx.exists())
		{
			Scanner pScanner = new Scanner(pFileIdx);
			while(pScanner.hasNext())
			{
				String strLine = pScanner.nextLine();
				
				String pSplit[] = strLine.split("\t");
				
				if(pSplit[0].equals(strGeneID))
				{					
					nFileOffset = Long.parseLong(pSplit[1]);
					break;
				}
			}
			pScanner.close();
		}
		
		if(nFileOffset == -1)
			return null;
		
		RandomAccessFile pFile = new RandomAccessFile(strFileData, "r");
		pFile.seek(nFileOffset);

		while(pFile.getFilePointer() < pFile.length())
		{
			SimpleSpliceScore hit = new SimpleSpliceScore();
			hit.ReadFromFile(pFile);
			
			if(hit.m_strGeneID.equals(gene.getGeneID()))
			{
				vcResults.add(hit);
			}
		}
		
		pFile.close();

		return vcResults;
	}

}
