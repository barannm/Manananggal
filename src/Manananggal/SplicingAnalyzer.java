package Manananggal;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.Scanner;
import java.util.TreeMap;
import java.util.Vector;

import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import BioKit.Exon;
import BioKit.ExonAndSpliceJunctionCoverage;
import BioKit.GTFParser;
import BioKit.Gene;
import BioKit.RandomAccessGFFReader;

public class SplicingAnalyzer
{
	static SplicingWebApp app;
	
	SplicingAnalyzer() throws IOException
	{
		app = new SplicingWebApp(true);
	}
	
	public static void main(String args[]) throws Exception
	{
		app = new SplicingWebApp(true);
		
		if(args.length >= 8)
		{
			System.out.println("searching for alternatively spliced exons");

			boolean bRunInDebugMode 		= false;
			boolean bSkipFirstAndLastExons 	= false;
			String strFirstGene 			= null;
			int	nThreads					= 4;
			
			for(String strParameter : args)
			{
				if(strParameter.equals("-debug")) 				bRunInDebugMode = true;
				if(strParameter.equals("-skipFirstAndLast")) 	bSkipFirstAndLastExons = true;
				if(strParameter.startsWith("-startAt"))			strFirstGene = strParameter.split("=")[1];
				if(strParameter.startsWith("-threads"))			nThreads = Integer.parseInt(strParameter.split("=")[1]);
			}
			
			app.SetParameters(Integer.parseInt(args[3]), Integer.parseInt(args[4]), Double.parseDouble(args[5]), Double.parseDouble(args[6]), 4, args[7].trim(), args[0].trim());
			app.RunAnalysis(args[0].trim(), args[1].trim(), args[2].trim(), strFirstGene, bSkipFirstAndLastExons, bRunInDebugMode, nThreads);
		}
		else if(args[0].trim().equals("merge") && args.length > 1)
		{
			String strProject = args[1];
			
			ProjectModel project = new ProjectModel();
			
			if(args.length == 6)
			{
				int		nMergingJunctions_minSamples				= Integer.parseInt(args[2]);
				int		nMergingJunctions_minSampleCoverage 		= Integer.parseInt(args[3]);
				double	fMergingJunctions_MinAverageCount			= Double.parseDouble(args[4]);
				int		nMergingJunctions_minCoverageSingleSample 	= Integer.parseInt(args[5]);
				
				project.Init(strProject, nMergingJunctions_minSamples, nMergingJunctions_minSampleCoverage, fMergingJunctions_MinAverageCount, nMergingJunctions_minCoverageSingleSample, true);
			}
			else
				project.Init(strProject, -1, -1, -1, -1, true);
		}
		else if(args.length == 2 && args[0].trim().equals("convert"))
		{
			System.out.println("converting hits to new format");
			app.ConvertHitList(args[1].trim());
		}
		else if(args.length == 3 && args[0].trim().equals("calculate_size_factors"))
		{
			System.out.println("calculating size factors");
			app.CalculateSizeFactors(args[1].trim(), args[2].trim());
		}
		else if((args.length == 4 || args.length == 5) && args[0].trim().equals("prepare_counts"))
		{
			String strFileReference 	= args[1].trim();
			String strFileIn 			= args[2].trim();
			String strOutputPrefix 		= args[3].trim();
			
			boolean bLenient = false;
			
			if(args.length == 5 && args[4].equals("-lenient"))
				bLenient = true;

			ExonAndSpliceJunctionCoverage.computeExonAndSpliceJunctionCoverageForBAMFileAndGTFTranscriptomeAnnotation(strFileReference, strFileIn, strOutputPrefix, true, false, bLenient);
		}
		else if( args.length == 4 && args[0].trim().equals("calculate_exon_inclusions"))
		{
			String strProject	= args[1].trim();
			String strFileExons = args[2].trim();
			String strFileGTF 	= args[3].trim();
			
			CalculateSingleSampleExonInclusionRate(strProject, strFileExons, strFileGTF);
		}
		else
		{
			System.out.println("invalid number of parameters specified");
		}
	}
	
	public static void CalculateSingleSampleExonInclusionRate(String strProject, String strFileExons, String strFileGTF) throws Exception
	{
		ProjectModel project = new ProjectModel();
		project.Init(strProject, -1, -1, -1, -1, false);
		
		TreeMap<String, String> mapBigWigFiles = project.GetBigWigFilesForSamples();
		
		Scanner pIn = new Scanner(new File(strFileExons));
		PrintWriter pOut = new PrintWriter(new File("sample_inclusion_rates.tsv"));

		RandomAccessGFFReader gtf_reader = new RandomAccessGFFReader(new File(strFileGTF), new File(strFileGTF + ".idx"));
		
		pOut.write("gene_id\tgene\tchrom\tstart\tend");
		for(String strSample : project.GetSamples())
			pOut.write("\t" + strSample);
		pOut.write("\n");

		while(pIn.hasNextLine())
		{
			String strLine = pIn.nextLine();
			
			// skip header
			if(strLine.startsWith("#"))
				continue;
			
			String pSplit[] = strLine.split("\t");
			
			String strGene			= pSplit[0];
			String strGeneSym		= pSplit[1];
			String strPos			= pSplit[2];
			String strRef			= strPos.split(":")[0];
			strPos					= strPos.split(":")[1];
			int nExStart			= Integer.parseInt(strPos.split("-")[0]);
			int nExEnd				= Integer.parseInt(strPos.split("-")[1]);
			
			pOut.write(strGene + "\t" + strGeneSym + "\t" + strRef + "\t" + nExStart + "\t" + nExEnd);
			
			System.out.println("processing: " + strGene + "\t" + strRef + "\t" + nExStart + "\t" + nExEnd);
			
			Gene gene = gtf_reader.ReadGene(strGene);
			// skip genes not contained in the reference file
			if(gene == null)
				continue;
			Vector<Exon> vcConstantExons = new Vector<Exon>();
			
			// identify constant gene regions
			Iterator<Exon> itExon = gene.getExonIterator(); 
			while(itExon.hasNext())
			{
				Exon ex = itExon.next();
				
				// skip the target exon
				if(ex.getCodingStart() >= nExStart && ex.getCodingStop() <= nExEnd)
					continue;
				
				int nIsoformsContainingExon = 0;
				for(String strIsoform : gene.getArrayOfGeneProductNames())
				{
					Exon[] pExons = gene.getSortedExonsForGeneProduct(strIsoform); 
					for(Exon ex2 : pExons)
					{
						if(ex.getExonID() == ex2.getExonID())
						{
							nIsoformsContainingExon += 1.0;
						}
					}
				}
				
				if(nIsoformsContainingExon > 5 || (nIsoformsContainingExon / gene.getArrayOfGeneProductNames().length) > 0.7)
				{
					// 'constant' exon detected
					vcConstantExons.add(ex);
				}
			}
			
			// get coverage for up to 100 bases of the upstream and downstream exon
			// and the alternatively spliced exon for each sample
			int nSamples = project.GetSamples().size();
			int nCurrentSample = 0;
			for(String strSample : project.GetSamples())
			{
				// open bigwig file
				String strFile = mapBigWigFiles.get(strSample);
				
				if(nCurrentSample % 3 == 0)
				{
					System.out.print("[");
					for(int k=0; k<50; k++)
					{
						if(nCurrentSample*100 / nSamples > k*2)
							System.out.print("*");
						else
							System.out.print(" ");
					}
					System.out.print("] " + nCurrentSample + " of " + nSamples + "\r");
				}
				
				BigWigReader reader = null;
				try
				{
					reader = new BigWigReader(strFile);
				}
				catch(IOException ex)
				{
					System.out.println(ex.getMessage());
					return;
				}

				BigWigIterator it = reader.getBigWigIterator(gene.getChromosome(), gene.getStart(), gene.getChromosome(), gene.getStop(), false);

				// fill array with values
				int pValues[] = new int[gene.getStop()-gene.getStart()+1];
				while(it.hasNext())
				{
					WigItem item = it.next();
					int nIdx = item.getEndBase() - gene.getStart();						
					int nValue = (int)item.getWigValue();
					
					pValues[nIdx] = nValue;
				}
				
				int nConstantValues			= 0;
				int nAltExonValues			= 0;
				double fConstantCoverage 	= 0.0;
				double fAltExonCoverage 	= 0.0;
				
				// get coverage for constant regions
				for(Exon ex : vcConstantExons)
				{
					for(int i=ex.getCodingStart(); i<ex.getCodingStop(); i++)
					{
						double fValue = pValues[i - gene.getStart()];
						if(fValue > 5.0)
						{
							fConstantCoverage += fValue;
							nConstantValues++;
						}
					}
				}
				
				// get alt. exon coverage
				for(int i=nExStart-1; i<nExEnd; i++)
				{
					fAltExonCoverage += pValues[i - gene.getStart()];
					nAltExonValues++;
				}
				
				fAltExonCoverage /= nAltExonValues;
				fConstantCoverage /= nConstantValues;
				
				double fInclusionRate = -1.0;
				if(fConstantCoverage > 10)
					fInclusionRate = fAltExonCoverage / fConstantCoverage;
				
				pOut.write("\t" + fInclusionRate);
				nCurrentSample++;
				
				reader.close();
			}
			System.out.println("");
			pOut.write("\n");
		}

		pIn.close();
		
		pOut.flush();
		pOut.close();
	}

	public static void CalculateSingleSampleExonInclusionRate2(String strProject, String strFileExons) throws IOException
	{
		ProjectModel project = new ProjectModel();
		project.Init(strProject, -1, -1, -1, -1, false);
		
		TreeMap<String, String> mapBigWigFiles = project.GetBigWigFilesForSamples();
		
		Scanner pIn = new Scanner(new File(strFileExons));
		PrintWriter pOut = new PrintWriter(new File("sample_inclusion_rates.tsv"));
		
		pOut.write("#gene\tchrom\tstart\tend");
		for(String strSample : project.GetSamples())
			pOut.write("\t" + strSample);
		pOut.write("\n");

		while(pIn.hasNextLine())
		{
			String strLine = pIn.nextLine();
			
			// skip header
			if(strLine.startsWith("#"))
				continue;
			
			String pSplit[] = strLine.split("\t");
			
			String strGene	= pSplit[0];
			String strRef	= pSplit[1];
			int nExStart	= Integer.parseInt(pSplit[2]);
			int nExEnd		= Integer.parseInt(pSplit[3]);
			int nUpStreamStart		= Integer.parseInt(pSplit[4]);
			int nUpStreamEnd		= Integer.parseInt(pSplit[5]);
			int nDownStreamStart	= Integer.parseInt(pSplit[6]);
			int nDownStreamEnd		= Integer.parseInt(pSplit[7]);
			
			pOut.write(strGene + "\t" + strRef + "\t" + nExStart + "\t" + nExEnd);
			
			System.out.println("processing: " + strGene + "\t" + strRef + "\t" + nExStart + "\t" + nExEnd + "\t" + nUpStreamStart + "-" + nDownStreamEnd);
			
			// get coverage for up to 100 bases of the upstream and downstream exon
			// and the alternatively spliced exon for each sample
			for(String strSample : project.GetSamples())
			{
				// open bigwig file
				String strFile = mapBigWigFiles.get(strSample);
				
				BigWigReader reader = null;
				
				try
				{
					reader = new BigWigReader(strFile);
				}
				catch(IOException ex)
				{
					System.out.println(ex.getMessage());
					return;
				}

				BigWigIterator it = reader.getBigWigIterator(strRef, nUpStreamStart, strRef, nDownStreamEnd, false);

				// fill array with values
				int pValues[] = new int[nDownStreamEnd-nUpStreamStart+1];
				while(it.hasNext())
				{
					WigItem item = it.next();
					int nIdx = item.getEndBase() - nUpStreamStart;						
					int nValue = (int)item.getWigValue();
					
					pValues[nIdx] = nValue;
				}
				
				double fUpStreamCoverage 	= 0.0;
				double fDownStreamCoverage 	= 0.0;
				double fExonCoverage 		= 0.0;
				
				// get upstream coverage			
				int nUpstreamLength = nUpStreamEnd - nUpStreamStart +1;
				if(nUpstreamLength > 100)
				{
					for(int i=nUpStreamEnd-100; i<nUpStreamEnd; i++)
						fUpStreamCoverage += pValues[i - nUpStreamStart];
					
					fUpStreamCoverage /= 100.0;
				}
				else
				{
					for(int i=nUpStreamStart-1; i<nUpStreamEnd; i++)
						fUpStreamCoverage += pValues[i - nUpStreamStart];
					
					fUpStreamCoverage /= nUpstreamLength;
				}
				
				// get downstream coverage			
				int nDownstreamLength = nDownStreamEnd - nDownStreamStart +1;
				if(nDownstreamLength > 100)
				{					
					for(int i=nDownStreamStart-1; i<nDownStreamStart+100; i++)
						fDownStreamCoverage += pValues[i - nUpStreamStart];
					
					fDownStreamCoverage /= 100.0;
				}
				else
				{
					for(int i=nDownStreamStart-1; i<nDownStreamEnd; i++)
						fDownStreamCoverage += pValues[i - nUpStreamStart];
					
					fDownStreamCoverage /= nDownstreamLength;
				}
				
				// get exon coverage
				for(int i=nExStart-1; i<nExEnd; i++)
					fExonCoverage += pValues[i - nUpStreamStart];
				
				fExonCoverage /= (nExEnd - nExStart + 1);
				
				double fInclusionRate = -1.0;
				if((fUpStreamCoverage + fDownStreamCoverage)*0.5 > 20)
					fInclusionRate = fExonCoverage / ((fUpStreamCoverage + fDownStreamCoverage)*0.5);
				
				pOut.write("\t" + fInclusionRate);
			}
			
			pOut.write("\n");
		}

		pIn.close();
		
		pOut.flush();
		pOut.close();
	}
}