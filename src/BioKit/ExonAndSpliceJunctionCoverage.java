package BioKit;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.Scanner;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

public class ExonAndSpliceJunctionCoverage 
{
	public static void main(String[] args) throws Exception
	{
		if(args.length != 3 && args.length != 4)
		{
			System.out.println("USAGE: ExonAndSpliceJunctionCoverage <gtfFile> <bamFile> <output prefix> [-multi]");
			System.out.println("Specifying the optional parameter -multi results in separate output files for"
					+ "splice junctions and exon read counts.");
			System.out.println("ExonAndSpliceJunctionCoverage -convertToDEXSeqInput <input> <output>");
			System.out.println("This method converts the standard exon count file (from split output) to DEXSeq compatible input.");
			System.out.println("ExonAndSpliceJunctionCoverage -convertFeatureCountsToDEXSeq <input_file> <output_file>");
			System.out.println("This method converts featureCounts exon output to DEXSeq exon counts");
			return;
		}
		
		if(args.length == 3 && args[0].equals("-convertToDEXSeqInput"))
		{
			exonCountsToDEXSeqInput(args[1], args[2]);
			return;
		}
		
		if(args.length == 3 && args[0].equals("-convertFeatureCountsToDEXSeq"))
		{
			convertFeatureCountsToDEXSeq(args[1], args[2]);
			return;
		}
		
		boolean bMultiFileOutput		= false;
		boolean bAdditionalFileOutput	= false;
		for(String strArg : args)
		{
			if(strArg.equals("-multi"))
				bMultiFileOutput = true;
			
			if(strArg.equals("-additional"))
			{
				bAdditionalFileOutput = true;
				bMultiFileOutput	  = true;
			}
		}
		
		if(args.length == 4 || args.length == 3)
			computeExonAndSpliceJunctionCoverageForBAMFileAndGTFTranscriptomeAnnotation(args[0], args[1], args[2], bMultiFileOutput, bAdditionalFileOutput, false);
	}
	
	public static void convertFeatureCountsToDEXSeq(String strFileIn, String strFileOut) throws FileNotFoundException
	{
		Scanner pScanner = new Scanner(new File(strFileIn));
		PrintWriter pWriter = new PrintWriter(new File(strFileOut));
		
		String strLastGeneID = "?";
		int nExonIdx = 1;
		while(pScanner.hasNextLine())
		{
			String strLine = pScanner.nextLine();
			
			// ignore comment lines
			if(strLine.startsWith("#") || strLine.startsWith("Geneid"))
				continue;
			
			String split[] = strLine.split("\\s+");
			String strGeneID = split[0];
			
			if(strGeneID.equals(strLastGeneID))
			{
				nExonIdx += 1;
			}
			else
			{
				strLastGeneID = strGeneID;
				nExonIdx = 1;
			}
			
			String strExonIdx = "???";
			if(nExonIdx < 10)				
				strExonIdx = ":E00" + Integer.toString(nExonIdx);
			else if(nExonIdx < 100)
				strExonIdx = ":E0" + Integer.toString(nExonIdx);
			else
				strExonIdx = ":E" + Integer.toString(nExonIdx);
			
			pWriter.print(strGeneID + strExonIdx + "\t" + split[6] + "\n");
		}
		
		pScanner.close();
		pWriter.close();
	}
	
	//##############################################################################################################
	//     bMultiFileOutput true/false        = generate separate files for junction and exon read counts
	//     bAdditionalFileOutput true/false	  = additionally generate intron and paired-end junction read counts
	//##############################################################################################################
	public static void computeExonAndSpliceJunctionCoverageForBAMFileAndGTFTranscriptomeAnnotation(String gtfFile, String bamFile, String outputPrefix, boolean bMultiFileOutput, boolean bAdditionalFileOutput, boolean bLenient) throws Exception
	{		
		HashMap<String, GenomeReadMappingReport> chromosomeToMappingReport = new HashMap<String, GenomeReadMappingReport>();
		
		System.out.println("using GTF file: " + gtfFile);
		System.out.println("using bam file: " + bamFile);
		
		//#############################
		//    create output file(s)
		//#############################
		PrintWriter pWriterExons 	 = null;	// this writer will be used for combined output as well
		PrintWriter pWriterJunctions = null;	// junction count only output file
		PrintWriter pWriterIntrons	 = null;	// intron count output file
		PrintWriter pWriterPaired 	 = null;	// junction evidence by paired-reads
		
		FileWriter pFileWriterExons 	= null;
		FileWriter pFileWriterJunctions = null;
		FileWriter pFileWriterIntrons	= null;
		FileWriter pFileWriterPaired	= null;
		
		if(bMultiFileOutput)
		{
			System.out.println("Generating separate files for junction and exon read counts.");
			
			pFileWriterExons	 = new FileWriter(outputPrefix+".exon_cnts.tsv");
			pFileWriterJunctions = new FileWriter(outputPrefix+".junction_cnts.tsv");
			
			pWriterExons	 = new PrintWriter(pFileWriterExons);
			pWriterJunctions = new PrintWriter(pFileWriterJunctions);
			
			if(bAdditionalFileOutput)
			{
				pFileWriterIntrons = new FileWriter(outputPrefix+".intron_cnts.tsv");
				pFileWriterPaired  = new FileWriter(outputPrefix+".paired_cnts.tsv");
				
				pWriterIntrons = new PrintWriter(pFileWriterIntrons);
				pWriterPaired  = new PrintWriter(pFileWriterPaired);
			}
		}
		else
		{
			pFileWriterExons	= new FileWriter(outputPrefix+".combined_cnts.tsv");
			pWriterExons	 	= new PrintWriter(pFileWriterExons);
		}
		
		//############################
		//        read genes
		//############################

		GTFParser parser = new GTFParser(gtfFile);

		TreeMap<String, Vector<Gene>> mapGenes = new TreeMap<String, Vector<Gene>>();
		
		GTFGene gtfGene;
		while((gtfGene = parser.nextGene()) != null)
		{
			String strChrom = gtfGene.getChromosome();
			if(mapGenes.containsKey(strChrom))
			{
				Vector<Gene> vcGenes = mapGenes.get(strChrom);
				Gene g = gtfGene.createGene();
				
				//System.out.println(gtfGene.toString());
				
				vcGenes.add(g);
				mapGenes.put(strChrom, vcGenes);
			}
			else
			{
				Vector<Gene> vcGenes = new Vector<Gene>();
				Gene g = gtfGene.createGene();
				
				//System.out.println(g.toShortString());
				
				vcGenes.add(g);
				mapGenes.put(strChrom, vcGenes);
			}
		}
		
		//#############################
		//    process bam/sam files
		//#############################
		SAMFileReader reader = new SAMFileReader(new File(bamFile));
		
		if(bLenient)
			reader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
		
		int lineCounter = 0;
		
		String lastChromosome = null;
		
		long startTime = System.currentTimeMillis();
		
		try
		{
			for(SAMRecord record : reader)
			{
				if(lineCounter % 1000000 == 0 && lineCounter != 0)
				{
					long stopTime = System.currentTimeMillis();
					System.out.println("Process line: " + lineCounter + " took " + (stopTime-startTime) + " milliseconds...");
					startTime = stopTime;
					
					System.out.println("Number of Chromosomes: " + chromosomeToMappingReport.size());
					
					for(String chr : chromosomeToMappingReport.keySet())
					{
						GenomeReadMappingReport report = chromosomeToMappingReport.get(chr);
						
						System.out.println("\t" + chr + ": #reads: " + report.getNumberOfReads() + ", #covered positions: " + report.getNumberOfCoveredPositions() + ", #junctions: " + report.getNumberOfSpliceJunctionsWithCoverage());
					}
				}
				
				// new chromosome started
				if(lastChromosome != null && !record.getReferenceName().equals(lastChromosome))
				{
					System.out.println("\tFinished chromosome: " + lastChromosome + ": #reads: " + chromosomeToMappingReport.get(lastChromosome).getNumberOfReads() + ", #covered positions: " + chromosomeToMappingReport.get(lastChromosome).getNumberOfCoveredPositions() + ", #junctions: " + chromosomeToMappingReport.get(lastChromosome).getNumberOfSpliceJunctionsWithCoverage());
					
					processChromosome(mapGenes, chromosomeToMappingReport.get(lastChromosome), lastChromosome, 
							pWriterExons, pWriterJunctions, pWriterIntrons, pWriterPaired,
							bMultiFileOutput, bAdditionalFileOutput);
					
					chromosomeToMappingReport.remove(lastChromosome);
				}
				
				lineCounter++;
				
				if(!chromosomeToMappingReport.containsKey(record.getReferenceName()))
				{
					System.out.println("current chromosome: " + record.getReadName());
					chromosomeToMappingReport.put(record.getReferenceName(), new GenomeReadMappingReport(record.getReferenceName(), -1, -1));
				}
				
				lastChromosome = record.getReferenceName();
				
				if(bAdditionalFileOutput)
					processSAMString(record, chromosomeToMappingReport.get(record.getReferenceName()), mapGenes);
				else
					processSAMString(record, chromosomeToMappingReport.get(record.getReferenceName()));
			}
		}
		catch(Exception ex)
		{
			ex.printStackTrace();
		}
		reader.close();
		
		if(lastChromosome != null)
		{
			System.out.println("\tFinished chromosome: " + lastChromosome + ": #reads: " + chromosomeToMappingReport.get(lastChromosome).getNumberOfReads() + ", #covered positions: " + chromosomeToMappingReport.get(lastChromosome).getNumberOfCoveredPositions() + ", #junctions: " + chromosomeToMappingReport.get(lastChromosome).getNumberOfSpliceJunctionsWithCoverage());
			
			processChromosome(mapGenes, chromosomeToMappingReport.get(lastChromosome), lastChromosome, 
					pWriterExons, pWriterJunctions, pWriterIntrons, pWriterPaired,
					bMultiFileOutput, bAdditionalFileOutput);
			
			chromosomeToMappingReport.remove(lastChromosome);
		}
		
		pWriterExons.flush();
		pWriterExons.close();
		
		if(bMultiFileOutput)
		{
			pWriterJunctions.flush();
			pWriterJunctions.close();

			if(bAdditionalFileOutput)
			{
				pWriterIntrons.flush();
				pWriterIntrons.close();
				
				pWriterPaired.flush();
				pWriterPaired.close();
			}
		}
	}
	
	private static void processChromosome(TreeMap<String, Vector<Gene>> mapGenes, GenomeReadMappingReport report, String chromosome, 
			PrintWriter pWriterExons, PrintWriter pWriterJunctions, PrintWriter pWriterIntrons, PrintWriter pWriterPaired,
			boolean bMultiFileOutput, boolean bAdditionalFileOutput)
	{
		//############################
		//        count reads
		//############################
		HashMap<String, TreeSet<String>> knownSpliceJunctionsPerChromosome = new HashMap<String, TreeSet<String>>();
		
		if(!mapGenes.containsKey(chromosome))
			return;
		
		for(Gene gene : mapGenes.get(chromosome))
		{			
			System.out.println("Gene: " + gene.getGeneID() + ", " + gene.getChromosome());
				
			if(!knownSpliceJunctionsPerChromosome.containsKey(gene.getChromosome()))
				knownSpliceJunctionsPerChromosome.put(gene.getChromosome(), new TreeSet<String>());
				
			if(report == null)
			{
				System.out.println("Missing report for chromosome: " + gene.getChromosome() + " ignore gene " + gene.getGeneID());
				continue;
			}
			
			int genePositionsCoveredCounter = 0;
			int totalExonicPositionsCounter = 0;
			int spliceJunctionsDetectedCounter = 0;
				
			//###################################
			//     generate exon read counts
			//###################################
			Iterator<Exon> exonIterator = gene.getExonIterator();
			
			while(exonIterator.hasNext())
			{
				Exon e = exonIterator.next();
				
				int sum = 0;
				int positionCounter = 0;
			
				for(int i=e.getGenomicStart(); i<=e.getGenomicStop(); i++)
				{
					sum += report.getReadCountAtPosition(i);
					positionCounter++;
					
					totalExonicPositionsCounter++;
					
					if(report.getReadCountAtPosition(i) > 0)
						genePositionsCoveredCounter++;
				}
				
				//System.out.println("Exon e " + e.getGenomicStart() + " " + genePositionsCoveredCounter + " " + totalExonicPositionsCounter);
				
				pWriterExons.println("Exon\t" + gene.getGeneID() + "\t" + gene.isPlusStrand() + "\t"+(e.getGenomicStart() + "-" + e.getGenomicStop()) + "\t"+  gene.getChromosome() + "\t" + e.getGenomicStart() + "\t" + e.getGenomicStop() + "\t" + positionCounter + "\t" + sum + "\t" + ((double)sum/positionCounter));
			}
				
			//###################################
			//     generate intron read counts
			//###################################
			if(bAdditionalFileOutput)
			{
				Iterator<Intron> intronIterator = gene.getIntronIterator();
				while(intronIterator.hasNext())
				{
					Intron intron = intronIterator.next();
					
					int sum = 0;
					int positionCounter = 0;
				
					for(int i=intron.getGenomicStart(); i<=intron.getGenomicStop(); i++)
					{
						sum += report.getReadCountAtPosition(i);
						positionCounter++;
					}
					
					pWriterIntrons.println("Intron\t" + gene.getGeneID() + "\t" + gene.isPlusStrand() + "\t"+(intron.getGenomicStart() + "-" + intron.getGenomicStop()) + "\t"+  gene.getChromosome() + "\t" + intron.getGenomicStart() + "\t" + intron.getGenomicStop() + "\t" + positionCounter + "\t" + sum + "\t" + ((double)sum/positionCounter));
				}
			}
			
			//###########################################
			//    generate known junction read counts
			//###########################################
			TreeSet<String> junctions = gene.getSpliceJunctionInformation();
			knownSpliceJunctionsPerChromosome.get(gene.getChromosome()).addAll(junctions);
			
			for(String junction : junctions)
			{
				int count = report.getReadCountForSpliceJunction(junction);
				
				if(count > 0)
					spliceJunctionsDetectedCounter++;
				
				if(bMultiFileOutput)
				{
					pWriterJunctions.println("Junction\t" + gene.getGeneID() + "\t" + gene.isPlusStrand() + "\t"+ (junction) + "\t"+  gene.getChromosome() + "\t" + junction.split("-")[0] + "\t" + junction.split("-")[1] + "\t" + 0 + "\t" + count + "\t" + count);
				}
				else
				{
					pWriterExons.println("Junction\t" + gene.getGeneID() + "\t" + gene.isPlusStrand() + "\t"+ (junction) + "\t"+  gene.getChromosome() + "\t" + junction.split("-")[0] + "\t" + junction.split("-")[1] + "\t" + 0 + "\t" + count + "\t" + count);
				}
			}
			
			//System.out.println("Gene\t" + gene.getGeneID() + "\t" + gene.isPlusStrand() + "\t" + totalExonicPositionsCounter + "\t" 
			//		+ genePositionsCoveredCounter + "\t" + (genePositionsCoveredCounter/(double)totalExonicPositionsCounter) + "\t" 
			//		+ junctions.size() + "\t" + spliceJunctionsDetectedCounter + "\t" + (spliceJunctionsDetectedCounter/(double)junctions.size()));
			
			pWriterExons.println("Gene\t" + gene.getGeneID() + "\t" + gene.isPlusStrand() + "\t" + totalExonicPositionsCounter + "\t" 
					+ genePositionsCoveredCounter + "\t" + (genePositionsCoveredCounter/(double)totalExonicPositionsCounter) + "\t" 
					+ junctions.size() + "\t" + spliceJunctionsDetectedCounter + "\t" + (spliceJunctionsDetectedCounter/(double)junctions.size()));
		}
		
		//##########################################l
		//    generate novel junction read counts
		//##########################################l
		TreeSet<String> knownJunctions = knownSpliceJunctionsPerChromosome.get(chromosome);
		
		if(report == null)
			return;
		
		if(knownJunctions == null)
			return;
		
		for(String junction : report.getSpliceJunctionNames())
		{
			if(!knownJunctions.contains(junction))
			{
				if(bMultiFileOutput)
				{
					pWriterJunctions.println("Novel_Junction\tunknown\tunknown\t"+ (junction) + "\t"+  chromosome + "\t" + junction.split("-")[0] + "\t" + junction.split("-")[1] + "\t" + 0 + "\t" + report.getReadCountForSpliceJunction(junction)+ "\t" + report.getReadCountForSpliceJunction(junction));
				}
				else
				{
					pWriterExons.println("Novel_Junction\tunknown\tunknown\t"+ (junction) + "\t"+  chromosome + "\t" + junction.split("-")[0] + "\t" + junction.split("-")[1] + "\t" + 0 + "\t" + report.getReadCountForSpliceJunction(junction)+ "\t" + report.getReadCountForSpliceJunction(junction));
				}
			}
		}
		
		//#######################################
		//    generate exon connection counts
		//#######################################
		HashMap<String, Integer> mapConnections = report.getAllReadCountsForConnectedExons();
			
		for(Entry<String, Integer> e : mapConnections.entrySet())
		{
			String[] split = e.getKey().split("_");
				
			pWriterPaired.println("connectedExons\tunknown\tunknown\t"  + chromosome + "\t" + split[0] + "\t" + split[1] + "\t" + split[2] + "\t" + split[3] + "\t" + e.getValue().toString());
			e.getValue();
		}
	}
	
	private static void processSAMString(SAMRecord samRecord, GenomeReadMappingReport report, TreeMap<String, Vector<Gene>> mapGenes)
	{
		try
		{
			String readID = samRecord.getReadName().split("#|\\s+")[0];
			samRecord.setReadName(readID);
			
			//samRecord.setReferenceName(samRecord.getReferenceName().replaceAll("chr", ""));
			
			report.setReadLength(samRecord.getReadLength());
			
			Cigar cigar = new Cigar(samRecord.getCigarString());
			
			int currentPosition = samRecord.getAlignmentStart();
			
			report.increaseNumberOfReads();
			
			//System.out.println("Current Position: " + currentPosition);
			
			boolean bIsSplitRead = false;
			
			for(int i=0; i<cigar.getNumberOfGroups(); i++)
			{
				// increase counts in corresponding range
				if(cigar.getSymbolOfGroup(i) == Cigar.ANY_MATCH)
				{
					//System.out.println("Increase counts in range: " + currentPosition + " - " + (currentPosition+cigar.getLengthOfGroup(i)-1));
					report.increaseCountsInRange(currentPosition, currentPosition+cigar.getLengthOfGroup(i)-1);
					//report.increaseCountsInRange(currentPosition, currentPosition+cigar.getLengthOfGroup(i));
				}
				// if there are larger skips (corresponding to introns in mRNAseq reads), add a junction to dataset
				else if(cigar.getSymbolOfGroup(i) == Cigar.SKIP_IN_REFERNCE)
				{
					//System.out.println(samString);
					//System.out.println("Add Junction: " + (currentPosition) + " - " + (currentPosition+cigar.getLengthOfGroup(i)));
					report.addSpliceJunction(currentPosition-1, currentPosition+cigar.getLengthOfGroup(i));
					//report.addSpliceJunction(currentPosition, currentPosition+cigar.getLengthOfGroup(i));
					bIsSplitRead = true;
				}
				// skip insertions
				else if(cigar.getSymbolOfGroup(i) == Cigar.DELETION_IN_REFERENCE)
					continue;

				currentPosition += cigar.getLengthOfGroup(i);
				//System.out.println("Updated position: " + currentPosition);
			}
			
			// use 'unsplit' reads to identify connected exons, restrict to first read in pair
			if(!bIsSplitRead && !samRecord.getSecondOfPairFlag())
			{
				Vector<Gene> vcGenes = mapGenes.get(samRecord.getReferenceName());
				
				int nFirst  = Math.min(samRecord.getAlignmentStart(), samRecord.getMateAlignmentStart());
				int nSecond = Math.max(samRecord.getAlignmentStart(), samRecord.getMateAlignmentStart());

				Exon exFirst	= null;
				Exon exSecond	= null;
				for(Gene gene : vcGenes)
				{
					if(gene.getStart() <= nSecond && gene.getStop() >= nFirst)
					{
						Iterator<Exon> it = gene.getExonIterator();
						while(it.hasNext())
						{
							Exon ex = it.next();
							
							if(ex.getCodingStart() <= nFirst && ex.getCodingStop() >= nFirst)
							{
								exFirst = ex;
								
								if(exSecond != null)
									report.addConnectedExons(exFirst.getCodingStart(), exFirst.getCodingStop(), exSecond.getCodingStart(), exSecond.getCodingStop());
							}
							else if(ex.getCodingStart() <= nSecond && ex.getCodingStop() >= nSecond)
							{
								exSecond = ex;
								
								if(exFirst != null)
									report.addConnectedExons(exFirst.getCodingStart(), exFirst.getCodingStop(), exSecond.getCodingStart(), exSecond.getCodingStop());
							}
						}
					}
				}
			}
		}
		catch(Exception ex)
		{
			// do nothing, this line could not be analyzed...
		}
	}
	
	
	// slower, but also counts reads in introns and how often exons are connected by paired reads
	private static void processSAMString(SAMRecord samRecord, GenomeReadMappingReport report)
	{
		try
		{
			String readID = samRecord.getReadName().split("#|\\s+")[0];
			samRecord.setReadName(readID);
			
			//samRecord.setReferenceName(samRecord.getReferenceName().replaceAll("chr", ""));
			
			report.setReadLength(samRecord.getReadLength());
			
			Cigar cigar = new Cigar(samRecord.getCigarString());
			
			int currentPosition = samRecord.getAlignmentStart();
			
			report.increaseNumberOfReads();
			
			//System.out.println("Current Position: " + currentPosition);
			
			for(int i=0; i<cigar.getNumberOfGroups(); i++)
			{
				// increase counts in corresponding range
				if(cigar.getSymbolOfGroup(i) == Cigar.ANY_MATCH)
				{
					//System.out.println("Increase counts in range: " + currentPosition + " - " + (currentPosition+cigar.getLengthOfGroup(i)-1));
					report.increaseCountsInRange(currentPosition, currentPosition+cigar.getLengthOfGroup(i)-1);
					//report.increaseCountsInRange(currentPosition, currentPosition+cigar.getLengthOfGroup(i));
				}
				// if there are larger skips (corresponding to introns in mRNAseq reads), add a junction to dataset
				else if(cigar.getSymbolOfGroup(i) == Cigar.SKIP_IN_REFERNCE)
				{
					//System.out.println(samString);
					//System.out.println("Add Junction: " + (currentPosition) + " - " + (currentPosition+cigar.getLengthOfGroup(i)));
					report.addSpliceJunction(currentPosition-1, currentPosition+cigar.getLengthOfGroup(i));
					//report.addSpliceJunction(currentPosition, currentPosition+cigar.getLengthOfGroup(i));
				}
				// skip insertions
				else if(cigar.getSymbolOfGroup(i) == Cigar.DELETION_IN_REFERENCE)
					continue;

				currentPosition += cigar.getLengthOfGroup(i);
				//System.out.println("Updated position: " + currentPosition);
			}
		}
		catch(Exception ex)
		{
			// do nothing, this line could not be analyzed...
		}
	}

	private static void exonCountsToDEXSeqInput(String strInput, String strOutput)
	{
		System.out.println("Input file: " + strInput);
		System.out.println("Output file: " + strOutput);
		
		Scanner pIn = null;
		PrintWriter pOut = null;
		try
		{
			pIn = new Scanner(new FileReader(strInput));
			pOut = new PrintWriter(new FileWriter(strOutput));
		}
		catch(FileNotFoundException e)
		{
			System.out.println("failed to open input file\n");
			return;
		}
		catch(IOException e)
		{
			System.out.println("failed to open output file\n");
			return;
		}
		
		int nExonID = 0;
		String strLastGene = "?";
		while(pIn.hasNextLine())
		{

			String strLine = pIn.nextLine();
			String[] split = strLine.split("\\s");
			
			if(split.length != 10)
			{
				System.out.println("failed to process line: " + strLine);
				continue;
			}
		
			if(!strLastGene.equals(split[1]))
			{
				nExonID = 0;
				strLastGene = split[1];
			}
			else
			{
				nExonID += 1;
			}
			
			String strID = "?";
			if(nExonID < 10)
				strID = split[1] + ":E00" + nExonID;
			else if(nExonID < 100)
				strID = split[1] + ":E0" + nExonID;
			else
				strID = split[1] + ":E" + nExonID;
			
			pOut.write(strID + "\t" + split[8] + "\n");
		}
		
		pIn.close();
		pOut.flush();
		pOut.close();
	}
}
