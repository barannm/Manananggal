package BioKit;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;
import java.util.zip.GZIPInputStream;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;

public class MMSeqUtils 
{
	public static void main(String[] args) throws Exception
	{
		if(args[0].equals("-hitsComparison"))
		{
			HashMap<String, TreeSet<String>> hitsOne = new HashMap<String, TreeSet<String>>();
			HashMap<String, TreeSet<String>> hitsTwo = new HashMap<String, TreeSet<String>>();
			
			BufferedReader reader = new BufferedReader(new FileReader(args[1]));
			
			String id = null;
			TreeSet<String> hits = null;
			
			int hitCounter = 0;
			
			String line = null;
			
			while((line = reader.readLine()) != null)
			{
				if(line.startsWith(">"))
				{
					if(id != null)
					{
						hitsOne.put(id.split("#|\\s+|/|F3|F5-RNA")[0], hits);
					}
					
					id = line;
					hits = new TreeSet<String>();
				}
				else if(id != null)
				{
					hits.add(line);
					hitCounter++;
				}
			}
			reader.close();
			System.out.println("IDs in file one: " + hitsOne.size() + ", hits: " + hitCounter);

			reader = new BufferedReader(new FileReader(args[2]));
			
			id = null;
			hits = null;
			line = null;
			
			hitCounter = 0;
			
			while((line = reader.readLine()) != null)
			{
				if(line.startsWith(">"))
				{
					if(id != null)
					{
						hitsTwo.put(id.split("#|\\s+|/|F3|F5-RNA")[0], hits);
					}
					
					id = line;
					hits = new TreeSet<String>();
				}
				else if(id != null)
				{
					hits.add(line);
					hitCounter++;
				}
			}
			reader.close();
			System.out.println("IDs in file two: " + hitsTwo.size() + ", hits: " + hitCounter);
			
			TreeSet<String> nonred= new TreeSet<String>();
			nonred.addAll(hitsOne.keySet());
			nonred.addAll(hitsTwo.keySet());
			
			System.out.println("Nonred ids: " + nonred.size());
			
			int identicalCounter = 0;
			int hitsOneIsSubsetCounter = 0;
			int hitsTwoIsSubsetCounter = 0;
			
			int hitsOneSmallerThanTwo = 0;
			int hitsTwoSmallerThanOne = 0;
			
			for(String header : nonred)
			{
				if(hitsOne.containsKey(header) && hitsTwo.containsKey(header))
				{
					TreeSet<String> nonredHits = new TreeSet<String>();
					nonredHits.addAll(hitsOne.get(header));
					nonredHits.addAll(hitsTwo.get(header));
					
					if(nonredHits.size() != hitsOne.get(header).size() || nonredHits.size() != hitsTwo.get(header).size())
					{
						System.out.println("Difference: " + header + ": " + hitsOne.get(header).size() + " vs. " +hitsTwo.get(header).size() + ", nonred: " + nonredHits.size());
					}
					
					if(nonredHits.size() == hitsOne.get(header).size() && nonredHits.size() == hitsTwo.get(header).size())
					{
						identicalCounter++;
					}
					
					if(nonredHits.size() == hitsOne.get(header).size() && hitsOne.get(header).size() > hitsTwo.get(header).size())
					{
						hitsTwoIsSubsetCounter++;
					}
					
					if(nonredHits.size() == hitsTwo.get(header).size() && hitsOne.get(header).size() < hitsTwo.get(header).size())
					{
						hitsOneIsSubsetCounter++;
					}
					
					if(hitsOne.get(header).size() > hitsTwo.get(header).size())
					{
						hitsTwoSmallerThanOne++;
					}
					
					if(hitsOne.get(header).size() < hitsTwo.get(header).size())
					{
						hitsOneSmallerThanTwo++;
					}
				}
				else if(hitsOne.containsKey(header))
				{
					System.out.println("Difference: " + header + ": " + hitsOne.get(header).size() + " vs. 0");
				}
				else if(hitsTwo.containsKey(header))
				{
					System.out.println("Difference: " + header + ": 0 vs. " + hitsTwo.get(header).size());
				}
			}
			
			System.out.println("Identical hit sets: " + identicalCounter + " out of " + nonred.size());
			System.out.println("Hits One subset of Two: " + hitsOneIsSubsetCounter);
			System.out.println("Hits Two subset of One: " + hitsTwoIsSubsetCounter);
			System.out.println("Hits One smaller than Two: " + hitsOneSmallerThanTwo);
			System.out.println("Hits Two smaller than One: " + hitsTwoSmallerThanOne);
			
		}
		else if(args[0].equals("-identifyDifferentialSplicing"))
		{
			String mmseqPredictionOne = null;
			String mmseqPredictionTwo = null;
			String gtfFile = null;
			
			for(int i=1; i<args.length; i++)
			{
				if(args[i].equals("-1"))
				{
					mmseqPredictionOne = args[i+1];
				}
				else if(args[i].equals("-2"))
				{
					mmseqPredictionTwo = args[i+1];
				}
				else if(args[i].equals("-g"))
				{
					gtfFile = args[i+1];
				}
			}
			
			identifyDifferentialSplicing(gtfFile, mmseqPredictionOne, mmseqPredictionTwo);
		}
		else if(args[0].equals("-generateHits"))
		{
			String samBamInputFile = null;
			String gtfFile = null;
			String outputFile = null;
			int readLength = 50;
			int minimalDistanceBetweenPairedEnds = 0;
			int maximalDistanceBetweenPairedEnds = 150;
			
			for(int i=1; i<args.length; i++)
			{
				if(args[i].equals("-s"))
				{
					samBamInputFile = args[i+1];
				}
				else if(args[i].equals("-g"))
				{
					gtfFile = args[i+1];
				}
				else if(args[i].equals("-l"))
				{
					readLength = Integer.parseInt(args[i+1]);
				}
				else if(args[i].equals("-o"))
				{
					outputFile = args[i+1];
				}
				else if(args[i].equals("-minpe"))
				{
					minimalDistanceBetweenPairedEnds = Integer.parseInt(args[i+1]);
				}
				else if(args[i].equals("-maxpe"))
				{
					maximalDistanceBetweenPairedEnds = Integer.parseInt(args[i+1]);
				}
			}
			
			generateHitsFile(gtfFile, samBamInputFile, outputFile, readLength, maximalDistanceBetweenPairedEnds, minimalDistanceBetweenPairedEnds);
		}
	}
	
	public static void generateHitsFile(String gtfFile, String samBamFile, String hitsFile, int readLength, int maximalPairedEndDistance, 
			int minimalPairedEndDistance) throws Exception
	{
		PrintWriter writer = new PrintWriter(new FileWriter(hitsFile));
	
		// writer mmseq header information
		GTFParser gtfParser = new GTFParser(gtfFile);
		
		String geneIsoformsInfo = "";
		String transcriptMetaData = "";
		
		System.out.println("Parse GTF file at " + gtfFile + " to obtain transcript structures");
		
		GTFGene gtfGene = null;
		
		while((gtfGene = gtfParser.nextGene()) != null)
		{
			Gene gene = gtfGene.createGene();
			
			geneIsoformsInfo += "@GeneIsoforms\t" + gene.getGeneID();
			
			for(String transcript : gene.getArrayOfGeneProductNames())
			{
				geneIsoformsInfo += "\t" + transcript;
				transcriptMetaData += "@TranscriptMetaData\t" + transcript + "\t" + gene.getTranscriptCDNALength(transcript) + "\t" + gene.getTranscriptCDNALength(transcript) + "\n";
			}
			geneIsoformsInfo += "\n";
		}
		gtfParser.closeReader();
		
		writer.print(transcriptMetaData);
		writer.print(geneIsoformsInfo);

		// now read BAM / SAM file and map to transcripts
		System.out.println("Process SAM file: " + samBamFile);
		
		// read mapping files (genomic) and generate read counts and RPKM values for all genes
		//BufferedReader reader = new BufferedReader(new FileReader(f));
		SAMFileReader reader = new SAMFileReader(new File(samBamFile));
		reader.setValidationStringency(ValidationStringency.SILENT);
		int lineCounter = 0;
		
		ReadMapping mapping = null;
		
		HashMap<String, TreeSet<String>> readIDToMappedTranscripts = new HashMap<String, TreeSet<String>>();
		HashMap<String, ReadMapping> readIDToMapping = new HashMap<String, ReadMapping>();
		HashMap<String, Integer> transcriptToUniqueReads = new HashMap<String, Integer>();
		HashMap<String, Integer> transcriptToCompatibleReads = new HashMap<String, Integer>();
		
		TreeSet<String> mappedTranscripts = null;

		System.out.println("Initialize mapping compatibility evaluator...");
		ReadMappingEvaluatorLight evaluator = new ReadMappingEvaluatorLight(gtfFile, readLength);
		int numOfMappings = 0;
		int numOfMappingsWithoutCompatibleTranscripts = 0;
		int numOfMappingsWithCompatibleTranscripts = 0;
		int numOfMappingsWithUniqueTranscripts = 0;
		
		for(SAMRecord record : reader)
		{
			if(lineCounter % 1000000 == 0)
			{
				System.out.println("Processed " + lineCounter + " lines... ReadID to mapped transcripts: " + readIDToMappedTranscripts.size() + ", ReadID to mapping: " + readIDToMapping.size());
			}
			
			lineCounter++;
			
			// if SAM format
			mapping = AlignmentFileParser.convertPicardSAMRecord(record);
			String readID = mapping.getReadID().split("#|\\s+|/|F3|F5-RNA")[0];
			
			if(mapping == null || mapping.isUnmapped())
				continue;
			
			numOfMappings++;
			
			//System.out.println(">" + mapping.getReadID());
			mappedTranscripts = evaluator.mappingToCompatibleTranscripts(mapping);
			
			if (mapping.pairedEndRead())
			{
				String hashKeyCurRead = mapping.getReadID().split("#|\\s+|/|F3|F5-RNA")[0] + ":" + mapping.getStartPosition() 
						+ "-" + mapping.getStartPositionOfPairedEnd();
				String hashKeyMateRead = mapping.getReadID().split("#|\\s+|/|F3|F5-RNA")[0] + ":" + mapping.getStartPositionOfPairedEnd() 
						+ "-" + mapping.getStartPosition();
				
				if(!readIDToMappedTranscripts.containsKey(hashKeyMateRead))
				{
					readIDToMappedTranscripts.put(hashKeyCurRead, mappedTranscripts);
					readIDToMapping.put(hashKeyCurRead, mapping);
				}
				else
				{
					TreeSet<String> consensusMapping = new TreeSet<String>();
					
					for(String t : mappedTranscripts)
					{
						if(readIDToMappedTranscripts.get(hashKeyMateRead).contains(t))
						{
							int startInTranscriptReadOne = evaluator.getTranscriptPositionForGenomicPosition(t, readIDToMapping.get(hashKeyMateRead).getStartPosition());
							int startInTranscriptReadTwo = evaluator.getTranscriptPositionForGenomicPosition(t, mapping.getStartPosition());
							int distance = Math.abs(startInTranscriptReadOne-startInTranscriptReadTwo);
							
							if(distance > minimalPairedEndDistance && distance <= maximalPairedEndDistance)
								consensusMapping.add(t);
						}
					}
					
					if(consensusMapping.size() != 0)
					{
						numOfMappingsWithCompatibleTranscripts++;
						writer.println(">"+readID);
						
						for(String c : consensusMapping)
						{
							writer.println(c);
							
							if(transcriptToCompatibleReads.containsKey(c))
							{
								transcriptToCompatibleReads.put(c, transcriptToCompatibleReads.get(c)+1);
							}
							else
							{
								transcriptToCompatibleReads.put(c, 1);
							}
						}
					}
					else numOfMappingsWithoutCompatibleTranscripts++;
					
					if(consensusMapping.size() == 1)
					{
						if(transcriptToUniqueReads.containsKey(consensusMapping.first()))
						{
							transcriptToUniqueReads.put(consensusMapping.first(), transcriptToUniqueReads.get(consensusMapping.first())+1);
						}
						else
						{
							transcriptToUniqueReads.put(consensusMapping.first(), 1);
						}
						numOfMappingsWithUniqueTranscripts++;
					}
					
					readIDToMappedTranscripts.remove(hashKeyMateRead);
					readIDToMapping.remove(hashKeyMateRead);
				}
			}
			else
			{
				if(mappedTranscripts.size() != 0)
				{
					numOfMappingsWithCompatibleTranscripts++;
					writer.println(">"+ readID);
					
					for(String c : mappedTranscripts)
					{
						writer.println(c);
						
						if(transcriptToCompatibleReads.containsKey(c))
						{
							transcriptToCompatibleReads.put(c, transcriptToCompatibleReads.get(c)+1);
						}
						else
						{
							transcriptToCompatibleReads.put(c, 1);
						}
					}
				} else numOfMappingsWithoutCompatibleTranscripts++;
				
				if(mappedTranscripts.size() == 1)
				{
					if(transcriptToUniqueReads.containsKey(mappedTranscripts.first()))
					{
						transcriptToUniqueReads.put(mappedTranscripts.first(), transcriptToUniqueReads.get(mappedTranscripts.first())+1);
					}
					else
					{
						transcriptToUniqueReads.put(mappedTranscripts.first(), 1);
					}
					numOfMappingsWithUniqueTranscripts++;
				}
			}
		}
		
		reader.close();
		writer.flush();
		writer.close();
		
		TreeSet<String> allTranscripts = new TreeSet<String>();
		
		allTranscripts.addAll(transcriptToCompatibleReads.keySet());
		allTranscripts.addAll(transcriptToUniqueReads.keySet());
		
		System.out.println("Read statistics");
		
		for(String t : allTranscripts)
		{
			int unique = 0;
			int compatible = 0;
			
			if(transcriptToCompatibleReads.containsKey(t))
				compatible = transcriptToCompatibleReads.get(t);
			
			if(transcriptToUniqueReads.containsKey(t))
				unique = transcriptToUniqueReads.get(t);
			
			System.out.println(t + "\t" + compatible + "\t" + unique);
		}
		
		System.out.println("Mappings: " + numOfMappings);
		System.out.println("Mappings withOUT compatible transcripts: " + numOfMappingsWithoutCompatibleTranscripts);
		System.out.println("Mappings with compatible transcripts: " + numOfMappingsWithCompatibleTranscripts);
		System.out.println("Mappings with unique transcripts: " + numOfMappingsWithUniqueTranscripts);
	}
	
	public static void generateHitsFileUsingGTFMapperSorted(String sortedGTFFile, String gtfFile, String samBamFile, String hitsFile) throws Exception
	{
		GTFMapperSorted gtfMapperSorted = new GTFMapperSorted(sortedGTFFile);
		MappedElementProcessingSAMBAMMMSeq mappedElementStoring = new MappedElementProcessingSAMBAMMMSeq(hitsFile, gtfFile);
		SortedSAMBAMReadProcessor readProcessor = new SortedSAMBAMReadProcessor(samBamFile, gtfMapperSorted, 
				mappedElementStoring, true, false, false, false, false);
		readProcessor.run();
		mappedElementStoring.finishProcessing();
	}
	
	public static void generateHitsFileForExonAndJunctionData(String gtfFile, String exonFile, String junctionFile, 
		String hitsFile, int readLength) throws Exception
	{
		System.out.println("Initialize mapping compatibility evaluator...");
		ReadMappingEvaluatorLight evaluator = new ReadMappingEvaluatorLight(gtfFile, readLength);
		
		PrintWriter writer = new PrintWriter(new FileWriter(hitsFile));
	
		// writer mmseq header information
		HashMap<String, TreeSet<String>> geneToTranscripts = new HashMap<String, TreeSet<String>>();
		HashMap<String, Integer> transcriptToLength = new HashMap<String, Integer>();
		
		GTFParser gtfParser = new GTFParser(gtfFile);
		
		System.out.println("Parse GTF file at " + gtfFile + " to obtain transcript structures");
		
		GTFGene gtfGene = null;
		
		while((gtfGene = gtfParser.nextGene()) != null)
		{
			Gene gene = gtfGene.createGene();
			
			geneToTranscripts.put(gtfGene.getId(), new TreeSet<String>());
			
			for(String transcript : gene.getArrayOfGeneProductNames())
			{
				geneToTranscripts.get(gene.getGeneID()).add(transcript);
				
				transcriptToLength.put(transcript, gene.getTranscriptCDNALength(transcript));
			}
		}
		
		for(String transcript : transcriptToLength.keySet())
		{
			writer.println("@TranscriptMetaData\t" + transcript + "\t" + transcriptToLength.get(transcript));
		}
		
		for(String gene : geneToTranscripts.keySet())
		{
			writer.print("@GeneIsoforms\t");
			writer.print(gene);
			
			for(String transcript : geneToTranscripts.get(gene))
				writer.print("\t" + transcript);
			
			writer.println();
		}

		// now read BAM / SAM file and map to transcripts
		//System.out.println("Process Exon file: "  + exonFile);
		
		BufferedReader reader = new BufferedReader(new FileReader(exonFile));
		String line = null;
		int lineCounter = 0;
		
		HashMap<String, Integer> transcriptToUniqueReads = new HashMap<String, Integer>();
		HashMap<String, Integer> transcriptToCompatibleReads = new HashMap<String, Integer>();
		
		TreeSet<String> mappedTranscripts = null;
		
		String[] split = null;
		String reference = null;
		int start = -1;
		int stop = -1;
		int numberOfReads = -1;
		
		while((line = reader.readLine()) != null)
		{
			if(lineCounter % 10 == 0)
			{
				System.out.println("Processed " + lineCounter + " lines...");
			}
			
			lineCounter++;
			
			System.out.println("Exon Line: " + line);
			
			split = line.split("\t");
			reference = split[0];
			start = Integer.parseInt(split[1]);
			stop = Integer.parseInt(split[2]);
			// how to normalize for very strange values in some of the files here???
			numberOfReads = (int)(Integer.parseInt(split[4]));
			
			if(numberOfReads > 0)
			{
				mappedTranscripts = evaluator.exonToCompatibleTranscripts(reference.replaceAll("chr", ""), start, stop);
				
				String lineID = line.replaceAll("\t", "-");
				
				for(int i=0; i<numberOfReads; i++)
				{
					writer.println(">"+lineID + "_" + i);
					
					for(String tr : mappedTranscripts)
					{
						writer.println(tr);
						
						if(transcriptToCompatibleReads.containsKey(tr))
						{
							transcriptToCompatibleReads.put(tr, transcriptToCompatibleReads.get(tr)+1);
						}
						else
						{
							transcriptToCompatibleReads.put(tr, 1);
						}
					}
					
					if(mappedTranscripts.size() == 1)
					{
						if(transcriptToUniqueReads.containsKey(mappedTranscripts.first()))
						{
							transcriptToUniqueReads.put(mappedTranscripts.first(), transcriptToUniqueReads.get(mappedTranscripts.first())+1);
						}
						else
						{
							transcriptToUniqueReads.put(mappedTranscripts.first(), 1);
						}
					}
				}
			}
		}
		reader.close();
		
		//System.out.println("Process Junctions file: "  + junctionFile);
		reader = new BufferedReader(new FileReader(junctionFile));
		line = null;
		lineCounter = 0;
		
		mappedTranscripts = null;
		
		split = null;
		reference = null;
		start = -1;
		stop = -1;
		numberOfReads = -1;
		
		while((line = reader.readLine()) != null)
		{
			if(lineCounter % 10 == 0)
			{
				System.out.println("Processed " + lineCounter + " lines...");
			}
			
			//System.out.println("Junction Line: " + line);
			
			lineCounter++;
			
			split = line.split("\t");
			reference = split[0];
			start = Integer.parseInt(split[1]);
			stop = Integer.parseInt(split[2]);
			numberOfReads = Integer.parseInt(split[4]);
			
			if(numberOfReads > 0)
			{
				mappedTranscripts = evaluator.junctionToCompatibleTranscripts(reference.replaceAll("chr", ""), start, stop);
				
				String lineID = line.replaceAll("\t", "-");
				
				for(int i=0; i<numberOfReads; i++)
				{
					writer.println(">" + lineID + "_" + i);
					
					for(String tr : mappedTranscripts)
					{
						writer.println(tr);
						
						if(transcriptToCompatibleReads.containsKey(tr))
						{
							transcriptToCompatibleReads.put(tr, transcriptToCompatibleReads.get(tr)+1);
						}
						else
						{
							transcriptToCompatibleReads.put(tr, 1);
						}
					}
					
					if(mappedTranscripts.size() == 1)
					{
						if(transcriptToUniqueReads.containsKey(mappedTranscripts.first()))
						{
							transcriptToUniqueReads.put(mappedTranscripts.first(), transcriptToUniqueReads.get(mappedTranscripts.first())+1);
						}
						else
						{
							transcriptToUniqueReads.put(mappedTranscripts.first(), 1);
						}
					}
				}
			}
		}
		
		writer.flush();
		writer.close();
		reader.close();
		
		TreeSet<String> allTranscripts = new TreeSet<String>();
		
		allTranscripts.addAll(transcriptToCompatibleReads.keySet());
		allTranscripts.addAll(transcriptToUniqueReads.keySet());
		
		System.out.println("Read statistics");
		
		for(String t : allTranscripts)
		{
			int unique = 0;
			int compatible = 0;
			
			if(transcriptToCompatibleReads.containsKey(t))
				compatible = transcriptToCompatibleReads.get(t);
			
			if(transcriptToUniqueReads.containsKey(t))
				unique = transcriptToUniqueReads.get(t);
			
			System.out.println(t + "\t" + compatible + "\t" + unique);
		}
	}
			
	public static void generateHitsFileUsingEJCCM(String gtfFile, String ejccmFile, String hitsFile, 
			boolean transcriptsMustContainExonMatchingExonicRegionExactly) throws Exception
	{
		PrintWriter writer = new PrintWriter(new FileWriter(hitsFile));
		
		// writer mmseq header information
		GTFParser gtfParser = new GTFParser(gtfFile);
		
		String geneIsoformsInfo = "";
		String transcriptMetaData = "";
		
		System.out.println("Parse GTF file " + gtfFile + " to obtain transcript structures");
		
		GTFGene gtfGene = null;
		
		while((gtfGene = gtfParser.nextGene()) != null)
		{
			Gene gene = gtfGene.createGene();
			
			geneIsoformsInfo += "@GeneIsoforms\t" + gene.getGeneID();
			
			for(String transcript : gene.getArrayOfGeneProductNames())
			{
				geneIsoformsInfo += "\t" + transcript;
				transcriptMetaData += "@TranscriptMetaData\t" + transcript + "\t" + gene.getTranscriptCDNALength(transcript) + "\n";
			}
			geneIsoformsInfo += "\n";
		}
		gtfParser.closeReader();
		
		writer.print(transcriptMetaData);
		writer.print(geneIsoformsInfo);

		System.out.println("\nProcess EJCCM file: " + ejccmFile);
		
		int combinationsMappingSingleGeneCounter = 0;
		// read ejccm files and generate read counts and RPKM values for all genes
		BufferedReader reader = null;
		GZIPInputStream gzip = null;
		if (ejccmFile.endsWith(".gz")) 
		{
			gzip = new GZIPInputStream(new FileInputStream(ejccmFile));
			reader = new BufferedReader(new InputStreamReader(gzip));
		}
		else reader = new BufferedReader(new FileReader(ejccmFile));
	
		ReadMappingEvaluatorLight evaluator = new ReadMappingEvaluatorLight(gtfFile, 1);
		
		String curLine = null;
		int lineCounter = 0;				
		while ((curLine = reader.readLine()) != null)
		{
			if(lineCounter % 100000 == 0)
			{
				System.out.println("Processed " + lineCounter + " lines...");
				System.out.println("\tCombinations mapping unique gene: " + combinationsMappingSingleGeneCounter);
			}			
			lineCounter++;

			// the lines starting with a hash symbol
			// contain general informatin that can be ignore here.
			if (curLine.startsWith("#")) continue;
			
			String[] curLineContent = curLine.split("\\s+");
			String combinationID = curLineContent[SAMBAMCompression.COL_EJCCM_COMB_ID - 1];
			int counts = Integer.parseInt(curLineContent[SAMBAMCompression.COL_EJCCM_COUNTS - 1]);
			
			TreeSet<String> compatibleTranscripts = getTranscriptsCompatibleWithCombination(combinationID, 
					evaluator, transcriptsMustContainExonMatchingExonicRegionExactly);
			
			if (compatibleTranscripts.isEmpty()) continue;
			
			String compatibleTranscriptsString = "";
			for (String transcript : compatibleTranscripts)
				compatibleTranscriptsString += transcript + "\n";
			
			for (int i = 0; i < counts; i++)
			{
				writer.println(">" + combinationID + ":" + i);
				writer.print(compatibleTranscriptsString);
			}
		}
		
		writer.flush();
		writer.close();
		gzip.close();
		
		System.out.println("DONE");
	}
	
	private static TreeSet<String> getTranscriptsCompatibleWithCombination(String combinationID, ReadMappingEvaluatorLight evaluator,
			boolean transcriptsMustContainExonMatchingExonicRegionExactly)
	{
		String[] features = combinationID.split(SAMBAMCompression.ELEMENT_ID_SEPARATOR);
		
		HashMap<String, Integer> transcriptsToNumMappedFeatures = new HashMap<String, Integer>();
		TreeMap<String, Vector<Integer>> transcriptToBorderTranscriptPositionFeatures = new TreeMap<String, Vector<Integer>>();
				
		for (String featureID : features)
		{
			String[] featureIDParts = featureID.split(",");
			// If the featureID does not contain any ',' then
			// the featureID describes an exon.
			// If the feature ID describes a junction then the 
			// junction start position and stop position are
			// separated by a comma s.t. the array has length 2.
			if (featureIDParts.length == 1)
			{
				featureIDParts = featureID.split(":");
				String seqNameExons = featureIDParts[0];
				String[] exons = null;
				/* Depending on the compression method the IDs can have
				 * different formats. In the case of multi-mapping reads
				 * either the exons the read mapped to are listed explicitly:
				 * 	chr1:(startExon1-endExon1/startExon2-endExon2/startExon3-endExon3)
				 * or only the minimum start and the maximum end are specified:
				 * 	chr1:minimumStart-maximumEnd
				 * In the first case mapped transcripts must contain an exon
				 * exactly matching one the regions whereas in the second they
				 * only need an overlapping exon.
				 */
				if (!featureIDParts[1].startsWith("\\d+"))
				{
					exons = featureIDParts[1].replaceAll("\\(", "").replaceAll("\\)", "").split("\\/"); 
				}
				else
				{
					exons = new String[1];
					exons[0] = featureIDParts[1];
				}
				for (String exon : exons)
				{
					String[] exonStartStop = exon.split("-"); 
					int exonStart = Integer.parseInt(exonStartStop[0]);
					int exonStop = Integer.parseInt(exonStartStop[1]);
					
					for (String transcript : evaluator.mappingToCompatibleTranscripts(exonStart, exonStop, seqNameExons,
							transcriptsMustContainExonMatchingExonicRegionExactly))
					{
						int exonStartInTranscript = evaluator.getTranscriptPositionForGenomicPosition(transcript, exonStart);
						int exonStopInTranscript = evaluator.getTranscriptPositionForGenomicPosition(transcript, exonStop);
						if (!transcriptsToNumMappedFeatures.containsKey(transcript)) transcriptsToNumMappedFeatures.put(transcript, 0);
						transcriptsToNumMappedFeatures.put(transcript, transcriptsToNumMappedFeatures.get(transcript) + 1);
						
						if (!transcriptToBorderTranscriptPositionFeatures.containsKey(transcript))
						{
							transcriptToBorderTranscriptPositionFeatures.put(transcript, new Vector<Integer>());
							// The first position of the first feature is not relevant since only the distance 
							// between the features are regarded
							transcriptToBorderTranscriptPositionFeatures.get(transcript).add(exonStopInTranscript);
						}
						else 
						{
							// If a feature is inbetween both borders are needed
							transcriptToBorderTranscriptPositionFeatures.get(transcript).add(exonStartInTranscript);
							transcriptToBorderTranscriptPositionFeatures.get(transcript).add(exonStopInTranscript);
						}
					}
				}
			}
			else if (featureIDParts.length == 2)
			{
				String junctionSeqName = null;
				String junction = null;
				Integer junctionStart = null;
				Integer junctionStop = null;
				for (String junctionInfo : featureIDParts)
				{
					String[] junctionInfoParts = junctionInfo.split(":");
					if (junctionSeqName == null) junctionSeqName = junctionInfoParts[0];
					if (junction == null)
					{
						junction = junctionInfoParts[1];
						junctionStart = Integer.parseInt(junctionInfoParts[1]);
					}
					else
					{
						junction += "-" + junctionInfoParts[1];
						junctionStop = Integer.parseInt(junctionInfoParts[1]);
					}
				}
								
				TreeSet<String> mappedTranscripts = evaluator.mappingToCompatibleTranscripts(junction, junctionSeqName);
				if (mappedTranscripts == null) continue;
				
				for (String transcript : mappedTranscripts)
				{
					int junctionStartInTranscript = evaluator.getTranscriptPositionForGenomicPosition(transcript, junctionStart);
					int junctionStopInTranscript = evaluator.getTranscriptPositionForGenomicPosition(transcript, junctionStop);
					
					if (!transcriptsToNumMappedFeatures.containsKey(transcript)) transcriptsToNumMappedFeatures.put(transcript, 0);
					transcriptsToNumMappedFeatures.put(transcript, transcriptsToNumMappedFeatures.get(transcript) + 1);
										
					if (!transcriptToBorderTranscriptPositionFeatures.containsKey(transcript))
					{
						transcriptToBorderTranscriptPositionFeatures.put(transcript, new Vector<Integer>());
						// The first position of the first feature is not relevant since only the distance 
						// between the features are regarded
						transcriptToBorderTranscriptPositionFeatures.get(transcript).add(junctionStopInTranscript);
					}
					else 
					{
						// If a feature is inbetween both borders are needed
						transcriptToBorderTranscriptPositionFeatures.get(transcript).add(junctionStartInTranscript);
						transcriptToBorderTranscriptPositionFeatures.get(transcript).add(junctionStopInTranscript);
					}
				}
			}	
			else
			{
				System.err.println("Cannot read feature ID : " + featureID);
				System.exit(0);
			}
		}
		
		/* Write transcripts mapped to combination to file.
		 * Originally the hits-file for MMSeq contains for each
		 * read the transcripts that could be mapped to that read.
		 * 
		 * The compressed sam/bam file basically counts how many reads
		 * mapped to each combination. In order to translate that 
		 * information to the hits-file format the combination with
		 * the mapped transcripts is written as often to the hits-file
		 * as the counts specify. The combination ID is then extended
		 * with a number indicating the repeat. 
		 */
		
		TreeSet<String> transcriptsMappedToCombination = new TreeSet<String>();
		for (String transcript : transcriptsToNumMappedFeatures.keySet())
		{
			if (transcriptsToNumMappedFeatures.get(transcript) == features.length)
			{
//				boolean approximateMinimumFragmentLengthAcceptable = true;
//				int approximateMinimumFragmentLength = 0;
//				
//				Vector<Integer> featureBordersCurTranscript = transcriptToBorderTranscriptPositionFeatures.get(transcript);
//				// Ignore last entry
//				featureBordersCurTranscript.remove(featureBordersCurTranscript.size() - 1);
//				
//				for (int i = 0; i < featureBordersCurTranscript.size(); i += 2)
//				{
//					int stopCurFeature = featureBordersCurTranscript.get(i);
//					int startNextFeature = featureBordersCurTranscript.get(i + 1);
//					approximateMinimumFragmentLength += (startNextFeature - stopCurFeature + 1);
//					
//					if (approximateMinimumFragmentLength > meanFragmentLength + readLength)
//					{
//						approximateMinimumFragmentLengthAcceptable = false;
//						break;
//					}
//				}
//				
//				if (!approximateMinimumFragmentLengthAcceptable) continue;
				
				transcriptsMappedToCombination.add(transcript);
			}
		}		
		return transcriptsMappedToCombination;
	}
		
	public static void identifyDifferentialSplicing(String gtfFile, String mmseqPredictionOne, String mmseqPredictionTwo) throws Exception
	{
		DecimalFormat df = new DecimalFormat( "0.00" ); 
		
		HashMap<String, Vector<String>> genesToTranscripts = new HashMap<String, Vector<String>>();
		
		GTFParser gtfParser = new GTFParser(gtfFile);
		
		System.out.println("Parse GTF file at " + gtfFile + " to obtain transcript structures");
		
		GTFGene gtfGene = null;
		
		int counter = 0;
		
		while((gtfGene = gtfParser.nextGene()) != null)
		{
			Gene gene = gtfGene.createGene();
			
			Vector<String> transcripts = new Vector<String>();
			
			for(String transcript : gene.getArrayOfGeneProductNames())
				transcripts.add(transcript);
			
			counter += transcripts.size();
			
			genesToTranscripts.put(gene.getGeneID(), transcripts);
		}
		
		gtfParser.closeReader();
		
		System.out.println("Identified " + genesToTranscripts.size() + " genes and " + counter + " transcripts...");
		
		HashMap<String, Double> expressionValuesPredictionOne = new HashMap<String, Double>();
		HashMap<String, Double> expressionValuesPredictionTwo = new HashMap<String, Double>();
		
		BufferedReader reader = new BufferedReader(new FileReader(mmseqPredictionOne));
		String line = null;
		
		while((line = reader.readLine()) != null)
		{
			try
			{
				String[] split = line.split("\\s+");
			
				expressionValuesPredictionOne.put(split[0], Math.exp(Double.parseDouble(split[2])));
			}
			catch(Exception ex)
			{
				System.err.println("Could not read line: " + line);
			}
		}
		reader.close();
		reader = new BufferedReader(new FileReader(mmseqPredictionTwo));
		line = null;
		
		while((line = reader.readLine()) != null)
		{
			try
			{
				String[] split = line.split("\\s+");
			
				expressionValuesPredictionTwo.put(split[0], Math.exp(Double.parseDouble(split[2])));
			}
			catch(Exception ex)
			{
				System.err.println("Could not read line: " + line);
			}
		}
		reader.close();
		// now analyze the two expression measurements
		for(String gene: genesToTranscripts.keySet())
		{
			double totalExpressionPredictionOne = 0;
			double totalExpressionPredictionTwo = 0;
			
			for(String transcript : genesToTranscripts.get(gene))
			{
				if(expressionValuesPredictionOne.containsKey(transcript))
					totalExpressionPredictionOne += expressionValuesPredictionOne.get(transcript);
				
				if(expressionValuesPredictionTwo.containsKey(transcript))
					totalExpressionPredictionTwo += expressionValuesPredictionTwo.get(transcript);
			}
			
			boolean headerPrinted = false;
			
			for(String transcript : genesToTranscripts.get(gene))
			{
				double expressionPredictionOne = 0;
				double expressionPredictionTwo = 0;
				
				if(expressionValuesPredictionOne.containsKey(transcript))
					expressionPredictionOne = expressionValuesPredictionOne.get(transcript);
				
				if(expressionValuesPredictionTwo.containsKey(transcript))
					expressionPredictionTwo = expressionValuesPredictionTwo.get(transcript);
				
				double relativeExpressionPredictionOne = (expressionPredictionOne/totalExpressionPredictionOne)*100;
				double relativeExpressionPredictionTwo = (expressionPredictionTwo/totalExpressionPredictionTwo)*100;
				
				if(totalExpressionPredictionOne > 5 && totalExpressionPredictionTwo > 5)
				{
					//if(Math.abs(relativeExpressionPredictionOne - relativeExpressionPredictionTwo) > 30)
					//{
						if(!headerPrinted)
						{
							System.out.println(">" + gene);
							headerPrinted = true;
						}
						
						System.out.println(transcript + "\t" + df.format(expressionPredictionOne) + "\t" + df.format(expressionPredictionTwo) + "\t"
							+ df.format(relativeExpressionPredictionOne) + "\t" + df.format(relativeExpressionPredictionTwo));
					//}
						
				}
			}
		}
	}

}
