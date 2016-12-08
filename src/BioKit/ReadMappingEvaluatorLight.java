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
 * @author Fabian Birzele
 */
package BioKit;

import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.Vector;

public class ReadMappingEvaluatorLight 
{
	// prepare junction data
	private HashMap<String, HashMap<String, Vector<Integer>>> chromosomeToExonJunctionPositionsToTranscriptIDs;
	private HashMap<String, HashMap<Integer, Integer>> chromosomeToGenomicPositionsToTranscriptGroupIDs;
	private HashMap<Integer, HashMap<Integer, Integer>> transcriptToGenomicPositionsToTranscriptPosition;
	private HashMap<Integer, int[]> transcriptToGenomicPositionsArray;
	private HashMap<Integer, int[]> transcriptToGenomicStartStop;
	private HashMap<Integer, String> idToTranscript;
	private HashMap<String, Integer> transcriptToId;
	private HashMap<Integer, HashSet<Integer>> idToTranscriptGroup;
	private HashMap<Integer, Vector<Integer>> transcriptToJunctions;
	private HashMap<String, Integer> junctionsToId;
	private HashMap<Integer, String> idToJunction;
	
	private int geneCounter;
	private int transcriptCounter;
	private int groupCounter;
	private int junctionCounter;
	private final int indexStart = 0;
	private final int indexStop = 1;
	private boolean sequenceNamesStartWithChr;
	
	public ReadMappingEvaluatorLight(String gtfFile, int readLength) throws Exception
	{
		this.sequenceNamesStartWithChr = false;
		initializeEvaluator(gtfFile, readLength);
	}
	
	public void initializeEvaluator(GTFGene gtfGene, int readLength) throws Exception
	{
		chromosomeToExonJunctionPositionsToTranscriptIDs = new HashMap<String, HashMap<String,Vector<Integer>>>();
		chromosomeToGenomicPositionsToTranscriptGroupIDs = new HashMap<String, HashMap<Integer,Integer>>();
		transcriptToGenomicPositionsToTranscriptPosition = new HashMap<Integer, HashMap<Integer,Integer>>();
		transcriptToGenomicPositionsArray = new HashMap<Integer, int[]>();
		transcriptToGenomicStartStop = new HashMap<Integer, int[]>();
		idToTranscript = new HashMap<Integer, String>();
		transcriptToId = new HashMap<String, Integer>();
		idToTranscriptGroup = new HashMap<Integer, HashSet<Integer>>();
		transcriptToJunctions = new HashMap<Integer, Vector<Integer>>();
		junctionsToId = new HashMap<String, Integer>();
		idToJunction = new HashMap<Integer, String>();
		
		processGTFGene(gtfGene, readLength);
	}
	
	public void initializeEvaluator(String gtfFile, int readLength) throws Exception
	{
		chromosomeToExonJunctionPositionsToTranscriptIDs = new HashMap<String, HashMap<String,Vector<Integer>>>();
		chromosomeToGenomicPositionsToTranscriptGroupIDs = new HashMap<String, HashMap<Integer,Integer>>();
		transcriptToGenomicPositionsToTranscriptPosition = new HashMap<Integer, HashMap<Integer,Integer>>();
		transcriptToGenomicPositionsArray = new HashMap<Integer, int[]>();
		transcriptToGenomicStartStop = new HashMap<Integer, int[]>();
		idToTranscript = new HashMap<Integer, String>();
		transcriptToId = new HashMap<String, Integer>();
		idToTranscriptGroup = new HashMap<Integer, HashSet<Integer>>();
		transcriptToJunctions = new HashMap<Integer, Vector<Integer>>();
		junctionsToId = new HashMap<String, Integer>();
		idToJunction = new HashMap<Integer, String>();
		
		// parse GTF file and generate HashMaps
		GTFParser gtfParser = new GTFParser(gtfFile);
		
		System.out.println("Parse GTF file at " + gtfFile + " to obtain transcript structures for evaluator");
		
		GTFGene gtfGene = null;
		
		while((gtfGene = gtfParser.nextGene()) != null)
		{
			processGTFGene(gtfGene, readLength);
		}
		
		//DEBUG
//		int groupID = chromosomeToGenomicPositionsToTranscriptGroupIDs.get("chr15").get(40509692);
//		System.out.println("chr5, position 40509692:");
//		for(Integer transcript : idToTranscriptGroup.get(groupID))
//		{
//			System.out.println(idToTranscript.get(transcript));
//		}
		//END DEBUG
	}
	
	public int getTranscriptPositionForGenomicPosition(String transcriptName, int genomicPosition)
	{
		if(!transcriptToId.containsKey(transcriptName))
			return -1;
		
		// lazy: generated only upon first request
		if(!transcriptToGenomicPositionsToTranscriptPosition.containsKey(transcriptToId.get(transcriptName)))
		{
			transcriptToGenomicPositionsToTranscriptPosition.put(transcriptToId.get(transcriptName), new HashMap<Integer, Integer>());
			
			int[] cDNAPositionToGenomicPosition = transcriptToGenomicPositionsArray.get(transcriptToId.get(transcriptName));
			
			for(int i=0; i<cDNAPositionToGenomicPosition.length; i++)
			{
				transcriptToGenomicPositionsToTranscriptPosition.get(transcriptToId.get(transcriptName)).put(cDNAPositionToGenomicPosition[i], i);
			}			
			transcriptToGenomicPositionsArray.remove(transcriptToId.get(transcriptName));
		}
		
		if(!transcriptToGenomicPositionsToTranscriptPosition.get(transcriptToId.get(transcriptName)).containsKey(genomicPosition))
			return -1;
		else
			return transcriptToGenomicPositionsToTranscriptPosition.get(transcriptToId.get(transcriptName)).get(genomicPosition);
	}
	
	private void processGTFGene(GTFGene gtfGene, int readLength)
	{
		Gene gene = gtfGene.createGene();
		
		if (gene.getChromosome().startsWith("chr")) this.sequenceNamesStartWithChr = true;
		
		geneCounter++;
		
		if(geneCounter % 1000 == 0)
			System.out.println("Processed " + geneCounter + " genes...");			
		
		String[] arrayOfGeneProductNames = gene.getArrayOfGeneProductNames();
		
		processPositionInformation(gene, arrayOfGeneProductNames, readLength);
		processJunctionInformation(gene, arrayOfGeneProductNames);
	}

	private void processPositionInformation(Gene gene, String[] arrayOfGeneProductNames, Integer readLength)
	{
		HashMap<Integer, HashSet<String>> positionToTranscripts = new HashMap<Integer, HashSet<String>>();
		
		for(String transcriptName : gene.getArrayOfGeneProductNames())
		{
			try
			{
				if(transcriptToId.containsKey(transcriptName))
				{
					System.out.println("Ignored duplicate transcript: " + transcriptName);
					continue;
				}
				
				if(!chromosomeToExonJunctionPositionsToTranscriptIDs.containsKey(gene.getChromosome()))
				{
					chromosomeToExonJunctionPositionsToTranscriptIDs.put(gene.getChromosome(), new HashMap<String, Vector<Integer>>());
					chromosomeToGenomicPositionsToTranscriptGroupIDs.put(gene.getChromosome(), new HashMap<Integer, Integer>());
				}
				
				int[] cDNAPositionToGenomicPosition = gene.getGenomicPositionsForGeneProduct(transcriptName);
				
				if(!gene.isPlusStrand())
				{
					cDNAPositionToGenomicPosition = reverseArray(cDNAPositionToGenomicPosition);
				}
				transcriptToGenomicPositionsArray.put(transcriptCounter, cDNAPositionToGenomicPosition);
				int[] genomicStartStop = {cDNAPositionToGenomicPosition[0], cDNAPositionToGenomicPosition[cDNAPositionToGenomicPosition.length - 1]};

				transcriptToGenomicStartStop.put(transcriptCounter, genomicStartStop);
			
				idToTranscript.put(transcriptCounter, transcriptName);
				transcriptToId.put(transcriptName, transcriptCounter);
								
				/* This loop needs to go until cDNAPositionToGenomicPosition.length-readLength + 1 and not
				 * cDNAPositionToGenomicPosition.length-readLength. The positionToTranscripts hash is used
				 * to generate the hash chromosomeToGenomicPositionsToTranscriptGroupIDs. The later is used
				 * to find transcripts a read can be mapped to. This is done by providing the first read 
				 * mapping position. In the case of letting the loop run until 
				 * cDNAPositionToGenomicPosition.length-readLength the first read position is not included
				 * for read mapping to end of a transcript. Thus the loop has to run one position further
				 * using cDNAPositionToGenomicPosition.length-readLength + 1. 
				 */
				for(int i=0; i<cDNAPositionToGenomicPosition.length-readLength + 1; i++)
				{
					// read would cross junction => ignore this position as it will be stored in the context of the junction information reads
					if(Math.abs(cDNAPositionToGenomicPosition[i] - cDNAPositionToGenomicPosition[Math.max(0, Math.min(i+readLength-1, cDNAPositionToGenomicPosition.length-1))]) > readLength) 
					{
						continue;
					}
					else
					{
						if(!positionToTranscripts.containsKey(cDNAPositionToGenomicPosition[i]))
							positionToTranscripts.put(cDNAPositionToGenomicPosition[i], new HashSet<String>());
						positionToTranscripts.get(cDNAPositionToGenomicPosition[i]).add(transcriptName);
						//DEBUG
//						if (cDNAPositionToGenomicPosition[i] == 40509692)
//							System.out.println(gene.getGeneID() + " - " + gene.getChromosome() + " 40509692 -> " + transcriptName);
						//END DEBUG
					}						
				}
				transcriptCounter++;
			}
			catch(Exception ex)
			{
				System.err.println("Could not process transcript " + transcriptName);
				ex.printStackTrace();
			}
		}
		
		HashMap<String, Integer> groupToId = new HashMap<String, Integer>();
		

		for(Integer position : positionToTranscripts.keySet())
		{
			if (positionToTranscripts.get(position).size() == 0) continue;
			
			String groups = "";
			HashSet<Integer> transcriptGroup = new HashSet<Integer>();
			
			if (chromosomeToGenomicPositionsToTranscriptGroupIDs.get(gene.getChromosome()).containsKey(position))
			{
				int groupID = chromosomeToGenomicPositionsToTranscriptGroupIDs.get(gene.getChromosome()).get(position);
				for (Integer transcriptID : idToTranscriptGroup.get(groupID))
				{
					transcriptGroup.add(transcriptID);
					groups += idToTranscript.get(transcriptID) + ",";
				}
			}
						
			for(String transcript : positionToTranscripts.get(position))
			{
				groups += transcript + ",";
				transcriptGroup.add(transcriptToId.get(transcript));
			}
										
			if(groupToId.containsKey(groups))
			{
				chromosomeToGenomicPositionsToTranscriptGroupIDs.get(gene.getChromosome()).put(position, groupToId.get(groups));
			}
			else
			{
				groupToId.put(groups, groupCounter);
				idToTranscriptGroup.put(groupCounter, transcriptGroup);
				groupCounter++;
				
				chromosomeToGenomicPositionsToTranscriptGroupIDs.get(gene.getChromosome()).put(position, groupToId.get(groups));
			}
			//DEBUG
//			if (position == 40509692)
//				System.out.println(gene.getGeneID() + " - " + gene.getChromosome() + " 40509692 -> " + groups + " ID: " + groupToId.get(groups));
			//END DEBUG
		}
	}
	
	private void processJunctionInformation(Gene gene, String[] arrayOfGeneProductNames)
	{
		for(String transcriptName : gene.getArrayOfGeneProductNames())
		{
			//DEBUG
			boolean debug = false;
//			if (transcriptName.equals("ENST00000391949_1"))
//			{
//				debug = true;
//			}
			String debugString = "";
			//END DEBUG
			int[] cDNAPositionToExonID = gene.getExonPositionsForGeneProduct(transcriptName);
			int[] cDNAPositionToGenomicPosition = gene.getGenomicPositionsForGeneProduct(transcriptName);
			
			if(!gene.isPlusStrand())
			{
				cDNAPositionToExonID = reverseArray(cDNAPositionToExonID);
				cDNAPositionToGenomicPosition = reverseArray(cDNAPositionToGenomicPosition);
			}
			
			int transcriptID = transcriptToId.get(transcriptName);
			HashMap<Integer, Integer> junctionToStartStop = new HashMap<Integer, Integer>();
			
			if (!transcriptToJunctions.containsKey(transcriptID))
			{
				transcriptToJunctions.put(transcriptID, new Vector<Integer>());
			}
			
			for(int i=0; i<cDNAPositionToExonID.length-1; i++)
			{
				// exon border found
				if(cDNAPositionToExonID[i] != cDNAPositionToExonID[i+1])
				{
					String junction = cDNAPositionToGenomicPosition[i] + "-" + cDNAPositionToGenomicPosition[i+1];
					int junctionStart = cDNAPositionToGenomicPosition[i];
					
//					if(!gene.isPlusStrand())
//					{
//						junction = cDNAPositionToGenomicPosition[i+1] + "-" + cDNAPositionToGenomicPosition[i];
//						junctionStart = cDNAPositionToGenomicPosition[i+1];
//					}
					//DEBUG
					debugString += "Found junction: " + junction + "\n";
					//END DEBUG
					Integer curJunctionID = null;
					if (!junctionsToId.containsKey(junction))
					{
						junctionCounter++;
						curJunctionID = junctionCounter;
						junctionsToId.put(junction, curJunctionID);
						idToJunction.put(curJunctionID, junction);
						junctionToStartStop.put(curJunctionID, junctionStart);
					}
					else
					{
						curJunctionID = junctionsToId.get(junction);
						junctionToStartStop.put(curJunctionID, junctionStart);
					}
					
					if(!chromosomeToExonJunctionPositionsToTranscriptIDs.get(gene.getChromosome()).containsKey(junction))
						chromosomeToExonJunctionPositionsToTranscriptIDs.get(gene.getChromosome()).put(junction, new Vector<Integer>(1,1));
					
					chromosomeToExonJunctionPositionsToTranscriptIDs.get(gene.getChromosome()).get(junction).add(transcriptID);
					
					/* For easier comparison of the junctions of two transcripts
					 * the junctions are stored sorted.
					 */
					int curIndex = 0;
					for (Integer curOtherJunctionID : transcriptToJunctions.get(transcriptID))
					{
						int otherJunctionStart = junctionToStartStop.get(curOtherJunctionID);
						if (otherJunctionStart > cDNAPositionToGenomicPosition[i]) break;
						curIndex++;
					}
					if (curIndex == transcriptToJunctions.get(transcriptID).size())
					{
						transcriptToJunctions.get(transcriptID).add(curJunctionID);
					}
					else
					{
						transcriptToJunctions.get(transcriptID).add(curIndex, curJunctionID);
					}
				}
			}
			//DEBUG
			if (debug) System.out.println(debugString);
			//END DEBUG
		}
	}
	
	public TreeSet<String> mappingToCompatibleTranscripts(ReadMapping mapping)
	{		
		if(mapping.isSplicedRead())
		{
			TreeSet<String> mappedTranscripts = new TreeSet<String>();
			Cigar cigar = mapping.getMappingCIGAR();
			String chromosome = mapping.getReference();
			if (this.sequenceNamesStartWithChr && !chromosome.startsWith("chr")) chromosome = "chr" + chromosome;
			int currentPosition = mapping.getStartPosition();
			
			HashMap<Integer, Integer> transcriptsToMappedJunctionsCount = new HashMap<Integer, Integer>();
			int junctionCount = 0;
			

			//DEBUG
//			System.out.println(mapping.getReadID());
//			System.out.println(mapping.getStartPosition() + "-" + mapping.getStopPosition());
//			System.out.println(mapping.getMappingCIGARString());
//			System.out.println("Found transcripts with junctions: ");
			//END DEBUG
			
			for(int i=0; i<cigar.getNumberOfGroups(); i++)
			{
				// if there are larger skips (corresponding to introns in mRNAseq reads), add a junction to dataset
				if(cigar.getSymbolOfGroup(i) == Cigar.SKIP_IN_REFERNCE)
				{
					String junction = (currentPosition-1) + "-" + (currentPosition+cigar.getLengthOfGroup(i));
					junctionCount++;
					Integer junctionID = junctionsToId.get(junction);
					if (!(junctionID == null))
					{
						if (!chromosomeToExonJunctionPositionsToTranscriptIDs.get(chromosome).containsKey(junction)) continue;
						//System.out.println("\tTranscripts containing junction: "+chromosomeToExonJunctionPositionsToTranscriptIDs.get(mapping.getReference()).get(junction).size());
						for(Integer transcript : chromosomeToExonJunctionPositionsToTranscriptIDs.get(chromosome).get(junction))
						{
							int[] genomicPositionStartAndStop = transcriptToGenomicStartStop.get(transcript);
							/* The read has to be mapped within the transcript. This means that the mapping has to start
							 * after the transcript and the mapping of the read has to end before
							 * the transcript ends. Otherwise the current transcript no candidate transcript for mapping.
							 */
							if (mapping.getStartPosition() < genomicPositionStartAndStop[indexStart] 
									|| mapping.getStartPosition() > genomicPositionStartAndStop[indexStop]
									|| mapping.getStopPosition() > genomicPositionStartAndStop[indexStop])
							{
								continue;
							}
							if (idToTranscript.get(transcript) == null) continue;
							if(!transcriptsToMappedJunctionsCount.containsKey(transcript))
								transcriptsToMappedJunctionsCount.put(transcript, 1);
							else
								transcriptsToMappedJunctionsCount.put(transcript, transcriptsToMappedJunctionsCount.get(transcript)+1);
						}
					}
				}
				
				currentPosition += cigar.getLengthOfGroup(i);
			}
			for(Integer transcript : transcriptsToMappedJunctionsCount.keySet())
			{		
				if(transcriptsToMappedJunctionsCount.get(transcript) == junctionCount)
				{
					mappedTranscripts.add(idToTranscript.get(transcript));
					//DEBUG
//					System.out.println(idToTranscript.get(transcript));
					//END DEBUG
				}	
			}
			return mappedTranscripts;
		}
		else
		{
			return getMappedTranscriptsForContinousRead(mapping);
		}		
	}
	
	private TreeSet<String> getMappedTranscriptsForContinousRead(ReadMapping mapping)
	{		
		TreeSet<String> mappedTranscripts = new TreeSet<String>();

		try
		{
			int mappingStartPosition = mapping.getStartPosition();
			int mappingEndPosition = mapping.getStopPosition();
			String chromosome = mapping.getReference();
			if (this.sequenceNamesStartWithChr && !chromosome.startsWith("chr")) chromosome = "chr" + chromosome;

//			//DEBUG
//			System.out.println(mapping.getReadID());
//			System.out.println(mapping.getStartPosition() + "-" + mapping.getStopPosition());
//			System.out.println(mapping.getMappingCIGARString());
//			System.out.println("Found transcripts for start position " + mapping.getStartPosition() + ":");
			//END DEBUG
			
			for(Integer transcript : idToTranscriptGroup.get(chromosomeToGenomicPositionsToTranscriptGroupIDs.get(chromosome).get(mapping.getStartPosition()))) 
			{
				//DEBUG
//				System.out.println(idToTranscript.get(transcript));
				//END DEBUG
				int[] genomicPositionStartAndStop = transcriptToGenomicStartStop.get(transcript);
				/* The read has to be mapped within the transcript. This means that the mapping has to start
				 * after the transcript begins but before it ends. And the mapping of the read has to end before
				 * the transcript ends.
				 */
				if (mappingStartPosition < genomicPositionStartAndStop[indexStart] 
						|| mappingStartPosition > genomicPositionStartAndStop[indexStop]
						|| mappingEndPosition > genomicPositionStartAndStop[indexStop])
				{
					continue;
				}
				/* If the read is within transcript boundaries there has to be checked if the transcript
				 * has junctions there which overlap the read mapping region.
				 * It can happen that the read was mapped to a transcript where intron retention occurred.
				 * Therefore for all candidate transcripts it has to be checked whether they have an
				 * exon-exon junction in the region where the read is supposed to map to.
				 */
				Vector<Integer> junctions = transcriptToJunctions.get(transcript);
				int junctionStart = 0;
				int junctionStop = 0;
				boolean junctionExists = false;
				for (Integer curJunction : junctions)
				{
					String junctionString = idToJunction.get(curJunction);
					//DEBUG
//					System.out.println("\t" + junctionString);
					//END DEBUG
					String[] splittedJunctionString = junctionString.split("-");
					junctionStart = Integer.parseInt(splittedJunctionString[0]);
					junctionStop = Integer.parseInt(splittedJunctionString[1]);
					/* A transcript is ignored:
					 * - if a junction lies within the mapped region
					 * - if the region start lies within a junction
					 * - if the region stop lies within a junction
					 */
					if ((mappingStartPosition <= junctionStart && mappingEndPosition > junctionStart)
							|| (mappingStartPosition > junctionStart && mappingStartPosition < junctionStop)
							|| (mappingEndPosition > junctionStart && mappingEndPosition < junctionStop))
					{
						junctionExists = true;
						break;
					}
				}	
				if (!junctionExists)
				{
					mappedTranscripts.add(idToTranscript.get(transcript));
					//DEBUG
//					System.out.println("--> ACCEPTED");
					//END DEBUG
				}
			}
		}
		catch(Exception ex)
		{
//			System.err.println("Could not identify transcript of mapping of continuous read:\n" + mapping.toString());
		}
		return mappedTranscripts;
	}
	
	public TreeSet<String> mappingToCompatibleTranscripts(String junction, String seqName)
	{
		TreeSet<String> mappedTranscripts = new TreeSet<String>();
		String chromosome = seqName;
		if (this.sequenceNamesStartWithChr && !chromosome.startsWith("chr")) chromosome = "chr" + chromosome;		
		if (!chromosomeToGenomicPositionsToTranscriptGroupIDs.containsKey(chromosome)) return null;
		if (!chromosomeToExonJunctionPositionsToTranscriptIDs.get(chromosome).containsKey(junction)) return null;
		
		//System.out.println("\tTranscripts containing junction: "+chromosomeToExonJunctionPositionsToTranscriptIDs.get(mapping.getReference()).get(junction).size());
		for(Integer transcriptID : chromosomeToExonJunctionPositionsToTranscriptIDs.get(chromosome).get(junction))
		{
			if (idToTranscript.get(transcriptID) == null) continue;
			mappedTranscripts.add(idToTranscript.get(transcriptID));
		}	
		return mappedTranscripts;
	}
	
	/**This method searches for transcripts in a given reference transcriptome which
	 * overlap with the specified region. If the parameter
	 * transcriptsMustContainExonMatchingSpecifiedRegionExactly is set to true only
	 * transcript where the specified region corresponds to an exon are returned.
	 * If the variable is set so false all transcripts for which there exist exons
	 * containing the specified region are returned. 
	 * 
	 * @param regionStart
	 * @param regionStop
	 * @param seqName
	 * @param transcriptsMustContainExonMatchingSpecifiedRegionExactly
	 * @return
	 */
	public TreeSet<String> mappingToCompatibleTranscripts(int regionStart, int regionStop, String seqName,
			boolean transcriptsMustContainExonMatchingSpecifiedRegionExactly)
	{
		TreeSet<String> mappedTranscripts = new TreeSet<String>();
		String chromosome = seqName;
		if (this.sequenceNamesStartWithChr && !chromosome.startsWith("chr")) chromosome = "chr" + chromosome;
		if (!chromosomeToGenomicPositionsToTranscriptGroupIDs.containsKey(chromosome)) return mappedTranscripts;
		if (!chromosomeToGenomicPositionsToTranscriptGroupIDs.get(chromosome).containsKey(regionStart))	return mappedTranscripts;
		/* Check all transcripts that overlap with the exonStart whether they contain that exon.
		 * Therefore for each transcript the genomic positions it covers are fetched and its
		 * junctions. Using the junctions the exon of the transcript are determined and compared
		 * to the specified boundaries. Depending on the boolean parameter 
		 * transcriptsMustContainExonMatchingSpecifiedRegionExactly either only transcripts that 
		 * exactly contain the specified exon are taken or transcripts which contain an exon lying within
		 * the specified boundaries.
		 */
		
		//DEBUG
		String debugString = "Region: " + seqName + ":" + regionStart + "-" + regionStop + "\n"; 
		debugString += "Found transcripts: \n";
		boolean debug = false;
		//END DEBUG
		for(Integer transcript : idToTranscriptGroup.get(chromosomeToGenomicPositionsToTranscriptGroupIDs.get(chromosome).get(regionStart))) 
		{
			//DEBUG
			debugString += idToTranscript.get(transcript) + "\n";
			debugString += "Known junctions: \n"; 
			//END DEBUG
			int[] genomicPositionStartAndStop = transcriptToGenomicStartStop.get(transcript);			
			Vector<Integer> junctions = transcriptToJunctions.get(transcript);
			Integer curExonStart = genomicPositionStartAndStop[0];
			Integer curExonStop = null;
			boolean transcriptContainsExon = false;
			
			for (Integer curJunction : junctions)
			{
				String junctionString = idToJunction.get(curJunction);
				// the start of a junction marks the end of an exon
				// and analogously the end of a junction marks the start
				// of an exon.
				String[] junctionStartStop = junctionString.split("-");
				curExonStop = Integer.parseInt(junctionStartStop[0]);
				//DEBUG
				debugString += "\tjunction: " + idToJunction.get(curJunction) + "\n";
				debugString += "\texon: " + curExonStart + "-" + curExonStop + "\n";
				debugString += "\texon start" + curExonStart + "\n";
				debugString += "\texon stop" + curExonStop + "\n";
				if (transcriptsMustContainExonMatchingSpecifiedRegionExactly
						&& regionStart == curExonStart && regionStop == curExonStop)
				{
					debugString += "CONTAINS REGION\n";
				}
				else if (transcriptsMustContainExonMatchingSpecifiedRegionExactly
						&& (regionStart != curExonStart || regionStop != curExonStop))
				{
					debugString += "DOES NOT CONTAIN REGION!\n";
				}
				//END DEBUG
				
				if (transcriptsMustContainExonMatchingSpecifiedRegionExactly && 
						regionStart == curExonStart && regionStop == curExonStop)
				{
					transcriptContainsExon = true;
					break;
				}
				else if (!transcriptsMustContainExonMatchingSpecifiedRegionExactly && 
						curExonStart <= regionStart && curExonStop >= regionStop)
				{
					transcriptContainsExon = true;
					//DEBUG
					debugString += "transcript contains region\n";
					//END DEBUG
					break;
				}	
				curExonStart = Integer.parseInt(junctionStartStop[1]);
			}	
			
			
			//Check the last exon of the transcript
			if (transcriptsMustContainExonMatchingSpecifiedRegionExactly && 
					regionStart == curExonStart 
					&& regionStop == genomicPositionStartAndStop[genomicPositionStartAndStop.length - 1])
			{
				transcriptContainsExon = true;
			}
			else if (!transcriptsMustContainExonMatchingSpecifiedRegionExactly && 
					curExonStart <= regionStart 
					&& genomicPositionStartAndStop[genomicPositionStartAndStop.length - 1] >= regionStop)
			{
				transcriptContainsExon = true;
				//DEBUG
				debugString += "Last transcript position: " + genomicPositionStartAndStop[genomicPositionStartAndStop.length - 1] + "\n";
				debugString += "transcript contains exon\n";
				//END DEBUG
			}					
		
			if (transcriptContainsExon) mappedTranscripts.add(idToTranscript.get(transcript));
			
		}
		
		//DEBUG
		if (debug) System.out.println(debugString);
		//END DEBUG
		
		return mappedTranscripts;
	}
	
	public TreeSet<String> exonToCompatibleTranscripts(String chromosome, int startPosition, int stopPosition)
	{
		//System.out.println("Identify compatible transcripts for exon: " + chromosome + ":" + startPosition + "-" + stopPosition );
		
		TreeSet<Integer> mappedTranscripts = new TreeSet<Integer>();
		
		try
		{
			for(int i=startPosition; i<=stopPosition; i++)
			{
				for(Integer transcript : idToTranscriptGroup.get(chromosomeToGenomicPositionsToTranscriptGroupIDs.get(chromosome).get(i)))
				{
					mappedTranscripts.add(transcript);
				}
			}
		}
		catch(Exception ex)
		{}
		
		TreeSet<String> results = new TreeSet<String>();
		
		for(Integer v : mappedTranscripts)
		{
			if(idToTranscript.containsKey(v))
				results.add(idToTranscript.get(v));
		}
		
		//System.out.println("Mapped " + mappedTranscripts.size() + " transcripts containing exon...");
		
		return results;
	}
	
	public TreeSet<String> junctionToCompatibleTranscripts(String chromosome, int start, int stop)
	{
		//System.out.println("Identify compatible transcripts for junction: " + chromosome + ":" + start + "-" + stop);
		
		TreeSet<String> mappedTranscripts = new TreeSet<String>();
		
		String junction = start + "-" + stop;
				
		try
		{
			//System.out.println("\tTranscripts containing junction: "+chromosomeToExonJunctionPositionsToTranscriptIDs.get(chromosome).get(junction).size());
			
			for(Integer transcript : chromosomeToExonJunctionPositionsToTranscriptIDs.get(chromosome).get(junction))
			{
				if(!idToTranscript.containsKey(transcript))
					continue; 
				
				mappedTranscripts.add(idToTranscript.get(transcript));
			}
		}
		catch(Exception ex)
		{
			//System.err.println("Could not identify junction: " + junction + " in chromosome " + chromosome +  " appearing in mapping: " + start + "-" + stop);
		}
		
		//System.out.println("Mapped " + mappedTranscripts.size() + " transcripts containing junction...");
		
		return mappedTranscripts;
	}

	private static int[] reverseArray(int[] toBeReversed)
	{
		int[] reversed = new int[toBeReversed.length];
		
		for(int i=0; i<toBeReversed.length; i++)
		{
			reversed[toBeReversed.length-1-i] = toBeReversed[i];
		}
		
		return reversed;
	}
	
	//DEBUG
	public int getIdForTranscript(String transcript)
	{
		return transcriptToId.get(transcript);
	}
	
	public String getTranscriptForID(int id)
	{
		return idToTranscript.get(id);
	}
	//END DEBUG
}
