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

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;
import java.util.concurrent.TimeUnit;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;
//DEBUG
import java.lang.management.*;
//END DEBUG

public class SortedSAMBAMReadProcessor 
{
	private String samBamSortedFile;
	private boolean pairedEndData;
	private boolean processMatesInAscendingStartPositionOrder;
	private boolean ignoreReadMappingsToDifferentChromosomes;
	private boolean ignoreReadMappingsToDifferentGenes;
	private boolean onlyUsePairedEndReadsWithBothMatesMapped;
	private GTFMapperSorted gtfMapperSorted;
	private MappedElementProcessingSAMBAM mappedElementProcessingSAMBAM;
	private boolean setCorrespondingExonsForJunctions;
	private boolean setCorrespondingTranscriptsForJunctions;
	
	//DEBUG
	private int maxNumOfMatesWaiting;
	private int maxNumOfMappedElementsStoredForWaitingMates;
	private int maxNumOfGenesStored;
	//END DEBUG
		
	public SortedSAMBAMReadProcessor(String samBamSortedFile, GTFMapperSorted gtfMapperSorted, 
			MappedElementProcessingSAMBAM mappedElementProcessingSAMBAM, boolean processMatesInAscendingStartPositionOrder,
			boolean ignoreReadMappingsToDifferentChromosomes, boolean ignoreReadMappingsToDifferentGenes,
			boolean onlyUsePairedEndReadsWithBothMatesMapped, boolean generateStrandSpecificCombIDs)
	{
		this.samBamSortedFile = samBamSortedFile;
		this.pairedEndData = false;
		this.processMatesInAscendingStartPositionOrder = processMatesInAscendingStartPositionOrder;
		this.ignoreReadMappingsToDifferentChromosomes = ignoreReadMappingsToDifferentChromosomes;
		this.ignoreReadMappingsToDifferentGenes = ignoreReadMappingsToDifferentGenes;
		this.onlyUsePairedEndReadsWithBothMatesMapped = onlyUsePairedEndReadsWithBothMatesMapped;
		this.mappedElementProcessingSAMBAM = mappedElementProcessingSAMBAM;
		this.gtfMapperSorted = gtfMapperSorted;
		this.setCorrespondingExonsForJunctions = true;
		this.setCorrespondingTranscriptsForJunctions = true;
		//DEBUG
		this.maxNumOfMatesWaiting = 0;
		this.maxNumOfMappedElementsStoredForWaitingMates = 0;
		this.maxNumOfGenesStored = 0;
		//END DEBUG
	}
	
	public void run() throws IOException
	{
		System.out.println("Read and process BAM file...");
		SAMFileReader reader = new SAMFileReader(new File(this.samBamSortedFile));
		reader.setValidationStringency(ValidationStringency.SILENT);
		int lineCounter = 0;
		
		Iterator<SAMRecord> samRecordIterator = reader.iterator();
//		SAMRecord record = null;
		
		HashMap<Integer, HashMap<String, HashSet<MappableElement>>> pairedEndReadsWithoutMate = 
				new HashMap<Integer, HashMap<String, HashSet<MappableElement>>>();
		HashMap<Integer, String> pairedEndReadsWithoutMateToReference = 
				new HashMap<Integer, String>();
		
//		//DEBUG
		HashMap<Integer, SAMRecord> pairedEndReadsWithoutMateToMapping = new HashMap<Integer, SAMRecord>();
		double timeStampMillis = System.currentTimeMillis();
		Calendar cal = Calendar.getInstance();
		Date timeStampCal = cal.getTime();
//		//END DEBUG
//		
		int startPositionLastRecord = -1;
		String referenceLastRecord = null;
				
		while(samRecordIterator.hasNext())
		{

			if(lineCounter % 100000 == 0)
			{
				SimpleDateFormat sdf = new SimpleDateFormat("HH:mm:ss");
				System.out.println("TIME\tLast time" + sdf.format(timeStampCal) );
				System.out.println("TIME\tCurrent time" + sdf.format(cal.getTime()) );
				
				double currentTime = System.currentTimeMillis();
				long minutesElapsed = TimeUnit.MILLISECONDS.toMinutes((long) (currentTime - timeStampMillis));
				System.out.println("TIME\tElapsed time processing last read: " + minutesElapsed + " minutes");
			}
			if(lineCounter % 1000000 == 0)
			{				
				System.out.println("Processed " + lineCounter + " lines... ");
				if (this.pairedEndData) System.out.println("Reads waiting for second mate:\t" + pairedEndReadsWithoutMate.keySet().size());
				System.out.println("Number of stored elements:\t" + mappedElementProcessingSAMBAM.getNumOfStoredElements());
				
				System.out.println("MEMORY\tCurrent Heap Memory Usage:\t" + (ManagementFactory.getMemoryMXBean().getHeapMemoryUsage().getUsed()/1024.0/1024.0) + " MB");
				System.out.println("MEMORY\tNon heap Memory Usage:\t" + (ManagementFactory.getMemoryMXBean().getNonHeapMemoryUsage().getUsed()/1024.0/1024.0) + " MB");
				for (MemoryPoolMXBean item : ManagementFactory.getMemoryPoolMXBeans())
				{
					System.out.println("MEMORY\tName:\t" + item.getName());
					System.out.println("MEMORY\tEstimated Current Usage:\t" + (item.getUsage().getUsed()/1024.0/1024.0) + " MB");
					System.out.println("MEMORY\tPeak Usage:\t" + (item.getPeakUsage().getUsed()/1024.0/1024.0) + " MB");
					System.out.println("MEMORY");
				}
				System.out.println("MEMORY\tMax genes stored: " + this.maxNumOfGenesStored);
				System.out.println("MEMORY\tMax mates waiting: " + this.maxNumOfMatesWaiting);
				System.out.println("MEMORY\tMax mapped elements of mates waiting: " + this.maxNumOfMappedElementsStoredForWaitingMates);
			}
			
			SAMRecord record = samRecordIterator.next();
			timeStampMillis = System.currentTimeMillis();
			timeStampCal = cal.getTime();
			lineCounter++;
			
			if (record == null) continue;
			
			// if SAM format
			String readID = record.getReadName().split("#|\\s+|/|F3|F5-RNA")[0];
			mappedElementProcessingSAMBAM.setReadID(readID);
			//DEBUG
			mappedElementProcessingSAMBAM.setReadID(readID + ":" + record.getReferenceName() + ":" + record.getAlignmentStart() + "-" + record.getAlignmentEnd());
			//END DEBUG
			
			if (pairedEndData == false && record.getReadPairedFlag()) pairedEndData = true;
			
			if (referenceLastRecord != null)
			{				
				assert(referenceLastRecord.equals(record.getReferenceName()) && record.getAlignmentStart() > startPositionLastRecord) : 
					"BAM file not sorted!! Start current mapping " + referenceLastRecord + ":" + record.getAlignmentStart() + 
					", start previous mapping " + record.getReferenceName() + ":" + startPositionLastRecord;
			}
			
			startPositionLastRecord = record.getAlignmentStart();
			referenceLastRecord = record.getReferenceName();
			
			//if(mapping == null || mapping.isUnmapped() || !mapping.allReadsProperlyAligned())	continue;			
			if(record.getReadUnmappedFlag() 
					|| (record.getReadPairedFlag() && record.getMateUnmappedFlag() && this.onlyUsePairedEndReadsWithBothMatesMapped)
					|| (record.getReadPairedFlag() && record.getReferenceName() != record.getMateReferenceName() && this.ignoreReadMappingsToDifferentChromosomes))
			{
				continue;
			}
						
			HashMap<String, HashSet<MappableElement>> mappedElements = new HashMap<String, HashSet<MappableElement>>();
			
			// Check if read is spliced.
			if (record.getCigarString().contains("N"))
			{
				processJunctionsCoveredByMapping(record, mappedElements);
			}
			else
			{
				processExonicRegionsCoveredByMapping(record.getReferenceName(), record.getAlignmentStart(), 
						record.getAlignmentEnd(), mappedElements);
			}
			
			// Mapping with an unmapped mate are treated as single-end reads
			if (record.getReadPairedFlag() && !record.getMateUnmappedFlag())
			{
				/* When the first mate was encountered the transcripts mapped to this read were
				 * stored using the read ID and the start position of its second mate, respectively
				 * the mate it is currently looked at.
				 */				
				int hashKeyCurRead = (readID + ":" + record.getAlignmentStart() 
						+ "-" + record.getMateAlignmentStart()).hashCode();
				int hashKeyMateRead = (readID + ":" + record.getMateAlignmentStart() 
						+ "-" + record.getAlignmentStart()).hashCode();
				
				
				if (pairedEndReadsWithoutMate.containsKey(hashKeyMateRead))
				{
					//DEBUG
					String debugString = "";
					if (mappedElementProcessingSAMBAM instanceof MappedElementProcessingDebug)
					{
						
						debugString = "\n---------------------------------------------------------------\n";
						debugString += "Process paired end read\n";
						debugString += "\n\nCurrent Mapping\n";
						debugString += "\t" + readID + "\t" + record.getReferenceName() 
								+ "\t" + record.getAlignmentStart() + "-" + record.getAlignmentEnd() + "\n";
						debugString += "\t" + record.getCigarString();
						
						SAMRecord mateRecord = pairedEndReadsWithoutMateToMapping.get(hashKeyMateRead);
						debugString += "\nMate Mapping\n";
						debugString += "\t" + readID + "\t" + mateRecord.getReferenceName() 
								+ "\t" + mateRecord.getAlignmentStart() + "-" + mateRecord.getAlignmentEnd() + "\n";
						debugString += "\t" + mateRecord.getCigarString();
	//					ReadMapping mappingMate = pairedEndReadsWithoutMateToMapping.get(hashKeyMateRead);
	//					debugString += "\nMapping mate" + "\n";
	//					debugString += "\t" + mappingMate.getReadID() + "\t" + mappingMate.getReference() 
	//							+ "\t" + mappingMate.getStartPosition() + "-" + mappingMate.getStopPosition() + "\n";
	//					debugString += "\t" + mappingMate.getMappingCIGARString() + "\n";
						((MappedElementProcessingDebug) mappedElementProcessingSAMBAM).debuggingMessage += debugString;
						System.out.println(debugString);
					}
					//END DEBUG
					if (this.processMatesInAscendingStartPositionOrder 
							&& (record.getAlignmentStart() > record.getMateAlignmentStart()))
					{
						processPaireEndReads(pairedEndReadsWithoutMate.get(hashKeyMateRead), mappedElements,
								pairedEndReadsWithoutMateToReference.get(hashKeyMateRead), record.getReferenceName(), debugString);
					}
					else
					{
						processPaireEndReads(mappedElements, pairedEndReadsWithoutMate.get(hashKeyMateRead),
								record.getReferenceName(), pairedEndReadsWithoutMateToReference.get(hashKeyMateRead), debugString);
					}
					
					//Mate of paired end was found s.t. it can be removed from hash collecting
					//reads without mates.
					/* It can occur that reads map to the exact same location, especially in the
					 * case of simulated data. In order to assign paired-end reads correctly in
					 * this case for these reads the occurrence is stored.
					 */
					pairedEndReadsWithoutMate.remove(hashKeyMateRead);
					pairedEndReadsWithoutMateToReference.remove(hashKeyMateRead);
//					//DEBUG
					pairedEndReadsWithoutMateToMapping.remove(hashKeyMateRead);
//					//END DEBUG
				}
				//Waiting for second mate
				else
				{
					pairedEndReadsWithoutMate.put(hashKeyCurRead, mappedElements);
					pairedEndReadsWithoutMateToReference.put(hashKeyCurRead, record.getReferenceName());
					//DEBUG
					pairedEndReadsWithoutMateToMapping.put(hashKeyCurRead, record);
					//END DEBUG
				}
			}
			else
			{
				//DEBUGString debugString = "";
//				if (mappedElementProcessingSAMBAM instanceof MappedElementProcessingDebug)
//				{
//					String debugString = "\n---------------------------------------------------------------\n";
//					debugString += "Process single end read\n";
//					debugString += "\n\nCurrent Mapping\n";
//					debugString += "\t" + readID + "\t" + record.getReferenceName() 
//							+ "\t" + record.getAlignmentStart() + "-" + record.getAlignmentEnd() + "\n";
//					debugString += "\t" + record.getCigarString();
//					((MappedElementProcessingDebug) mappedElementProcessingSAMBAM).debuggingMessage += debugString;
//				}
				//END DEBUG
				processSingleEndReadsOrSingleMates(mappedElements, record.getReferenceName());
			}
			
			//DEBUG
			if (this.maxNumOfGenesStored < this.gtfMapperSorted.getCurGenes().size())
			{
				this.maxNumOfGenesStored = this.gtfMapperSorted.getCurGenes().size();
			}
			if (this.maxNumOfMappedElementsStoredForWaitingMates < pairedEndReadsWithoutMate.values().size())
			{
				this.maxNumOfMappedElementsStoredForWaitingMates = pairedEndReadsWithoutMate.values().size();
			}
			if (this.maxNumOfMatesWaiting < pairedEndReadsWithoutMate.keySet().size())
			{
				this.maxNumOfMatesWaiting = pairedEndReadsWithoutMate.keySet().size();
			}
			//END DEBUG
		}
		reader.close();		
		mappedElementProcessingSAMBAM.finishProcessing();
		
		//Take of paired-end reads where only one of the mates was mapped.
		if (!this.onlyUsePairedEndReadsWithBothMatesMapped)
		{
			for (int keySingleMateRead : pairedEndReadsWithoutMate.keySet())
			{
				processSingleEndReadsOrSingleMates(pairedEndReadsWithoutMate.get(keySingleMateRead),
						pairedEndReadsWithoutMateToReference.get(keySingleMateRead));
			}
		}
		
		System.out.println("Processed " + lineCounter + " lines... ");
		System.out.println("Reads waiting for second mate: " + pairedEndReadsWithoutMate.keySet().size());
		System.out.println("Number of stored elements:\t" + mappedElementProcessingSAMBAM.getNumOfStoredElements());
		
		//DEBUG
		System.out.println("MEMORY\tCurrent Heap Memory Usage:\t" + (ManagementFactory.getMemoryMXBean().getHeapMemoryUsage().getUsed()/1024.0/1024.0) + " MB");
		System.out.println("MEMORY\tNon heap Memory Usage:\t" + (ManagementFactory.getMemoryMXBean().getNonHeapMemoryUsage().getUsed()/1024.0/1024.0) + " MB");
		for (MemoryPoolMXBean item : ManagementFactory.getMemoryPoolMXBeans())
		{
			System.out.println("MEMORY\tName:\t" + item.getName());
			System.out.println("MEMORY\tEstimated Current Usage:\t" + (item.getUsage().getUsed()/1024.0/1024.0) + " MB");
			System.out.println("MEMORY\tPeak Usage:\t" + (item.getPeakUsage().getUsed()/1024.0/1024.0) + " MB");
			System.out.println("MEMORY");
		}
		System.out.println("MEMORY\tMax genes stored: " + this.maxNumOfGenesStored);
		System.out.println("MEMORY\tMax mates waiting: " + this.maxNumOfMatesWaiting);
		System.out.println("MEMORY\tMax mapped elements of mates waiting: " + this.maxNumOfMappedElementsStoredForWaitingMates);
		
		//END DEBUG
		 
	}
	
	private void processExonicRegionsCoveredByMapping(String seqName,
			int regionStart, int regionStop,
			HashMap<String, HashSet<MappableElement>> genesToCoveredElements) throws IOException
	{			
		
		HashMap<String, Vector<MappableElement>> geneToMappedGeneticRegions = 
				this.gtfMapperSorted.mapPositionToGeneticElements(seqName, regionStart, regionStop, true, false);
			
		for (String geneID : geneToMappedGeneticRegions.keySet())
		{
			for (MappableElement curGeneticRegion : geneToMappedGeneticRegions.get(geneID))
			{
				if (!genesToCoveredElements.containsKey(geneID))
					genesToCoveredElements.put(geneID, new HashSet<MappableElement>());
				
				genesToCoveredElements.get(geneID).add(curGeneticRegion);
			}
		}	
	}
	
	private void processJunctionsCoveredByMapping(SAMRecord record,
			HashMap<String, HashSet<MappableElement>> genesToCoveredElements) throws IOException
	{				
		Cigar cigar = new Cigar(record.getCigarString());
		/* In the case of a spliced read getting all exons between
		 * mapping start and end is inaccurate since between these
		 * position there can exist other exons that are not covered
		 * by a read due to being spliced out in the transcript the
		 * read originates from. Therefore the cigar string is 
		 * iterated for getting the precise positions mapped by the
		 * read.
		 */
		int currentPosition = record.getAlignmentStart();
		int curJunctionStart = -1;
		int curJunctionEnd;
		
		int numberOfJunctions = 0;
		HashMap<String, Integer> geneToNumOfJunctionsFound = new HashMap<String, Integer>();

		for (int i = 0; i < cigar.getNumberOfGroups(); i++)
		{
			if (cigar.getSymbolOfGroup(i) == Cigar.SKIP_IN_REFERNCE)
			{
				curJunctionStart = currentPosition - 1;
				curJunctionEnd = currentPosition + cigar.getLengthOfGroup(i);
				Junction curJunction = new Junction(record.getReferenceName(), curJunctionStart, curJunctionEnd);
				numberOfJunctions++;
						
				//Find out genes containing this junction.
				HashMap<String, HashSet<MappableElement>> genesToOverlappingExonicRegionsJunctionStart = new HashMap<String, HashSet<MappableElement>>();
				HashMap<String, HashSet<MappableElement>> genesToOverlappingExonicRegionsJunctionStop = new HashMap<String, HashSet<MappableElement>>();
								
				processExonicRegionsCoveredByMapping(record.getReferenceName(), curJunctionStart, curJunctionStart,  
						genesToOverlappingExonicRegionsJunctionStart);					
				processExonicRegionsCoveredByMapping(record.getReferenceName(), curJunctionEnd, curJunctionEnd, 
						genesToOverlappingExonicRegionsJunctionStop);
				
				HashMap<String, HashSet<String>> geneToTranscriptsWithExonWithinJunction = new HashMap<String, HashSet<String>>();
				
				for (String gene : genesToOverlappingExonicRegionsJunctionStart.keySet())
				{
					
					if (!geneToTranscriptsWithExonWithinJunction.containsKey(gene))
					{
						geneToTranscriptsWithExonWithinJunction.put(gene, 
								getTranscriptsWithExonicRegionWithinJunction(gene, curJunctionStart, curJunctionEnd));
					}
				}
				for (String gene : genesToOverlappingExonicRegionsJunctionStop.keySet())
				{
					if (!geneToTranscriptsWithExonWithinJunction.containsKey(gene))
					{
						geneToTranscriptsWithExonWithinJunction.put(gene, 
								getTranscriptsWithExonicRegionWithinJunction(gene, curJunctionStart, curJunctionEnd));
					}
				}
				
				//Identify exons that define the junction
				for (String gene : genesToOverlappingExonicRegionsJunctionStart.keySet())
				{
					if (!genesToOverlappingExonicRegionsJunctionStop.containsKey(gene)) continue;
					
					/* At first the exons that overlap the junction start or junction
					 * stop are identified. For these exons it is not given that they
					 * end at the junction start or start at the exon end
					 */
					HashSet<MappableElement> overlappingExonicRegionsJunctionStart = 
							genesToOverlappingExonicRegionsJunctionStart.get(gene);
					HashSet<MappableElement> overlappingExonicRegionsJunctionStop =
							genesToOverlappingExonicRegionsJunctionStop.get(gene);
					
					HashSet<MappableElement> exonicRegionsBeforeJunction = new HashSet<MappableElement>();
					HashSet<MappableElement> exonicRegionsAfterJunction = new HashSet<MappableElement>();
						
					// Only used if transcripts for junction should be set.
					HashSet<String> transcriptsWithExonicRegionWithinJunction = new HashSet<String>();
					
					if (this.setCorrespondingTranscriptsForJunctions)
					{				
						transcriptsWithExonicRegionWithinJunction = getTranscriptsWithExonicRegionWithinJunction(gene, 
								curJunctionStart, curJunctionEnd);
					}
					
					HashSet<String> transcriptsOverlappingJunctionStart = new HashSet<String>();
					/* From the exons overlapping the junction start and stop select those
					 * that either end at the junction start or start at the junction stop.
					 * Thus find out the exons defining the junction. 
					 */	
					for (MappableElement exonicRegionOverlappingJunctionStart : overlappingExonicRegionsJunctionStart)
					{	
						// Check whether the current exon can be part of the junction.
						// This is the case if the exon ends at the junctionStart.
						if (exonicRegionOverlappingJunctionStart.getCodingStop() != curJunctionStart) continue;
						exonicRegionsBeforeJunction.add(exonicRegionOverlappingJunctionStart);

						if (this.setCorrespondingExonsForJunctions)
						{
							for (Exon e : exonicRegionOverlappingJunctionStart.getMappedExons())
							{							
								if (e.getCodingStop() != curJunctionStart) continue;
								curJunction.addMappedExon(e);
								
								if (this.setCorrespondingTranscriptsForJunctions)
								{
									transcriptsOverlappingJunctionStart.addAll(e.getMappedTranscripts());
								}
							}
						}
					}
					for (MappableElement exonicRegionOverlappingJunctionStop : overlappingExonicRegionsJunctionStop)
					{
						// Check whether the current exon can be part of the junction.
						// This is the case if the exon start at the junction end.
						if (exonicRegionOverlappingJunctionStop.getCodingStart() != curJunctionEnd) continue;
						exonicRegionsAfterJunction.add(exonicRegionOverlappingJunctionStop);

						if (this.setCorrespondingExonsForJunctions)
						{
							for (Exon e : exonicRegionOverlappingJunctionStop.getMappedExons())
							{		
								if (e.getCodingStart() != curJunctionEnd) continue;
								curJunction.addMappedExon(e);
								
								if (this.setCorrespondingTranscriptsForJunctions)
								{
									for (String transcript : e.getMappedTranscripts())
									{
										if (transcriptsOverlappingJunctionStart.contains(transcript)
												&& !transcriptsWithExonicRegionWithinJunction.contains(transcript))
										{
											curJunction.addMappedTranscript(transcript);
										}
									}
								}
							}
						}
					}
					
					// If either of the hashes is empty then the junction does not exist for
					// the current gene.
					if (exonicRegionsBeforeJunction.isEmpty() || exonicRegionsAfterJunction.isEmpty()) continue;
										
					if (!genesToCoveredElements.containsKey(gene))
					{
						genesToCoveredElements.put(gene, new HashSet<MappableElement>());
					}
					
					if (!geneToNumOfJunctionsFound.containsKey(gene))
					{
						geneToNumOfJunctionsFound.put(gene, 0);
					}
					geneToNumOfJunctionsFound.put(gene, geneToNumOfJunctionsFound.get(gene) + 1);
					
					genesToCoveredElements.get(gene).addAll(exonicRegionsBeforeJunction);
					genesToCoveredElements.get(gene).addAll(exonicRegionsAfterJunction);
					genesToCoveredElements.get(gene).add(curJunction);
				}				

				//DEBUG
//				if (record.getReadName().equals("ERR030880.32808262"))
//				{
//					System.out.println(debugString);
//				}
				//END DEBUG
				
			}
			
			if (cigar.getSymbolOfGroup(i) == 'I') continue;			
			currentPosition += cigar.getLengthOfGroup(i);
		}
		
		for (String gene : geneToNumOfJunctionsFound.keySet())
		{
			if (geneToNumOfJunctionsFound.get(gene) != numberOfJunctions)
			{
				genesToCoveredElements.remove(gene);
			}
		}
		//DEBUG
//		System.out.println(debugString);
		//END DEBUG
	}
	
	private HashSet<String> getTranscriptsWithExonicRegionWithinJunction(String gene, int junctionStart, int junctionStop)
	{
		HashSet<String> transcriptsWithExonWithinJunction = new HashSet<String>();
		Vector<GeneticRegion> exonicRegionsForGene = this.gtfMapperSorted.getExonicRegionsForGene(gene);
		for (GeneticRegion exonicRegion : exonicRegionsForGene)
		{
			if (exonicRegion.getCodingStart() > junctionStart
					&& exonicRegion.getCodingStop() < junctionStop)
			{
				for (Exon e : exonicRegion.getMappedExons())
				{
					if (e.getCodingStart() > junctionStart
							&& e.getCodingStop() < junctionStop)
					{
						transcriptsWithExonWithinJunction.addAll(e.getMappedTranscripts());
					}
				}
			}
		}
		return transcriptsWithExonWithinJunction;
	}
	
	private void processPaireEndReads(HashMap<String, HashSet<MappableElement>> genesToCoveredElementsFirstMate, 
			HashMap<String, HashSet<MappableElement>> genesToCoveredElementsSecondMate, 
			String referenceFirstMate, String referenceSecondMate, String debugString) throws IOException
	{		
		if (this.ignoreReadMappingsToDifferentChromosomes &&
				!referenceFirstMate.equals(referenceSecondMate))
		{
			return;
		}
						
		if (this.onlyUsePairedEndReadsWithBothMatesMapped &&
				genesToCoveredElementsFirstMate.isEmpty())
		{
			return;
		}
		
		if (this.onlyUsePairedEndReadsWithBothMatesMapped &&
				genesToCoveredElementsSecondMate.isEmpty())
		{
			return;
		}
				
		if (genesToCoveredElementsFirstMate.isEmpty() 
				&& genesToCoveredElementsSecondMate.isEmpty()) return;
		
		HashSet<String> genesMappedByBothReads = new HashSet<String>();
		for (String gene : genesToCoveredElementsFirstMate.keySet())
		{
			if (genesToCoveredElementsSecondMate.containsKey(gene))
				genesMappedByBothReads.add(gene);
		}		
		//DEBUG
		debugString += "\nElements mapped first mate\n";
		for (String gene : genesToCoveredElementsFirstMate.keySet())
		{
			debugString += "\t" + gene + "\n";
			for (MappableElement element : genesToCoveredElementsFirstMate.get(gene))
			{
				debugString += "\t\t" + element.getCodingStart() + "-" + element.getCodingStop() + "\n";
			}
		}
		
		debugString += "Elements mapped second mate" + "\n";
		for (String gene : genesToCoveredElementsSecondMate.keySet())
		{
			debugString += "\t" + gene + "\n";
			for (MappableElement element : genesToCoveredElementsSecondMate.get(gene))
			{
				debugString += "\t\t" + element.getCodingStart() + "-" + element.getCodingStop() + "\n";
			}
		}
		
		debugString += "Genes mapped by both:" + "\n";
		for (String gene : genesMappedByBothReads)
		{
			debugString += "\t" + gene + "\n";
		}		
		
//		System.out.println(debugString);
		//END DEBUG
			
		this.mappedElementProcessingSAMBAM.storeProcessingResultsPairedEnd(genesMappedByBothReads, genesToCoveredElementsFirstMate, genesToCoveredElementsSecondMate, 
				referenceFirstMate, referenceSecondMate, this.ignoreReadMappingsToDifferentGenes);		
	}
	
	private void processSingleEndReadsOrSingleMates(HashMap<String, HashSet<MappableElement>> genesToCoveredElements,
			String reference) throws IOException
	{				
		if (genesToCoveredElements.isEmpty()) return;
		
		this.mappedElementProcessingSAMBAM.storeProcessingResultsSingleEnd(genesToCoveredElements, reference);
	}
	
	public boolean isPairedEndData()
	{
		return this.pairedEndData;
	}

	public void setSetCorrespondingExonsForJunctions(
			boolean setCorrespondingExonsForJunctions) 
	{
		this.setCorrespondingExonsForJunctions = setCorrespondingExonsForJunctions;
	}

	public void setSetCorrespondingTranscriptsForJunctions(
			boolean setCorrespondingTranscriptsForJunctions) 
	{
		this.setCorrespondingTranscriptsForJunctions = setCorrespondingTranscriptsForJunctions;
	}	
}
