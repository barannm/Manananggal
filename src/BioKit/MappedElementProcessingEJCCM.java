package BioKit;

import java.util.Vector;

public interface MappedElementProcessingEJCCM 
{
	/* Processes the genes or transcripts mapped
	 * by a combination ID.
	 */
	public void processElements(Vector<String> elementsMappedByCombID, String combID, int occurrence);
	
	public void finishProcessing();
}
