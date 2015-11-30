package Manananggal;

import java.io.IOException;

import org.broad.igv.bbfile.BBFileReader;

public class BigWigReader extends BBFileReader
{
	public BigWigReader(String path) throws IOException
	{
		super(path);
	}
	
	public void close() throws IOException
	{
		getBBFis().close();
	}

}
