package BioKit;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

public class StreamGobbler extends Thread
{
	private InputStream is;
	private BufferedReader br;
	//private String type;
	private boolean collectOutput;
	
	private StringBuffer text;

	public StreamGobbler(InputStream is, String type, boolean collectOutput)
	{
		this.is = is;
		//this.type = type;
		text = new StringBuffer();
		this.collectOutput = collectOutput;
	}

	public void run()
	{
		try
		{
			InputStreamReader isr = new InputStreamReader(is);
			br = new BufferedReader(isr);
			String line = null;
			while ((line = br.readLine()) != null)
			{
				if(collectOutput)
					text.append(line + "\n");
			}
		}
		catch (IOException ioe)
		{
			//ioe.printStackTrace();
		}
	}
	
	public String getMessages()
	{
		return text.toString();
	}
	
	public void kill()
	{
		try
		{
			br.close();
		}
		catch(Exception ex)
		{
			// do nothing
		}
	}
}