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