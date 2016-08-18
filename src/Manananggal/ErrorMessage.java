package Manananggal;

import org.zkoss.zul.Messagebox;

/**
 *    Helper class used to output errors. If a GUI is available, the error
 *    will also be shown in a popup window.
 */
public class ErrorMessage
{
	ErrorMessage()
	{
	}
	
	static void ShowError(Object... args)
	{
		String strOut = "";
		for(Object item : args)
		{
			strOut += item;
		}
		strOut += "\n";
		
		// output to console
		System.out.println(strOut);
		
		// try to show a popup window
		try
		{
			Messagebox.show(strOut);
		}
		catch(Exception e)
		{
			// fails when not using a GUI
		}
	}
}
