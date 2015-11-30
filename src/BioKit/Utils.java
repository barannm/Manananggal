package BioKit;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.net.MalformedURLException;
import java.net.URL;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.ResultSetMetaData;
import java.sql.Statement;
import java.util.Collection;
import java.util.HashMap;
import java.util.Properties;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;
import java.util.zip.GZIPInputStream;

/**
 * This class contains some methods useful for many non-project specific tasks
 * 
 * @author Fabian Birzele
 */
public class Utils
{ 
	/**
	 * Method copies a file from one direction to the other by reading the file line by line 
	 * and writing another file in the target directory with the same content.
	 * 
	 * @param from
	 * @param to
	 * @return
	 */
	public static boolean copyPerStream(File from, File to)
	{
		try
		{
			BufferedReader reader = new BufferedReader(new FileReader(from));
			PrintWriter writer = new PrintWriter(new FileWriter(to));
			
			String currentLine = null;
			
			while((currentLine = reader.readLine()) != null)
			{
				writer.println(currentLine);
			}
			
			reader.close();
			
			writer.flush();
			writer.close();
			
			return true;
		}
		catch(IOException ex)
		{
			return false;
		}
	}
	
	/**
	 * Method can be used to save an instance as serialzed object to the specified
	 * output directory. The name of the file will be ID + ".ser"
	 * 
	 * @param outputDirectory
	 * @param object
	 * @param id
	 * @throws IOException
	 */
	public static void serializeObject(String outputDirectory, String id, Object object) throws IOException
	{
		FileOutputStream out = new FileOutputStream(outputDirectory + File.separatorChar + id + ".ser");
		ObjectOutputStream s = new ObjectOutputStream(out);
		s.writeObject(object);
		s.flush();
		s.close();
	}
	
	/**
	 * Method can be used to save an instance as serialzed object. The name of the file is specified in the given parameter
	 * 
	 * @param outputDirectory
	 * @param object
	 * @param id
	 * @throws IOException
	 */
	public static void serializeObject(String fileName, Object object) throws IOException
	{
		FileOutputStream out = new FileOutputStream(fileName);
		ObjectOutputStream s = new ObjectOutputStream(out);
		s.writeObject(object);
		s.flush();
		s.close();
	}
	
	/**
	 * Method can be used to read an instance from a serialized object file. 
	 * 
	 * @param objectFile
	 * @return
	 * @throws IOException
	 */
	public static Object readSerializedObject(File objectFile) throws IOException
	{
		FileInputStream in = new FileInputStream(objectFile);
		ObjectInputStream s = new ObjectInputStream(in);
		
		try
		{
			Object object = s.readObject();
			s.close();
			return object;
		}
		catch(ClassNotFoundException ex)
		{
			s.close();
			return null;
		}
	}
	
	/**
	 * Tries 100x to delete the specific file
	 * 
	 * @param file
	 */
	public static void removeFile(String file)
	{
		File f = new File(file);
		
		int counter=0;
		
		while(!f.delete())
		{
			counter++;
			
			if(counter > 100)
				break;
		}
	}
	
	/**
	 * This method uses wget to obain a file from the internet.
	 * 
	 * @param url
	 * @param fileToBeWritten
	 * @throws Exception
	 */
	public static void downloadDataFromInternet(String url, String fileToBeWritten, String proxyUser, String proxy) throws Exception
	{
		writeAndExecuteShellScriptOnTheFly("#!/bin/bash\ncurl --speed-time 7200 --speed-limit 0 --proxy-ntlm --proxy-user "
			+proxyUser+" --proxy "+proxy+" --url '"+url+"' -o " +fileToBeWritten);
	}
	
	/**
	 * Uses gunzip to extract a gz file
	 * @param file
	 * @throws Exception
	 */
	public static void extractGZIPFile(String file) throws Exception
	{
		writeAndExecuteShellScriptOnTheFly("gunzip " + file);
	}
	
	/**
	 * Method allows to concatenate a set of files
	 * 
	 * @param directory
	 * @param outputFile
	 * @param filter
	 * @throws Exception
	 */
	public static void catFiles(String directory, String outputFile, String filter) throws Exception
	{
		File[] files = (new File(directory)).listFiles();
		
		for(File f : files)
		{
			// ignore files not satisfying filter criterion
			if(!f.getAbsolutePath().contains(filter))
			{
				System.out.println("Ignore file: " + f.getAbsolutePath() + " does not match filter: " + filter);
				continue;
			}
			
			System.out.println("Concatenate " + f.getAbsolutePath() + " to file " + outputFile);
			
			BufferedReader reader = new BufferedReader(new FileReader(f));
			String line = null;
			
			PrintWriter writer = new PrintWriter(new FileWriter(outputFile, true));
			
			while((line = reader.readLine()) != null)
			{
				writer.println(line);
			}
			reader.close();
			
			writer.flush();
			writer.close();
		}
	}
	
	/**
	 * Method to write and execute a shell script on the fly
	 * 
	 * @param text
	 * @throws Exception
	 */
	public static void writeAndExecuteShellScriptOnTheFly(String text) throws Exception
	{
		File dummy = new File("dummyScript" +(int)(Math.random()*10000) + ".sh");
		
		PrintWriter writer = new PrintWriter(new FileWriter(dummy));
		writer.println(text);
		writer.flush();
		writer.close();
		
		Utils.runProcess("chmod 700 " + dummy.getAbsolutePath(), true, true, true);
		Utils.runProcess(dummy.getAbsolutePath(), true, true, true);
		
		dummy.delete();
	}
	
	/**
	 * Method to list files in a directory in an array
	 * 
	 * @param directory
	 * @param filter
	 * @param recurse
	 * @return
	 */
	public static File[] listFilesAsArray(File directory, String filter, boolean recurse)
	{
		Collection<File> files = listFiles(directory, filter, recurse);
		File[] arr = new File[files.size()];
		return files.toArray(arr);
	}

	/**
	 * Recursive function to list all files in a directory and all subdirectories 
	 * @param directory
	 * @param filter
	 * @param recurse
	 * @return
	 */
	public static Collection<File> listFiles(File directory, String filter, boolean recurse)
	{
		// List of files / directories
		Vector<File> files = new Vector<File>();
	
		// Get files / directories in the directory
		File[] entries = directory.listFiles();
		
		// Go over entries
		for (File entry : entries)
		{
			// If there is no filter or the filter accepts the 
			// file / directory, add it to the list
			if (entry.getAbsolutePath().contains(filter))
			{
				files.add(entry);
			}
			
			// If the file is a directory and the recurse flag
			// is set, recurse into the directory
			if (recurse && entry.isDirectory())
			{
				files.addAll(listFiles(entry, filter, recurse));
			}
		}
		
		// Return collection of files
		return files;		
	}

	
	/**
	 * Method runs a process with the specified process call. The output of the process is returned as a String[].
	 * If the process caused errors, null is returned in order to indicate that the process was not successful. 
	 * In that case, detailed information about the error type is written to the Standard Error stream.
	 * 
	 * @param processCall
	 * @param returnNullIfErrors TODO
	 * @param outputErrorMessages TODO
	 * @return
	 */
	public static String[] runProcess(String processCall, boolean printOutputToStandardOut, boolean collectOutput, boolean returnNullIfErrors)
	{
		try
		{
			Runtime rt = Runtime.getRuntime();
			Process proc = rt.exec(processCall);
			
			// any error message?
		    StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR", collectOutput);            
		            
		    // any output?
		    StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT", collectOutput);
		                
		    // kick them off
		    errorGobbler.start();
		    outputGobbler.start();
		        
		    proc.waitFor();
		    
		    //System.out.println("Process finished...");
		        
		    String errorText = errorGobbler.getMessages();
		    String[] lines = errorText.split("\\n");
	
			// print out system out information
		    boolean successful = true;
		        
			for(int i=0; i<lines.length; i++)
			{
				if(!lines[i].trim().equals(""))
				{
					if(printOutputToStandardOut)
						System.err.println("\nERROR Utils.runProcess(): " + lines[i]);
					
					successful = false;
				}
			}
			
			// return null if call has caused errors
			if ( !successful ) {
				outputGobbler.kill();
				errorGobbler.kill();
				
				//System.out.println("Killed threads...");
				
				// close all streams
				proc.getErrorStream().close();
				proc.getInputStream().close();
				proc.getOutputStream().close();
				
				//System.out.println("Closed streams...");
				
				proc.destroy();
				
				//System.out.println("Destroyed process...");

				if(returnNullIfErrors)
					return null;				
			}
			
			String outputText = outputGobbler.getMessages();
			lines = outputText.split("\\n");
			
			for(int i=0; i<lines.length; i++)
			{
				if(!lines[i].trim().equals(""))
				{
					if(printOutputToStandardOut)
						System.out.println("\nOUTPUT Utils.runProcess(): " + lines[i]);
				}
			}
			
			outputGobbler.kill();
			errorGobbler.kill();
			
			// close all streams
			proc.getErrorStream().close();
			proc.getInputStream().close();
			proc.getOutputStream().close();
			proc.destroy();
			
			return lines;
		}
		catch(Exception ex)
		{
			ex.printStackTrace();
			System.out.println("Could not run " + processCall + " due to error " + ex.toString());
		}
		
		return null;
	}
	
	/**
	 * Method runs a process with the specified process call. The output of the process is returned as a String[].
	 * If the process caused errors, null is returned in order to indicate that the process was not successful. 
	 * In that case, detailed information about the error type is written to the Standard Error stream.
	 * 
	 * @param processCall
	 * @param returnNullIfErrors TODO
	 * @param shellSpecification TODO
	 * @param outputErrorMessages TODO
	 * @return
	 */
	public static String[] runProcessInShell(String processCall, boolean outputMessages, boolean returnNullIfErrors, String shellSpecification)
	{
		try
		{
			String shell = "/bin/csh";
			
			if(shellSpecification != null)
				shell = shellSpecification;
			
			String[] cmd = new String[]{shell, "-c", processCall};
			
			Runtime rt = Runtime.getRuntime();
			Process proc = rt.exec(cmd);
			
			// any error message?
		    StreamGobbler errorGobbler = new StreamGobbler(proc.getErrorStream(), "ERROR", outputMessages);            
		            
		    // any output?
		    StreamGobbler outputGobbler = new StreamGobbler(proc.getInputStream(), "OUTPUT", outputMessages);
		                
		    // kick them off
		    errorGobbler.start();
		    outputGobbler.start();
		        
		    proc.waitFor();
		    
		    //System.out.println("Process finished...");
		        
		    String errorText = errorGobbler.getMessages();
		    String[] lines = errorText.split("\\n");
	
			// print out system out information
		    boolean successful = true;
		        
			for(int i=0; i<lines.length; i++)
			{
				if(!lines[i].trim().equals(""))
				{
					if(outputMessages)
						System.err.println("\nERROR Utils.runProcess(): " + lines[i]);
					
					successful = false;
				}
			}
			
			// return null if call has caused errors
			if ( !successful ) {
				outputGobbler.kill();
				errorGobbler.kill();
				
				//System.out.println("Killed threads...");
				
				// close all streams
				proc.getErrorStream().close();
				proc.getInputStream().close();
				proc.getOutputStream().close();
				
				//System.out.println("Closed streams...");
				
				proc.destroy();
				
				//System.out.println("Destroyed process...");

				if(returnNullIfErrors)
					return null;				
			}
			
			String outputText = outputGobbler.getMessages();
			lines = outputText.split("\\n");
			
			for(int i=0; i<lines.length; i++)
			{
				if(!lines[i].trim().equals(""))
				{
					if(outputMessages)
						System.err.println("\nOUTPUT Utils.runProcess(): " + lines[i]);
				}
			}
	
			outputGobbler.kill();
			errorGobbler.kill();
			
			// close all streams
			proc.getErrorStream().close();
			proc.getInputStream().close();
			proc.getOutputStream().close();
			proc.destroy();
			
			return lines;
		}
		catch(Exception ex)
		{
			ex.printStackTrace();
			System.out.println("Could not run " + processCall + " due to error " + ex.toString());
		}
		
		return null;
	}
	
	/**
	 * Method retrieves the content of an internet page at the specified url from the internet.
	 * All textual content of the page is returned as a string instance.
	 * 
	 * @param url
	 * @return
	 */
	public static String getWebpageFromInternet(String url)
	{
		StringBuffer entryText = new StringBuffer();
		
		 try 
		 {
			 URL u = new URL(url);
			 InputStream in = u.openStream();
			  
			int b;
			
			 while ((b = in.read()) != -1) 
			 {
				 entryText.append((char)(b));
			 }
		 }
		catch (MalformedURLException e) 
		{
			System.err.println(e);
		}
		catch (IOException e) 
		{
			System.err.println(e);
		}

		return entryText.toString();
	}
	
	/**
	 * Method retrieves the content of an internet page at the specified url from the internet.
	 * All textual content of the page is written to the specified file.
	 * 
	 * @param url
	 * @return
	 */
	public static void getWebpageFromInternetAndWriteToFile(String url, String outputFile)
	{
		try 
		{
			PrintWriter writer = new PrintWriter(new FileWriter(outputFile));
			
			URL u = new URL(url.replaceAll("\\s+", "%20"));
			
			InputStream in = u.openStream();
			
			System.out.println("Opened stream...");
			  
			int b;
			
			boolean bitWritten = false;
			
			while ((b = in.read()) != -1) 
			{
				if(!bitWritten)
				{
					bitWritten = true;
					System.out.println("Started writing...");
					
				}
				writer.print((char)(b));
			}
			
			writer.flush();
			writer.close();
		}
		catch (MalformedURLException e) 
		{
			e.printStackTrace();
			System.err.println(e);
		}
		catch (IOException e) 
		{
			e.printStackTrace();
			System.err.println(e);
		}
	}
	
	/**
	 * Method converts an array of strings to a string with elements separated by the specified separator
	 * @param array
	 * @param separator
	 * @return
	 */
	public static String stringArrayToString(String[] array, String separator)
	{
		StringBuffer buffer = new StringBuffer();
		
		for(int i=0; i<array.length; i++)
		{
			buffer.append(array[i]);
			
			if(i<array.length-1)
				buffer.append(separator);
		}
		
		return buffer.toString();
	} 
	
	/**
	 * Method converts an array of ints to a string with elements separated by the specified separator
	 * 
	 * @param array
	 * @param separator
	 * @return
	 */
	public static String intArrayToString(int[] array, String separator)
	{
		StringBuffer buffer = new StringBuffer();
		
		for(int i=0; i<array.length; i++)
		{
			buffer.append(array[i]);
			
			if(i<array.length-1)
				buffer.append(separator);
		}
		
		return buffer.toString();
	}

	/**
	 * Method adds leading zeros to a number such that each number has digits digits ;-)
	 * The number is returned as a String.
	 * 
	 * @param digits
	 * @param number
	 * @return
	 */
	public static String addLeadingZeroes(int digits, int number)
	{
		StringBuffer string = new StringBuffer();
		int additionalDigits;
		
		if(number == 0)
			additionalDigits = digits -1 ;
		else
			additionalDigits = digits - 1 - (int)Math.floor(Math.log10(number));
	
		for(int i=0;i<additionalDigits;i++)
		{
			string.append("0");
		}
		string.append(number);
		return string.toString();
	}
	
	/**
	 * This little class extracts numeric values from a tab or space separated
	 * file. The values are written to standard out.
	 * 
	 * @param args
	 */
	public static void extractNumericValues(String[] args) throws Exception
	{
		int[] values = new int[args.length-2];
		
		File file = new File(args[1]);
		
		for(int i=2; i<args.length; i++)
		{
			values[i-2]=Integer.parseInt(args[i]);
		}
		
		BufferedReader reader = new BufferedReader(new FileReader(file));
		
		String line = null;
		
		while((line = reader.readLine()) != null)
		{
			String[] splitArray = line.split("\\s+");
			
			StringBuffer output = new StringBuffer();
			
			for(int i=0; i<values.length; i++)
			{
				output.append(splitArray[values[i]].replace(',', ' '));
				output.append("\t");
			}
			
			System.out.println(output.toString());
		}
		reader.close();
	}
	
	public static void rowsToColumns(String file, int valueColumn, int columnHeaderColumn, 
		String rowNameColumnsString, String splitRegExp, String idSplitRegExp) throws Exception
	{
		Vector<Integer> rowNameColumns = new Vector<Integer>();
		
		String[] rowNameCols = rowNameColumnsString.split(",");
		
		for(String n : rowNameCols)
		{
			rowNameColumns.add(Integer.parseInt(n));
		}
		
		TreeSet<String> columnHeaders = new TreeSet<String>();
		
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line = null;
		
		HashMap<String, HashMap<String, String>> rowsToColumnsToValues = new HashMap<String, HashMap<String,String>>();
		
		while((line = reader.readLine()) != null)
		{
			String[] split = line.split(splitRegExp);
			
			StringBuffer rowIDBuffer = new StringBuffer();
			
			for(int col : rowNameColumns)
			{
				rowIDBuffer.append(split[col] + " ");
			}
			
			String rowID = rowIDBuffer.toString();
			
			String colID = split[columnHeaderColumn].split(idSplitRegExp)[0];
			
			columnHeaders.add(colID);
			
			String value = split[valueColumn];
			
			if(!rowsToColumnsToValues.containsKey(rowID))
				rowsToColumnsToValues.put(rowID, new HashMap<String, String>());
			
			rowsToColumnsToValues.get(rowID).put(colID, value);
		}
		reader.close();
		
		System.out.print("Row ID");
		
		for(String c : columnHeaders)
		{
			System.out.print("\t"+c);
		}
		
		System.out.println();
		
		for(String r : rowsToColumnsToValues.keySet())
		{
			System.out.print(r);
			
			HashMap<String, String> values = rowsToColumnsToValues.get(r);
			
			for(String c : columnHeaders)
			{
				if(values.containsKey(c))
					System.out.print("\t"+values.get(c));
				else
					System.out.print("\t");
			}
			
			System.out.println();
		}
		
	}
	
	/**
	 * This method counts unique values in one column of a a tab or space separated file.
	 * 
	 * @param file
	 * @param column
	 * @throws Exception
	 */
	public static void uniqueColumnCounter(String file, int column, String splitRegExp) throws Exception
	{
		TreeSet<String> uniqueThings = new TreeSet<String>();
		
		BufferedReader reader = new BufferedReader(new FileReader(file));
		
		String line = null;
		
		while((line = reader.readLine()) != null)
		{
			try
			{
				String[] split = line.split(splitRegExp);
				uniqueThings.add(split[column]);
			}
			catch(Exception ex)
			{
				System.out.println("ERROR reading line: " + line + " " + ex.toString());;
			}
		}
		reader.close();
		
		System.out.println("File contains " + uniqueThings.size() + " unique entries in column " + column);
	}
	
	/**
	 * This method counts unique values in one column of a a tab or space separated file and outputs them afterwards
	 * 
	 * @param file
	 * @param column
	 * @throws Exception
	 */
	public static void uniqueColumnCounterWithOutput(String file, int column, String splitCriterion, boolean completeLine) throws Exception
	{
		TreeMap<String, String> uniqueThings = new TreeMap<String, String>();
		
		BufferedReader reader = new BufferedReader(new FileReader(file));
		
		String line = null;
		
		while((line = reader.readLine()) != null)
		{
			try
			{
				String[] split = line.split(splitCriterion);
				uniqueThings.put(split[column], line);
			}
			catch(Exception ex)
			{
				System.out.println("ERROR reading line: " + line + " " + ex.toString());;
			}
		}
		reader.close();
		
		if(!completeLine)
		{
			for(String s : uniqueThings.keySet())
			{
				System.out.println(s);
			}
		}
		else
		{
			for(String s : uniqueThings.keySet())
			{
				System.out.println(uniqueThings.get(s));
			}
		}
	}
    
    public static void printHistogram(String file, int column, String splitter, String columnSplitter, int columnValuePosition) throws Exception
    {
    	TreeMap<String, Integer> histogram = new TreeMap<String, Integer>();
    	
    	BufferedReader reader = new BufferedReader(new FileReader(file));
		
		String line = null;
		
		while((line = reader.readLine()) != null)
		{
			try
			{
				String[] split = line.split(splitter);
				
				String key = split[column].trim();
				
				if(columnSplitter != null)
					key = key.split(columnSplitter)[columnValuePosition];
				
				if(histogram.containsKey(key))
					histogram.put(key, histogram.get(key) + 1);
				else 
					histogram.put(key, 1);
			}
			catch(Exception ex)
			{
				System.out.println("ERROR reading line: " + line + " " + ex.toString());;
			}
		}
		reader.close();
		
		System.out.println("Histogram:");
		System.out.println("Class\tFrequency");
		
		for(String k : histogram.keySet())
		{
			System.out.println(k + "\t" + histogram.get(k));
		}
    }
    
    public static void printHistogramTwoColumns(String file, int columnOne, int columnTwo) throws Exception
    {
    	TreeSet<Integer> classes = new TreeSet<Integer>();
    	TreeMap<Integer, Integer> histogram = new TreeMap<Integer, Integer>();
    	TreeMap<Integer, Integer> histogram2 = new TreeMap<Integer, Integer>();
    	
    	BufferedReader reader = new BufferedReader(new FileReader(file));
		
		String line = null;
		
		while((line = reader.readLine()) != null)
		{
			try
			{
				String[] split = line.split("\\s+");
				
				int key = Integer.parseInt(split[columnOne]);
				int key2 = Integer.parseInt(split[columnTwo]);
				
				classes.add(key);
				classes.add(key2);
				
				if(histogram.containsKey(key))
					histogram.put(key, histogram.get(key) + 1);
				else 
					histogram.put(key, 1);
				

				if(histogram2.containsKey(key2))
					histogram2.put(key2, histogram2.get(key2) + 1);
				else 
					histogram2.put(key2, 1);
				
				
			}
			catch(Exception ex)
			{
				System.out.println("ERROR reading line: " + line + " " + ex.toString());;
			}
		}
		reader.close();
		
		System.out.println("Histogram:");
		System.out.println("Class\tF1\tF2");
		
		int totalCounts = 0;
		
		for(Integer k : classes)
		{
			int value1 = 0;
			
			if(histogram.containsKey(k))
				value1 = histogram.get(k);
			
			int value2 = 0;
			
			if(histogram2.containsKey(k))
				value2 = histogram2.get(k);
			
			totalCounts += value1;
			
			System.out.println(k + "\t" + value1 + "\t" + value2);
		}
		
		System.out.println("Total number of values in F1: " + totalCounts);
    }
    
    public static void printHistogramOfNumericValues(String file, int column, int bins, String splitCrit, int maxValue) throws Exception
    {
    	Vector<Double> values = new Vector<Double>();
    	
    	BufferedReader reader = new BufferedReader(new FileReader(file));
    	
    	String line = null;
    	
    	double max = 0;
    	
    	while((line = reader.readLine()) != null)
    	{
    		String[] split = line.split(splitCrit);
    		
    		try
    		{
    			double v = Double.parseDouble(split[column].trim());
    			values.add(v);
    			
    			if(max < v)
    				max = v;
    		}
    		catch(Exception ex)
    		{
    			System.out.println("Can not parse line: " + line);
    		}
    	}
    	reader.close();
    	
    	int[] counterArray = new int[bins];
    	
    	for(int i=0; i<counterArray.length; i++)
    		counterArray[i] = 0;
    	
    	double step = max / bins;
    	
    	if(maxValue != -1)
    		step = (maxValue / (double)bins);
    	
    	for(double d : values)
    	{
    		int bin = (int) Math.floor(d / step);
    		
    		if(bin >= bins)
    			bin = bins-1;
    		
    		counterArray[bin]++;
    	}
    	
    	for(int i=0; i<counterArray.length; i++)
    		System.out.println((i*step) + "\t" + counterArray[i]);
    }
    
    public static void sortLinesWithRespectToColumn(String file, int sorterColumn, String splitCrit) throws Exception
    {
    	TreeSet<NumericSorterContainer<String>> sorter = new TreeSet<NumericSorterContainer<String>>();
    	
    	BufferedReader reader = new BufferedReader(new FileReader(file));
    	
    	String line = null;
    	
    	while((line = reader.readLine()) != null)
    	{
    		String[] split = line.split(splitCrit);
    		
    		try
    		{
    			sorter.add(new NumericSorterContainer<String>(new Double(Double.parseDouble(split[sorterColumn].trim())), line, true));
    		}
    		catch(Exception ex)
    		{
    			System.out.println("Can not parse line: " + line);
    		}
    	}
    	
    	reader.close();
    	
    	for(NumericSorterContainer<String> c : sorter)
    		System.out.println(c.getContent());
    }
    
    public static void sortLinesWithRespectToColumnNonNumeric(String file, int sorterColumn, String splitCrit) throws Exception
    {
    	TreeMap<String, Vector<String>> sorter = new TreeMap<String, Vector<String>>();
    	
    	BufferedReader reader = new BufferedReader(new FileReader(file));
    	
    	String line = null;
    	
    	while((line = reader.readLine()) != null)
    	{
    		String[] split = line.split(splitCrit);
    		
    		try
    		{
    			String id = split[sorterColumn].trim();
    			
    			if(!sorter.containsKey(id))
    				sorter.put(id, new Vector<String>());
    			
    			sorter.get(id).add(line);
    		}
    		catch(Exception ex)
    		{
    			System.out.println("Can not parse line: " + line);
    		}
    	}
    	reader.close();
    	
    	for(String id : sorter.keySet())
    	{
    		for(String l : sorter.get(id))
    			System.out.println(l);
    	}
    }
    
    public static void extractLinesWhereValueInColumnMatchesCriterion(String file, int matcherColumn, String valueRegExp, String splitCrit) throws Exception
    {
    	BufferedReader reader = new BufferedReader(new FileReader(file));
    	
    	String line = null;
    	
    	while((line = reader.readLine()) != null)
    	{
    		String[] split = line.split(splitCrit);
    		
    		try
    		{
    			if(split[matcherColumn].matches(valueRegExp))
    				System.out.println(line);
    		}
    		catch(Exception ex)
    		{
    			System.out.println("Can not parse line: " + line);
    		}
    	}
    	reader.close();
    }
    
    public static void filterLines(String file, int filterColumn, String filterFile, int filterIDColumn, String splitCrit, boolean outputLinesWhichDoNotContainFilterIDs) throws Exception
    {
    	TreeSet<String> filter = new TreeSet<String>();
    	BufferedReader reader = new BufferedReader(new FileReader(filterFile));
    	
    	String line = null;
    	
    	while((line = reader.readLine()) != null)
    	{
    		String[] split = line.split(splitCrit);
    		
    		filter.add(split[filterIDColumn].trim());
    		
    		//System.out.println("Filter ID: " + split[filterIDColumn]);
    	}
    	reader.close();
    	
    	reader = new BufferedReader(new FileReader(file));
    	
    	line = null;
    	
    	while((line = reader.readLine()) != null)
    	{
    		try
    		{
	    		String[] split = line.split(splitCrit);
	    		
	    		if(!outputLinesWhichDoNotContainFilterIDs && filter.contains(split[filterColumn]))
	    			System.out.println(line);
	    		else if(outputLinesWhichDoNotContainFilterIDs && !filter.contains(split[filterColumn]))
	    			System.out.println(line);
    		}
    		catch(Exception ex)
    		{
    			System.err.println("Could not read line " + line);
    		}
    			
    	}
    	reader.close();
    }
    
    public static void extractLinesWhereValueInColumndDoesNotMatchCriterion(String file, int matcherColumn, String valueRegExp, String splitCrit) throws Exception
    {
    	BufferedReader reader = new BufferedReader(new FileReader(file));
    	
    	String line = null;
    	
    	while((line = reader.readLine()) != null)
    	{
    		String[] split = line.split(splitCrit);
    		
    		if(!split[matcherColumn].matches(valueRegExp))
    			System.out.println(line);
    	}
    	reader.close();
    }
    
    public static void extractLinesWhereNumericValueInColumnIsSmallerEqualOrGreaterThanUserValue(String file, int matcherColumn, double value, String splitCrit, boolean smaller, boolean equal, boolean greater) throws Exception
    {
    	BufferedReader reader = new BufferedReader(new FileReader(file));
    	
    	String line = null;
    	
    	while((line = reader.readLine()) != null)
    	{
    		try
    		{
	    		String[] split = line.split(splitCrit);
	    		
	    		double v = Double.parseDouble(split[matcherColumn]);
	    		
	    		if(smaller && v < value)
	    			System.out.println(line);
	    		else if(equal && v == value)
	    			System.out.println(line);
	    		else if(greater && v > value)
	    			System.out.println(line);
    		}
    		catch(Exception ex)
    		{
    			System.err.println("Could not parse line: " + line);
    		}
    	}
    	reader.close();
    }
    
    public static void extractLinesWhereNumericValueInColumnIsWithinInterval(String file, int matcherColumn, double lowerBound, double upperBound, String splitCrit, boolean searchForWithinInterval) throws Exception
    {
    	BufferedReader reader = new BufferedReader(new FileReader(file));
    	
    	String line = null;
    	
    	while((line = reader.readLine()) != null)
    	{
    		String[] split = line.split(splitCrit);
    		
    		try
    		{
    			double v = Double.parseDouble(split[matcherColumn]);
    		
	    		if(searchForWithinInterval && v <= upperBound && v >= lowerBound)
	    		{
	    			System.out.println(line);
	    		}
	    		else if(!searchForWithinInterval && (v > upperBound || v < lowerBound))
	    		{
	    			System.out.println(line);
	    		}
    		}
    		catch(Exception ex)
    		{
    			// do nothing, just ignore line.
    		}
    	}
    	reader.close();
    }
    
    /**
	 * This method counts how many values in the corresponding column of file 1 are unique in File 1 and not contained in File 2.
	 * 
	 * @param file
	 * @param column
	 * @throws Exception
	 */
	public static void nonUniqueColumnCounter(String fileOne, String fileTwo, int columnOne, int columnTwo) throws Exception
	{
		TreeSet<String> uniqueThingsOne = new TreeSet<String>();
		
		BufferedReader reader = new BufferedReader(new FileReader(fileOne));
		
		String line = null;
		
		while((line = reader.readLine()) != null)
		{
			try
			{
				String[] split = line.split("\\s+");
				uniqueThingsOne.add(split[columnOne]);
			}
			catch(Exception ex)
			{
				System.out.println("ERROR reading line: " + line + " " + ex.toString());;
			}
		}
		reader.close();
		
		System.out.println("File "+fileOne+" contains " + uniqueThingsOne.size() + " unique entries in column " + columnOne);
		
		TreeSet<String> uniqueThingsTwo = new TreeSet<String>();
		
		reader = new BufferedReader(new FileReader(fileTwo));
		
		line = null;
		
		while((line = reader.readLine()) != null)
		{
			try
			{
				String[] split = line.split("\\s+");
				uniqueThingsTwo.add(split[columnTwo]);
			}
			catch(Exception ex)
			{
				System.out.println("ERROR reading line: " + line + " " + ex.toString());;
			}
		}
		reader.close();
		
		System.out.println("File "+fileTwo+" contains " + uniqueThingsTwo.size() + " unique entries in column " + columnTwo);
		
		int counter = 0;
		
		TreeSet<String> inBoth = new TreeSet<String>();
		
		for(String s : uniqueThingsOne)
		{
			if(!uniqueThingsTwo.contains(s))
				counter++;
			else
				inBoth.add(s);
		}
		
		System.out.println("Number of values in file " + fileOne + " not contained in " + fileTwo + ": " + counter);
		
		counter = 0;
		
		for(String s : uniqueThingsTwo)
		{
			if(!uniqueThingsOne.contains(s))
				counter++;
			else
				inBoth.add(s);
		}
		
		System.out.println("Number of values in file " + fileTwo + " not contained in " + fileOne + ": " + counter);
		
		System.out.println("Values contained in both files: " + inBoth.size());
	}
	
	public static void joinTwoFiles(String fileOne, String fileTwo, int idColumnFileOne, int idColumnFileTwo, boolean inBothOnly, String splitter) throws Exception
	{
		HashMap<String, String> idToValueFileOne = new HashMap<String, String>();
		HashMap<String, String> idToValueFileTwo = new HashMap<String, String>();
		
		TreeSet<String> uniqueIDs = new TreeSet<String>();
		
		BufferedReader reader = new BufferedReader(new FileReader(fileOne));
		String line = null;
		
		while((line = reader.readLine()) !=null)
		{
			try
			{
				String[] split = line.split(splitter);
				
				StringBuffer b = new StringBuffer();
				
				for(int i=0; i<split.length; i++)
				{
					if(i == idColumnFileOne)
						continue;
					
					b.append(split[i] + "\t");
				}
				
				idToValueFileOne.put(split[idColumnFileOne].trim(), b.toString());
				uniqueIDs.add(split[idColumnFileOne].trim());
			}
			catch(Exception ex)
			{
				//System.out.println("Problem with line: " + line);
			}
		}
		reader.close();
		
		reader = new BufferedReader(new FileReader(fileTwo));
		line = null;
		
		while((line = reader.readLine()) !=null)
		{
			try
			{
				String[] split = line.split(splitter);
				
				StringBuffer b = new StringBuffer();
				
				for(int i=0; i<split.length; i++)
				{
					if(i == idColumnFileTwo)
						continue;
					
					b.append(split[i] + "\t");
				}
				
				idToValueFileTwo.put(split[idColumnFileTwo].trim(), b.toString());
				uniqueIDs.add(split[idColumnFileTwo].trim());
			}
			catch(Exception ex)
			{
				//System.out.println("Problem with line: " + line);
			}
		}
		reader.close();
		
		for(String id : uniqueIDs)
		{
			String valueOne = "unknown";
			String valueTwo = "unknown";
			
			// test if values are integer values only => default value is 0
			try
			{
				if(idToValueFileTwo.containsKey(id))
					Double.parseDouble(idToValueFileTwo.get(id));
				
				if(idToValueFileOne.containsKey(id))
					Double.parseDouble(idToValueFileOne.get(id));
				
				valueOne = "0";
				valueTwo = "0";
			}
			catch(NumberFormatException ex){}
			
			boolean notContainedInOne = false;
			
			if(idToValueFileOne.containsKey(id))
			{
				valueOne = idToValueFileOne.get(id);
			}
			else
			{
				notContainedInOne = true;
			}
			
			if(idToValueFileTwo.containsKey(id))
			{
				valueTwo = idToValueFileTwo.get(id);
			}
			else
			{
				notContainedInOne = true;
			}
			
			if(!inBothOnly)
				System.out.println(id + "\t" + valueOne + "\t" + valueTwo);
			else if(inBothOnly && !notContainedInOne)
				System.out.println(id + "\t" + valueOne + "\t" + valueTwo);
				
		}
		
	}
	
	public static void leftJoinTwoFiles(String fileOne, String fileTwo, int idColumnFileOne, int idColumnFileTwo, String splitter) throws Exception
	{
		HashMap<String, String> idToValueFileOne = new HashMap<String, String>();
		HashMap<String, String> idToValueFileTwo = new HashMap<String, String>();
		
		TreeSet<String> uniqueIDs = new TreeSet<String>();
		
		BufferedReader reader = new BufferedReader(new FileReader(fileOne));
		String line = null;
		
		while((line = reader.readLine()) !=null)
		{
			try
			{
				String[] split = line.split(splitter);
				
				StringBuffer b = new StringBuffer();
				
				for(int i=0; i<split.length; i++)
				{
					if(i == idColumnFileOne)
						continue;
					
					b.append(split[i] + "\t");
				}
				
				idToValueFileOne.put(split[idColumnFileOne].trim(), b.toString());
				uniqueIDs.add(split[idColumnFileOne].trim());
			}
			catch(Exception ex)
			{
				//System.out.println("Problem with line: " + line);
			}
		}
		reader.close();
		
		reader = new BufferedReader(new FileReader(fileTwo));
		line = null;
		
		while((line = reader.readLine()) !=null)
		{
			try
			{
				String[] split = line.split(splitter);
				
				StringBuffer b = new StringBuffer();
				
				for(int i=0; i<split.length; i++)
				{
					if(i == idColumnFileTwo)
						continue;
					
					b.append(split[i] + "\t");
				}
				
				idToValueFileTwo.put(split[idColumnFileTwo].trim(), b.toString());
			}
			catch(Exception ex)
			{
				//System.out.println("Problem with line: " + line);
			}
		}
		reader.close();
		
		// only ids in first file are used => left join
		for(String id : uniqueIDs)
		{
			String valueOne = "unknown";
			String valueTwo = "unknown";
			
			// test if values are integer values only => default value is 0
			try
			{
				if(idToValueFileTwo.containsKey(id))
					Double.parseDouble(idToValueFileTwo.get(id));
				
				if(idToValueFileOne.containsKey(id))
					Double.parseDouble(idToValueFileOne.get(id));
				
				valueOne = "0";
				valueTwo = "0";
			}
			catch(NumberFormatException ex){}
			
			if(idToValueFileOne.containsKey(id))
			{
				valueOne = idToValueFileOne.get(id);
			}
			
			if(idToValueFileTwo.containsKey(id))
			{
				valueTwo = idToValueFileTwo.get(id);
			}
			
			System.out.println(id + "\t" + valueOne + "\t" + valueTwo);
		}
		
	}
	
	public static void sumInformationForIDs(String file, int idColumn, String summarizedSplitter, String splitter) throws Exception
	{
		HashMap<String, HashMap<Integer, TreeSet<String>>> informationForID = new HashMap<String, HashMap<Integer,TreeSet<String>>>();
		
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line = null;
		
		int maxCols = 0;
		
		while((line = reader.readLine()) !=null)
		{
			String[] split = line.split(splitter);
			
			if(split.length > maxCols)
				maxCols = split.length;
			
			String id = split[idColumn];
			
			for(int i=0; i<split.length; i++)
			{
				if(i == idColumn)
					continue;
				
				String information = split[i];
				
				if(informationForID.containsKey(id) && informationForID.get(id).containsKey(i))
				{
					informationForID.get(id).get(i).add(information);
				}
				else if(informationForID.containsKey(id) && !informationForID.get(id).containsKey(i))
				{
					TreeSet<String> s = new TreeSet<String>();
					s.add(information);
					informationForID.get(id).put(i, s);
				}
				else
				{
					HashMap<Integer, TreeSet<String>> h = new HashMap<Integer, TreeSet<String>>();
					TreeSet<String> s = new TreeSet<String>();
					s.add(information);
					h.put(i, s);
					informationForID.put(id, h);		
				}
			}
		}
		reader.close();
		
		for(String id : informationForID.keySet())
		{
			System.out.print(id + "\t");
			
			for(int i=0; i<maxCols; i++)
			{
				if(i == idColumn)
					continue;
				
				TreeSet<String> s = informationForID.get(id).get(i);
				
				if(s != null)
				{
					boolean first = true;
					
					for(String info : s)
					{
						if(first)
						{
							System.out.print(info);
							first = false;
						}
						else
						{
							System.out.print(summarizedSplitter + info);
						}
					}
				}
				else
				{
					System.out.print("\t");
				}
				
				System.out.print("\t");
			}
			
			System.out.println();
		}
	}
	
	public static void joinTwoFilesOneToN(String fileOne, String fileTwo, int idColumnFileOne, int idColumnFileTwo, boolean inBothOnly, String splitter) throws Exception
	{
		HashMap<String, String> idToValueFileOne = new HashMap<String, String>();
		
		BufferedReader reader = new BufferedReader(new FileReader(fileOne));
		String line = null;
		
		while((line = reader.readLine()) !=null)
		{
			try
			{
				String[] split = line.split(splitter);
				
				StringBuffer b = new StringBuffer();
				
				for(int i=0; i<split.length; i++)
				{
					if(i == idColumnFileOne)
						continue;
					
					b.append(split[i] + "\t");
				}
				
				idToValueFileOne.put(split[idColumnFileOne], b.toString());
			}
			catch(Exception ex)
			{
				//System.out.println("Problem with line: " + line);
			}
		}
		reader.close();
		
		reader = new BufferedReader(new FileReader(fileTwo));
		line = null;
		
		while((line = reader.readLine()) !=null)
		{
			try
			{
				String[] split = line.split(splitter);
				
				StringBuffer b = new StringBuffer();
				
				for(int i=0; i<split.length; i++)
				{
					if(i == idColumnFileTwo)
						continue;
					
					b.append(split[i] + "\t");
				}
				
				String id = split[idColumnFileTwo];
				
				String valueOne = "unknown";
				String valueTwo = "unknown";
				
				// test if values are integer values only => default value is 0
				try
				{
					Double.parseDouble(b.toString());
					
					if(idToValueFileOne.containsKey(id))
						Double.parseDouble(idToValueFileOne.get(id));
					
					valueOne = "0";
					valueTwo = "0";
				}
				catch(NumberFormatException ex){}
				
				boolean notContainedInOne = false;
				
				if(idToValueFileOne.containsKey(id))
				{
					valueOne = idToValueFileOne.get(id);
				}
				else
				{
					notContainedInOne = true;
				}
				
				valueTwo = b.toString();
				
				if(!inBothOnly)
					System.out.println(id + "\t" + valueOne + "\t" + valueTwo);
				else if(inBothOnly && !notContainedInOne)
					System.out.println(id + "\t" + valueOne + "\t" + valueTwo);
					
			}
			catch(Exception ex)
			{
				//System.out.println("Problem with line: " + line);
			}
		}
		reader.close();
	}
	
	public static void joinMultipleFiles(boolean inAllOnly, String[] args) throws Exception
	{
		HashMap<String, HashMap<String, String>> fileToIDToValueInFile = new HashMap<String, HashMap<String,String>>();
		
		TreeSet<String> uniqueIDs = new TreeSet<String>();
		
		boolean numericColumnsOnly = true;
		
		for(int i=2; i<args.length; i+=2)
		{
			String file = args[i];
			int column = Integer.parseInt(args[i+1]);
			
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String line = null;
			
			HashMap<String, String> map = new HashMap<String, String>();
			
			while((line = reader.readLine()) !=null)
			{
				String[] split = line.split("\\s+");
				
				StringBuffer b = new StringBuffer();
				
				boolean first = true;
				
				for(int j=0; j<split.length; j++)
				{
					if(j == column)
						continue;
					
					if(!first)
						b.append("|");
					
					if(numericColumnsOnly)
					{
						try
						{
							Double.parseDouble(split[j]);
						}
						catch(Exception ex)
						{
							numericColumnsOnly = false;
						}
					}
					
					b.append(split[j]);
					
					first = false;
				}
				
				map.put(split[column], b.toString());
				uniqueIDs.add(split[column]);
			}
			
			fileToIDToValueInFile.put(file, map);
			
			reader.close();
		}
		
		// print header
		System.out.print("ID\t");
		for(String file : fileToIDToValueInFile.keySet())
			System.out.print(file + "\t");
		System.out.println();
		
		for(String id : uniqueIDs)
		{
			boolean availableInAll = true;
			
			// test if id is available in all
			for(HashMap<String, String> map : fileToIDToValueInFile.values())
			{
				if(!map.containsKey(id))
				{
					availableInAll = false;
					break;
				}
			}
			
			if(inAllOnly && availableInAll)
			{
				System.out.print(id + "\t");
				
				for(String file : fileToIDToValueInFile.keySet())
				{
					if(fileToIDToValueInFile.get(file).containsKey(id))
					{
						System.out.print(fileToIDToValueInFile.get(file).get(id) +"\t");
					}
					else
					{
						if(numericColumnsOnly)
							System.out.print("0\t");
						else
							System.out.print("unknown\t");
					}
				}
				
				System.out.println();
			}
			else if(!inAllOnly)
			{
				System.out.print(id + "\t");
				
				for(String file : fileToIDToValueInFile.keySet())
				{
					if(fileToIDToValueInFile.get(file).containsKey(id))
					{
						System.out.print(fileToIDToValueInFile.get(file).get(id) +"\t");
					}
					else
					{
						if(numericColumnsOnly)
							System.out.print("0\t");
						else
							System.out.print("unknown\t");
					}
				}
				
				System.out.println();
			}
		}
	}
	
	public static void joinMultipleFilesWithExtractedColumns(boolean inAllOnly, String[] args) throws Exception
	{
		HashMap<String, HashMap<String, String>> fileToIDToValueInFile = new HashMap<String, HashMap<String,String>>();
		
		TreeSet<String> uniqueIDs = new TreeSet<String>();
		
		boolean numericColumnsOnly = true;
		
		for(int i=2; i<args.length; i+=3)
		{
			String file = args[i];
			String[] idColumnsString = args[i+1].split(",");
			Vector<Integer> idColumns = new Vector<Integer>();
		
			for(String n : idColumnsString)
				idColumns.add(Integer.parseInt(n));
			
			String[] extractColumnsString = args[i+2].split(",");
			Vector<Integer> extractColumns = new Vector<Integer>();
			
			for(String n : extractColumnsString)
				extractColumns.add(Integer.parseInt(n));

			
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String line = null;
			
			HashMap<String, String> map = new HashMap<String, String>();
			
			while((line = reader.readLine()) !=null)
			{
				try
				{
					String[] split = line.split("\\s+");
					
					StringBuffer b = new StringBuffer();
					
					boolean first = true;
					
					for(int exc : extractColumns)
					{
						if(!first)
							b.append("|");
						
						if(numericColumnsOnly)
						{
							try
							{
								Double.parseDouble(split[exc]);
							}
							catch(Exception ex)
							{
								numericColumnsOnly = false;
							}
						}
						
						b.append(split[exc]);
						
						first = false;
					}
					
					StringBuffer combinedID = new StringBuffer();
					
					for(int idc : idColumns)
						combinedID.append(split[idc] + "|");
					
					map.put(combinedID.toString(), b.toString());
					uniqueIDs.add(combinedID.toString());
				}
				catch(Exception ex)
				{
					System.err.println("Problem parsing line (ignored...): " + line);
				}
			}
			
			fileToIDToValueInFile.put(file, map);
			
			reader.close();
		}
		
		// print header
		System.out.print("ID\t");
		for(String file : fileToIDToValueInFile.keySet())
			System.out.print(file + "\t");
		System.out.println();
		
		for(String id : uniqueIDs)
		{
			boolean availableInAll = true;
			
			// test if id is available in all
			for(HashMap<String, String> map : fileToIDToValueInFile.values())
			{
				if(!map.containsKey(id))
				{
					availableInAll = false;
					break;
				}
			}
			
			if(inAllOnly && availableInAll)
			{
				System.out.print(id + "\t");
				
				for(String file : fileToIDToValueInFile.keySet())
				{
					if(fileToIDToValueInFile.get(file).containsKey(id))
					{
						System.out.print(fileToIDToValueInFile.get(file).get(id) +"\t");
					}
					else
					{
						if(numericColumnsOnly)
							System.out.print("0\t");
						else
							System.out.print("unknown\t");
					}
				}
				
				System.out.println();
			}
			else if(!inAllOnly)
			{
				System.out.print(id + "\t");
				
				for(String file : fileToIDToValueInFile.keySet())
				{
					if(fileToIDToValueInFile.get(file).containsKey(id))
					{
						System.out.print(fileToIDToValueInFile.get(file).get(id) +"\t");
					}
					else
					{
						if(numericColumnsOnly)
							System.out.print("0\t");
						else
							System.out.print("unknown\t");
					}
				}
				
				System.out.println();
			}
		}
	}
	
	public static void joinMultipleFilesInDirectoryWithExtractedColumns(boolean inAllOnly, String idColumnString, String extractColumnsRawString, 
		String directory, String fileExtension, String splitArgument) throws Exception
	{
		HashMap<String, HashMap<String, String>> fileToIDToValueInFile = new HashMap<String, HashMap<String,String>>();
		
		TreeSet<String> uniqueIDs = new TreeSet<String>();
		
		boolean numericColumnsOnly = true;
		
		File[] files = (new File(directory)).listFiles();
		
		for(File f : files)
		{
			if(!f.getAbsolutePath().endsWith(fileExtension))
				continue;
			
			int idColumn = Integer.parseInt(idColumnString);
			String[] extractColumnsString = extractColumnsRawString.split(",");
			
			Vector<Integer> extractColumns = new Vector<Integer>();
			
			for(String n : extractColumnsString)
			{
				extractColumns.add(Integer.parseInt(n));
			}
			
			BufferedReader reader = new BufferedReader(new FileReader(f));
			String line = null;
			
			HashMap<String, String> map = new HashMap<String, String>();
			
			while((line = reader.readLine()) !=null)
			{
				try
				{
					String[] split = line.split(splitArgument);
					
					StringBuffer b = new StringBuffer();
					
					boolean first = true;
					
					for(int exc : extractColumns)
					{
						if(!first)
							b.append("|");
						
						if(numericColumnsOnly)
						{
							try
							{
								Double.parseDouble(split[exc]);
							}
							catch(Exception ex)
							{
								numericColumnsOnly = false;
							}
						}
						
						b.append(split[exc]);
						
						first = false;
					}
					
					map.put(split[idColumn], b.toString());
					uniqueIDs.add(split[idColumn]);
				}
				catch(Exception ex)
				{
					System.out.println("Could not parse line: " + line);
				}
			}
			
			fileToIDToValueInFile.put(f.getName(), map);
			
			reader.close();
		}
		
		// print header
		System.out.print("ID\t");
		for(String file : fileToIDToValueInFile.keySet())
			System.out.print(file + "\t");
		System.out.println();
		
		for(String id : uniqueIDs)
		{
			boolean availableInAll = true;
			
			// test if id is available in all
			for(HashMap<String, String> map : fileToIDToValueInFile.values())
			{
				if(!map.containsKey(id))
				{
					availableInAll = false;
					break;
				}
			}
			
			if(inAllOnly && availableInAll)
			{
				System.out.print(id + "\t");
				
				for(String file : fileToIDToValueInFile.keySet())
				{
					if(fileToIDToValueInFile.get(file).containsKey(id))
					{
						String value = fileToIDToValueInFile.get(file).get(id).trim();
						
						if(value.equals(id))
							System.out.print("known\t");
						else
							System.out.print(fileToIDToValueInFile.get(file).get(id) +"\t");
					}
					else
					{
						if(numericColumnsOnly)
							System.out.print("0\t");
						else
							System.out.print("unknown\t");
					}
				}
				
				System.out.println();
			}
			else if(!inAllOnly)
			{
				System.out.print(id + "\t");
				
				for(String file : fileToIDToValueInFile.keySet())
				{
					if(fileToIDToValueInFile.get(file).containsKey(id))
					{
						String value = fileToIDToValueInFile.get(file).get(id).trim();
						
						if(value.equals(id))
							System.out.print("known\t");
						else
							System.out.print(fileToIDToValueInFile.get(file).get(id) +"\t");
					}
					else
					{
						if(numericColumnsOnly)
							System.out.print("0\t");
						else
							System.out.print("unknown\t");
					}
				}
				
				System.out.println();
			}
		}
	}
	
	public static void extractColumns(String file, String[] cols, String splitter, int offset, boolean outputFasta) throws Exception
	{
		int[] columns = new int[cols.length-offset];
		
		int maxCol = -1;
		
		for(int i=offset; i<cols.length; i++)
		{
			columns[i-offset] = Integer.parseInt(cols[i]);
			
			if(Integer.parseInt(cols[i]) > maxCol)
				maxCol = Integer.parseInt(cols[i]);
		}
		
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line = null;
		
		int counter = 0;
		
		while((line = reader.readLine()) !=null)
		{
			try
			{
				String[] split = line.split(splitter);
				
				if(outputFasta)
					System.out.println("> " + counter);
				
				for(int i : columns)
					System.out.print(split[i] + "\t");
				
				counter++;
				
				System.out.println();
			}
			catch(Exception ex)
			{}
		}
		reader.close();
	}
	
	public static void extractColumnsByColumnHeaderNames(String file, String fileWithColumnHeaders, String splitter) throws Exception
	{
		TreeSet<Integer> columns = new TreeSet<Integer>();
		Vector<String> columnHeaderNames = new Vector<String>();
		
		BufferedReader reader = new BufferedReader(new FileReader(fileWithColumnHeaders));
		String line = null;
		
		while((line = reader.readLine()) != null)
		{
			columnHeaderNames.add(line.trim());
		}
		reader.close();
		
		reader = new BufferedReader(new FileReader(file));
		line = null;
		
		boolean headerRead = false;
		
		while((line = reader.readLine()) !=null)
		{
			try
			{
				String[] split = line.split(splitter);
				
				if(!headerRead)
				{
					int columnCounter = 0;
					
					StringBuffer b = new StringBuffer();
					
					for(String s : split)
					{
						for(String header : columnHeaderNames)
						{
							if(s.contains(header))
							{
								columns.add(columnCounter);
								
								b.append(header +"\t");
							}
						}
						
						columnCounter++;
					}
					
					System.out.println(b.toString());
					
					headerRead = true;
					
					
				}
				else
				{				
					for(int i : columns)
						System.out.print(split[i] + "\t");
					
					System.out.println();
				}
			}
			catch(Exception ex)
			{}
		}
		reader.close();
	}
	
	public static void addIDColumn(String file, String[] cols, String splitter, int offset) throws Exception
	{
		int[] columns = new int[cols.length-offset];
		
		int maxCol = -1;
		
		for(int i=offset; i<cols.length; i++)
		{
			columns[i-offset] = Integer.parseInt(cols[i]);
			
			if(Integer.parseInt(cols[i]) > maxCol)
				maxCol = Integer.parseInt(cols[i]);
		}
		
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line = null;
		
		while((line = reader.readLine()) !=null)
		{
			try
			{
				String[] split = line.split(splitter);
				
				TreeSet<String> idComponents = new TreeSet<String>();
				
				for(int i : columns)
					idComponents.add(split[i]);
				
				StringBuffer b = new StringBuffer();
				
				boolean first = true;
				
				for(String str : idComponents)
				{
					if(!first)
						b.append("-");
					
					b.append(str );
					
					first = false;
				}
				
				System.out.println(b.toString() + "\t" + line);
			}
			catch(Exception ex)
			{}
		}
		reader.close();
	}
	
	public static void mergeColumns(String file, int columnOne, int columnTwo, String mergeSeparator, String splitter) throws Exception
	{
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line = null;
		
		while((line = reader.readLine()) !=null)
		{
			String[] split = line.split(splitter);
			
			try
			{
				StringBuffer b = new StringBuffer();
				
				b.append(split[columnOne] + mergeSeparator + split[columnTwo] + "\t");
				
				for(int i=0; i<split.length; i++)
				{
					b.append(split[i] + "\t");
				}
				
				System.out.println(b.toString());
			}
			catch(Exception ex)
			{
				// do nothing
			}
		}
		reader.close();
	}
	
	public static void splitColumns(String file, int columnOne, String splitSeparator, String splitter) throws Exception
	{
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line = null;
		
		while((line = reader.readLine()) !=null)
		{
			String[] split = line.split(splitter);
			
			StringBuffer b = new StringBuffer();
			
			for(int i=0; i<split.length; i++)
			{
				if(i != columnOne)
				{
					b.append(split[i] + "\t");
				}
				else
				{
					String[] split2 = split[i].split(splitSeparator);
					
					for(String s : split2)
						b.append(s + "\t");
				}
				
			}
			
			System.out.println(b.toString());
		}
		reader.close();
	}
	
	public static void append(String file, int column, String addOnToColumn, String splitter) throws Exception
	{
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line = null;
		
		while((line = reader.readLine()) !=null)
		{
			String[] split = line.split(splitter);
			
			StringBuffer b = new StringBuffer();
			
			for(int i=0; i<split.length; i++)
			{
				if(i == column)
				{
					b.append(addOnToColumn + split[i] + "\t");
				}
				else
				{
					b.append(split[i] + "\t");
				}
				
			}
			
			System.out.println(b.toString());
		}
		reader.close();
	}
	
	public static void splitAndSummarize(String file, int idColumn, String idColumnSpliter, int valueColumn, String overallSplitter, String separator) throws Exception
	{
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line = null;
		
		HashMap<String, TreeSet<String>> idToValues = new HashMap<String, TreeSet<String>>();
		
		while((line = reader.readLine()) !=null)
		{
			// split total line
			String[] split = line.split(overallSplitter);
			
			// value
			String value = split[valueColumn].trim();
			
			// all ids
			String[] ids = split[idColumn].split(idColumnSpliter);
			
			for(String id : ids)
			{
				String tid = id.trim();
				
				if(idToValues.containsKey(tid))
				{
					idToValues.get(tid).add(value);
				}
				else
				{
					TreeSet<String> set = new TreeSet<String>();
					set.add(value);
					idToValues.put(tid, set);
				}
			}
		}
		reader.close();
		
		for(String id : idToValues.keySet())
		{
			StringBuffer b = new StringBuffer();
			
			boolean first = true;
			
			for(String value : idToValues.get(id))
			{
				if(!first)
					b.append(separator);
				
				first = false;
				
				b.append(value);
			}
			
			System.out.println(id + "\t" + b.toString());
		}
	}
	
	/**
	 * This little class extracts numeric values from a tab or space separated
	 * file. The values are written to standard out.
	 * 
	 * @param args
	 */
	public static void signumForColumns(String[] args) throws Exception
	{
		TreeSet<Integer> columns = new TreeSet<Integer>();
		
		File file = new File(args[1]);
		
		for(int i=2; i<args.length; i++)
		{
			columns.add(Integer.parseInt(args[i]));
		}
		
		BufferedReader reader = new BufferedReader(new FileReader(file));
		
		String line = null;
		
		while((line = reader.readLine()) != null)
		{
			String[] split = line.split("\\s+");
			
			StringBuffer b = new StringBuffer();
			
			for(int i=0; i<split.length; i++)
			{
				if(!columns.contains(i))
					b.append(split[i] + "\t");
				else
					b.append((Double.parseDouble(split[i])*-1) + "\t");
			}
			
			System.out.println(b.toString());
		}
		reader.close();
	}
	
	/**
	 * This little class extracts numeric values from a tab or space separated
	 * file. The values are written to standard out.
	 * 
	 * @param args
	 */
	public static void log2ForColumns(String[] args) throws Exception
	{
		TreeSet<Integer> columns = new TreeSet<Integer>();
		
		File file = new File(args[1]);
		
		for(int i=2; i<args.length; i++)
		{
			columns.add(Integer.parseInt(args[i]));
		}
		
		BufferedReader reader = new BufferedReader(new FileReader(file));
		
		String line = null;
		
		while((line = reader.readLine()) != null)
		{
			String[] split = line.split("\\s+");
			
			StringBuffer b = new StringBuffer();
			
			for(int i=0; i<split.length; i++)
			{
				try
				{
					if(!columns.contains(i))
						b.append(split[i] + "\t");
					else
						b.append((Math.log((Double.parseDouble(split[i])))/Math.log(2)) + "\t");
				}
				catch(Exception ex)
				{
					b.append(split[i] + "\t");
				}
			}
			
			System.out.println(b.toString());
		}
		reader.close();
	}
	
	public static void signumForColumn(String file, int columnOne) throws Exception
	{
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line = null;
		
		while((line = reader.readLine()) !=null)
		{
			String[] split = line.split("\\s+");
			
			StringBuffer b = new StringBuffer();
			
			for(int i=0; i<split.length; i++)
			{
				if(i != columnOne)
					b.append(split[i] + "\t");
				else
					b.append((Double.parseDouble(split[i])*-1) + "\t");
			}
			
			System.out.println(b.toString());
		}
		reader.close();
	}
	
	public static void sumColumn(String file, int columnOne, String splitter) throws Exception
	{
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line = null;
		
		double sum = 0;
		double counter = 0;
		
		while((line = reader.readLine()) !=null)
		{
			String[] split = line.split(splitter);
			
			for(int i=0; i<split.length; i++)
			{
				if(i == columnOne)
				{
					try
					{
						sum += Double.parseDouble(split[i]);
						counter++;
					}
					catch(Exception ex)
					{
						
					}
				}
			}
		}
		reader.close();
		
		System.out.println("Sum " + sum + ", Average: " + (sum / counter));
	}
      
	
	public static void computeVennDiagramForMFiles(String splitter, boolean output, String[] args) throws Exception
	{
		HashMap<String, Integer> fileToColumn = new HashMap<String, Integer>(); 
		TreeMap<String, String> fileToID = new TreeMap<String, String>();
		
		for(int i=3; i<args.length; i+=3)
		{
			String file = args[i];
			int column = Integer.parseInt(args[i+1]);
			String id = args[i+2];
			
			fileToColumn.put(file, column);
			fileToID.put(file, id);
		}
		
		TreeMap<String, TreeSet<String>> fileToContent = new TreeMap<String, TreeSet<String>>();
		
		TreeSet<String> allIDs = new TreeSet<String>();
		
		for(String file : fileToColumn.keySet())
		{
			int column = fileToColumn.get(file);
			
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String line = null;
			
			TreeSet<String> idsInFile = new TreeSet<String>();
			
			while((line = reader.readLine()) !=null)
			{
				try
				{
					String[] split = line.split(splitter);
				
					allIDs.add(split[column].trim());
					idsInFile.add(split[column].trim());
				}
				catch(Exception ex)
				{
					System.err.println("Could not pare line " + line + " in file " + file);
				}
			}
			
			fileToContent.put(file, idsInFile);
			reader.close();
		}
		
		TreeMap<String, Integer> combinationCounter = new TreeMap<String, Integer>();
		
		for(String id : allIDs)
		{
			StringBuffer b = new StringBuffer();
			
			for(String file : fileToContent.keySet())
			{
				if(fileToContent.get(file).contains(id))
					b.append(fileToID.get(file) + "|");
			}
			
			if(combinationCounter.containsKey(b.toString()))
				combinationCounter.put(b.toString(), combinationCounter.get(b.toString()) + 1);
			else
				combinationCounter.put(b.toString(), 1);
			
			if(output)
				System.out.println(id + "\tGroup:" + b.toString() + ":Group");
		}
		
		for(String key : combinationCounter.keySet())
		{
			System.out.println("Group:" + key + "\t" + combinationCounter.get(key));
		}
		
		System.out.println("All ids: " + allIDs.size());
		
		for(String key : fileToContent.keySet())
			System.out.println(fileToID.get(key) + " contains " + fileToContent.get(key).size() + " ids...");
	}
	
	public static void computeVennDiagramForMFilesCombinedID(String splitter, boolean output, String[] args) throws Exception
	{
		HashMap<String, Vector<Integer>> fileToColumns = new HashMap<String, Vector<Integer>>(); 
		HashMap<String, String> fileToID = new HashMap<String, String>();
		
		for(int i=3; i<args.length; i+=3)
		{
			String file = args[i];
			
			Vector<Integer> columns = new Vector<Integer>();
			
			for(String s : args[i+1].split(","))
				columns.add(Integer.parseInt(s));
			
			String id = args[i+2];
			
			fileToColumns.put(file, columns);
			fileToID.put(file, id);
		}
		
		TreeMap<String, TreeSet<String>> fileToContent = new TreeMap<String, TreeSet<String>>();
		
		TreeSet<String> allIDs = new TreeSet<String>();
		
		for(String file : fileToColumns.keySet())
		{
			Vector<Integer> columns = fileToColumns.get(file);
			
			BufferedReader reader = null;
			
			if(file.endsWith(".gz"))
			{
				System.out.println("File is zipped with gzip. Will use corresponding reader...");
				reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
			}
			else
				reader = new BufferedReader(new FileReader(file));
			
			String line = null;
			
			TreeSet<String> idsInFile = new TreeSet<String>();
			
			while((line = reader.readLine()) !=null)
			{
				try
				{
					String[] split = line.split(splitter);
					
					StringBuffer idBuffer = new StringBuffer();
					
					for(int i : columns)
						idBuffer.append("-" + split[i]);
				
					allIDs.add(idBuffer.toString());
					idsInFile.add(idBuffer.toString());
				}
				catch(Exception ex)
				{
					System.err.println("Could not pare line " + line + " in file " + file);
				}
			}
			reader.close();
			fileToContent.put(file, idsInFile);
		}
		
		TreeMap<String, Integer> combinationCounter = new TreeMap<String, Integer>();
		
		for(String id : allIDs)
		{
			StringBuffer b = new StringBuffer();
			
			for(String file : fileToContent.keySet())
			{
				if(fileToContent.get(file).contains(id))
					b.append(fileToID.get(file) + "|");
			}
			
			if(combinationCounter.containsKey(b.toString()))
				combinationCounter.put(b.toString(), combinationCounter.get(b.toString()) + 1);
			else
				combinationCounter.put(b.toString(), 1);
			
			if(output)
				System.out.println(id + "\tGroup:" + b.toString() + ":Group");
		}
		
		for(String key : combinationCounter.keySet())
		{
			System.out.println("Group:" + key + "\t" + combinationCounter.get(key));
		}
		
		System.out.println("All ids: " + allIDs.size());
		
		for(String key : fileToContent.keySet())
			System.out.println(fileToID.get(key) + " contains " + fileToContent.get(key).size() + " ids...");
	}
	
	public static void dumpMySQL(String statement, String config, boolean fasta) throws Exception
	{
		Properties properties = new Properties();
		properties.load(new FileInputStream(config));
		String host = properties.getProperty("host");
		String user = properties.getProperty("user");
		String db = properties.getProperty("db");
		String passwd = properties.getProperty("password");
		
		// Load the JDBC driver
		String driverName = "org.gjt.mm.mysql.Driver";
		Class.forName(driverName);

		// Create a connection to the database
		String url = "jdbc:mysql://" + host + "/" + db; // a JDBC url
		Connection connection = DriverManager.getConnection(url, user, passwd);
		
		Statement stat = connection.createStatement();
		ResultSet set = stat.executeQuery(statement);
		
		ResultSetMetaData data = set.getMetaData();
		
		while(set.next())
		{
			if(fasta)
			{
				System.out.print(">"+set.getObject(1) +"\n");
				System.out.print(set.getObject(2) +"\n");
				System.out.println();
			}
			else
			{
				for(int i=1; i<=data.getColumnCount(); i++)
				{
					System.out.print(set.getObject(i) +"\t");
				}
				System.out.println();
			}
		}
	}
	
	public static void head(String file, int numberOfLines) throws Exception
	{
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line = null;
		
		int counter = 0;
		
		while((line = reader.readLine()) != null)
		{
			if(counter > numberOfLines)
				break;
			
			System.out.println(line);
			counter++;
		}
		reader.close();
	}
	
	public static void replace(String file, String query, String replacement) throws Exception
	{
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line = null;
				
		while((line = reader.readLine()) != null)
		{
			System.out.println(line.replaceAll(query, replacement));
		}
		reader.close();
	}
	
	public static void replaceValuesInColumn(String file, int column, String columnSplitter, String columnInternalSplitter, String query, String replacement) throws Exception
	{
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line = null;
				
		while((line = reader.readLine()) != null)
		{
			String[] split = line.split(columnSplitter);
			
			String[] columnSplit = split[column].split(columnInternalSplitter);
			
			TreeSet<String> sortedColumnValues = new TreeSet<String>();
			
			StringBuffer b = new StringBuffer();
			
			for(int i=0; i<columnSplit.length; i++)
			{
				String result = null;
				
				if(columnSplit[i].contains(query))
					result = replacement;
				else
					result = columnSplit[i];
				
				sortedColumnValues.add(result + ",");
			}
			
			for(String s : sortedColumnValues)
				b.append(s);
			
			StringBuffer b2 = new StringBuffer();
			
			for(int i=0;i<split.length; i++)
			{
				if(i!= column)
					b2.append(split[i] + "\t");
				else
					b2.append(b.toString() + "\t");
			}
			
			System.out.println(b2.toString());
		}
		reader.close();
	}
	
	public static void replaceMulti(String file, String fromToFile) throws Exception
	{
		HashMap<String, String> fromToMap = new HashMap<String, String>();
		
		BufferedReader reader = new BufferedReader(new FileReader(fromToFile));
		String line = null;
		
		while((line = reader.readLine()) != null)
		{
			String[] split = line.split("\t");
			
			if(split.length >= 2)
				fromToMap.put(split[0], split[1]);
		}
		
		reader.close();
		
		reader = new BufferedReader(new FileReader(file));
		line = null;
				
		while((line = reader.readLine()) != null)
		{
			boolean nothingChanged = true;
			
			for(String f : fromToMap.keySet())
			{
				String newLine = line.replaceAll(f, fromToMap.get(f));
				
				if(!line.equals(newLine))
				{
					System.out.println(newLine);
					nothingChanged = false;
				}
			}
			
			if(nothingChanged)
				System.out.println(line);
		}
	}
	
	public static void tab2FASTA(String file, int sequenceCol, int idCol) throws Exception
	{
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line = null;
				
		while((line = reader.readLine()) != null)
		{
			String[] split = line.split("\t");
			
			System.out.println(">" + split[idCol]);
			System.out.println(split[sequenceCol]);
		}
		reader.close();
	}
	
	public static void maplot(String fileOne, String fileTwo, int idColumnOne, int idColumnTwo, int intensityColumnOne, 
		int intensityColumnTwo, String splitter) throws Exception
	{
		TreeSet<String> allGenes = new TreeSet<String>();
		HashMap<String, Double> log2ValuesDatasetOne = new HashMap<String, Double>();
		HashMap<String, Double> log2ValuesDatasetTwo = new HashMap<String, Double>();
		
		BufferedReader reader = new BufferedReader(new FileReader(fileOne));
		String line = null;
		
		while((line = reader.readLine()) != null)
		{
			String[] split = line.split(splitter);
			
			allGenes.add(split[idColumnOne]);
			log2ValuesDatasetOne.put(split[idColumnOne], Math.log(Double.parseDouble(split[intensityColumnOne]))/Math.log(2));
		}
		reader.close();
		
		reader = new BufferedReader(new FileReader(fileTwo));
		line = null;
		
		while((line = reader.readLine()) != null)
		{
			String[] split = line.split(splitter);
			
			allGenes.add(split[idColumnTwo]);
			log2ValuesDatasetTwo.put(split[idColumnTwo], Math.log(Double.parseDouble(split[intensityColumnTwo]))/Math.log(2));
		}
		reader.close();
		
		System.err.println("All genes: " + allGenes.size());
		System.err.println("Dataset one: " + log2ValuesDatasetOne.size());
		System.err.println("Dataset two: " + log2ValuesDatasetTwo.size());
		
		for(String gene : allGenes)
		{
			if(log2ValuesDatasetOne.containsKey(gene) && log2ValuesDatasetTwo.containsKey(gene))
			{
				double m = log2ValuesDatasetOne.get(gene) - log2ValuesDatasetTwo.get(gene);
				double a = 0.5 * (log2ValuesDatasetOne.get(gene) + log2ValuesDatasetTwo.get(gene));
				
				System.out.println(gene + "\t" + log2ValuesDatasetOne.get(gene) + "\t" + log2ValuesDatasetTwo.get(gene) + "\t" + m + "\t" + a);
			}
		}
	}
	
	public static void extractSequencesFromFASTA(String file, String headerContainsInformaion) throws Exception
	{
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line = null;
		
		boolean printSequence = false;
		
		while((line = reader.readLine()) != null)
		{
			if(line.startsWith(">"))
			{
				if(line.contains(headerContainsInformaion))
				{
					printSequence = true;
					System.out.println(line);
				}
				else
				{
					printSequence = false;
				}
			}
			else
			{
				if(printSequence)
					System.out.println(line);
			}
		}
		reader.close();
	}
	
	public static void computeStandardDeviationsAndMeanForValuesInRows(String file, int idColumn, String splitter) throws Exception
	{
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line = null;
		
		while((line = reader.readLine()) != null)
		{
			String[] split = line.split(splitter);
			
			String id = split[idColumn];
			Vector<Double> values = new Vector<Double>();
			
			for(int i=0; i<split.length; i++)
			{
				if(i==idColumn)
					continue;
				
				try
				{
					values.add(Double.parseDouble(split[i]));
				}
				catch(Exception ex)
				{
					// do nothing, just ignore value
				}
			}
			
			// compute standard deviation and mean value
			double mean = 0;
			
			for(double value : values)
				mean += value;
			
			mean /= values.size();
			
			double standardDeviation = 0;
			
			for(double value : values)
				standardDeviation += Math.pow((value - mean), 2);
			
			standardDeviation /= values.size();
			standardDeviation = Math.sqrt(standardDeviation);
			
			System.out.println(id + "\t" + mean + "\t" + standardDeviation + "\t" + values.size());
			
		}
		reader.close();
	}
	
	public static void main(String[] args) throws Exception
	{
		if(args.length == 0)
		{
			System.out.println("-------------- WELCOME TO BIOKIT UTILS ------------- ");
			System.out.println();
			System.out.println("Extractraction columns:");
			System.out.println("USAGE to extract numeric values: java -cp . biokit/util/Utils -e File {columns}");
			System.out.println("USAGE to extract columns from a file: java -cp . biokit/util/Utils -ec File {columns}");
			System.out.println("USAGE to extract columns from a file with your own split argument: java -cp . biokit/util/Utils -ecs File splitter {columns}");
			System.out.println("USAGE to extract columns from a file, output in fasta format: java -cp . biokit/util/Utils -ecf File splitter {columns}");
			System.out.println("USAGE to extract columns from a file with column names specified in a file: java -cp . biokit/util/Utils -ecl File FileContainingColumnNamesOneIDPerLine splitter");
			System.out.println();
			System.out.println("Column manipulation:");
			System.out.println("USAGE to merge two columns: java -cp . biokit/util/Utils -m File col1 col2 mergeSeparator splitter");
			System.out.println("USAGE to split two columns: java -cp . biokit/util/Utils -um File col splitterForTextInSpecificColumn splitterForAllColumns");
			System.out.println("USAGE append a value to each entry in a column: java -cp . biokit/util/Utils -append file column valueToBeAdded splitter");
			System.out.println("USAGE to replace a query string by a replacement in a file: java -cp . biokit/util/Utils -replace file query target");
			System.out.println("USAGE to replace a set of query strings in a file by new values in a file: java -cp . biokit/util/Utils -replaceMulti file fromToMappingFile");
			System.out.println("USAGE to add an id to a column based file that is composed of multiple columns of the file and sorts components lexicographically (separated by -): java -cp . biokit/util/Utils -id File splitter {columns}");
			System.out.println();
			System.out.println("Extract lines");
			System.out.println("USAGE to extract the first n lines from a file: java -cp . biokit/util/Utils -head File numberOfLines");
			System.out.println("USAGE to extract all lines which contain a specified value (reg. exp.) in the given colum: java -cp . biokit/util/Utils -l File columnMatching matchinValueRegExp splitter");
			System.out.println("USAGE to extract all lines which DO NOT contain a specified value (reg. exp.) in the given colum: java -cp . biokit/util/Utils -ln File columnMatching matchinValueRegExp splitter");
			System.out.println("USAGE to extract all lines for which the numeric value in the specifie column is larger, smaller or equal (user specified) than a user specified value: java -cp . biokit/util/Utils -lseg File columnMatching numericValue splitter smaller(boolean) equal(boolean) greater(boolean)");
			System.out.println("USAGE to extract all lines for which the numeric value in the specifie column is within a user specified interval: java -cp . biokit/util/Utils -linter File columnMatching lowerBound upperBound splitter searchForWithinInterval");
			System.out.println();
			System.out.println("Line filtering");
			System.out.println("USAGE to filter lines which contain an id specified in another file: java -cp . biokit/util/Utils -fl fileToBeFiltered ColumnWhichContainsFilterID fileContainingIDsUsedForFiltering ColumnContainingFilterIDs splitter onlyLinesWhichDoNOTContainFilterID");
			System.out.println();
			System.out.println("Summarization");
			System.out.println("USAGE to summarize data in a file for a specified id column: java -cp . biokit/util/Utils -sumCols file idColumn splitterForMergedValues splitter");
			System.out.println("USAGE to summarize data in a file for an id column containing multiple ids: java -cp . biokit/util/Utils -splitSum file idColumn idSplit valueColumn overallSplitter valueSeparator");
			System.out.println("USAGE to summarize rows to columns and rows: -rowsToColumns file valueColumn headerColumn rowNameColumnsCommaSeparated, splitArgument");
			System.out.println();
			System.out.println("Unique values:");
			System.out.println("USAGE to count unique values in a file: java -cp . biokit/util/Utils -c File column");
			System.out.println("USAGE to count unique values in a file with your own split argument: java -cp . biokit/util/Utils -cs File column split");
			System.out.println("USAGE to count unique values in a file not contained in another file: java -cp . biokit/util/Utils -cn FileOne columnOne FileTwo columnTwo");
			System.out.println();
			System.out.println("Output:");
			System.out.println("USAGE to output unique values in a file: java -cp . biokit/util/Utils -o File column");
			System.out.println("USAGE to output unique values in a file with user specified split argument: java -cp . biokit/util/Utils -os File column splitter completeLineFlag");
			System.out.println("USAGE to make a tab delimited file to a FASTA file: java -cp . biokit/util/Utils -tab2FASTA inFile idCol sequenceCol");
			System.out.println();
			System.out.println("Histograms:");
			System.out.println("USAGE to output a histogram on frequencies of classes in a column: java -cp . biokit/util/Utils -h File column");
			System.out.println("USAGE to output a histogram on frequencies of classes in two columns: java -cp . biokit/util/Utils -h2 File columnOne columnTwo");
			System.out.println("USAGE to output a histogram on frequencies of values in a column: java -cp . biokit/util/Utils -hn File column numberOfBins splitter [maxValue]");
			System.out.println();
			System.out.println("Venn diagrams:");
			System.out.println("USAGE to compute a venn diagram of ids in two files: java -cp . biokit/util/Utils -venn2 fileA fileB columnA columnB labelA labelB title splitter outputBoolean");
			System.out.println("USAGE to compute a venn diagram of ids in three files: java -cp . biokit/util/Utils -venn3 fileA fileB fileC columnA columnB columnC labelA labelB labelC title splitter outputBoolean");
			System.out.println("USAGE to compute overlap of ids in m files: java -cp . biokit/util/Utils -vennm splitter outputBoolean [file column id file2 column2 id2 ...]");
			System.out.println();
			System.out.println("Join files:");
			System.out.println("USAGE to join two files: java -cp . biokit/util/Utils -j FileOne FileTwo inBothOnly [idColumnFileOne idColumnFileTwo]");
			System.out.println("USAGE to left join two files (all ids in first are joined with ids in second file (if available): java -cp . biokit/util/Utils -jl FileOne FileTwo idColumnFileOne idColumnFileTwo splitter");
			System.out.println("USAGE to join two files where the second file may contain n entries which map to one id in file one: java -cp . biokit/util/Utils -j1n FileOne FileTwo inBothOnly [idColumnFileOne idColumnFileTwo splitter]");
			System.out.println("USAGE to join multiple files: java -cp . biokit/util/Utils -jm inAllOnly fileOne columnOne ... fileN columnN");
			System.out.println("USAGE to join multiple files and extract specific columns: java -cp . biokit/util/Utils -jmec inAllOnly fileOne idColumnOne columnsForJoinOneCommaSeparated ... fileN idColumnN columnsForJoinNCommaSeparated");
			System.out.println("USAGE to join multiple files with specified extension in a directory and extract specific columns: java -cp . biokit/util/Utils -jmecd inAllOnly idColumn columnsForJoinCommaSeparated directory fileExtension");
			System.out.println();
			System.out.println("Sorting:");
			System.out.println("USAGE to sort lines in a file with respect to a value in a specified column (lexicographical order): java -cp . biokit/util/Utils -s File columnWithNumericValue splitter");
			System.out.println("USAGE to sort lines in a file with respect to a numeric value in a specified column: java -cp . biokit/util/Utils -sn File columnWithNumericValue splitter");
			System.out.println();
			System.out.println("Math:");
			System.out.println("USAGE to apply signum to one column: java -cp . biokit/util/Utils -signum File col1 [col2...colN]");
			System.out.println("USAGE to apply log2 to a set of columns: java -cp . biokit/util/Utils -log2 File col1 [col2...colN]");
			System.out.println("USAGE to compute MA plot columns: java -cp . biokit/util/Utils -maplot FileOne FileTwo idColumnOne idColumnTwo intensityColumnOne intensityColumnTwo splitter");
			System.out.println("USAGE to compute the sum of values in a column: java -cp . biokit/util/Utils -sum File col splitter");
			System.out.println("USAGE to compute mean and sd for values in a row: java -cp . biokit/util/Utils -meansd File idCol splitter");
			System.out.println();
			System.out.println("MySQL");
			System.out.println("USAGE to retrieve data from mysql: java -cp . biokit/util/Utils -mysql \"statement replace all \" by \\\" \" config fastaOutput(boolean)");
			System.out.println("FASTA");
			System.out.println("USAGE to extract columns from a file, output in fasta format: java -cp . biokit/util/Utils -ecf File splitter {columns}");
			System.out.println("USAGE to extract sequences that fit header information from FASTA file: java -cp . biokit/util/Utils -ef File headerContainsString");
			System.out.println("----------------- HAVE FUN ----------------");
			
			

			System.exit(0);
		}
		
		if(args[0].equals("-e"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to extract numeric values: java -cp . biokit/util/Utils -e File {columns}");
				System.exit(0);
			}
			
			Utils.extractNumericValues(args);
		}
		else if(args[0].equals("-ec"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to extract columns from a file: java -cp . biokit/util/Utils -ec File {columns}");
				System.exit(0);
			}
			
			Utils.extractColumns(args[1], args, "\\s+", 2, false);
		}
		else if(args[0].equals("-ecs"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to extract columns from a file with your own split argument: java -cp . biokit/util/Utils -ecs File splitter {columns}");
				System.exit(0);
			}
			
			Utils.extractColumns(args[1], args, args[2], 3, false);
		}
		else if(args[0].equals("-ecl"))
		{
			if(args.length < 4)
			{
				System.out.println("USAGE to extract columns from a file with column names specified in a file: java -cp . biokit/util/Utils -ecl File FileContainingColumnNamesOneIDPerLine splitter {columns}");
				System.exit(0);
			}
			
			Utils.extractColumnsByColumnHeaderNames(args[1], args[2], args[3]);
		}
		else if(args[0].equals("-ecf"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to extract columns from a file, output in fasta format: java -cp . biokit/util/Utils -ecf File splitter {columns}");
				System.exit(0);
			}
			
			Utils.extractColumns(args[1], args, "\\s+", 2, true);
		}
		else if(args[0].equals("-id"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to add an id to a column based file that is composed of multiple columns of the file and sorts components lexicographically: java -cp . biokit/util/Utils -id File splitter {columns}");
				System.exit(0);
			}
			
			Utils.addIDColumn(args[1], args, args[2], 3);
		}
		else if(args[0].equals("-head"))
		{
			if(args.length < 3)
			{
				System.out.println("USAGE to extract the first n lines from a file: java -cp . biokit/util/Utils -head File numberOfLines");
				System.exit(0);
			}
			
			Utils.head(args[1], Integer.parseInt(args[2]));
		}
		else if(args[0].equals("-c"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to count unique values in a file: java -cp . biokit/util/Utils -c File column");
				System.exit(0);
			}
			
			Utils.uniqueColumnCounter(args[1], Integer.parseInt(args[2]), "\\s+");
		}
		else if(args[0].equals("-cs"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to count unique values in a file with your own split argument: java -cp . biokit/util/Utils -cs File column split");
				System.exit(0);
			}
			
			Utils.uniqueColumnCounter(args[1], Integer.parseInt(args[2]), args[3]);
		}
		else if(args[0].equals("-cn"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to count unique values in a file not contained in another file: java -cp . biokit/util/Utils -cn FileOne columnOne FileTwo columnTwo");
				System.exit(0);
			}
			
			Utils.nonUniqueColumnCounter(args[1], args[2], Integer.parseInt(args[3]), Integer.parseInt(args[4]));
		}
		else if(args[0].equals("-copy"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to copy files in a list to a new directory: java -cp . biokit/util/Utils -copy ListOfFiles outputDirectory fileIDColumn [additionalSubdirectoryColumn]");
				System.exit(0);
			}
			
			BufferedReader reader = new BufferedReader(new FileReader(args[1]));
			String line = null;
			
			int fileColumn = Integer.parseInt(args[3]);
			int additionalSubdirectoryColumn = -1;
			
			if(args.length == 5)
				additionalSubdirectoryColumn = Integer.parseInt(args[4]);
			
			while((line = reader.readLine()) != null)
			{
				String[] split = line.split("\\s+");
				
				File f = new File(split[fileColumn]);
				
				String subdirectory = "";
				
				if(additionalSubdirectoryColumn != -1)
				{
					subdirectory = split[additionalSubdirectoryColumn] + File.separatorChar;
					
					File subdir = new File(args[2] + File.separatorChar + subdirectory);
					
					if(!subdir.exists())
					{
						System.out.println("Create subdirectory " + subdir.getAbsolutePath());
						subdir.mkdir();
					}
				}
				
				if(f.exists())
				{
					File toFile = new File(args[2] + File.separatorChar + subdirectory + f.getName());
					
					System.out.println("Copy: " + f.getAbsolutePath() + " to " + toFile.getAbsolutePath());
					
					copyPerStream(f, toFile);
				}
			}
			reader.close();
		}
		else if(args[0].equals("-m"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to merge two columns: java -cp . biokit/util/Utils -m File col1 col2 mergeSeparator splitter");
				System.exit(0);
			}
			
			Utils.mergeColumns(args[1], Integer.parseInt(args[2]), Integer.parseInt(args[3]), args[4], args[5]);
		}
		else if(args[0].equals("-um"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to split two columns: java -cp . biokit/util/Utils -um File col splitterForTextInSpecificColumn splitterForAllColumns");
				System.exit(0);
			}
			
			Utils.splitColumns(args[1], Integer.parseInt(args[2]), args[3], args[4]);
		}
		else if(args[0].equals("-o"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to output unique values in a file: java -cp . biokit/util/Utils -o File column completeLineFlag");
				System.exit(0);
			}
			
			Utils.uniqueColumnCounterWithOutput(args[1], Integer.parseInt(args[2]), "\\s+", Boolean.parseBoolean(args[3]));
		}
		else if(args[0].equals("-os"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to output unique values in a file with user specified split argument: java -cp . biokit/util/Utils -os File column splitter completeLineFlag");
				System.exit(0);
			}
			
			Utils.uniqueColumnCounterWithOutput(args[1], Integer.parseInt(args[2]), args[3], Boolean.parseBoolean(args[4]));
		}
		else if(args[0].equals("-h"))
		{
			if(args.length < 4)
			{
				System.out.println("USAGE to output a histogram on frequencies of classes in a column: java -cp . biokit/util/Utils -h File column splitter [columnSplitter columnValuePositionAfterSplit]");
				System.exit(0);
			}
			
			if(args.length==4)
				Utils.printHistogram(args[1], Integer.parseInt(args[2]), args[3], null, -1);
			else if(args.length==6)
				Utils.printHistogram(args[1], Integer.parseInt(args[2]), args[3], args[4], Integer.parseInt(args[5]));
		}
		else if(args[0].equals("-h2"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to output a histogram on frequencies of classes in two columns: java -cp . biokit/util/Utils -h2 File columnOne columnTwo");
				System.exit(0);
			}
			
			Utils.printHistogramTwoColumns(args[1], Integer.parseInt(args[2]), Integer.parseInt(args[3]));
		}
		else if(args[0].equals("-hn"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to output a histogram on frequencies of values in a column: java -cp . biokit/util/Utils -hn File column numberOfBins splitter [maxValue]");
				System.exit(0);
			}
			
			if(args.length == 6)
				Utils.printHistogramOfNumericValues(args[1], Integer.parseInt(args[2]), Integer.parseInt(args[3]), args[4], Integer.parseInt(args[5]));
			else
				Utils.printHistogramOfNumericValues(args[1], Integer.parseInt(args[2]), Integer.parseInt(args[3]), args[4], -1);
		}
		else if(args[0].equals("-j") && args.length<6)
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to join two files: java -cp . biokit/util/Utils -j FileOne FileTwo inBothOnly [idColumnFileOne idColumnFileTwo splitter]");
				System.exit(0);
			}
			
			String splitter = "\\s+";
			
			Utils.joinTwoFiles(args[1], args[2], 0, 0, Boolean.parseBoolean(args[3]), splitter);
		}
		else if(args[0].equals("-j"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to join two files: java -cp . biokit/util/Utils -j FileOne FileTwo inBothOnly [idColumnFileOne idColumnFileTwo splitter]");
				System.exit(0);
			}
			
			String splitter = "\\s+";
			
			if(args.length == 7)
				splitter = args[6];
			
			Utils.joinTwoFiles(args[1], args[2], Integer.parseInt(args[4]), Integer.parseInt(args[5]), Boolean.parseBoolean(args[3]), splitter);
		}
		else if(args[0].equals("-jl"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to left join two files (all ids in first are joined with ids in second file (if available): java -cp . biokit/util/Utils -jl FileOne FileTwo idColumnFileOne idColumnFileTwo splitter");
				System.exit(0);
			}
			
			Utils.leftJoinTwoFiles(args[1], args[2], Integer.parseInt(args[3]), Integer.parseInt(args[4]), args[5]);
		}
		else if(args[0].equals("-j1n"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to join two files where the second file may contain n entries which map to one id in file one: java -cp . biokit/util/Utils -j1n FileOne FileTwo inBothOnly [idColumnFileOne idColumnFileTwo splitter]");
				System.exit(0);
			}
			
			String splitter = "\\s+";
			
			if(args.length == 7)
				splitter = args[6];
				
			
			Utils.joinTwoFilesOneToN(args[1], args[2], Integer.parseInt(args[4]), Integer.parseInt(args[5]), Boolean.parseBoolean(args[3]), splitter);
		}
		else if(args[0].equals("-jm"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to join multiple files: java -cp . biokit/util/Utils -jm inAllOnly fileOne columnOne ... fileN columnN");
				System.exit(0);
			}
			
			Utils.joinMultipleFiles(Boolean.parseBoolean(args[1]), args);
		}
		else if(args[0].equals("-jmec"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to join multiple files and extract specific columns: java -cp . biokit/util/Utils -jmec inAllOnly fileOne idColumnsOneCommaSeparatedForCombinedIDs columnsForJoinOneCommaSeparated ... fileN idColumnsNCommaSeparatedForCombinedIDs columnsForJoinNCommaSeparated");
				
				System.exit(0);
			}
			
			Utils.joinMultipleFilesWithExtractedColumns(Boolean.parseBoolean(args[1]), args);
		}
		else if(args[0].equals("-jmecd"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to join multiple files with specified extension in a directory and extract specific columns: java -cp . biokit/util/Utils -jmecd inAllOnly idColumn columnsForJoinCommaSeparated directory fileExtension splitArgument");
				
				System.exit(0);
			}
			
			Utils.joinMultipleFilesInDirectoryWithExtractedColumns(Boolean.parseBoolean(args[1]), args[2], args[3], args[4], args[5], args[6]);
		}
		else if(args[0].equals("-sn"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to sort lines in a file with respect to a numeric value in a specified column: java -cp . biokit/util/Utils -sn File columnWithNumericValue splitter");
				System.exit(0);
			}
			
			Utils.sortLinesWithRespectToColumn(args[1], Integer.parseInt(args[2]), args[3]);
		}
		else if(args[0].equals("-s"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to sort lines in a file with respect to a specified column: java -cp . biokit/util/Utils -s File columnWithValue splitter");
				System.exit(0);
			}
			
			Utils.sortLinesWithRespectToColumnNonNumeric(args[1], Integer.parseInt(args[2]), args[3]);
		}
		else if(args[0].equals("-l"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to extract all lines which contain a specified value (reg. exp.) in the given colum: java -cp . biokit/util/Utils -l File columnMatching matchinValueRegExp splitter");
				System.exit(0);
			}
			
			Utils.extractLinesWhereValueInColumnMatchesCriterion(args[1], Integer.parseInt(args[2]), args[3], args[4]);
		}
		else if(args[0].equals("-ln"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to extract all lines which DO NOT contain a specified value (reg. exp.) in the given colum: java -cp . biokit/util/Utils -ln File columnMatching matchinValueRegExp splitter");
				System.exit(0);
			}
			
			Utils.extractLinesWhereValueInColumndDoesNotMatchCriterion(args[1], Integer.parseInt(args[2]), args[3], args[4]);
		}
		else if(args[0].equals("-lseg"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to extract all lines for which the numeric value in the specifie column is larger, smaller or equal (user specified) than a user specified value: java -cp . biokit/util/Utils -lseg File columnMatching numericValue splitter smaller(boolean) equal(boolean) greater(boolean)");
				System.exit(0);
			}
			
			Utils.extractLinesWhereNumericValueInColumnIsSmallerEqualOrGreaterThanUserValue(args[1], Integer.parseInt(args[2]), Double.parseDouble(args[3]), 
				args[4], Boolean.parseBoolean(args[5]), Boolean.parseBoolean(args[6]), Boolean.parseBoolean(args[7]));
		}
				
		else if(args[0].equals("-linter"))
		{	
			if(args.length < 2)
			{
				System.out.println("USAGE to extract all lines for which the numeric value in the specifie column is within a user specified interval: java -cp . biokit/util/Utils -linter File columnMatching lowerBound upperBound splitter matchEntriesWhichAreWithinInterval");
				System.exit(0);
			}
			
			Utils.extractLinesWhereNumericValueInColumnIsWithinInterval(args[1], Integer.parseInt(args[2]), Double.parseDouble(args[3]), 
					Double.parseDouble(args[4]), args[5], Boolean.parseBoolean(args[6]));
		}
		else if(args[0].equals("-signum"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to apply signum to one column: java -cp . biokit/util/Utils -signum File col1 [col2...colN]");
				System.exit(0);
			}
			
			Utils.signumForColumns(args);
		}
		else if(args[0].equals("-log2"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to apply log2 to a set of columns: java -cp . biokit/util/Utils -log2 File col1 [col2...colN]");
				System.exit(0);
			}
			
			Utils.log2ForColumns(args);
		}
		else if(args[0].equals("-vennm"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to compute overlap of ids in m files: java -cp . biokit/util/Utils -vennm splitter outputBoolean [file column id file2 column2 id2 ...]");
				System.exit(0);
			}
			
			Utils.computeVennDiagramForMFiles(args[1], Boolean.parseBoolean(args[2]), args);
		}
		else if(args[0].equals("-vennmcid"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to compute overlap of ids in m files: java -cp . biokit/util/Utils -vennmcid splitter outputBoolean [file columnsCommaSparated id file2 columns2CommaSeparated id2 ...]");
				System.exit(0);
			}
			
			Utils.computeVennDiagramForMFilesCombinedID(args[1], Boolean.parseBoolean(args[2]), args);
		}
		else if(args[0].equals("-fl"))
		{
			if(args.length < 2)
			{
				System.out.println("USAGE to filter lines which contain an id specified in another file: java -cp . biokit/util/Utils -fl fileToBeFiltered ColumnWhichContainsFilterID fileContainingIDsUsedForFiltering ColumnContainingFilterIDs splitter onlyLinesWhichDoNOTContainFilterID");
				System.exit(0);
			}
			
			Utils.filterLines(args[1], Integer.parseInt(args[2]), args[3], Integer.parseInt(args[4]), args[5], Boolean.parseBoolean(args[6]));
		}
		else if(args[0].equals("-mysql"))
		{
			if(args.length < 3)
			{
				System.out.println("USAGE to retrieve data from mysql: java -cp . biokit/util/Utils -mysql \"statement replace all \" by \\\" \" config fastaOutput(boolean)");
				System.exit(0);
			}
			Utils.dumpMySQL(args[1], args[2], Boolean.parseBoolean(args[3]));
		}
		else if(args[0].equals("-sumCols"))
		{
			if(args.length != 5)
			{
				System.out.println("USAGE to summarize data in a file for a specified id column: java -cp . biokit/util/Utils -sumCols file idColumn splitterForMergedValues splitter");
				System.exit(0);
			}
			
			Utils.sumInformationForIDs(args[1], Integer.parseInt(args[2]), args[3], args[4]);
		}
		else if(args[0].equals("-append"))
		{
			if(args.length != 5)
			{
				System.out.println("USAGE append a value to each entry in a column: java -cp . biokit/util/Utils -append file column valueToBeAdded splitter");
				System.exit(0);
			}
			
			Utils.append(args[1], Integer.parseInt(args[2]), args[3], args[4]);
		}
		else if(args[0].equals("-splitSum"))
		{
			if(args.length != 7)
			{
				System.out.println("USAGE to summarize data in a file for an id column containing multiple ids: java -cp . biokit/util/Utils -splitSum file idColumn idSplit valueColumn overallSplitter valueSeparator");
				System.exit(0);
			}
			
			Utils.splitAndSummarize(args[1], Integer.parseInt(args[2]), args[3], Integer.parseInt(args[4]), args[5], args[6]);
		}
		else if(args[0].equals("-replace"))
		{
			if(args.length != 4)
			{
				System.out.println("USAGE to replace a query string by a replacement in a file: java -cp . biokit/util/Utils -replace file query target");
				System.exit(0);
			}
			
			Utils.replace(args[1], args[2], args[3]);
		}
		else if(args[0].equals("-replaceMulti"))
		{
			if(args.length != 3)
			{
				System.out.println("USAGE to replace a set of query strings in a file by new values in a file: java -cp . biokit/util/Utils -replaceMulti file fromToMappingFile");
				System.exit(0);
			}
			
			Utils.replaceMulti(args[1], args[2]);
		}
		else if(args[0].equals("-replaceComplicated"))
		{
			if(args.length != 7)
			{
				System.out.println("USAGE to replace a string in a column that needs to be splitted first: java -cp . biokit/util/Utils -replaceComplicated file column columnSplitter columnInternalSplitter query replacement");
				System.exit(0);
			}
			
			Utils.replaceValuesInColumn(args[1], Integer.parseInt(args[2]), args[3], args[4], args[5], args[6]);
		}
		else if(args[0].equals("-tab2FASTA"))
		{
			if(args.length != 4)
			{
				System.out.println("USAGE to make a tab delimited file to a FASTA file: java -cp . biokit/util/Utils -tab2FASTA inFile idCol sequenceCol");
				System.exit(0);
			}
			
			Utils.tab2FASTA(args[1], Integer.parseInt(args[3]), Integer.parseInt(args[2]));
		}
		else if(args[0].equals("-ef"))
		{
			if(args.length != 3)
			{
				System.out.println("USAGE to extract sequences that fit header information from FASTA file: java -cp . biokit/util/Utils -ef File headerContainsString");
				System.exit(0);
			}
			
			Utils.extractSequencesFromFASTA(args[1], args[2]);
		}
		else if(args[0].equals("-maplot"))
		{
			if(args.length != 7)
			{
				System.out.println("USAGE to compute MA plot columns: java -cp . biokit/util/Utils -maplot FileOne FileTwo idColumnOne idColumnTwo intensityColumnOne intensityColumnTwo splitter");
			}
			
			Utils.maplot(args[1], args[2], Integer.parseInt(args[3]), Integer.parseInt(args[4]), Integer.parseInt(args[5]), Integer.parseInt(args[6]), args[7]);
		}
		else if(args[0].equals("-sum"))
		{
			if(args.length != 4)
				System.out.println("USAGE to compute the sum of values in a column: java -cp . biokit/util/Utils -sum File col splitter");
			
			Utils.sumColumn(args[1], Integer.parseInt(args[2]), args[3]);
		}
		else if(args[0].equals("-rowsToColumns"))
		{
			if(args.length != 5)
				System.out.println("USAGE to summarize rows to columns and rows: -rowsToColumns file valueColumn headerColumn rowNameColumnsCommaSeparated, splitArgument idSplitArgument");
			
			Utils.rowsToColumns(args[1], Integer.parseInt(args[2]), Integer.parseInt(args[3]), args[4], args[5], args[6]);
		}
		else if(args[0].equals("-meansd"))
		{
			if(args.length != 4)
				System.out.println("USAGE to compute mean and sd for values in a row: java -cp . biokit/util/Utils -meansd File idCol splitter");
			
			Utils.computeStandardDeviationsAndMeanForValuesInRows(args[1], Integer.parseInt(args[2]), args[3]);
		}
		
		
	}
}