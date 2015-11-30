package BioKit;

public class NumericSorterContainer<E> implements Comparable<NumericSorterContainer<E>>
{
	private double value;
	private E content;
	boolean ignoreEqual;
	
	public NumericSorterContainer(double value, E content, boolean ignoreEqual)
	{
		this.value = value;
		this.content = content;
		this.ignoreEqual = ignoreEqual;
	}
	
	public E getContent()
	{
		return content;
	}
	
	public void setContent(E content)
	{
		this.content = content;
	}
	
	public double getValue()
	{
		return value;
	}

	public int compareTo(NumericSorterContainer<E> o)
	{
		if(value > o.value)
			return 1;
		else if(value == o.value && !ignoreEqual)
			return 0;
		else
			return -1;
	}	
}
