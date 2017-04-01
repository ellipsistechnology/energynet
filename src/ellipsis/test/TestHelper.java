package ellipsis.test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import com.mls.util.Timer;

import ellipsis.util.Pair;


public class TestHelper
{
	public static void assertAssertsOn()
	{
		boolean assertEnabled = false;
		assert assertEnabled = true;
		if(!assertEnabled)
			throw new RuntimeException("asserts are NOT enabled.");
	}

	public static boolean equals(double d1, double d2, double accuracy)
	{
		return Math.abs(d1 - d2) < accuracy;
	}
	
	public boolean continueOnFail = false;
	protected boolean passing = true;
	public boolean verbose = false;
	public int iterations = -1;

	public void parseArgs(String[] args)
	{
		iterations = -1;
		for (int i = 0; i < args.length; i++)
		{
			if(args[i].startsWith("-m"))
			{
				String n = args[i].substring(2);
				iterations = Integer.valueOf(n);
			}
			else if(args[i].equals("-v"))
			{
				verbose = true;
			}
			else if(args[i].equals("-c"))
			{
				continueOnFail = true;
			}
		}
	}
	
	public void assertTrue(boolean test, String error)
	{
		if(continueOnFail)
		{
			if(!test)
			{
				System.err.println("Assertion failed: "+error);
				passing = false;
			}
		}
		else
		{
			assert test : error;
		}
	}
	
	protected void startTest()
	{
		passing = true;
		Timer.getTimer("test").start();
	}
	
	protected void endTest(String testName)
	{
		System.out.print(testName+" "+(passing?"passed":"failed"));
		long time = Timer.getTimer("test").stop();
		System.out.println(" in "+time+"ms");
	}
	
	
	//// IO ////
	
	protected void save(Object o, String fileName)
	{
		try
		{
			ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(new File(fileName)));
			out.writeObject(o);
			out.close();
		} 
		catch (IOException e)
		{
			throw new RuntimeException(e);
		}
	}
	
	protected Object read(String fileName)
	{
		try
		{
			File file = new File(fileName);
			if(!file.exists())
				return null;
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(file));
			Object o = in.readObject();
			in.close();
			return o;
		}
		catch (IOException | ClassNotFoundException e)
		{
			throw new RuntimeException(e);
		}
		
	}

	@SuppressWarnings("rawtypes")
	protected List<Pair<RealVector, Double>> convertFromSerializable(List<Pair> list)
	{
		if(list == null)
			return null;
		
		List<Pair<RealVector, Double>> newList = new ArrayList<>();
		for (Pair pair : list)
		{
			double[] key = (double[]) pair.getKey();
			Double value = (Double) pair.getValue();
			newList.add(new Pair<RealVector, Double>(key == null ? null : new ArrayRealVector(key), value));
		}
		return newList;
	}

	@SuppressWarnings({ "unchecked", "rawtypes" })
	protected Object convertToSerializable(List<Pair<RealVector, Double>> list)
	{
		List<Pair> newList = new ArrayList<>();
		for (Pair<RealVector, Double> pair : list)
		{
			RealVector k = pair.getKey();
			double[] key = k == null ? null : k.toArray();
			Double value = pair.getValue();
			newList.add(new Pair(key, value));
		}
		return newList;
	}
}
