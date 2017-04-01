package ellipsis.util;

import static java.lang.Math.abs;

public class ArrayHelper
{
	/**
	 * Maximum element by magnitude.
	 * @param x
	 * @return
	 */
	public static double maxMagnitude(double[] x)
	{
		double max = 0;
		for (int i = 0; i < x.length; i++)
		{
			double abs = abs(x[i]);
			if(abs > max)
				max = abs;
		}
		return max;
	}
	public static double minMagnitude(double[] x)
	{
		double min = Double.MAX_VALUE;
		for (int i = 0; i < x.length; i++)
		{
			double abs = abs(x[i]);
			if(abs < min)
				min = abs;
		}
		return min;
	}

	public static double[] subtract(double[] x1, double[] x2)
	{
		double[] ret = new double[x1.length];
		for (int i = 0; i < ret.length; i++)
		{
			ret[i] = x1[i] - x2[i];
		}
		return ret;
	}

	public static double[] add(double[] x1, double[] x2)
	{
		double[] ret = new double[x1.length];
		for (int i = 0; i < ret.length; i++)
		{
			ret[i] = x1[i] + x2[i];
		}
		return ret;
	}

	public static double[] mul(double a, double[] x)
	{
		double[] ret = new double[x.length];
		for (int i = 0; i < ret.length; i++)
		{
			ret[i] = a*x[i];
		}
		return ret;
	}
	
	public static double[] ebeDiv(double[] x1, double[] x2)
	{
		double[] ret = new double[x1.length];
		for (int i = 0; i < ret.length; i++)
		{
			ret[i] = x1[i]/x2[i];
		}
		return ret;
	}

	public static double[] array(int length, double preset)
	{
		double[] data = new double[length];
		for (int i = 0; i < data.length; i++)
		{
			data[i] = preset;
		}
		return data;
	}
	
	public static double[] array(double... ds)
	{
	    return ds;
	}
}