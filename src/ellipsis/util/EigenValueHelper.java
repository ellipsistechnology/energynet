package ellipsis.util;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;

public class EigenValueHelper
{
	public static void main(String[] args) throws InterruptedException
	{
//		int count = 1000;
//		long[] time = new long[count];
//		long[] nextTime = new long[count];
//		for(int i = 0; i < count; ++i)
//		{
//			time[i] = System.nanoTime();
//			nextTime[i] = System.nanoTime();
//			Thread.sleep(1);
//		}
//		
//		for(int i = 0; i < count; ++i)
//		{
//			System.out.println(nextTime[i]-time[i]);
//		}
		
		//
//		double d = 67;
//		double i = 0.03/365;
//		int T = 20*365;
//		double x0 = 0;
//		for(int j = 0; j < T-1; ++j)
//		{
//			x0 += Math.pow(1+i,j);
//		}
//		x0 *= d/Math.pow(1+i,T);
//		
//		System.out.println(x0);
//		
//		//
//		double xt = 0;
////		int count = 0;
//		for(int j = T; j >= 0; --j)
//		{
////			System.out.println(xt);
//			xt = (xt+d)/(1+i);
////			if(count%100 == 0)
////				try
////				{
////					System.in.read();
////				} catch (IOException e)
////				{
////					// TODO Auto-generated catch block
////					e.printStackTrace();
////				}
////			++count;
//		}
//		System.out.println(xt);
		
		//
		// Lambda (sensitivities):
		double[][] data = new double[][]{
				{0.001288998,	4.62E-04,	7.74E-04,	7.79E-04,	7.79E-04,	0.001117052,	0.001116541,	0.001117313,	0.001116594,	0.001116945,	-0.010291368,	-6.34E-05,	-0.004165141,	-0.004161934,	-0.004161792,	-0.008253614,	-0.008253928,	-0.008253502,	-0.008253896,	-0.008253716},
				{4.58E-04,	4.62E-04,	4.63E-04,	4.63E-04,	4.63E-04,	4.61E-04,	4.61E-04,	4.61E-04,	4.61E-04,	4.61E-04,	                        -5.27E-05,	-5.78E-05,	-5.43E-05,	-5.43E-05,	-5.41E-05,	-5.23E-05,	-5.27E-05,	-5.23E-05,	-5.27E-05,	-5.24E-05},
				{7.97E-04,	4.62E-04,	0.004623717,	8.02E-04,	8.02E-04,	8.00E-04,	7.99E-04,	8.00E-04,	7.99E-04,	8.00E-04,	                    -0.00415576,	-6.09E-05,	-0.00945217,	-0.004157359,	-0.004157208,	-0.004155399,	-0.004155761,	-0.004155342,	-0.004155738,	-0.004155481},
				{7.92E-04,	4.62E-04,	7.93E-04,	0.001411464,	0.001411722,	7.96E-04,	7.95E-04,	7.96E-04,	7.95E-04,	7.95E-04,	-0.004156669,	-6.14E-05,	-0.004161458,	-0.010334764,	-0.010334616,	-0.004156312,	-0.00415667,	-0.004156255,	-0.004156647,	-0.004156394},
				{7.91E-04,	4.62E-04,	7.91E-04,	0.00140811,	0.00177886,	7.94E-04,	7.94E-04,	7.94E-04,	7.94E-04,	7.94E-04,	                        -0.004156945,	-6.15E-05,	-0.004161729,	-0.010335236,	-0.014046099,	-0.004156589,	-0.004156945,	-0.004156531,	-0.004156923,	-0.00415667},
				{0.001112295,	4.62E-04,	7.73E-04,	7.78E-04,	7.78E-04,	0.001903936,	0.001114589,	0.001115358,	0.001114641,	0.001482314,	-0.008254214,	-6.35E-05,	-0.004165339,	-0.004162132,	-0.00416199,	-0.017287629,	-0.008254215,	-0.008253791,	-0.008254184,	-0.011946875},
				{0.001118229,	4.62E-04,	7.76E-04,	7.81E-04,	7.81E-04,	0.001121075,	0.001120558,	0.001121336,	0.001120611,	0.001120967,	-0.008253336,	-6.32E-05,	-0.004164734,	-0.004161527,	-0.004161384,	-0.00825302,	-0.008253337,	-0.008252908,	-0.008253305,	-0.008253122},
				{0.001112163,	4.62E-04,	7.73E-04,	7.78E-04,	7.78E-04,	0.001114965,	0.001114457,	0.001290277,	0.001114467,	0.001114857,	-0.008254234,	-6.35E-05,	-0.004165353,	-0.004162145,	-0.004162004,	-0.008253923,	-0.008254235,	-0.010565344,	-0.008311854,	-0.008254024},
				{0.001117998,	4.62E-04,	7.76E-04,	7.81E-04,	7.81E-04,	0.001120842,	0.001120326,	0.001121096,	0.001120378,	0.001120735,	-0.00825337,	-6.32E-05,	-0.004164757,	-0.00416155,	-0.004161407,	-0.008253054,	-0.008253371,	-0.008310596,	-0.00831099,	-0.008253156},
				{0.001113422,	4.62E-04,	7.74E-04,	7.78E-04,	7.79E-04,	0.00148411,	0.001115723,	0.001116493,	0.001115775,	0.001853152,	    -0.008254048,	-6.35E-05,	-0.004165224,	-0.004162017,	-0.004161875,	-0.011946526,	-0.008254048,	-0.008253623,	-0.008254017,	-0.015642056},
				
				{0.010182062,	5.75E-05,	0.004144131,	0.004144458,	0.004147347,	0.008207049,	0.008195174,	0.00820818,	0.008195557,	0.008206693,	0.001381581,	4.61E-04,	8.31E-04,	8.31E-04,	8.33E-04,	0.001213233,	0.001206102,	0.001213245,	0.001206334,	0.001212152},
				{5.69E-05,	5.78E-05,	5.84E-05,	5.79E-05,	5.80E-05,	5.75E-05,	5.74E-05,	5.76E-05,	5.74E-05,	5.75E-05,	4.64E-04,	4.63E-04,	4.64E-04,	4.63E-04,	4.63E-04,	4.64E-04,	4.64E-04,	4.64E-04,	4.64E-04,	4.64E-04},
				{0.004121287,	5.79E-05,	0.009464105,	0.004167549,	0.004170454,	0.004151253,	0.004145215,	0.004151822,	0.004145409,	0.004151075,	8.50E-04,	4.63E-04,	0.004666561,	8.36E-04,	8.37E-04,	8.54E-04,	8.50E-04,	8.54E-04,	8.50E-04,	8.53E-04},
				{0.004120304,	5.78E-05,	0.004166226,	0.010347101,	0.010354272,	0.004150263,	0.004144226,	0.004150832,	0.00414442,	0.004150085,	8.50E-04,	4.63E-04,	8.36E-04,	0.00145685,	0.001460964,	8.54E-04,	8.50E-04,	8.54E-04,	8.50E-04,	8.53E-04},
				{0.004121793,	5.79E-05,	0.004167731,	0.010350839,	0.014065248,	0.004151762,	0.004145724,	0.004152332,	0.004145918,	0.004151584,	8.50E-04,	4.63E-04,	8.36E-04,	0.001457376,	0.001832821,	8.54E-04,	8.50E-04,	8.54E-04,	8.50E-04,	8.53E-04},
				{0.008177435,	5.77E-05,	0.004159051,	0.004159379,	0.004162279,	0.017284436,	0.008224679,	0.008237732,	0.008225063,	0.011937843,	0.001210456,	4.62E-04,	8.34E-04,	8.34E-04,	8.36E-04,	0.002013322,	0.001210444,	0.001217613,	0.001210677,	0.001588124},
				{0.008171274,	5.77E-05,	0.004155918,	0.004156246,	0.004159143,	0.008230391,	0.008218482,	0.008231525,	0.008218867,	0.008230034,	0.001209544,	4.62E-04,	8.34E-04,	8.33E-04,	8.35E-04,	0.001216683,	0.001209532,	0.001216696,	0.001209765,	0.001215599},
				{0.008177953,	5.78E-05,	0.004159315,	0.004159643,	0.004162542,	0.008237118,	0.0082252,	0.010554246,	0.00828335,	0.008236761,	0.001210532,	4.62E-04,	8.35E-04,	8.34E-04,	8.36E-04,	0.001217678,	0.001210521,	0.001394796,	0.001210754,	0.001216593},
				{0.008171465,	5.77E-05,	0.004156015,	0.004156343,	0.00415924,	0.008230583,	0.008218675,	0.008289526,	0.008276779,	0.008230227,	0.001209572,	4.62E-04,	8.34E-04,	8.33E-04,	8.35E-04,	0.001216712,	0.001209561,	0.001216774,	0.001209794,	0.001215628},
				{0.008177231,	5.77E-05,	0.004158947,	0.004159276,	0.004162175,	0.011938065,	0.008224474,	0.008237526,	0.008224859,	0.015636471,	0.001210426,	4.62E-04,	8.34E-04,	8.34E-04,	8.36E-04,	0.001589656,	0.001210414,	0.001217583,	0.001210647,	0.001958325}
		};
		RealMatrix L = new Array2DRowRealMatrix(data);
		
		testEigenValues(L);
		
//		RealMatrix M = new Array2DRowRealMatrix(new double[][] {
//				{0.001288998, -0.010291368},
//				{0.010182062, 0.001381581},
//		});
//		testEigenValues(M);
		
		
		//
//		double[][] data = new double[][]{
//				{0.9,0.1},
//				{0.0999,0.9}
//		};
//		RealMatrix m = new Array2DRowRealMatrix(data);
//		
//		for(int i = 0; i < 100; ++i)
//		{
//			printMatrix(m);
//			m = m.multiply(m);
//		}
	}

	public static void testEigenValues(RealMatrix L)
	{
		// Constraint matrix:
		RealMatrix N1 = inverse(L);
		System.out.println("N=Inverse L:");
		testStepSize(0.8, L, N1);

		RealMatrix N2 = MatrixUtils.createRealIdentityMatrix(L.getColumnDimension());
		System.out.println("\nN=I:");
		testStepSize(0.8, L, N2);
		
		RealMatrix N3 = MatrixUtils.createRealIdentityMatrix(L.getColumnDimension());
		for(int i = 0; i < L.getColumnDimension(); ++i)
			N3.setEntry(i, i, 1+0.1*(2*Math.random()-1));
		System.out.println("\nN~I:");
		testStepSize(0.8, L, N3);

		RealMatrix N4 = new Array2DRowRealMatrix(L.getRowDimension(), L.getColumnDimension());
		for(int i = 0; i < L.getRowDimension(); ++i)
			for(int j = 0; j < L.getColumnDimension(); ++j)
				N4.setEntry(i, j, 2*Math.random()-1);
		System.out.println("\nN=Rand:");
		testStepSize(0.8, L, N4);

		RealMatrix N5 = new Array2DRowRealMatrix(L.getRowDimension(), L.getColumnDimension());
		for(int i = 0; i < L.getRowDimension(); ++i)
			for(int j = 0; j < L.getColumnDimension(); ++j)
				N5.setEntry(i, j, 0.1);
		System.out.println("\nN=[1]:");
		testStepSize(0.8, L, N5);
	}

	public static void testStepSize(double gamma, RealMatrix L, RealMatrix N)
	{
		// Product:
		RealMatrix LN = L.multiply(N);
		
		// Eigen values:
		Complex[] values = eigenValues(LN);
		System.out.println("LN eigen values:");
		printVector(values);
		
		// Max eigen value:
		double maxValue = maxValue(values);
		System.out.println("Max value = "+maxValue);
		
		// Divide N by max eigen value:
		N = N.scalarMultiply(gamma/maxValue);
		LN = L.multiply(N);
		
		// New eigen values:
		values = eigenValues(LN);
		System.out.println("New LN eigen values:");
		printVector(values);
		maxValue = maxValue(values);
		System.out.println("New max value = "+maxValue);
		
		// Sequence:
		RealMatrix LNj = LN;
		for(int j = 0; j < 10; ++j)
		{
//			double max = maxValue(eigenValues(LNj));
			System.out.println(j+","+LNj.getFrobeniusNorm());
			LNj = LNj.multiply(LN);
		}
		
		// Show LN:
//		System.out.println("");
		printMatrix(LN);
	}

	@SuppressWarnings("unused")
	private static double maxValue(double[] values1)
	{
		double max = 0;
		for (int i = 0; i < values1.length; i++)
		{
			if(Math.abs(values1[i]) > max)
				max = Math.abs(values1[i]);
		}
		return max;
	}

	public static double maxValue(Complex[] values1)
	{
		double max = 0;
		for (int i = 0; i < values1.length; i++)
		{
			if(values1[i].abs() > max)
				max = values1[i].abs();
		}
		return max;
	}

	public static Complex[] eigenValues(RealMatrix M)
	{
		EigenDecomposition decomp = new EigenDecomposition(M);
		double[] r = decomp.getRealEigenvalues();
		double[] i = decomp.getImagEigenvalues();
		Complex[] c = new Complex[r.length];
		for (int j = 0; j < c.length; j++)
		{
			c[j] = new Complex(r[j], i[j]);
		}
		return c;
		
//		EigenvalueDecomposition decomp = new EigenvalueDecomposition(apacheToColt(M));
//		DoubleMatrix1D values = decomp.getRealEigenvalues();
//		
//		double[] vs = new double[values.size()];
//		for (int i = 0; i < vs.length; i++)
//		{
//			vs[i] = values.get(i);
//		}
//		return vs;
	}

	public static void printVector(Object[] values)
	{
		for (int i = 0; i < values.length; i++)
		{
			if(values[i] instanceof Complex && ((Complex)values[i]).getImaginary() == 0)
				System.out.print(((Complex)values[i]).getReal());
			else
				System.out.print(values[i]);
			System.out.print(',');
		}
		System.out.println();
	}

	@SuppressWarnings("unused")
	private static DoubleMatrix2D apacheToColt(RealMatrix m)
	{
		return new DenseDoubleMatrix2D(m.getData());
	}

	private static RealMatrix inverse(RealMatrix M)
	{
		SingularValueDecomposition decomp = new SingularValueDecomposition(M);
		return decomp.getSolver().getInverse();
	}

	private static void printMatrix(RealMatrix m)
	{
		printMatrix(m, 1e-12);
	}
	private static void printMatrix(RealMatrix m, double minMagnitude)
	{
		StringBuffer sb = new StringBuffer();
		sb.append('[');
		for(int i = 0; i < m.getRowDimension(); ++i)
		{
			sb.append('[');
			for(int j = 0; j < m.getColumnDimension(); ++j)
			{
				double m_ij = m.getEntry(i, j);
				if(Math.abs(m_ij) < minMagnitude)
					sb.append('0');
				else
					sb.append(m_ij);
				sb.append(',');
			}
			sb.append("],\n");
		}
		sb.append(']');
		
		System.out.println(sb);
	}
}