package ellipsis.util;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.linear.SingularValueDecomposition;

public class MatrixHelper
{
	public static RealMatrix invert(RealMatrix m)
	{
		DecompositionSolver solver = new LUDecomposition(m).getSolver();
        RealMatrix inv = solver.getInverse();
        return inv;
	}
	
	public static int rank(RealMatrix m)
	{
	    return new SingularValueDecomposition(m).getRank(); 
	}

	public static RealVector solve(RealMatrix m, RealVector v)
	{
		try
		{
			DecompositionSolver solver = new LUDecomposition(m).getSolver();
			return solver.solve(v);
		}
		catch(SingularMatrixException sme)
		{
			System.err.println("Singular matrix given to MatrixHelper.solve().");
			return null;
		}
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
	}
	
	public static RealMatrix matrix(double[]... rows)
	{
	    return new Array2DRowRealMatrix(rows, false);
	}
	
	public static void main(String[] args)
	{
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
		
		Complex[] l = eigenValues(L);
		for (int i = 0; i < l.length; i++)
		{
			System.out.println(l[i]);
		}
	}
}