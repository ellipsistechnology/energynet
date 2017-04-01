package ellipsis.energy.calculation;

import static ellipsis.util.Sum.sum;
import static java.lang.Math.cos;
import static java.lang.Math.max;
import static java.lang.Math.sin;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealVector;

import ellipsis.util.Sum.Adder;

public class NewtonRaphsonPowerFlowSolver implements PowerFlowSolver 
{
    // Inputs:
	private Complex[] power;
	private boolean[] isSlackSource;
	private AdmittanceMatrix admittances;
	
	// Outputs:
	private Array2DRowRealMatrix jacobian;
    private double error = 0;
	
	public NewtonRaphsonPowerFlowSolver(Complex[] power, boolean[] isSlackSource, AdmittanceMatrix admittances)
	{
		this.power = power;
		this.isSlackSource = isSlackSource; // TODO need to handle more than one slack bus
		this.admittances = admittances;
	}
	
	@Override
	public double getError() 
	{
		return error;
	}

	@Override
	public Complex[] solve(Complex[] v) 
	{
		// Calculate delta P and delta Q:
		ArrayRealVector deltaPQ = calculatePowerDeltas(v);
		
		// Calculate Jacobian Matrix:
		jacobian = calculateJacobianMatrix(v);
		
//test(jacobian, v);

		// Solve for delta v:
		DecompositionSolver solver = new LUDecomposition(jacobian).getSolver();
		RealVector deltaV = solver.solve(deltaPQ);

		// Update voltages:
		int offset = 0;
		int size = v.length-1;
		for (int i = 0; i < v.length; i++) 
		{
			if(isSlackSource[i])
			{
				--offset;
			}
			else
			{
				double v_i_abs = v[i].abs();
				double v_i_arg = v[i].getArgument();
				int i2 = i+offset;
				v[i] = ComplexUtils.polar2Complex(
						v_i_abs + deltaV.getEntry(size + i2), 
						v_i_arg + deltaV.getEntry(i2));
			}
		}
		
		return v;
	}

//boolean test = false;
//public Map<String, Integer> busNumbers;
//private void test(RealMatrix jacobian, Complex[] v)
//{
//    if(!test)
//        return;
//    
//    System.out.print("dv/dv");
//    int size = v.length-1;
//    String busNames[] = new String[v.length];
//    for (String name : busNumbers.keySet())
//    {
//        busNames[busNumbers.get(name)] = name;
//                 
//    }
//    for (int i = 0; i < busNames.length; i++)
//    {
//        if(!isSlackSource[i])
//        {
//            System.out.print(',');
//            System.out.print(busNames[i]);
//        }
//    }
//    System.out.println();
//    
//    DecompositionSolver solver = new LUDecomposition(jacobian).getSolver();
//    RealMatrix inv = solver.getInverse();
//    int iOffset = 0;
//    for (int i = 0; i < v.length; i++)
//    {
//        if(isSlackSource[i])
//        {
//            --iOffset;
//            continue;
//        }
//
//        System.out.print(busNames[i]);
//
//        int jOffset = 0;
//        for (int j = 0; j < v.length; j++)
//        {
//            if(isSlackSource[j])
//            {
//                --jOffset;
//                continue;
//            }
//            
//            double ratio;
//            if(i == j)
//            {
//                ratio = 2;
//            }
//            else
//            {
//                int i2 = i+iOffset;
//                int j2 = j+jOffset;
//                
//                double dv_idP_i = inv.getEntry(i2+size, i2);
//                double dP_idV_j = jacobian.getEntry(i2, j2+size);
//                
//                double dv_idP_j = inv.getEntry(i2+size, j2);
//                double dP_jdV_j = jacobian.getEntry(j2, j2+size);
//                
//                double dv_idQ_i = inv.getEntry(i2+size, i2+size);
//                double dQ_idv_j = jacobian.getEntry(i2+size, j2+size);
//                
//                double dv_idQ_j = inv.getEntry(i2+size, j2+size);
//                double dQ_jdv_j = jacobian.getEntry(j2+size, j2+size);
//                
//                double dv_jdP_j = inv.getEntry(j2+size, j2);
//                double dv_jdQ_j = inv.getEntry(j2+size, j2+size);
//                
//                ratio = /*dv_idP_j/dv_jdP_j +*/ dv_idQ_j/*/dv_jdQ_j*/;//dv_idP_i*dP_idV_j + dv_idP_j*dP_jdV_j + dv_idQ_i*dQ_idv_j + dv_idQ_j*dQ_jdv_j;
//            }
//
//            System.out.print(',');
//            System.out.print(ratio);
//        }
//
//        System.out.println();
//    }
//}

    private Array2DRowRealMatrix calculateJacobianMatrix(Complex[] v) {
		int size = v.length-1; // without slack bus
		double[][] aJacobian = new double[size*2][size*2];
		int iOffset = 0; // Slack bus offsets
		
		for (int i = 0; i < v.length; i++) 
		{
			if(isSlackSource[i])
			{
				--iOffset;
				continue;
			}
			
			int jOffset = 0; // Slack bus offsets
			
			for (int j = 0; j < v.length; j++) 
			{
				if(isSlackSource[j])
				{
					--jOffset;
					continue;
				}

				int j2 = j+jOffset;
				int i2 = i+iOffset;
				
				Complex Y_ij = admittances.get(i, j);
				double Y_ij_abs = Y_ij.abs();
				double Y_ij_arg = Y_ij.getArgument();
				double v_i_abs = v[i].abs();
				double v_i_arg = v[i].getArgument();
				double v_j_abs = v[j].abs();
				double v_j_arg = v[j].getArgument();
				
				// dP/dd (top-left):
				aJacobian[i2][j2] = i == j ? 
						dP_idd_i(i, Y_ij_abs, Y_ij_arg, v_i_abs, v_i_arg, v) : 
						dP_idd_j(Y_ij_abs, Y_ij_arg, v_i_abs, v_i_arg, v_j_abs, v_j_arg);
						
				// dP/dv (top-right):
				aJacobian[i2][j2 + size] = i == j ?
						dP_idv_i(i, Y_ij_abs, Y_ij_arg, v_i_abs, v_i_arg, v) : 
						dP_idv_j(Y_ij_abs, Y_ij_arg, v_i_abs, v_i_arg, v_j_arg);
				
				// dQ/dd (bottom-left):
				aJacobian[i2 + size][j2] = i == j ? 
						dQ_idd_i(i, v_i_abs, v_i_arg, v) : 
						dQ_idd_j(Y_ij_abs, Y_ij_arg, v_i_abs, v_i_arg, v_j_abs, v_j_arg);
				
				// dQ/dv (bottom-right):
				aJacobian[i2 + size][j2 + size] = i == j ? 
						dQ_idv_i(i, Y_ij_abs, Y_ij_arg, v_i_abs, v_i_arg, v) : 
						dQ_idv_j(Y_ij_abs, Y_ij_arg, v_i_abs, v_i_arg, v_j_abs, v_j_arg);
			}
		}
		Array2DRowRealMatrix jacobian = new Array2DRowRealMatrix(aJacobian, false);
		return jacobian;
	}

	private double dQ_idv_j(double Y_ij_abs, double Y_ij_arg, double v_i_abs, double v_i_arg, double v_j_abs, double v_j_arg) 
	{
		return -v_i_abs*Y_ij_abs*sin(Y_ij_arg - v_i_arg + v_j_arg);
	}

	private double dQ_idv_i(final int i, double Y_ii_abs, double Y_ii_arg, double v_i_abs, final double v_i_arg, final Complex[] v) 
	{
		return -2*v_i_abs*Y_ii_abs*sin(Y_ii_arg) - sum(0, v.length-1, new Adder(){ public double add(int j) {
				if(i == j) return 0;
				Complex Y_ij = admittances.get(i, j);
				double Y_ij_abs = Y_ij.abs();
				double Y_ij_arg = Y_ij.getArgument();
				double v_j_abs = v[j].abs();
				double v_j_arg = v[j].getArgument();
				return v_j_abs*Y_ij_abs*sin(Y_ij_arg - v_i_arg + v_j_arg);
			}});
	}

	private double dQ_idd_j(double Y_ij_abs, double Y_ij_arg, double v_i_abs, double v_i_arg, double v_j_abs, double v_j_arg) 
	{
		return -v_i_abs*v_j_abs*Y_ij_abs*cos(Y_ij_arg - v_i_arg + v_j_arg);
	}

	private double dQ_idd_i(final int i, final double v_i_abs, final double v_i_arg, final Complex[] v) 
	{
		return sum(0, v.length-1, new Adder() { public double add(int j) {
			if(i == j) return 0;
			Complex Y_ij = admittances.get(i, j);
			double Y_ij_abs = Y_ij.abs();
			double Y_ij_arg = Y_ij.getArgument();
			double v_j_abs = v[j].abs();
			double v_j_arg = v[j].getArgument();
			return v_i_abs*v_j_abs*Y_ij_abs*cos(Y_ij_arg - v_i_arg + v_j_arg);
		}});
	}

	private double dP_idv_j(double Y_ij_abs, double Y_ij_arg, double v_i_abs, double v_i_arg, double v_j_arg) 
	{
		return v_i_abs*Y_ij_abs*cos(Y_ij_arg - v_i_arg + v_j_arg);
	}

	private double dP_idv_i(final int i, double Y_ii_abs, double Y_ii_arg, double v_i_abs, final double v_i_arg, final Complex[] v) 
	{
		return 2*v_i_abs*Y_ii_abs*cos(Y_ii_arg) + sum(0, v.length-1, new Adder(){ public double add(int j) {
			if(i == j) return 0;
			Complex Y_ij = admittances.get(i, j);
			double Y_ij_abs = Y_ij.abs();
			double Y_ij_arg = Y_ij.getArgument();
			double v_j_abs = v[j].abs();
			double v_j_arg = v[j].getArgument();
			return v_j_abs*Y_ij_abs*cos(Y_ij_arg - v_i_arg + v_j_arg);
		}});
	}

	private double dP_idd_j(double Y_ij_abs, double Y_ij_arg, double v_i_abs, double v_i_arg, double v_j_abs, double v_j_arg) 
	{
		return -v_i_abs*v_j_abs*Y_ij_abs*sin(Y_ij_arg - v_i_arg + v_j_arg);
	}

	private double dP_idd_i(final int i, double Y_ij_abs, double Y_ij_arg, final double v_i_abs, final double v_i_arg, final Complex[] v) 
	{
		return sum(0, v.length-1, new Adder() { public double add(int j) {
			if(i == j) return 0;
			Complex Y_ij = admittances.get(i, j);
			double Y_ij_abs = Y_ij.abs();
			double Y_ij_arg = Y_ij.getArgument();
			double v_j_abs = v[j].abs();
			double v_j_arg = v[j].getArgument();
			return v_i_abs*v_j_abs*Y_ij_abs*sin(Y_ij_arg - v_i_arg + v_j_arg);
		}});
	}

	private ArrayRealVector calculatePowerDeltas(Complex[] v) 
	{
		int size = v.length-1;
		double deltaPQ[] = new double[size*2];
		error = 0; // error will be set below
		
		int offset = 0; // Offset to skip slack bus
		for (int i = 0; i < v.length; i++) 
		{
			if(!isSlackSource[i])
			{
				double sumP = 0;
				double sumQ = 0;

				double v_i_abs = v[i].abs();
				double v_i_arg = v[i].getArgument();

				for (int j = 0; j < v.length; j++) 
				{
					Complex Y_ij = admittances.get(i, j);
					double Y_ij_abs = Y_ij.abs();
					double Y_ij_arg = Y_ij.getArgument();
					
					double v_j_arg = v[j].getArgument();
					double v_j_abs = v[j].abs();
					
					double vvY_abs = v_j_abs*Y_ij_abs;
                    double vvY_arg = Y_ij_arg + v_j_arg - v_i_arg;

                    sumP += vvY_abs * cos(vvY_arg);
					sumQ += vvY_abs * sin(vvY_arg);
				}

				double P = v_i_abs*sumP;
				double Q = -v_i_abs*sumQ;
				
				int i2 = i+offset;
				deltaPQ[i2] = power[i].getReal() - P;
				deltaPQ[size + i2] = power[i].getImaginary() - Q;

				// Update error:
				error = max(Math.abs(deltaPQ[i2]), 
						max(Math.abs(deltaPQ[size + i2]), 
							error));
			}
			else
			{
				--offset;
			}
		}
		
		return new ArrayRealVector(deltaPQ, false);
	}

    public Array2DRowRealMatrix getJacobian()
    {
        return jacobian;
    }
}
