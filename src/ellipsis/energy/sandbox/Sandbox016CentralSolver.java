package ellipsis.energy.sandbox;

import static ellipsis.util.Sum.sum;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import ellipsis.energy.calculation.AnalysisResults;

public class Sandbox016CentralSolver extends Sandbox016
{
	public static void main(String[] args)
	{
		new Sandbox016CentralSolver().run();
	}
	
//	@Override
//	protected void initGrid()
//	{
//		super.initGrid();
//		
//		grid.getSource("DG1").setPowerOutput(0.0,  0.0);
//		grid.getSource("DG3").setPowerOutput(0.0,  0.0);
//		grid.getSource("DG8").setPowerOutput(0.0,  0.0);
//		grid.getSource("DG14").setPowerOutput(0.0,  0.0);
//		grid.getSource("DG18").setPowerOutput(0.0,  0.0);
//		grid.getSource("DG26").setPowerOutput(0.0,  0.0);
//		grid.getSource("DG29").setPowerOutput(0.0,  0.0);
//	}
	
//	@Override
//	protected void initGrid()
//	{
//		// Copied from Sandbox008:
//		grid = 
//            Grid.grid().
//            
//                Bus(SLACK_BUS).
//                
//                    SlackSource(SLACK_SOURCE, 1.00*BASE_VOLTAGE, 0, 0, 0).
//                    
//                    Line(LINE_1_2, 1, 0.02*BASE_IMPEDANCE, 0.04*BASE_IMPEDANCE).
//                        Bus(BUS_2).
//                            Load(LOAD_2).
//                            Line(LINE_2_3, BUS_3, 1, 0.0125*BASE_IMPEDANCE, 0.025*BASE_IMPEDANCE).
//                            Line(LINE_2_4, 1.5, 0.0125*BASE_IMPEDANCE, 0.025*BASE_IMPEDANCE).
//	                        	Bus(BUS_4).
//	                        		DistributedSource(DG_4).
//	                        	terminate().
//	                        terminate().
//                        terminate().
//                    terminate().
//                    
//                    Line(LINE_1_3, 1, 0.01*BASE_IMPEDANCE, 0.03*BASE_IMPEDANCE).
//                        Bus(BUS_3).
//                        	DistributedSource(DG_3).
//                        	Line(LINE_3_5, 1.5, 0.0125*BASE_IMPEDANCE, 0.025*BASE_IMPEDANCE).
//	                        	Bus(BUS_5).
//	                        		Load(LOAD_5).
//	                        	terminate().
//	                        terminate().
//                        terminate().
//                    terminate().
//                    
//                terminate().
//                
//            grid();
//		
//		grid.setBaseVoltage(BASE_VOLTAGE);
//		grid.setBasePower(BASE_POWER);
//		
//		grid.getLoad(LOAD_2).setLoad(new Complex(150e6, 50e6));
////		grid.getLoad(LOAD_5).setLoad(new Complex(80e6, 25e6));
//		
//        grid.getSource(DG_3).setPmax(650e6);
//        grid.getSource(DG_3).setPowerOutput(100e6, 0);
////        grid.getSource(DG_4).setPmax(650e6);
////        grid.getSource(DG_4).setPowerOutput(100e6, 0);
//	}
	
	@Override
	protected void run()
	{
//		START_WITH_TRUE_VOTLAGES = true;
		AnalysisResults results = init();
//e.setEntry(1, 1.0); // FIXME DO NOT COMMIT
//f.setEntry(0, 0.0); // FIXME DO NOT COMMIT
		
		solve(); // Centralized
//		loop(); // Distributed
		
        testResults(results);
		finish();
	}

	private void solve()
	{
		AugmentedLagrangeSolver solver = new AugmentedLagrangeSolver();
		solver.setCostFunction(this::cost);
		solver.setCostGradient(this::gradC);
		solver.setConstraintFunction(this::g);
		solver.setConstraintJacobian(this::gradG);
		solver.setProjection(this::projectX);
		solver.setBeta(1.1);
		solver.setIterations(100);
		solver.setGradDecIterations(10);
		solver.setGradDecTargetGradient(1.0);
		
		RealVector initialX = p.append(q).append(e).append(f);
		RealVector x = solver.solve(initialX);

		int length = x.getDimension()/4;
		p = x.getSubVector(0, length);
		q = x.getSubVector(length, length);
		e = x.getSubVector(2*length, length);
		f = x.getSubVector(3*length, length);
//		System.out.println("Solution:\nindex,p,q,e,f");
//		for(int i = 0; i < length; ++i)
//			System.out.println(i+","+p.getEntry(i)+","+q.getEntry(i)+","+e.getEntry(i)+","+f.getEntry(i)+",");
	}
	
	public double cost(RealVector x)
	{
		int length = x.getDimension()/4;
		RealVector p = x.getSubVector(0, length);
		RealVector q = x.getSubVector(length, length);
		return cost(p, q);
	}
	
	public RealVector gradC(RealVector x)
	{
		int length = x.getDimension()/4;
		RealVector p = x.getSubVector(0, length);
		RealVector q = x.getSubVector(length, length);
		RealVector grad = new ArrayRealVector(4*length, 0.0);
		for(int i = 0; i < length; ++i)
		{
			if(generatorMask.getEntry(i) == 1)
			{
				double gradC_p_i = gradC_p(p, q, i);
				double gradC_q_i = gradC_q(p, q, i);
				grad.setEntry(i, gradC_p_i);
				grad.setEntry(i+length, gradC_q_i);
			}
		}
		return grad;
	}
	
	public RealVector g(RealVector x)
	{
		int length = x.getDimension()/4;
		RealVector p = x.getSubVector(0, length);
		RealVector q = x.getSubVector(length, length);
		RealVector e = x.getSubVector(2*length, length);
		RealVector f = x.getSubVector(3*length, length);
		RealVector gp = gp(p, q, e, f);
		RealVector gq = gq(p, q, e, f);
		return gp.append(gq);
	}
	
	public RealMatrix gradG(RealVector x)
	{
		RealMatrix grad_eGp = gradG(x, this::gradGpi_dej);
		RealMatrix grad_fGp = gradG(x, this::gradGpi_dfj);
		RealMatrix grad_eGq = gradG(x, this::gradGqi_dej);
		RealMatrix grad_fGq = gradG(x, this::gradGqi_dfj);
		
		int dimension = x.getDimension()/4;
		RealMatrix grad = new Array2DRowRealMatrix(4*dimension, 2*dimension);
		
		// [   -I          0    ]
		// [    0         -I    ]
		// [ grad_eGp  grad_eGq ]
		// [ grad_fGp  grad_fGq ]
		for(int i = 0; i < dimension; ++i)
		{
			for(int j = 0; j < dimension; ++j)
			{
				// -I (generators only):
				if(i == j && generatorMask.getEntry(i) == 1)
				{
					grad.setEntry(i, j, -1.0);
					grad.setEntry(i+dimension, j+dimension, -1.0);
				}
				
				grad.setEntry(2*dimension+i, j, grad_eGp.getEntry(i, j));
				grad.setEntry(2*dimension+i, j+dimension, grad_eGq.getEntry(i, j));

				grad.setEntry(3*dimension+i, j, grad_fGp.getEntry(i, j));
				grad.setEntry(3*dimension+i, j+dimension, grad_fGq.getEntry(i, j));
			}
		}
		
		return grad;
	}
    
	protected static interface Grad
	{
		/**
		 * Gradient of constraint i W.R.T. index j.
		 */
		double value(RealVector e, RealVector f, int i, int j);
	}
    protected RealMatrix gradG(RealVector x, Grad grad)
    {
		int length = x.getDimension()/4;
		RealVector e = x.getSubVector(2*length, length);
		RealVector f = x.getSubVector(3*length, length);
    	int dimension = e.getDimension();
		RealMatrix jacobian = new Array2DRowRealMatrix(dimension, dimension);
		for(int i = 0; i < jacobian.getRowDimension(); ++i)
		{
			if(i == slackIndex)
				continue;
			
			for(int j = 0; j < jacobian.getColumnDimension(); ++j)
			{
				if(j == slackIndex)
					continue;
				
				jacobian.setEntry(i, j, grad.value(e, f, j, i));
			}
		}
		return jacobian;
    }
	
	private double gradGpi_dej(RealVector e, RealVector f, int i, int j)
	{
		double e_i = e.getEntry(i);
		double f_i = f.getEntry(i);
		
		if(i != j)
			return e_i*G(i, j) + f_i*B(i, j);
		
		int dimension = p.getDimension();
		return sum(n -> e.getEntry(n)*G(i, n) - f.getEntry(n)*B(i, n), dimension, i) + 2*G(i, i)*e_i;
	}
	
	private double gradGpi_dfj(RealVector e, RealVector f, int i, int j)
	{
		double e_i = e.getEntry(i);
		double f_i = f.getEntry(i);
		
		if(i != j)
			return f_i*G(i, j) - e_i*B(i, j);
		
		int dimension = p.getDimension();
		return sum(n -> f.getEntry(n)*G(i, n) + e.getEntry(n)*B(i, n), dimension, i) + 2*G(i, i)*f_i;
	}
	
	private double gradGqi_dej(RealVector e, RealVector f, int i, int j)
	{
		double e_i = e.getEntry(i);
		double f_i = f.getEntry(i);
		
		if(i != j)
			return f_i*G(i, j) - e_i*B(i, j);
		
		int dimension = p.getDimension();
		return sum(n -> -f.getEntry(n)*G(i, n) - e.getEntry(n)*B(i, n), dimension, i) - 2*B(i, i)*e_i;
	}
	
	private double gradGqi_dfj(RealVector e, RealVector f, int i, int j)
	{
		double e_i = e.getEntry(i);
		double f_i = f.getEntry(i);
		
		if(i != j)
			return -e_i*G(i, j) - f_i*B(i, j);
		
		int dimension = p.getDimension();
		return sum(n -> e.getEntry(n)*G(i, n) - f.getEntry(n)*B(i, n), dimension, i) - 2*B(i, i)*f_i;
	}

	public RealVector projectX(RealVector x)
	{
		return x; // TODO
	}
}
