package ellipsis.energy.sandbox;

import static ellipsis.energy.sandbox.GridConstants.BASE_IMPEDANCE;
import static ellipsis.energy.sandbox.GridConstants.BASE_POWER;
import static ellipsis.energy.sandbox.GridConstants.BASE_VOLTAGE;
import static ellipsis.energy.sandbox.GridConstants.BUS_2;
import static ellipsis.energy.sandbox.GridConstants.BUS_3;
import static ellipsis.energy.sandbox.GridConstants.BUS_4;
import static ellipsis.energy.sandbox.GridConstants.BUS_5;
import static ellipsis.energy.sandbox.GridConstants.DG_3;
import static ellipsis.energy.sandbox.GridConstants.DG_4;
import static ellipsis.energy.sandbox.GridConstants.LINE_1_2;
import static ellipsis.energy.sandbox.GridConstants.LINE_1_3;
import static ellipsis.energy.sandbox.GridConstants.LINE_2_3;
import static ellipsis.energy.sandbox.GridConstants.LINE_2_4;
import static ellipsis.energy.sandbox.GridConstants.LINE_3_5;
import static ellipsis.energy.sandbox.GridConstants.LOAD_2;
import static ellipsis.energy.sandbox.GridConstants.LOAD_5;
import static ellipsis.energy.sandbox.GridConstants.SLACK_BUS;
import static ellipsis.energy.sandbox.GridConstants.SLACK_SOURCE;
import static ellipsis.energy.sandbox.Sandbox004.logBusses;
import static ellipsis.energy.sandbox.Sandbox004.maxPower;
import static ellipsis.energy.sandbox.Sandbox004.p0;
import static ellipsis.energy.sandbox.Sandbox004.q0;
import static ellipsis.energy.sandbox.Sandbox008.debugLogNames;

import java.util.Map;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.calculation.LoadFlowAnalyser;
import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.Grid;

/**
 * Test gradient method optimisation of OPF problem
 * with real and imaginary split (rather than the usual
 * abs/arg split for voltage).
 * Vector matrix form.
 * @author bmillar
 *
 */
public class Sandbox012
{
	public static void main(String[] args)
	{
		new Sandbox012().run();
	}
	
	protected double gamma = 0.01;
	protected RealVector p_max;
	protected RealMatrix generatorMask;
	
	protected RealVector p, q, e, f, l_p, l_q;
	protected RealMatrix G, B;
	private int slackIndex;
	private Grid grid;
	
	protected void run()
	{
		init();
		loop();
	}
	
	protected void init()
	{
		initGrid();
		
		// Analyse grid to get parameters:
		LoadFlowAnalyser lfa = new LoadFlowAnalyser(grid);
        lfa.setBasePower(grid.getBasePower());
        lfa.setBaseVoltage(grid.getBaseVoltage());
        AnalysisResults results = lfa.analyse();
        
        logBusses(results);
        
        // Initialise variables:
        initVariables(results);
        
        // Debug:
		debugHeader(results);
	}
	
	private void initGrid()
	{
		// Copied from Sandbox008:
		grid = 
            Grid.grid().
            
                Bus(SLACK_BUS).
                
                    SlackSource(SLACK_SOURCE, 1.00*BASE_VOLTAGE, 0, 0, 0).
                    
                    Line(LINE_1_2, 1, 0.02*BASE_IMPEDANCE, 0.04*BASE_IMPEDANCE).
                        Bus(BUS_2).
                            Load(LOAD_2).
                            Line(LINE_2_3, BUS_3, 1, 0.0125*BASE_IMPEDANCE, 0.025*BASE_IMPEDANCE).
                            Line(LINE_2_4, 1.5, 0.0125*BASE_IMPEDANCE, 0.025*BASE_IMPEDANCE).
	                        	Bus(BUS_4).
	                        		DistributedSource(DG_4).
	                        	terminate().
	                        terminate().
                        terminate().
                    terminate().
                    
                    Line(LINE_1_3, 1, 0.01*BASE_IMPEDANCE, 0.03*BASE_IMPEDANCE).
                        Bus(BUS_3).
                        	DistributedSource(DG_3).
                        	Line(LINE_3_5, 1.5, 0.0125*BASE_IMPEDANCE, 0.025*BASE_IMPEDANCE).
	                        	Bus(BUS_5).
	                        		Load(LOAD_5).
	                        	terminate().
	                        terminate().
                        terminate().
                    terminate().
                    
                terminate().
                
            grid();
		
		grid.setBaseVoltage(BASE_VOLTAGE);
		grid.setBasePower(BASE_POWER);
		
		grid.getLoad(LOAD_2).setLoad(new Complex(150e6, 50e6));
		grid.getLoad(LOAD_5).setLoad(new Complex(80e6, 25e6));
		
        grid.getSource(DG_3).setPmax(650e6);
        grid.getSource(DG_3).setPowerOutput(100e6, 0);
        grid.getSource(DG_4).setPmax(650e6);
        grid.getSource(DG_4).setPowerOutput(100e6, 0);
	}

	protected void initVariables(AnalysisResults results)
	{
        slackIndex = results.getSlackIndex();
        
		// Get admittances:
        FieldMatrix<Complex> Y = results.getAdmittanceMatrix().Y;
        int dimension = Y.getRowDimension();
        B = new Array2DRowRealMatrix(dimension, dimension);
        G = new Array2DRowRealMatrix(dimension, dimension);
        for(int i = 0; i < dimension; ++i)
        {
        	for(int j = 0; j < dimension; ++j)
        	{
        		Complex Y_ij = Y.getEntry(i, j);
				G.setEntry(i, j, Y_ij.getReal());
				B.setEntry(i, j, Y_ij.getImaginary());
        	}
        }
        
        // Prepare generator mask:
        generatorMask = new Array2DRowRealMatrix(dimension, dimension);
        for (String name : results.getBusNumbers().keySet())
		{
			int index = results.getBusNumbers().get(name);
			Bus bus = grid.getBus(name);
			if(!bus.getGeneratedPower().equals(Complex.ZERO))
				generatorMask.setEntry(index, index, 1);
		}
        
        // Init power vector:
        p = p0(results);
        q = q0(results);
        p.setEntry(slackIndex, 0);
        q.setEntry(slackIndex, 0);
        
        p_max = maxPower(results, grid);
        
        // Init voltage vectors:
        e = new ArrayRealVector(dimension, 1);
        f = new ArrayRealVector(dimension, 0);
        
        // Init Lagrange multipliers:
        l_p = new ArrayRealVector(dimension, 0);
        l_q = new ArrayRealVector(dimension, 0);
	}
	
	public void loop()
	{
		int K = 1000;
		for(int k = 0; k < K; ++k)
		{
			// Get diagonal matrix versions of lambda:
			RealMatrix L_p = diag(l_p);
			RealMatrix L_q = diag(l_q);
			
			// e := (I + 2*gamma(Lambda_p*G - Lambda_q*B))e
			int dimension = e.getDimension();
			RealMatrix I = MatrixUtils.createRealIdentityMatrix(dimension);
			e = I.add( L_p.multiply(G).subtract(L_q.multiply(B)).scalarMultiply(2*gamma) ).operate(e);
			
			// f := (I + 2*gamma(Lambda_p*G - Lambda_q*B)f
			f = I.add( L_p.multiply(G).subtract(L_q.multiply(B)).scalarMultiply(2*gamma) ).operate(f);
			
			// p := p + gamma*(grad_p.f(x) - lambda_p)
			RealVector grad_p = gradF_p(p, q).subtract(l_p).mapMultiply(gamma);
			grad_p = generatorMask.operate(grad_p);
//			p = p.add( grad_p ); FIXME
			
			// q := q + gamma*(grad_q.f(x) - lambda_q)
			RealVector grad_q = gradF_q(p, q).subtract(l_q).mapMultiply(gamma);
			grad_q = generatorMask.operate(grad_q);
//			q = q.add( grad_q ); FIXME
			
			// lambda_p := lambda_p - gamma*(diag(e)Ge + diag(f)Gf + diag(f)Be - diag(e)Bf - p)
			l_p = l_p.subtract(
						(   diag(e).multiply(G).operate(e)
							.add(diag(f).multiply(G).operate(f))
							.add(diag(f).multiply(B).operate(e))
							.subtract(diag(e).multiply(B).operate(f))
							.subtract(p)
						).mapMultiply(gamma)
				  );
			
			// lambda_q := lambda_p - gamma*(diag(f)Ge - diag(e)Gf - diag(e)Be - diag(f)Bf - q)
			l_q = l_q.subtract(
					(   diag(f).multiply(G).operate(e)
						.subtract(diag(e).multiply(G).operate(f))
						.subtract(diag(e).multiply(B).operate(e))
						.subtract(diag(f).multiply(B).operate(f))
						.subtract(q)
					).mapMultiply(gamma)
			  );
			
			debug(k);
		}
	}

	private RealVector gradF_p(RealVector p, RealVector q)
	{
		return p.subtract(p_max).mapMultiply(-2);
	}

	private RealVector gradF_q(RealVector p, RealVector q)
	{
		return q.mapMultiply(-2);
	}

	private RealMatrix diag(RealVector v)
	{
		int dimension = v.getDimension();
		RealMatrix M = new Array2DRowRealMatrix(dimension, dimension);
		for (int i = 0; i < dimension; i++)
		{
			M.setEntry(i, i, v.getEntry(i));
		}
		return M;
	}
	
	
	//// Debug ////

//	protected double pSlack()
//	{
//		double sum = 0;
//		int dimension = p.getDimension();
//		for(int j = 0; j < dimension; ++j)
//		{
//			double theta_0j = Y.getEntry(slackIndex, j).getArgument();
//			double delta_j = varg.getEntry(j);
//			double cos = Math.cos(theta_0j+delta_j);
//			sum += vabs.getEntry(j)*Y.getEntry(slackIndex, j).abs()*cos;
//		}
//		return sum;
//	}
//	
//	protected double p_i(int i)
//	{
//		double sum = 0;
//		int dimension = p.getDimension();
//		for(int j = 0; j < dimension; ++j)
//		{
//			double theta_ij = Y.getEntry(i, j).getArgument();
//			double delta_j = varg.getEntry(j);
//			double delta_i = varg.getEntry(i);
//			double cos = cos(theta_ij-delta_i+delta_j);
//			sum += vabs.getEntry(i)*vabs.getEntry(j)*Y.getEntry(i, j).abs()*cos;
//		}
//		return sum;
//	}
//
//	protected double qSlack()
//	{
//		double sum = 0;
//		int dimension = p.getDimension();
//		for(int j = 0; j < dimension; ++j)
//		{
//			double theta_0j = Y.getEntry(slackIndex, j).getArgument();
//			double delta_j = varg.getEntry(j);
//			double sin = Math.sin(theta_0j+delta_j);
//			sum += vabs.getEntry(j)*Y.getEntry(slackIndex, j).abs()*sin;
//		}
//		return -sum;
//	}
//	
//	protected double q_i(int i)
//	{
//		double sum = 0;
//		int dimension = q.getDimension();
//		for(int j = 0; j < dimension; ++j)
//		{
//			double theta_ij = Y.getEntry(i, j).getArgument();
//			double delta_j = varg.getEntry(j);
//			double delta_i = varg.getEntry(i);
//			double sin = sin(theta_ij-delta_i+delta_j);
//			sum += vabs.getEntry(i)*vabs.getEntry(j)*Y.getEntry(i, j).abs()*sin;
//		}
//		return -sum;
//	}

	protected void debugHeader(AnalysisResults results)
	{
		Map<String, Integer> numbers = results.getBusNumbers();
		System.out.print("-,");
		
		int dimension = p.getDimension();
		debugLogNames(numbers, dimension, "p");
		debugLogNames(numbers, dimension, "q");
		debugLogNames(numbers, dimension, "e");
		debugLogNames(numbers, dimension, "f");
		debugLogNames(numbers, dimension, "l_p");
		debugLogNames(numbers, dimension, "l_q");
//		System.out.print("p_slack,q_slack,");
//		System.out.print("U(pq),L(pql)");
		System.out.println();
	}
	
	protected void debug(int k)
	{
		System.out.print(k);
		System.out.print(',');
		int dimension = p.getDimension();
		
		// Power:
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(p.getEntry(i)));
			System.out.print(',');
		}
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(q.getEntry(i)));
			System.out.print(',');
		}
		
		// Voltage:
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(e.getEntry(i)));
			System.out.print(',');
		}
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(f.getEntry(i)));
			System.out.print(',');
		}
		
		// Lambda:
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(l_p.getEntry(i)));
			System.out.print(',');
		}
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(l_q.getEntry(i)));
			System.out.print(',');
		}
		
		// Slack power:
//		System.out.print(format(pSlack()));
//		System.out.print(',');
//		System.out.print(format(qSlack()));
//		System.out.print(',');
		
		// Utility:
//		System.out.print(format(utility()));
//		System.out.print(',');
//		System.out.print(format(lagrangian()));
//		System.out.print(',');
		
		System.out.println();
	}
	
	public static String format(double d)
	{
		return ""+d;//String.format("%.4f", d);
	}
	
//	public static double TODO output utility and lagrandian, and perhaps slack power
}