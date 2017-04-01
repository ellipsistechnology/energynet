package ellipsis.energy.sandbox;

import static ellipsis.energy.sandbox.Sandbox004.generators;
import static ellipsis.energy.sandbox.Sandbox004.logBusses;
import static ellipsis.energy.sandbox.Sandbox004.maxPower;
import static ellipsis.energy.sandbox.Sandbox004.minPower;
import static ellipsis.energy.sandbox.Sandbox004.p0;
import static ellipsis.energy.sandbox.Sandbox004.q0;
import static ellipsis.energy.sandbox.GridConstants.*;
//import static ellipsis.energy.sandbox.Sandbox004.vArg0;
import static ellipsis.energy.sandbox.Sandbox007.analyse;
import static java.lang.Math.cos;
import static java.lang.Math.sin;

import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.linear.RealVector;

import com.mls.util.Util;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.calculation.LoadFlowAnalyser;
import ellipsis.energy.grid.Grid;

/**
 * Min/max Lagrange gradient method for OPF.
 * @author bmillar
 *
 */
public class Sandbox008
{	
    //// Variables ////
    
    // Parameters:
    protected Grid grid;
    protected FieldMatrix<Complex> Y;
    protected int slackIndex;
    protected Set<Integer> G; // The set of generator busses.
    protected RealVector s_min;
    protected RealVector s_max;
//	protected RealVector s;
	
	// Variables:
    protected RealVector p;
    protected RealVector q;
    protected RealVector vabs;
    protected RealVector varg;
    protected RealVector lambdaP;
    protected RealVector lambdaQ;
	
	// Step sizes:
    protected double pStepSize = 0.1;
    protected double qStepSize = 0.1;
    protected double vabsStepSize = 0.012; // 0.012 seems best
    protected double vargStepSize = 0.012;
    protected double lambdaPStepSize = 0.012;
    protected double lambdaQStepSize = 0.012;
	

	public static void main(String[] args)
	{
		new Sandbox008().run();
	}

	protected void run()
	{
		AnalysisResults results = init();
		loop();
		debugResults(results);
//		testCase();
	}
	
	protected void testCase()
	{
		grid.getSource(DG_3).setPowerOutput(6.486*BASE_POWER, 0.0092*BASE_POWER);
		grid.getSource(DG_4).setPowerOutput(6.467*BASE_POWER, 0.02*BASE_POWER);
        
        analyse(grid);
	}
    
    
    //// Setup ////
	
	protected AnalysisResults init()
	{
		initGrid();
		
		// Analyse grid to get parameters:
		LoadFlowAnalyser lfa = new LoadFlowAnalyser(grid);
        lfa.setBasePower(grid.getBasePower());
        lfa.setBaseVoltage(grid.getBaseVoltage());
        AnalysisResults results = lfa.analyse();
        slackIndex = results.getSlackIndex();
        
        logBusses(results);
        
        // Initialise variables:
        initVariables(results);
        
        debugResults(results);
        
        // Debug:
		debugHeader(results);
		
		return results;
	}

	protected void initVariables(AnalysisResults results)
	{
		// Set G:
        G = generators(grid, results.getBusNumbers());
        
        // Get Y:
        Y = results.getAdmittanceMatrix().Y;
//        logSensitivities(Y);
     
        // Init power vector:
        p = p0(results);
        q = q0(results);
        p.setEntry(slackIndex, 0);
        q.setEntry(slackIndex, 0);
        int dimension = p.getDimension();
        
        // Init voltage vectors:
        vabs = new ArrayRealVector(dimension, 1);
        varg = new ArrayRealVector(dimension, 0);
//        vabs = vAbs0(results);
//        varg = vArg0(results);
        
        // Set p_min and p_max:
        s_max = maxPower(results, grid);
        s_min = minPower(results, grid);
//        s = s_max.mapMultiply(0.8); // power output is set to 80% of DG power rating
        
        // Init lambdas:
		lambdaP = new ArrayRealVector(dimension, 0.0); // FIXME these need some more thought; their initial value is likely to dictate convergence and/or optimality.
		lambdaQ = new ArrayRealVector(dimension, 0.0);
	}

	public void initGrid()
	{
		// Copied from Sandbox007:
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
	
	
	//// Optimisation ////

//	private void testCase()
//	{
//		
//	}

	public void loop()
	{
		int K = 1000;
		for(int k = 0; k < K; ++k)
		{
if(k == 999)
	Util.nullop();
			optimiseStep();
			constraintStep();
			debug(k);
		}
	}

	private void optimiseStep()
	{
		// Gradient accent:
		RealVector gradP = gradP();
		RealVector gradQ = gradQ();
		RealVector gradVabs = gradVabs();
		RealVector gradVarg = gradVarg();
		RealVector nextP = p.add(gradP.mapMultiply(pStepSize));
		RealVector nextQ = q.add(gradQ.mapMultiply(qStepSize));
		RealVector nextVabs = vabs.add(gradVabs.mapMultiply(vabsStepSize));
		RealVector nextVarg = varg.add(gradVarg.mapMultiply(vargStepSize));
		
		// Bound p & q:
//		FIXME clampPQ(nextP, nextQ);
		
		p = nextP;
		q = nextQ;
		vabs = nextVabs;
		varg = nextVarg;
	}

	protected void clampPQ(RealVector nextP, RealVector nextQ)
	{
		for (int i = 0; i < nextP.getDimension(); i++)
		{
			if(!G.contains(i))
				continue;
			
			double s_min_i = s_min.getEntry(i);
			double s_max_i = s_max.getEntry(i);
			double p_i = nextP.getEntry(i);
			if(p_i < 0)
				p_i = 0;
			double q_i = nextQ.getEntry(i);
			Complex s_i = new Complex(p_i, q_i);
			
			double s_i_abs = s_i.abs();
			if(s_i_abs < s_min_i)
				s_i = s_i.multiply(s_min_i/s_i_abs);
			else if(s_i_abs > s_max_i)
				s_i = s_i.multiply(s_max_i/s_i_abs);
			
			nextP.setEntry(i, s_i.getReal());
			nextQ.setEntry(i, s_i.getImaginary());
		}
	}

	private RealVector gradP()
	{
		int dimension = p.getDimension();
		RealVector grad = new ArrayRealVector(dimension);
		for(int i = 0; i < dimension; ++i)
			grad.setEntry(i, gradP_i(i));
		
		return grad;
	}

	private double gradP_i(int i)
	{
		if(!G.contains(i))
			return 0;
		return gradF_p_i(i) - lambdaP.getEntry(i);
	}

	/**
	 * Utility function gradient W.R.T. p_i.
	 * @param i
	 * @return
	 */
	private double gradF_p_i(int i)
	{
		double p_i = p.getEntry(i);
		double p_max_i = s_max.getEntry(i);
		return -2*(p_i-p_max_i);
	}

	private RealVector gradQ()
	{
		int dimension = q.getDimension();
		RealVector grad = new ArrayRealVector(dimension);
		for(int i = 0; i < dimension; ++i)
			grad.setEntry(i, gradQ_i(i));
		
		return grad;
	}

	private double gradQ_i(int i)
	{
		if(!G.contains(i))
			return 0;
		return gradF_q_i(i) - lambdaQ.getEntry(i);
	}

	/**
	 * Utility function gradient W.R.T. q_i.
	 * @param i
	 * @return
	 */
	private double gradF_q_i(int i)
	{
		return -2*q.getEntry(i);
	}

	private RealVector gradVabs()
	{
		int dimension = vabs.getDimension();
		RealVector grad = new ArrayRealVector(dimension);
		for(int i = 0; i < dimension; ++i)
			grad.setEntry(i, gradVabs_i(i));
		
		return grad;
	}

	private double gradVabs_i(int i)
	{
		if(i == slackIndex)
			return 0;
		return  2*(1-vabs.getEntry(i))
				// lambdaP.getEntry(slackIndex)*gradVabs_i_p0(i) - lambdaQ.getEntry(slackIndex)*gradVabs_i_q0(i)
				+ lambdaP.getEntry(i)*gradVabs_i_pi(i) + gradVabs_i_pj(i)
				- lambdaQ.getEntry(i)*gradVabs_i_qi(i) - gradVabs_i_qj(i);
	}

//	private double gradVabs_i_p0(int i)
//	{
//		Complex Y_0i = Y.getEntry(slackIndex, i);
//		double varg_i = varg.getEntry(i);
//		return Y_0i.abs()*cos(Y_0i.getArgument() + varg_i);
//	}
//
//	private double gradVabs_i_q0(int i)
//	{
//		Complex Y_0i = Y.getEntry(slackIndex, i);
//		double varg_i = varg.getEntry(i);
//		return Y_0i.abs()*sin(Y_0i.getArgument() + varg_i);
//	}

	private double gradVabs_i_pi(int i)
	{
		double varg_i = varg.getEntry(i);
		double vabs_i = vabs.getEntry(i);
		Complex Y_ii = Y.getEntry(i, i);
		
		double sum = 0;
		int dimension = vabs.getDimension();
		for(int j = 0; j < dimension; ++j)
		{
			if(i == j)
				continue;
			double vabs_j = vabs.getEntry(j);
			double varg_j = varg.getEntry(j);
			Complex Y_ij = Y.getEntry(i, j);
			sum += vabs_j*Y_ij.abs()*cos(Y_ij.getArgument() - varg_i + varg_j);
		}
		return sum + 2*vabs_i*Y_ii.abs()*cos(Y_ii.getArgument());
	}

	private double gradVabs_i_pj(int i)
	{
		double varg_i = varg.getEntry(i);
		double sum = 0;
		int dimension = vabs.getDimension();
		for(int j = 0; j < dimension; ++j)
		{
			if(i == j || i == slackIndex)
				continue;
			double vabs_j = vabs.getEntry(j);
			double varg_j = varg.getEntry(j);
			Complex Y_ji = Y.getEntry(j, i);
			sum += lambdaP.getEntry(j)*vabs_j*Y_ji.abs()*cos(Y_ji.getArgument() - varg_j + varg_i);
		}
		
		return sum;
	}

	private double gradVabs_i_qi(int i)
	{
		double varg_i = varg.getEntry(i);
		double vabs_i = vabs.getEntry(i);
		Complex Y_ii = Y.getEntry(i, i);
		
		double sum = 0;
		int dimension = vabs.getDimension();
		for(int j = 0; j < dimension; ++j)
		{
			if(i == j)
				continue;
			double vabs_j = vabs.getEntry(j);
			double varg_j = varg.getEntry(j);
			Complex Y_ij = Y.getEntry(i, j);
			sum += vabs_j*Y_ij.abs()*sin(Y_ij.getArgument() - varg_i + varg_j);
		}
		return sum + 2*vabs_i*Y_ii.abs()*sin(Y_ii.getArgument());
	}

	private double gradVabs_i_qj(int i)
	{
		double varg_i = varg.getEntry(i);
		double sum = 0;
		int dimension = vabs.getDimension();
		for(int j = 0; j < dimension; ++j)
		{
			if(i == j || i == slackIndex)
				continue;
			double vabs_j = vabs.getEntry(j);
			double varg_j = varg.getEntry(j);
			Complex Y_ji = Y.getEntry(j, i);
			sum += lambdaQ.getEntry(j)*vabs_j*Y_ji.abs()*sin(Y_ji.getArgument() - varg_j + varg_i);
		}
		
		return sum;
	}

	private RealVector gradVarg()
	{
		int dimension = vabs.getDimension();
		RealVector grad = new ArrayRealVector(dimension);
		for(int i = 0; i < dimension; ++i)
			grad.setEntry(i, gradVarg_i(i));
		
		return grad;
	}

	private double gradVarg_i(int i)
	{
		if(i == slackIndex)
			return 0;
		return  -2*(varg.getEntry(i))
				// - lambdaP.getEntry(slackIndex)*gradVarg_i_p0(i) - lambdaQ.getEntry(slackIndex)*gradVarg_i_q0(i)
				+ lambdaP.getEntry(i)*gradVarg_i_pi(i) - gradVarg_i_pj(i) 
				+ lambdaQ.getEntry(i)*gradVarg_i_qi(i) - gradVarg_i_qj(i);
	}

//	private double gradVarg_i_p0(int i)
//	{
//		double vabs_i = vabs.getEntry(i);
//		Complex Y_0i = Y.getEntry(slackIndex, i);
//		double varg_i = varg.getEntry(i);
//		return vabs_i*Y_0i.abs()*sin(Y_0i.getArgument() + varg_i);
//	}
//
//	private double gradVarg_i_q0(int i)
//	{
//		double vabs_i = vabs.getEntry(i);
//		Complex Y_0i = Y.getEntry(slackIndex, i);
//		double varg_i = varg.getEntry(i);
//		return vabs_i*Y_0i.abs()*cos(Y_0i.getArgument() + varg_i);
//	}

	private double gradVarg_i_pi(int i)
	{
		double vabs_i = vabs.getEntry(i);
		double varg_i = varg.getEntry(i);
				
		double sum = 0;
		int dimension = varg.getDimension();
		for(int j = 0; j < dimension; ++j)
		{
			if(j == i)
				continue;
			double vabs_j = vabs.getEntry(j);
			double varg_j = varg.getEntry(j);
			Complex Y_ij = Y.getEntry(i, j);
			sum += vabs_i*vabs_j*Y_ij.abs()*sin(Y_ij.getArgument() - varg_i + varg_j);
		}
		
		return sum;
	}

	private double gradVarg_i_pj(int i)
	{
		double vabs_i = vabs.getEntry(i);
		double varg_i = varg.getEntry(i);
				
		double sum = 0;
		int dimension = varg.getDimension();
		for(int j = 0; j < dimension; ++j)
		{
			if(j == i)
				continue;
			double vabs_j = vabs.getEntry(j);
			double varg_j = varg.getEntry(j);
			Complex Y_ji = Y.getEntry(j, i);
			sum += lambdaP.getEntry(j)*vabs_i*vabs_j*Y_ji.abs()*sin(Y_ji.getArgument() - varg_j + varg_i);
		}
		
		return sum;
	}

	private double gradVarg_i_qi(int i)
	{
		double vabs_i = vabs.getEntry(i);
		double varg_i = varg.getEntry(i);
				
		double sum = 0;
		int dimension = varg.getDimension();
		for(int j = 0; j < dimension; ++j)
		{
			if(j == i)
				continue;
			double vabs_j = vabs.getEntry(j);
			double varg_j = varg.getEntry(j);
			Complex Y_ij = Y.getEntry(i, j);
			sum += vabs_i*vabs_j*Y_ij.abs()*cos(Y_ij.getArgument() - varg_i + varg_j);
		}
		
		return sum;
	}

	private double gradVarg_i_qj(int i)
	{
		double vabs_i = vabs.getEntry(i);
		double varg_i = varg.getEntry(i);
				
		double sum = 0;
		int dimension = varg.getDimension();
		for(int j = 0; j < dimension; ++j)
		{
			if(j == i)
				continue;
			double vabs_j = vabs.getEntry(j);
			double varg_j = varg.getEntry(j);
			Complex Y_ji = Y.getEntry(j, i);
			sum += lambdaQ.getEntry(j)*vabs_i*vabs_j*Y_ji.abs()*cos(Y_ji.getArgument() - varg_j + varg_i);
		}
		
		return sum;
	}

	private void constraintStep()
	{
		// Gradient decent:
		RealVector gradLambdaP = gradLambdaP();
		RealVector gradLambdaQ = gradLambdaQ();
		RealVector nextLambdaP = lambdaP.subtract(gradLambdaP.mapMultiply(lambdaPStepSize));
		RealVector nextLambdaQ = lambdaQ.subtract(gradLambdaQ.mapMultiply(lambdaQStepSize));

		lambdaP = nextLambdaP;
		lambdaQ = nextLambdaQ;
	}
	
	private RealVector gradLambdaP()
	{
		int dimension = vabs.getDimension();
		RealVector grad = new ArrayRealVector(dimension);
		for(int i = 0; i < dimension; ++i)
			grad.setEntry(i, gradLambdaP_i(i));
		
		return grad;
	}

	private double gradLambdaP_i(int i)
	{
if(i == slackIndex)
	return 0; // FIXME Include slack bus here for power balancing.
		double vabs_i = vabs.getEntry(i);
		double varg_i = varg.getEntry(i);
		
		double sum = 0;
		int dimension = varg.getDimension();
		for(int j = 0; j < dimension; ++j)
		{
			double vabs_j = vabs.getEntry(j);
			double varg_j = varg.getEntry(j);
			Complex Y_ij = Y.getEntry(i, j);
			sum += vabs_i*vabs_j*Y_ij.abs()*cos(Y_ij.getArgument() - varg_i + varg_j);
		}
		
		double p_i = i == slackIndex ? 0 : p.getEntry(i);
		return sum - p_i;
	}

	private RealVector gradLambdaQ()
	{
		int dimension = vabs.getDimension();
		RealVector grad = new ArrayRealVector(dimension);
		for(int i = 0; i < dimension; ++i)
			grad.setEntry(i, gradLambdaQ_i(i));
		
		return grad;
	}

	private double gradLambdaQ_i(int i)
	{
if(i == slackIndex)
	return 0; // FIXME
		double vabs_i = vabs.getEntry(i);
		double varg_i = varg.getEntry(i);
		
		double sum = 0;
		int dimension = varg.getDimension();
		for(int j = 0; j < dimension; ++j)
		{
			double vabs_j = vabs.getEntry(j);
			double varg_j = varg.getEntry(j);
			Complex Y_ij = Y.getEntry(i, j);
			sum += vabs_i*vabs_j*Y_ij.abs()*sin(Y_ij.getArgument() - varg_i + varg_j);
		}
		
		double q_i = i == slackIndex ? 0 : q.getEntry(i);
		return -(sum + q_i);
	}
	
	
	//// DEBUG ////

	protected double pSlack()
	{
		double sum = 0;
		int dimension = p.getDimension();
		for(int j = 0; j < dimension; ++j)
		{
			double theta_0j = Y.getEntry(slackIndex, j).getArgument();
			double delta_j = varg.getEntry(j);
			double cos = Math.cos(theta_0j+delta_j);
			sum += vabs.getEntry(j)*Y.getEntry(slackIndex, j).abs()*cos;
		}
		return sum;
	}
	
	protected double p_i(int i)
	{
		double sum = 0;
		int dimension = p.getDimension();
		for(int j = 0; j < dimension; ++j)
		{
			double theta_ij = Y.getEntry(i, j).getArgument();
			double delta_j = varg.getEntry(j);
			double delta_i = varg.getEntry(i);
			double cos = cos(theta_ij-delta_i+delta_j);
			sum += vabs.getEntry(i)*vabs.getEntry(j)*Y.getEntry(i, j).abs()*cos;
		}
		return sum;
	}

	protected double qSlack()
	{
		double sum = 0;
		int dimension = p.getDimension();
		for(int j = 0; j < dimension; ++j)
		{
			double theta_0j = Y.getEntry(slackIndex, j).getArgument();
			double delta_j = varg.getEntry(j);
			double sin = Math.sin(theta_0j+delta_j);
			sum += vabs.getEntry(j)*Y.getEntry(slackIndex, j).abs()*sin;
		}
		return -sum;
	}
	
	protected double q_i(int i)
	{
		double sum = 0;
		int dimension = q.getDimension();
		for(int j = 0; j < dimension; ++j)
		{
			double theta_ij = Y.getEntry(i, j).getArgument();
			double delta_j = varg.getEntry(j);
			double delta_i = varg.getEntry(i);
			double sin = sin(theta_ij-delta_i+delta_j);
			sum += vabs.getEntry(i)*vabs.getEntry(j)*Y.getEntry(i, j).abs()*sin;
		}
		return -sum;
	}
	
	protected void debugResults(AnalysisResults results)
	{
		// Header:
		Map<String, Integer> numbers = results.getBusNumbers();
		int dimension = p.getDimension();
		debugLogNames(numbers, dimension, "p (calculated)");
		debugLogNames(numbers, dimension, "q (calculated)");
		System.out.println();
		
		// Power:
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(p_i(i)));
			System.out.print(',');
		}
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(q_i(i)));
			System.out.print(',');
		}
		System.out.println();
	}

	protected void debugHeader(AnalysisResults results)
	{
		Map<String, Integer> numbers = results.getBusNumbers();
		System.out.print("-,");
		
		int dimension = p.getDimension();
		debugLogNames(numbers, dimension, "p");
		debugLogNames(numbers, dimension, "q");
		debugLogNames(numbers, dimension, "|v|");
		debugLogNames(numbers, dimension, "/_v");
		debugLogNames(numbers, dimension, "lambdaP");
		debugLogNames(numbers, dimension, "lambdaQ");
		System.out.print("p_slack,q_slack,");
		System.out.print("U(pq),L(pql)");
		System.out.println();
	}

	public static void debugLogNames(Map<String, Integer> numbers, int dimension, String type)
	{
		for(int i = 0; i < dimension; ++i)
		{
			for (String name : numbers.keySet())
			{
				int index = numbers.get(name);
				if(index == i)
				{
					System.out.print(name);
					System.out.print('(');
					System.out.print(type);
					System.out.print("),");
					break;
				}
			}
		}
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
			System.out.print(format(vabs.getEntry(i)));
			System.out.print(',');
		}
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(varg.getEntry(i)));
			System.out.print(',');
		}
		
		// Lambda:
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(lambdaP.getEntry(i)));
			System.out.print(',');
		}
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(lambdaQ.getEntry(i)));
			System.out.print(',');
		}
		
		// Slack power:
		System.out.print(format(pSlack()));
		System.out.print(',');
		System.out.print(format(qSlack()));
		
		// Utility:
		System.out.print(',');
		System.out.print(format(utility()));
		System.out.print(',');
		System.out.print(format(lagrangian()));
		
		System.out.println();
	}
	
	private double utility()
	{
		double sum = 0;
		for (int i : G)
		{
			double p_i = p.getEntry(i);
			double q_i = q.getEntry(i);
			double s_max_i = s_max.getEntry(i);
			double pDiff = p_i-s_max_i;
			sum += -pDiff*pDiff - q_i*q_i;
		}
		return sum;
	}
	
	private double lagrangian()
	{
		double sumP = 0;
		double sumQ = 0;
		int dimension = p.getDimension();
		for(int i = 0; i < dimension; ++i)
		{
			double sum_j = 0;
			for(int j = 0; j < dimension; ++j)
			{
				Complex Y_ij = Y.getEntry(i, j);
				double vabs_i = vabs.getEntry(i);
				double vabs_j = vabs.getEntry(j);
				double varg_i = varg.getEntry(i);
				double varg_j = varg.getEntry(j);
				sum_j += vabs_i*vabs_j*Y_ij.abs()*cos(Y_ij.getArgument() - varg_i + varg_j);
			}
			double p_i = i == slackIndex ? 0 : p.getEntry(i);
			sumP += lambdaP.getEntry(i)*(sum_j - p_i);
		}
		for(int i = 0; i < dimension; ++i)
		{
			double sum_j = 0;
			for(int j = 0; j < dimension; ++j)
			{
				Complex Y_ij = Y.getEntry(i, j);
				double vabs_i = vabs.getEntry(i);
				double vabs_j = vabs.getEntry(j);
				double varg_i = varg.getEntry(i);
				double varg_j = varg.getEntry(j);
				sum_j += vabs_i*vabs_j*Y_ij.abs()*sin(Y_ij.getArgument() - varg_i + varg_j);
			}
			double q_i = i == slackIndex ? 0 : q.getEntry(i);
			sumP += lambdaQ.getEntry(i)*(sum_j + q_i);
		}
		return utility() + sumP - sumQ;
	}
	
	public static String format(double d)
	{
		return String.format("%.4f", d);
	}
}