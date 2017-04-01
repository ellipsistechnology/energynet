package ellipsis.energy.sandbox;

import static ellipsis.energy.test.GridlessADPIEEE13BusGridTest.makeDG;

import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import com.mls.util.Util;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.calculation.LoadFlowAnalyser;
import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.Grid;
import ellipsis.energy.grid.Source;
import ellipsis.energy.test.IEEE13BusGrid;

public class Sandbox004
{
	/**
	 * @param args
	 */
	public static void main(String[] args)
	{
		IEEE13BusGrid grid = new IEEE13BusGrid();

		grid.bus680.addChild(makeDG("DG680"));
		grid.bus675.addChild(makeDG("DG675"));		
		grid.bus684.addChild(makeDG("DG684"));
		grid.bus645.addChild(makeDG("DG645"));
		grid.bus611.addChild(makeDG("DG611"));
		grid.bus646.addChild(makeDG("DG646"));
		
		new Sandbox004().run(grid);
	}

	private int k;
	private static final int K = 1000;
	private Grid grid;
	protected void run(Grid grid)
	{
		this.grid = grid;
		init();
		for(k = 0; k < K; ++k)
		{
			step();
			debug(k);
		}
	}
	
	// Parameters:
	private FieldMatrix<Complex> Y;
	private RealMatrix sensitivities;
	private int slackIndex;
	private Set<Integer> G; // The set of generator busses.
	private RealVector s_min;
	private RealVector s_max;
//	private RealVector s;
	
	// Variables:
	private RealVector p;
	private RealVector q;
	private RealVector vAbs;
	private RealVector vArg;
	private double lambdaP;
	private double lambdaQ;
	private double[] lambdaS;
	
	// Step sizes:
	private double pStepSize = 0.1;
	private double qStepSize = 0.1;
	private double lambdaPStepSize = 0.1;
	private double lambdaQStepSize = 0.1;
//	private double lambdaSStepSize = 0.1;
	
	
	//// Initialisation ////
	
	protected void init()
	{
		// Analyse grid to get parameters and initial variables:
		LoadFlowAnalyser lfa = new LoadFlowAnalyser(grid);
        lfa.setBasePower(grid.getBasePower());
        lfa.setBaseVoltage(grid.getBaseVoltage());
        AnalysisResults results = lfa.analyse();
        slackIndex = results.getSlackIndex();
        
        logBusses(results);
        
        // Get sensitivities:
        sensitivities = results.sensitivities();
        
        // Get Y:
        Y = results.getAdmittanceMatrix().Y;
        logSensitivities(Y);
        
        // Set G:
        G = generators(grid, results.getBusNumbers());
        
        // Init power vector:
        p = p0(results);
        q = q0(results);
        
        // Set p_min and p_max:
        s_max = maxPower(results, grid);
        s_min = minPower(results, grid);
//        s = s_max.mapMultiply(0.8); // power output is set to 80% of DG power rating
        
        // Init voltage vectors:
        vAbs = vAbs0(results);
        vArg = vArg0(results);
        
        // Init lambdas:
        lambdaP = 0; // FIXME these need some more thought; their initial value is likely to dictate convergence and/or optimality.
        lambdaQ = 0;
        int dimension = p.getDimension();
		lambdaS = new double[dimension];
        for(int i = 0; i < dimension; ++i)
        	lambdaS[i] = 0;
        
        // Debug:
		debugHeader(results);
	}

	public static RealVector maxPower(AnalysisResults results, Grid grid)
	{
		Map<String, Integer> numbers = results.getBusNumbers();
		RealVector powers = new ArrayRealVector(numbers.size());
		for (String name : numbers.keySet())
		{
			int index = numbers.get(name);
			Bus bus = grid.getBus(name);
			Collection<Source> sources = bus.getSources();
			if(sources.size() > 0)
			{
				Source s = sources.iterator().next();
				powers.setEntry(index, s.getPmax()/grid.getBasePower());
			}
		}
		return powers;
	}

	public static RealVector minPower(AnalysisResults results, Grid grid)
	{
		Map<String, Integer> numbers = results.getBusNumbers();
		RealVector powers = new ArrayRealVector(numbers.size());
		for (String name : numbers.keySet())
		{
			int index = numbers.get(name);
			Bus bus = grid.getBus(name);
			Collection<Source> sources = bus.getSources();
			if(sources.size() > 0)
			{
				Source s = sources.iterator().next();
				powers.setEntry(index, s.getPmin()/grid.getBasePower());
			}
		}
		return powers;
	}

	public static Set<Integer> generators(Grid grid, Map<String, Integer> numbers)
	{
		Set<Integer> generatorBusses = new HashSet<>();
		for (String name : numbers.keySet())
		{
			int i = numbers.get(name);
			Bus bus = grid.getBus(name);
			if(!bus.getGeneratedPower().equals(Complex.ZERO))
				generatorBusses.add(i);
		}
		return generatorBusses;
	}

	public static RealVector vAbs0(AnalysisResults results)
	{
		Map<String, Integer> numbers = results.getBusNumbers();
		double[] v0 = new double[numbers.size()];
		
		for (String name : numbers.keySet())
		{
			int i = numbers.get(name);
			v0[i] = results.getBusVoltage(name).abs();
		}
		
		return new ArrayRealVector(v0);
	}
	
	public static RealVector vArg0(AnalysisResults results)
	{
		Map<String, Integer> numbers = results.getBusNumbers();
		double[] v0 = new double[numbers.size()];
		
		for (String name : numbers.keySet())
		{
			int i = numbers.get(name);
			v0[i] = results.getBusVoltage(name).getArgument();
		}
		
		return new ArrayRealVector(v0);
	}

	public static RealVector p0(AnalysisResults results)
	{
		Map<String, Integer> numbers = results.getBusNumbers();
		double[] p0 = new double[numbers.size()];
		
		for (String name : numbers.keySet())
		{
			int i = numbers.get(name);
			p0[i] = results.getBusPower(name).getReal();
		}
		
		return new ArrayRealVector(p0);
	}

	public static RealVector q0(AnalysisResults results)
	{
		Map<String, Integer> numbers = results.getBusNumbers();
		double[] q0 = new double[numbers.size()];
		
		for (String name : numbers.keySet())
		{
			int i = numbers.get(name);
			q0[i] = results.getBusPower(name).getImaginary();
		}
		
		return new ArrayRealVector(q0);
	}
	
	
	//// Iteration ////

	protected void step()
	{
if(k == 400)
	Util.nullop();
for(int x = 0; x < 100; ++x) // FIXME
		nextPQ();
		nextLambda();
	}

	protected void nextV(RealVector nextP, RealVector nextQ)
	{
		int dimension = vAbs.getDimension();
		for(int i = 0; i < dimension; ++i)
		{
			if(i == slackIndex)
				continue;
			vAbs.setEntry(i, nextVAbs(i, dimension, nextP, nextQ));
			vArg.setEntry(i, nextVArg(i, dimension, nextP, nextQ));
		}
	}

	protected double nextVAbs(int i, int dimension, RealVector nextP, RealVector nextQ)
	{
		// Slack-free i:
		int i_noSlack = i > slackIndex ? i-1 : i;
		
		double sum = 0;
		for(int j = 0; j < dimension; ++j)
		{
			if(j == slackIndex)
				continue;
			
			// Slack-free j:
			int j_noSlack = j > slackIndex ? j-1 : j;
			
			double nextP_j = nextP.getEntry(j);
			double nextQ_j = nextQ.getEntry(j);
			double p_j = p.getEntry(j);
			double q_j = q.getEntry(j);
			double LambdaVAbsP_ij = sensitivities.getEntry(i_noSlack+dimension-1, j_noSlack);
			double LambdaVAbsQ_ij = sensitivities.getEntry(i_noSlack+dimension-1, j_noSlack+dimension-1);
			sum += LambdaVAbsP_ij*(nextP_j-p_j) + LambdaVAbsQ_ij*(nextQ_j-q_j);
		}
		return vAbs.getEntry(i) + sum;
	}

	protected double nextVArg(int i, int dimension, RealVector nextP, RealVector nextQ)
	{
		// Slack-free i:
		int i_noSlack = i > slackIndex ? i-1 : i;
		
		double sum = 0;
		for(int j = 0; j < dimension; ++j)
		{
			if(j == slackIndex)
				continue;
			
			// Slack-free j:
			int j_noSlack = j > slackIndex ? j-1 : j;
			
			double nextP_j = nextP.getEntry(j);
			double nextQ_j = nextQ.getEntry(j);
			double p_j = p.getEntry(j);
			double q_j = q.getEntry(j);
			double LambdaVArgP_ij = sensitivities.getEntry(i_noSlack, j_noSlack);
			double LambdaVArgQ_ij = sensitivities.getEntry(i_noSlack, j_noSlack+dimension-1);
			sum += LambdaVArgP_ij*(nextP_j-p_j) + LambdaVArgQ_ij*(nextQ_j-q_j);
		}
		return vArg.getEntry(i) + sum;
	}

	/**
	 * Gradient decent step to maximise augmented cost,
	 * utilising the current estimate of lambda.
	 */
	protected void nextPQ()
	{
		// p = p + eta*gradP:
		RealVector nextP = p.add(gradP().mapMultiply(pStepSize));
		RealVector nextQ = q.add(gradQ().mapMultiply(qStepSize));
		
		// p_min <= p <= p_max:
		clampPQ(nextP, nextQ);
		
		nextV(nextP, nextQ);
		p = nextP;
		q = nextQ;
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

	protected RealVector gradP()
	{
		int dimension = p.getDimension();
		RealVector grad = new ArrayRealVector(dimension);
		for(int i = 0; i < dimension; ++i)
			grad.setEntry(i, gradP_i(i));
		
		return grad;
	}

	protected RealVector gradQ()
	{
		int dimension = q.getDimension();
		RealVector grad = new ArrayRealVector(dimension);
		for(int i = 0; i < dimension; ++i)
			grad.setEntry(i, gradQ_i(i));
		
		return grad;
	}

	protected double gradP_i(int i)
	{
		// If i is not a DG bus, return 0:
		if(!G.contains(i))
			return 0;
		
		return gradUtilityP(i) + lambdaP*pConstraintP(i) - lambdaQ*qConstraintP(i);// FIXME + sConstraintP(i);
	}

	public double sConstraintP(int i)
	{
		return 2*lambdaS[i]*p.getEntry(i);
	}

	protected double gradQ_i(int i)
	{
		// If i is not a DG bus, return 0:
		if(!G.contains(i))
			return 0;
		
		return gradUtilityQ(i) + lambdaP*pConstraintQ(i) - lambdaQ*qConstraintQ(i);// FIXME + sConstraintQ(i);
	}

	public double sConstraintQ(int i)
	{
		return 2*lambdaS[i]*q.getEntry(i);
	}

	/**
	 * Quadratic utility function: -(p-p_max)^2
	 * @param i
	 * @return
	 */
	protected double gradUtilityP(int i)
	{
		double p_i = p.getEntry(i);
		double p_max_i = s_max.getEntry(i);
		return -2*(p_i-p_max_i);
	}

	/**
	 * Quadratic utility function: -q^2
	 * @param i
	 * @return
	 */
	protected double gradUtilityQ(int i)
	{
		double q_i = q.getEntry(i);
		return -2*(q_i-s_max.getEntry(i));
	}

	protected double pConstraintP(int i)
	{
		double sum = 0;
		int dimension = p.getDimension();
		for(int j = 0; j < dimension; ++j)
		{
			// Skip j = slack:
			if(j == slackIndex)
				continue;
			
			sum += pConstraintP_j(i, j, dimension);
		}
		return sum;
	}

	protected double pConstraintP_j(int i, int j, int dimension)
	{
		// Slack-free j:
		int j_noSlack = j > slackIndex ? j-1 : j;
		int i_noSlack = i > slackIndex ? i-1 : i;
		
		// Arg terms:
		double theta_0j = Y.getEntry(slackIndex, j).getArgument();
		double delta_j = vArg.getEntry(j);
		double LambdaVargP_ji = sensitivities.getEntry(j_noSlack, i_noSlack);
		double cos = Math.cos(theta_0j+delta_j);
		double sin = Math.sin(theta_0j+delta_j);

		// Abs terms:
		double Y_0j = Y.getEntry(slackIndex, j).abs();
		double LambdaVabsP_ji = sensitivities.getEntry(j_noSlack+dimension-1, i_noSlack); // -1 because dimension includes slack
		double v_j = vAbs.getEntry(j);
		
		return Y_0j*(LambdaVabsP_ji*cos - v_j*LambdaVargP_ji*sin);
	}

	protected double pConstraintQ(int i)
	{
		double sum = 0;
		int dimension = p.getDimension();
		for(int j = 0; j < dimension; ++j)
		{
			// Skip j = slack:
			if(j == slackIndex)
				continue;
			
			sum += pConstraintQ_j(i, j, dimension);
		}
		return sum;
	}

	protected double pConstraintQ_j(int i, int j, int dimension)
	{
		// Slack-free j:
		int j_noSlack = j > slackIndex ? j-1 : j;
		int i_noSlack = i > slackIndex ? i-1 : i;
		
		// Arg terms:
		double theta_0j = Y.getEntry(slackIndex, j).getArgument();
		double delta_j = vArg.getEntry(j);
		double LambdaVargQ_ji = sensitivities.getEntry(j_noSlack, i_noSlack+dimension-1);
		double cos = Math.cos(theta_0j+delta_j);
		double sin = Math.sin(theta_0j+delta_j);

		// Abs terms:
		double Y_0j = Y.getEntry(slackIndex, j).abs();
		double LambdaVabsQ_ji = sensitivities.getEntry(j_noSlack+dimension-1, i_noSlack+dimension-1); // -1 because dimension includes slack
		double v_j = vAbs.getEntry(j);
		
		return Y_0j*(LambdaVabsQ_ji*cos - v_j*LambdaVargQ_ji*sin);
	}

	protected double qConstraintP(int i)
	{
		double sum = 0;
		int dimension = p.getDimension();
		for(int j = 0; j < dimension; ++j)
		{
			// Skip j = slack:
			if(j == slackIndex)
				continue;

			sum += qConstraintP_j(i, j, dimension);
		}
		return sum;
	}

	protected double qConstraintP_j(int i, int j, int dimension)
	{
		// Slack-free j:
		int j_noSlack = j > slackIndex ? j-1 : j;
		int i_noSlack = i > slackIndex ? i-1 : i;
		
		// Arg terms:
		double theta_0j = Y.getEntry(slackIndex, j).getArgument();
		double delta_j = vArg.getEntry(j);
		double LambdaVargP_ji = sensitivities.getEntry(j_noSlack, i_noSlack);
		double cos = Math.cos(theta_0j+delta_j);
		double sin = Math.sin(theta_0j+delta_j);

		// Abs terms:
		double Y_0j = Y.getEntry(slackIndex, j).abs();
		double LambdaVabsP_ji = sensitivities.getEntry(j_noSlack+dimension-1, i_noSlack); // -1 because dimension includes slack
		double v_j = vAbs.getEntry(j);
		
		return Y_0j*(LambdaVabsP_ji*sin + v_j*LambdaVargP_ji*cos);
	}

	protected double qConstraintQ(int i)
	{
		double sum = 0;
		int dimension = p.getDimension();
		for(int j = 0; j < dimension; ++j)
		{
			// Skip j = slack:
			if(j == slackIndex)
				continue;

			sum += qConstraintQ_j(i, j, dimension);
		}
		return sum;
	}

	protected double qConstraintQ_j(int i, int j, int dimension)
	{
		// Slack-free j:
		int j_noSlack = j > slackIndex ? j-1 : j;
		int i_noSlack = i > slackIndex ? i-1 : i;
		
		// Arg terms:
		double theta_0j = Y.getEntry(slackIndex, j).getArgument();
		double delta_j = vArg.getEntry(j);
		double LambdaVargQ_ji = sensitivities.getEntry(j_noSlack, i_noSlack+dimension-1);
		double cos = Math.cos(theta_0j+delta_j);
		double sin = Math.sin(theta_0j+delta_j);

		// Abs terms:
		double Y_0j = Y.getEntry(slackIndex, j).abs();
		double LambdaVabsQ_ji = sensitivities.getEntry(j_noSlack+dimension-1, i_noSlack+dimension-1); // -1 because dimension includes slack
		double v_j = vAbs.getEntry(j);
		
		return Y_0j*(LambdaVabsQ_ji*sin + v_j*LambdaVargQ_ji*cos);
	}

	/**
	 * Gradient decent step to minimise augmented cost,
	 * utilising the current estimate of p.
	 */
	protected void nextLambda()
	{
		lambdaP -= pSlack()*lambdaPStepSize;
		lambdaQ -= qSlack()*lambdaQStepSize;
		
		// FIXME Update each lambdaS:
//		int dimension = p.getDimension();
//		for(int i = 0; i < dimension; ++i)
//		{
//			double s_i = s.getEntry(i);
//			double q_i = q.getEntry(i);
//			double p_i = p.getEntry(i);
//			lambdaS[i] -= (p_i*p_i + q_i*q_i - s_i*s_i)*lambdaSStepSize;
//		}
	}

	protected double pSlack()
	{
		double sum = 0;
		int dimension = p.getDimension();
		for(int j = 0; j < dimension; ++j)
		{
			double theta_0j = Y.getEntry(slackIndex, j).getArgument();
			double delta_j = vArg.getEntry(j);
			double cos = Math.cos(theta_0j+delta_j);
			sum += vAbs.getEntry(j)*Y.getEntry(slackIndex, j).abs()*cos;
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
			double delta_j = vArg.getEntry(j);
			double sin = Math.sin(theta_0j+delta_j);
			sum += vAbs.getEntry(j)*Y.getEntry(slackIndex, j).abs()*sin;
		}
		return -sum;
	}
	
	
	//// Testing ////
	
	protected void debugHeader(AnalysisResults results)
	{
		Map<String, Integer> numbers = results.getBusNumbers();
		System.out.print("-,");
		
		int dimension = p.getDimension();
		debugLogNames(numbers, dimension, "p");
		debugLogNames(numbers, dimension, "q");
		debugLogNames(numbers, dimension, "|v|");
		debugLogNames(numbers, dimension, "/_v");
		System.out.print("p_slack,q_slack,");
		System.out.print("U(pq),f(pql)");
		System.out.println();
	}

	public void debugLogNames(Map<String, Integer> numbers, int dimension, String type)
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
			System.out.print(format(vAbs.getEntry(i)));
			System.out.print(',');
		}
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(vArg.getEntry(i)));
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
		return utility() + lambdaP*pSlack() - lambdaQ*(qSlack());
	}
	
	private String format(double d)
	{
		return String.format("%.4f", d);
	}

	protected void logSensitivities(FieldMatrix<Complex> Y)
	{
		System.out.println("Slack index = "+slackIndex);
		Complex Y_00 = Y.getEntry(slackIndex, slackIndex);
		System.out.println("Y_slack-slack = "+Y_00);
		System.out.println("Y_slack-slack = "+Y_00.abs()+"/_"+Y_00.getArgument());
	}

	public static void logBusses(AnalysisResults results)
	{
		System.out.println("index,bus,|v|,arg,p,q");
		Map<String, Integer> numbers = results.getBusNumbers();
		for (String name : numbers.keySet())
		{
			int index = numbers.get(name);
			System.out.print(index);
			System.out.print(',');
			System.out.print(name);
			System.out.print(',');
			Complex v = results.getBusVoltage(name);
			System.out.print(v.abs());
			System.out.print(',');
			System.out.print(v.getArgument());
			System.out.print(',');
			Complex s = results.getBusPower(name);
			System.out.print(s.getReal());
			System.out.print(',');
			System.out.print(s.getImaginary());
			System.out.println();
		}
	}
}