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
import static ellipsis.energy.sandbox.Sandbox004.maxPower;
import static ellipsis.energy.sandbox.Sandbox004.p0;
import static ellipsis.energy.sandbox.Sandbox004.q0;
import static ellipsis.energy.sandbox.Sandbox008.debugLogNames;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.FieldMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.calculation.LoadFlowAnalyser;
import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.Grid;
import ellipsis.energy.grid.Source;

/**
 * Distributed version of {@link Sandbox013}.
 * @author bmillar
 *
 */
public class Sandbox014
{
	public static void main(String[] args)
	{
		new Sandbox014().run();
	}
	
	
	//// Lambda expression interfaces ////
	
	public static interface IndexedFunction
	{
		double value(int i);
	}
	
	public static interface IndexedBooleanFunction
	{
		boolean value(int i);
	}
	
	public static interface NeighbourFunction
	{
		double value(Agent neighbour);
	}
	
	
	//// Global variables ////
	
	protected RealVector p, q, e, f, lp, lq;
	protected RealMatrix G, B;
	protected RealVector p_max;
	private int slackIndex;
	protected RealVector generatorMask;
	protected Set<Agent> agents;
	protected Grid grid;
	protected double c = 0.001;
	
	protected void run()
	{
		AnalysisResults results = init();
		loop();
		testResults(results);
	}
	
	protected void testResults(AnalysisResults results)
	{
		// Set values to results:
		Map<String, Integer> busNumbers = results.getBusNumbers();
		for (String busName : busNumbers.keySet())
		{
			Bus bus = grid.getBus(busName);
			Source dg = null;
			for (Source source : bus.getSources())
			{
				dg = source;
				break;
			}
			if(dg != null)
			{
				int i = busNumbers.get(busName);
				dg.setPowerOutput(BASE_POWER*p(i), BASE_POWER*q(i));
			}
		}
//        grid.getSource(DG_3).setPowerOutput(BASE_POWER*6.5, BASE_POWER*0);
//        grid.getSource(DG_4).setPowerOutput(BASE_POWER*6.5, BASE_POWER*0);
		
		// Analyse grid to get parameters:
		LoadFlowAnalyser lfa = new LoadFlowAnalyser(grid);
        lfa.setBasePower(grid.getBasePower());
        lfa.setBaseVoltage(grid.getBaseVoltage());
        AnalysisResults results2 = lfa.analyse();
        
        System.out.println("\nFinal Values:");
        System.out.println("index,bus,e,f,p,q");
        for (String busName : busNumbers.keySet())
		{
			int i = busNumbers.get(busName);
        	System.out.println(i+","+busName+','+format(e(i))+','+format(f(i))+','+format(p(i))+','+format(q(i)));
        }
        System.out.println("\nAnalysis Results:");
        logBusses(results2);
	}
	
	
	//// Agents ////
	
	public static class Agent
	{
		public double p, q, e, f; // primal variables
		double p_new, q_new, e_new, f_new;
		public double p_max;
		public double lp, lq; // dual variables
		double lp_new, lq_new;
		public double gp, gq; // derived values
		public double c = 0.001;
		public boolean isSlack;
		public boolean isGenerator;
		public int index;
		private Map<Agent, Complex> neighbours = new HashMap<>();
		public double G, B; // Self admittances.
		
		public Agent(int i)
		{
			this.index = i;
		}
		
		public void addNeighbour(Agent neighbour, Complex admittance)
		{
			neighbours.put(neighbour, admittance);
		}

		public void stepPrimal()
		{
			if(isSlack)
				return;
			
			// Find x step direction:
			double dp = isGenerator ? gradL_p() : 0;
			double dq = isGenerator ? gradL_q() : 0;
			double de = isSlack ? 0 : gradL_e();
			double df = isSlack ? 0 : gradL_f();
			
			// Check stopping criterion:
//						double gradNorm = dp.getL1Norm() + dq.getL1Norm() + de.getL1Norm() + df.getL1Norm();
//						if(gradNorm < 10e-6)
//							break;
			
			// Backtrack to find best step size:
			double t = 0.05; // Step size
			double alpha = 0.6; // Step update size.
			double beta = 0.5;
//			double L = lagrangian(p, q, e, f, lp, lq);
//			double L_new;
			double step = alpha*(dp*dp + dq*dq + de*de + df*df);
			double minImprovement;
//int x = 0;
//			boolean debugLoop = false;
//			do
//			{
				p_new = p - dp*t;
				q_new = q - dq*t;
				e_new = e - de*t;
				f_new = f - df*t;
				
//				L_new = lagrangian(p_new, q_new, e_new, f_new, lp, lq);
//				minImprovement = -t*step;
				
//System.out.println(t+","+L_new+","+minImprovement);
				
//				t = t*beta;
//if(x++ > 100)
//	break;
//			} while(L_new - L > minImprovement || debugLoop);
			
			// Update x:
//			p = p_new;
//			q = q_new;
//			e = e_new;
//			f = f_new;
			
			// Update cost function values:
//			gp = gp();
//			gq = gq();
				
//			stepDual();
		}

		/**
		 * The dual step is completely local.
		 */
		protected void stepDual()
		{//FIXME think i need to skip slack here
			if(isSlack)
				return;
			
			gp = gp();
			gq = gq();
			
			// Maximise (l := l - c*g(x), ref http://en.wikipedia.org/wiki/Augmented_Lagrangian_method
			//           and Bertsekas' book):
			double dlp = isSlack ? 0 : gp;
			lp = lp + dlp*c;
			double dlq = isSlack ? 0 : gq;
			lq = lq + dlq*c;
		}
		
		/**
		 * Synchronous update to new values.
		 * This is only called at the end of an iteration, thus ensuring that
		 * the iterations are executed synchronously by the agents; i.e. only
		 * values from the previous iteration are used to calculate the new values.
		 */
		public void updatePrimal()
		{
			if(isSlack)
				return;
			
			p = p_new;
			q = q_new;
			e = e_new;
			f = f_new;
		}
		
		private double gp()
		{
			double sum = 
					0;
					//sum(neighbours, agent -> e*agent.e*G(agent) + f*agent.f*G(agent) + f*agent.e*B(agent) - e*agent.f*B(agent));
			for (Agent agent : neighbours.keySet())
			{
				sum += e*agent.e*G(agent) + f*agent.f*G(agent) + f*agent.e*B(agent) - e*agent.f*B(agent);
			}
			sum += e*e*G + f*f*G; // self
			return 	sum	- p;
		}
		
		private double gq()
		{
			double sum1 = sum(neighbours, agent -> f*agent.e*G(agent) - e*agent.f*G(agent) - e*agent.e*B(agent) - f*agent.f*B(agent));
			double sum = sum1
					- e*e*B -f*f*B;
			return 
				sum
				- q;
		}

		protected double gradL_p()
		{
			return 	gradF_p() - lp

					// Penalty:
					-c*gp()
					;
		}
		
		/**
		 * Assuming F_i(p_i) = 0.5*(p_i - p_max)^2
		 * => dF_i/dp_i = p_i - p_max
		 * @return
		 */
		private double gradF_p()
		{
			return p - p_max;
		}

		protected double gradL_q()
		{
//System.out.println(gradF_q());
//System.out.println(lq);
//System.out.println(c*gq());
			return gradF_q() - lq

					// Penalty:
					-c*gq()
					;
		}
		
		/**
		 * Assuming F_i(q_i) = 0.5*q_i^2
		 * => dF_i/dp_i = q_i
		 * @param i
		 * @return
		 */
		private double gradF_q()
		{
			return q;
		}

		protected double gradL_e()
		{
			double constraintTerms = 
						 sum(neighbours, agent -> agent.lp*(agent.e*G(agent) + agent.f*B(agent)))
						+sum(neighbours, agent -> agent.lq*(agent.f*G(agent) - agent.e*B(agent)))
						+lp*(sum(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent))
						     +2*G*e)
						+lq*(sum(neighbours, agent -> -agent.f*G(agent) - agent.e*B(agent))
						     -2*B*e);

			double pt1 = c*sum(neighbours, agent -> agent.isSlack ? 0 : 
				agent.gp*(agent.e*G(agent) + agent.f*B(agent))
				+agent.gq*(agent.f*G(agent) - agent.e*B(agent)));
			double pt2 = c*gp*(sum(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent)) 
					           +2*e*G);
			double pt3 = c*gq*(sum(neighbours, agent -> -agent.f*G(agent) - agent.e*B(agent)) 
					           -2*e*B);
			double penaltyTerms = pt1
								+
								(isSlack ? 0 : 
								(
									+pt2
									+pt3
								));
			return
				constraintTerms
				+penaltyTerms
				;
		}

		protected double gradL_f()
		{
			double constraintTerms = 
						 sum(neighbours, agent -> agent.lp*(agent.f*G(agent) - agent.e*B(agent)))
						+sum(neighbours, agent -> agent.lq*(-agent.e*G(agent) - agent.f*B(agent)))
						+lp*(sum(neighbours, agent -> agent.f*G(agent) + agent.e*B(agent))
						     +2*G*f)
						+lq*(sum(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent))
						     -2*B*f);

			double pt1 = c*sum(neighbours, agent -> agent.isSlack ? 0 : 
							   agent.gp*(agent.f*G(agent) - agent.e*B(agent))
							   +agent.gq*(-agent.e*G(agent) - agent.f*B(agent)));
			double pt2 = c*gp*(sum(neighbours, agent -> agent.f*G(agent) + agent.e*B(agent)) 
				               +2*f*G);
			double pt3 = c*gq*(sum(neighbours, agent -> agent.e*G(agent) - agent.f*B(agent)) 
				               -2*f*B);
			double penaltyTerms = pt1
								+
								(isSlack ? 0 : 
								(
									+pt2
									+pt3
								));
			return
				constraintTerms
				+penaltyTerms
				;
		}

		private double G(Agent agent)
		{
			return neighbours.get(agent).getReal();
		}

		private double B(Agent agent)
		{
			return neighbours.get(agent).getImaginary();
		}

		public void init()
		{
			gp = gp();
			gq = gq();
		}
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
        
        logBusses(results);
        
        // Initialise variables:
        initVariables(results);
        
        // Create agents:
        initAgents();
        
        // Debug:
		debugHeader(results);
		
		return results;
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
        
        // Prepare generator mask (1 if generator, 0 otherwise):
        generatorMask = new ArrayRealVector(dimension);
        for (String name : results.getBusNumbers().keySet())
		{
			int index = results.getBusNumbers().get(name);
			Bus bus = grid.getBus(name);
			if(!bus.getGeneratedPower().equals(Complex.ZERO))
				generatorMask.setEntry(index, 1);
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
        lp = new ArrayRealVector(dimension, 0);
        lq = new ArrayRealVector(dimension, 0);
	}
	
	protected void initAgents()
	{
		// Make agents:
		agents = new HashSet<>();
		int dimension = p.getDimension();
		for(int i = 0; i < dimension; ++i)
		{
			Agent agent = new Agent(i);
			agent.p = p(i);
			agent.q = q(i);
			agent.e = e(i);
			agent.f = f(i);
			agent.p_max = p_max.getEntry(i);
			agent.lp = lp(i);
			agent.lq = lq(i);
			agent.isGenerator = generatorMask.getEntry(i) == 1;
			agent.isSlack = i == slackIndex;
			agents.add(agent);
		}
		
		// Setup Neighbours:
		for (Agent agent_i : agents)
		{
			int i = agent_i.index;
			for (Agent agent_j : agents)
			{
				int j = agent_j.index;
				Complex Y_ij = new Complex(G(i,j), B(i,j));
				if(!Y_ij.equals(Complex.ZERO))
				{
					if(i == j)
					{
						agent_i.G = Y_ij.getReal();
						agent_i.B = Y_ij.getImaginary();
					}
					else
					{
						agent_i.addNeighbour(agent_j, Y_ij);
					}
				}
			}
		}
		
		for (Agent agent : agents)
			agent.init();
	}
	
	
	//// Iterations ////

	public void loop()
	{
		final int K = 1000;
		for(int k = 0; k < K; ++k)
		{
			for (Agent agent : agents)
			{
				agent.stepPrimal();
			}
			
			// The following ensure parallel processing of primal:
			for (Agent agent : agents)
			{
				agent.updatePrimal();
			}
			
			for (Agent agent : agents)
			{
				agent.stepDual();
			}
			
			update();

			// Log:
			debug(k);
		}
	}
	
	/**
	 * Update global variables from agent values.
	 */
	private void update()
	{
		for (Agent agent : agents)
		{
			int i = agent.index;
			p.setEntry(i, agent.p);
			q.setEntry(i, agent.q);
			e.setEntry(i, agent.e);
			f.setEntry(i, agent.f);
			lp.setEntry(i, agent.lp);
			lq.setEntry(i, agent.lq);
		}
	}
	
	
	//// Convenience methods ////

	public static String format(double d)
	{
		return String.format("%.4f", d);
	}
	
	public static double sum(IndexedFunction f, int dimension, int exclude)
	{
		double sum = 0;
		for(int i = 0; i < dimension; ++i)
		{
			if(i != exclude)
			{
				sum += f.value(i);
			}
		}
		
		return sum;
	}
	
	public static double sum(IndexedFunction f, int dimension)
	{
		double sum = 0;
		for(int i = 0; i < dimension; ++i)
		{
			sum += f.value(i);
		}
		
		return sum;
	}
	
	public static double sum(Map<Agent, Complex> neighbours, NeighbourFunction f)
	{
		double sum = 0;
		for (Agent agent : neighbours.keySet())
		{
			double v = f.value(agent);
			sum += v;
		}
		
		return sum;
	}

	protected double p(int i)
	{
		return p.getEntry(i);
	}

	protected double q(int i)
	{
		return q.getEntry(i);
	}

	protected double lp(int i)
	{
		return lp.getEntry(i);
	}

	protected double lq(int i)
	{
		return lq.getEntry(i);
	}

	protected double B(int i, int j)
	{
		return B.getEntry(i, j);
	}

	protected double G(int i, int j)
	{
		return G.getEntry(i, j);
	}

	protected double e(int i)
	{
		return e.getEntry(i);
	}

	protected double f(int i)
	{
		return f.getEntry(i);
	}
	
	
	//// Debug ////

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
		debugLogNames(numbers, dimension, "p(error)");
		debugLogNames(numbers, dimension, "q(error)");
		debugLogNames(numbers, dimension, "gradL_p");
		debugLogNames(numbers, dimension, "gradL_q");
		debugLogNames(numbers, dimension, "gradL_e");
		debugLogNames(numbers, dimension, "gradL_f");
		System.out.print("C(pq),L(pql),t");
		System.out.println();
	}
	
	/**
	 * 
	 * @param k Iteration counter.
	 */
	protected void debug(int k)
	{
		System.out.print(k);
		System.out.print(',');
		int dimension = p.getDimension();
		
		// Power:
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(p(i)));
			System.out.print(',');
		}
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(q(i)));
			System.out.print(',');
		}
		
		// Voltage:
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(e(i)));
			System.out.print(',');
		}
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(f(i)));
			System.out.print(',');
		}
		
		// Lambda:
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(lp(i)));
			System.out.print(',');
		}
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(lq.getEntry(i)));
			System.out.print(',');
		}
		
		// Calculated Power Error:
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(/*lp(i)**/gp(i)));
			System.out.print(',');
		}
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(/*lq(i)**/gq(i)));
			System.out.print(',');
		}
		
		// Power derivatives:
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(gradL_p(i)));
			System.out.print(',');
		}
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(gradL_q(i)));
			System.out.print(',');
		}
		
		// Voltage derivatives: // FIXME these are different compared to Sandbox013
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(gradL_e(i)));
			System.out.print(',');
		}
		for(int i = 0; i < dimension; ++i)
		{
			System.out.print(format(gradL_f(i)));
			System.out.print(',');
		}
		
		// Slack power:
//		System.out.print(format(pSlack()));
//		System.out.print(',');
//		System.out.print(format(qSlack()));
//		System.out.print(',');
		
		// Costs:
		System.out.print(format(cost()));
		System.out.print(',');
		System.out.print(format(L()));
		System.out.print(',');

		// Step size:
		System.out.print(format(0));
		System.out.print(',');
		
		System.out.println();
	}

	protected double gradL_p(int i)
	{
		return 	gradF_p(i) - lp(i)

				// Penalty:
				-c*gp(i)
				;
	}
	
	/**
	 * Assuming F_i(p_i) = 0.5*(p_i - p_max)^2
	 * => dF_i/dp_i = p_i - p_max
	 * @param i
	 * @return
	 */
	private double gradF_p(int i)
	{
		return p(i) - p_max.getEntry(i);
	}

	protected double gradL_q(int i)
	{
		return gradF_q(i) - lq(i)

				// Penalty:
				-c*gq(i)
				;
	}
	
	/**
	 * Assuming F_i(q_i) = 0.5*q_i^2
	 * => dF_i/dp_i = q_i
	 * @param i
	 * @return
	 */
	private double gradF_q(int i)
	{
		return q(i);
	}
	
	private double gp(int i)
	{
		int dimension = p.getDimension();
		return 
			sum(n -> e(i)*e(n)*G(i,n) + f(i)*f(n)*G(i,n) + f(i)*e(n)*B(i,n) - e(i)*f(n)*B(i,n), dimension)
			- p(i);
	}
	
	private double gq(int i)
	{
		int dimension = q.getDimension();
		return 
			sum(n -> f(i)*e(n)*G(i,n) - e(i)*f(n)*G(i,n) - e(i)*e(n)*B(i,n) - f(i)*f(n)*B(i,n), dimension)
			- q(i);
	}

	protected double gradL_e(int j)
	{
		int dimension = p.getDimension();
		double constraintTerms = sum(i -> lp(i)*(e(i)*G(i,j) + f(i)*B(i,j)), dimension, j) // i != j
					+sum(i -> lq(i)*(f(i)*G(i,j) - e(i)*B(i,j)), dimension, j) // i != j
					+lp(j)*(sum(n -> e(n)*G(j,n) - f(n)*B(j,n), dimension, j) // n != j
					        +2*G(j,j)*e(j))
					+lq(j)*(sum(n -> -f(n)*G(j,n) - e(n)*B(j,n), dimension, j) // n != j
					        -2*B(j,j)*e(j));
		
		double pt1 = c*sum(i -> i == slackIndex ? 0 : gp(i)*(e(i)*G(i,j) + f(i)*B(i,j))
				   		+gq(i)*(f(i)*G(i,j) - e(i)*B(i,j)),
				   			dimension, j);
		double pt2 = +c*gp(j)*(sum(n -> e(n)*G(j,n) - f(n)*B(j,n), dimension, j) 
				+2*e(j)*G(j,j));
		double pt3 = c*gq(j)*(sum(n -> -f(n)*G(j,n) - e(n)*B(j,n), dimension, j) 
				-2*e(j)*B(j,j));
		double penaltyTerms = pt1
							+
							(j == slackIndex ? 0 : 
							(
								pt2
								+pt3
							));
		return
			constraintTerms
			+penaltyTerms
			;
	}

	protected double gradL_f(int j)
	{
		int dimension = p.getDimension();
		double constraintTerms = sum(i -> lp(i)*(f(i)*G(i,j) - e(i)*B(i,j)), dimension, j) // i != j
					+sum(i -> lq(i)*(-e(i)*G(i,j) - f(i)*B(i,j)), dimension, j) // i != j
					+lp(j)*(sum(n -> f(n)*G(j,n) + e(n)*B(j,n), dimension, j) // n != j
					        +2*G(j,j)*f(j))
					+lq(j)*(sum(n -> e(n)*G(j,n) - f(n)*B(j,n), dimension, j) // n != j
					        -2*B(j,j)*f(j));
		
		double penaltyTerms = c*sum(i -> i == slackIndex ? 0 : gp(i)*( f(i)*G(i,j) - e(i)*B(i,j))
								                			   +gq(i)*(-e(i)*G(i,j) - f(i)*B(i,j)),
								    dimension, j)
							+ 
							(j == slackIndex ? 0 : 
							(
								+c*gp(j)*(sum(n -> f(n)*G(j,n) + e(n)*B(j,n), dimension, j) 
										  +2*f(j)*G(j,j))
								+c*gq(j)*(sum(n -> e(n)*G(j,n) - f(n)*B(j,n), dimension, j) 
										  -2*f(j)*B(j,j))
							));
		
		return 
			constraintTerms
			+penaltyTerms
			;
	}
	
	/**
	 * Lagrangian.
	 * @return
	 */
	protected double L()
	{
		int dimension = p.getDimension();
		
		return 
				cost()
				+sum(i -> lp(i)*gp(i), dimension)
				+sum(i -> lq(i)*gq(i), dimension)
				+0.5*c*sum(i -> gp(i)*gp(i), dimension, slackIndex)
				+0.5*c*sum(i -> gq(i)*gq(i), dimension, slackIndex);
	}

	private double cost()
	{
		int dimension = p.getDimension();
		return 
				sum(i -> generatorMask.getEntry(i) == 1 ? cost(i) : 0, dimension);
//				sum(i -> (1-e(i))*(1-e(i)) + f(i)*f(i), dimension);
	}

	protected double cost(int i)
	{
		return cost_i(p(i), q(i), p_max.getEntry(i));
	}

	protected double cost_i(double p_i, double q_i, double p_max_i)
	{
		double d = p_i - p_max_i;
		return 
				0.5*d*d
				+0.5*q_i*q_i;
	}
	
	public static void logBusses(AnalysisResults results)
	{
		System.out.println("index,bus,e,f,p,q");
		Map<String, Integer> numbers = results.getBusNumbers();
		for (String name : numbers.keySet())
		{
			int index = numbers.get(name);
			System.out.print(index);
			System.out.print(',');
			System.out.print(name);
			System.out.print(',');
			Complex v = results.getBusVoltage(name);
			System.out.print(format(v.getReal()));
			System.out.print(',');
			System.out.print(format(v.getImaginary()));
			System.out.print(',');
			Complex s = results.getBusPower(name);
			System.out.print(format(s.getReal()));
			System.out.print(',');
			System.out.print(format(s.getImaginary()));
			System.out.println();
		}
	}
}