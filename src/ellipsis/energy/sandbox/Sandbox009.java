package ellipsis.energy.sandbox;

import static java.lang.Math.cos;
import static java.lang.Math.sin;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.complex.Complex;

import ellipsis.energy.calculation.AnalysisResults;

/**
 * Distributed version of {@link Sandbox008}.
 * @author bmillar
 *
 */
public class Sandbox009 extends Sandbox008
{
	public static void main(String[] args)
	{
		new Sandbox009().run();
	}
	
	
	//// Setup ////
	
	@Override
	protected AnalysisResults init()
	{
		AnalysisResults results = super.init();
		
		// Make agents:
		agents = new HashSet<>();
		int dimension = p.getDimension();
		for(int i = 0; i < dimension; ++i)
		{
			Agent agent = new Agent(i);
			agent.setP(p.getEntry(i));
			agent.setQ(q.getEntry(i));
			agent.setVabs(vabs.getEntry(i));
			agent.setVarg(varg.getEntry(i));
			agent.setSmin(s_min.getEntry(i));
			agent.setSmax(s_max.getEntry(i));
			agent.setLambdaP(lambdaP.getEntry(i));
			agent.setLambdaQ(lambdaQ.getEntry(i));
			agent.setIsGenerator(G.contains(i));
			agent.setIsSlack(i == slackIndex);
			agents.add(agent);
		}
		
		// Setup Neighbours:
		for (Agent agent_i : agents)
		{
			int i = agent_i.getIndex();
			for (Agent agent_j : agents)
			{
				int j = agent_j.getIndex();
				Complex Y_ij = Y.getEntry(i, j);
				if(!Y_ij.equals(Complex.ZERO))
				{
					agent_i.addNeighbour(agent_j, Y_ij);
				}}
		}
		
		return results;
	}
	
	
	//// Iterations ////
	
	@Override
	public void loop()
	{
		final int K = 1000;
		for(int k = 0; k < K; ++k)
		{
			for (Agent agent : agents)
			{
				agent.step();
			}
			
			update();
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
			int i = agent.getIndex();
			p.setEntry(i, agent.getP());
			q.setEntry(i, agent.getQ());
			vabs.setEntry(i, agent.getVabs());
			varg.setEntry(i, agent.getVarg());
			lambdaP.setEntry(i, agent.getLambdaP());
			lambdaQ.setEntry(i, agent.getLambdaQ());
		}
	}
	
	
	//// Agents ////


	protected Set<Agent> agents;

	// Step sizes:
    protected static double pStepSize = 0.01;
    protected static double qStepSize = 0.01;
    protected static double vabsStepSize = 0.005; // 0.012 seems best
    protected static double vargStepSize = 0.01;
    protected static double lambdaPStepSize = 0.006;
    protected static double lambdaQStepSize = 0.006;
    
	public static class Agent
	{
	    
	    // Parameters:
		private int index;
		private Map<Agent, Complex> neighbours = new HashMap<>(); // Neighbours contains this agent too.
		protected boolean isGenerator;
		protected boolean isSlack;
		protected double s_min;
	    protected double s_max;
		protected double v_min = 0.95;
	    protected double v_max = 1.05;
		
		// Variables:
	    protected double p;
	    protected double q;
	    protected double vabs;
	    protected double varg;
	    protected double lambdaP;
	    protected double lambdaQ;

		public Agent(int i)
		{
			this.setIndex(i);
		}

		public void step()
		{
			optimiseStep();
			constraintStep();
		}

		// Gradient accent:
		public static final int P = 0;
		public static final int Q = 1;
		public static final int ABS = 0;
		public static final int ARG = 1;
		private void optimiseStep()
		{
			// P & Q:
			double gradP = gradP_i();
			double gradQ = gradQ_i();
			double nextP = p + gradP*pStepSize;
			double nextQ = q + gradQ*qStepSize;
			
			double[] nextPQ = clampPQ(nextP, nextQ);
			
			p = nextPQ[P];
			q = nextPQ[Q];
			
			// Voltage abs and arg:
			double gradVabs = gradVabs_i();
			double gradVarg = gradVarg_i();
			double nextVabs = vabs + gradVabs*vabsStepSize;
			double nextVarg = varg + gradVarg*vargStepSize;
			
			double nextV[] = clampV(nextVabs, nextVarg);
			
			vabs = nextV[ABS];
			varg = nextV[ARG];
		}

		/**
		 * Only clamps absolute, not argument.
		 * @param nextVabs
		 * @param nextVarg
		 * @return
		 */
		protected double[] clampV(double nextVabs, double nextVarg)
		{
//			if(nextVabs < v_min)
//				nextVabs = v_min;
//			else if(nextVabs > v_max)
//				nextVabs = v_max;
			
			return new double[]{nextVabs, nextVarg};
		}

		protected double[] clampPQ(double nextP, double nextQ)
		{
			if(!getIsGenerator())
				return new double[]{nextP, nextQ};
			
			double s_min_i = s_min;
			double s_max_i = s_max;
			double p_i = nextP;
			if(p_i < 0)
				p_i = 0;
			double q_i = nextQ;
			Complex s_i = new Complex(p_i, q_i);
			
			double s_i_abs = s_i.abs();
			if(s_i_abs < s_min_i)
				s_i = s_i.multiply(s_min_i/s_i_abs);
			else if(s_i_abs > s_max_i)
				s_i = s_i.multiply(s_max_i/s_i_abs);
			
			return new double[]{s_i.getReal(), s_i.getImaginary()};
		}

		private double gradP_i()
		{
			if(!getIsGenerator())
				return 0;
			return gradF_p_i() - lambdaP;
		}
		
		private double gradF_p_i()
		{
			double p_i = p;
			double p_max_i = s_max;
			return -2*(p_i-p_max_i);
		}

		private double gradQ_i()
		{
			if(!getIsGenerator())
				return 0;
			return gradF_q_i() - lambdaQ;
		}
		
		private double gradF_q_i()
		{
			return -2*q;
		}

		private double gradVabs_i()
		{
			if(getIsSlack())
				return 0;
			return  2*(1-vabs)
					+ lambdaP*gradVabs_i_pi() + gradVabs_i_pj()
					- lambdaQ*gradVabs_i_qi() - gradVabs_i_qj();
		}

		private double gradVabs_i_pi()
		{
			double varg_i = varg;
			double vabs_i = vabs;
			Complex Y_ii = selfAdmittance();
			
			double sum = 0;
			for (Agent agent_j : neighbours.keySet())
			{
				int j = agent_j.getIndex();
				if(index == j)
					continue;
				double vabs_j = agent_j.getVabs();
				double varg_j = agent_j.getVarg();
				Complex Y_ij = neighbours.get(agent_j);
				sum += vabs_j*Y_ij.abs()*cos(Y_ij.getArgument() - varg_i + varg_j);
			}
			return sum + 2*vabs_i*Y_ii.abs()*cos(Y_ii.getArgument());
		}

		private double gradVabs_i_pj()
		{
			double varg_i = varg;
			double sum = 0;
			for (Agent agent_j : neighbours.keySet())
			{
				int j = agent_j.getIndex();
				if(index == j)// FIXME || getIsSlack())
					continue;
				double vabs_j = agent_j.getVabs();
				double varg_j = agent_j.getVarg();
				Complex Y_ji = agent_j.neighbours.get(this);
				double lambdaP_j = agent_j.getLambdaP();
				sum += lambdaP_j*vabs_j*Y_ji.abs()*cos(Y_ji.getArgument() - varg_j + varg_i);
			}
			
			return sum;
		}

		private double gradVabs_i_qi()
		{
			double varg_i = varg;
			double vabs_i = vabs;
			Complex Y_ii = selfAdmittance();
			
			double sum = 0;
			for (Agent agent_j : neighbours.keySet())
			{
				int j = agent_j.getIndex();
				if(index == j)
					continue;
				double vabs_j = agent_j.getVabs();
				double varg_j = agent_j.getVarg();
				Complex Y_ij = neighbours.get(agent_j);
				sum += vabs_j*Y_ij.abs()*sin(Y_ij.getArgument() - varg_i + varg_j);
			}
			return sum + 2*vabs_i*Y_ii.abs()*sin(Y_ii.getArgument());
		}

		private double gradVabs_i_qj()
		{
			double varg_i = varg;
			double sum = 0;
			for (Agent agent_j : neighbours.keySet())
			{
				int j = agent_j.getIndex();
				if(index == j)// FIXME || getIsSlack())
					continue;
				double vabs_j = agent_j.getVabs();
				double varg_j = agent_j.getVarg();
				Complex Y_ji = agent_j.neighbours.get(this);
				double lambdaQ_j = agent_j.getLambdaQ();
				sum += lambdaQ_j*vabs_j*Y_ji.abs()*sin(Y_ji.getArgument() - varg_j + varg_i);
			}
			
			return sum;
		}

		private double gradVarg_i()
		{
			if(getIsSlack())
				return 0;
			return  -2*(varg)
					+ lambdaP*gradVarg_i_pi() - gradVarg_i_pj() 
					+ lambdaQ*gradVarg_i_qi() - gradVarg_i_qj();
		}
		
		private double gradVarg_i_pi()
		{
			double vabs_i = vabs;
			double varg_i = varg;
					
			double sum = 0;
			for (Agent agent_j : neighbours.keySet())
			{
				int j = agent_j.getIndex();
				if(index == j)
					continue;
				double vabs_j = agent_j.getVabs();
				double varg_j = agent_j.getVarg();
				Complex Y_ij = neighbours.get(agent_j);
				sum += vabs_i*vabs_j*Y_ij.abs()*sin(Y_ij.getArgument() - varg_i + varg_j);
			}
			
			return sum;
		}

		private double gradVarg_i_pj()
		{
			double vabs_i = vabs;
			double varg_i = varg;
					
			double sum = 0;
			for (Agent agent_j : neighbours.keySet())
			{
				int j = agent_j.getIndex();
				if(index == j)
					continue;
				double vabs_j = agent_j.getVabs();
				double varg_j = agent_j.getVarg();
				Complex Y_ji = agent_j.neighbours.get(this);
				double lambdaP_j = agent_j.getLambdaP();
				sum += lambdaP_j*vabs_i*vabs_j*Y_ji.abs()*sin(Y_ji.getArgument() - varg_j + varg_i);
			}
			
			return sum;
		}

		private double gradVarg_i_qi()
		{
			double vabs_i = vabs;
			double varg_i = varg;
					
			double sum = 0;
			for (Agent agent_j : neighbours.keySet())
			{
				int j = agent_j.getIndex();
				if(index == j)
					continue;
				double vabs_j = agent_j.getVabs();
				double varg_j = agent_j.getVarg();
				Complex Y_ij = neighbours.get(agent_j);
				sum += vabs_i*vabs_j*Y_ij.abs()*cos(Y_ij.getArgument() - varg_i + varg_j);
			}
			
			return sum;
		}

		private double gradVarg_i_qj()
		{
			double vabs_i = vabs;
			double varg_i = varg;
					
			double sum = 0;
			for (Agent agent_j : neighbours.keySet())
			{
				int j = agent_j.getIndex();
				if(index == j)
					continue;
				double vabs_j = agent_j.getVabs();
				double varg_j = agent_j.getVarg();
				Complex Y_ji = agent_j.neighbours.get(this);
				double lambdaQ_j = agent_j.getLambdaQ();
				sum += lambdaQ_j*vabs_i*vabs_j*Y_ji.abs()*cos(Y_ji.getArgument() - varg_j + varg_i);
			}
			
			return sum;
		}

		private void constraintStep()
		{
			// Gradient decent:
			double gradLambdaP = gradLambdaP_i();
			double gradLambdaQ = gradLambdaQ_i();
			double nextLambdaP = lambdaP - gradLambdaP*lambdaPStepSize;
			double nextLambdaQ = lambdaQ - gradLambdaQ*lambdaQStepSize;

			lambdaP = nextLambdaP;
			lambdaQ = nextLambdaQ;
		}

		private double gradLambdaP_i()
		{
if(getIsSlack())
	return 0; // FIXME Include slack bus here for power balancing.
			double vabs_i = vabs;
			double varg_i = varg;
			
			double sum = 0;
			for (Agent agent_j : neighbours.keySet())
			{
				double vabs_j = agent_j.getVabs();
				double varg_j = agent_j.getVarg();
				Complex Y_ij = neighbours.get(agent_j);
				sum += vabs_i*vabs_j*Y_ij.abs()*cos(Y_ij.getArgument() - varg_i + varg_j);
			}
			
			double p_i = getIsSlack() ? 0 : p;
			return sum - p_i;
		}

		private double gradLambdaQ_i()
		{
if(getIsSlack())
	return 0; // FIXME
			double vabs_i = vabs;
			double varg_i = varg;
			
			double sum = 0;
			for (Agent agent_j : neighbours.keySet())
			{
				double vabs_j = agent_j.getVabs();
				double varg_j = agent_j.getVarg();
				Complex Y_ij = neighbours.get(agent_j);
				sum += vabs_i*vabs_j*Y_ij.abs()*sin(Y_ij.getArgument() - varg_i + varg_j);
			}
			
			double q_i = getIsSlack() ? 0 : q;
			return -(sum + q_i);
		}
		
		
		//// Agent Accessors ////

		protected Complex selfAdmittance()
		{
			return neighbours.get(this);
		}
		
		public void addNeighbour(Agent neighbour, Complex admittance)
		{
			neighbours.put(neighbour, admittance);
		}

		public int getIndex()
		{
			return index;
		}

		public void setIndex(int index)
		{
			this.index = index;
		}
		
	    public void setNeighbours(Map<Agent, Complex> neighbours)
		{
			this.neighbours = neighbours;
		}

		public void setIsGenerator(boolean isGenerator)
		{
			this.isGenerator = isGenerator;
		}
		
		public void setIsSlack(boolean isSlack)
		{
			this.isSlack = isSlack;
		}

		public void setSmin(double s_min)
		{
			this.s_min = s_min;
		}

		public void setSmax(double s_max)
		{
			this.s_max = s_max;
		}

		public void setP(double p)
		{
			this.p = p;
		}

		public void setQ(double q)
		{
			this.q = q;
		}

		public void setVabs(double vabs)
		{
			this.vabs = vabs;
		}

		public void setVarg(double varg)
		{
			this.varg = varg;
		}

		public void setLambdaP(double lambdaP)
		{
			this.lambdaP = lambdaP;
		}

		public void setLambdaQ(double lambdaQ)
		{
			this.lambdaQ = lambdaQ;
		}

		public Map<Agent, Complex> getNeighbours()
		{
			return neighbours;
		}

		public boolean getIsGenerator()
		{
			return isGenerator;
		}
		
		public boolean getIsSlack()
		{
			return isSlack;
		}

		public double getS_min()
		{
			return s_min;
		}

		public double getS_max()
		{
			return s_max;
		}

		public double getP()
		{
			return p;
		}

		public double getQ()
		{
			return q;
		}

		public double getVabs()
		{
			return vabs;
		}

		public double getVarg()
		{
			return varg;
		}

		public double getLambdaP()
		{
			return lambdaP;
		}

		public double getLambdaQ()
		{
			return lambdaQ;
		}
	}
}