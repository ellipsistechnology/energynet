package ellipsis.energy.sandbox.Sandbox021;

import static ellipsis.util.VectorHelper.vector;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import ellipsis.energy.sandbox.Sandbox020SmartBuilding;
import ellipsis.energy.sandbox.Sandbox021_AugLagSolver;
import ellipsis.energy.sandbox.Sandbox021_AugLagSolver.Problem;
import ellipsis.energy.sandbox.Sandbox021_AugLagSolver.ProblemConfiguration;
import ellipsis.energy.sandbox.Sandbox021_AugLagSolver.Solution;
import ellipsis.energy.sandbox.sandbox020.Agent;
import ellipsis.energy.sandbox.sandbox020.ConstantCurrentAgent;
import ellipsis.energy.sandbox.sandbox020.ConstantPowerAgent;
import ellipsis.energy.sandbox.sandbox020.VoltageControlledAgent;

public class HEMSProblem implements Problem
{	
	public ArrayList<Agent> agents = new ArrayList<>();
	
	@Override
	public int getEqualityConstraintDimension()
	{
		return agents.size()*2;
	}

	@Override
	public RealVector project(RealVector x)
	{
		RealVector projX = new ArrayRealVector(x.getDimension());
		int i = 0; // agent index
		for (Agent agent : agents)
		{
			double vplus = x.getEntry(i*2);
			double vminus = x.getEntry(i*2+1);
			
			if(agent.isGrounded())
	        {
	            if(vplus > agent.getVmax())
	            {
	                vplus = agent.getVmax();
	            }
	            else if(vplus < agent.getVmin())
	            {
	                vplus = agent.getVmin();
	            }
	        }
	        else
	        {
	            double v = vplus-vminus;
	            double vlimit = 0;
	            if(v > agent.getVmax())
	                vlimit = agent.getVmax();
	            else if(v < agent.getVmin())
	                vlimit = agent.getVmin();
	            if(vlimit != 0)
	            {
	                double delta = (vminus - vplus + vlimit)/2;
	                vplus += delta;
	                vminus -= delta;
	            }
	        }
			
			projX.setEntry(i*2, vplus);
			projX.setEntry(i*2+1, vminus);
			
			++i;
		}
		
		return projX;
//		return x;
	}

	@Override
	public double cost(RealVector x)
	{
		double x2 = 0;
		int i = 0;
		for (Agent agent : agents)
		{
			if(agent instanceof VoltageControlledAgent)
			{
				double vplus = x.getEntry(i*2);
				x2 += vplus*vplus;

				if(!agent.isGrounded())
				{
					double vminus = x.getEntry(i*2+1);
					x2 += vminus*vminus;
				}
			}
			++i;
		}
		return x2;
	}

	@Override
	public RealVector g(RealVector x)
	{
		RealVector g = new ArrayRealVector(x.getDimension());
		int i = 0; // agent index
		for (Agent agent : agents)
		{
			double vplus = x.getEntry(i*2);
			double vminus = x.getEntry(i*2+1);
			if(agent instanceof VoltageControlledAgent)
			{
				g.setEntry(i*2, gVC(x, agent, vplus, vminus));
			}
			else if(agent instanceof ConstantCurrentAgent)
			{
				g.setEntry(i*2, gCC1(x, (ConstantCurrentAgent)agent, vplus, vminus));
				g.setEntry(i*2+1, gCC2(x, (ConstantCurrentAgent)agent, vplus, vminus));
			}
			else if(agent instanceof ConstantPowerAgent)
			{
				g.setEntry(i*2, gCP1(x, (ConstantPowerAgent)agent, vplus, vminus));
				g.setEntry(i*2+1, gCP2(x, (ConstantPowerAgent)agent, vplus, vminus));
			}
			else
				throw new RuntimeException("Invalid Agent");
			++i;
		}
		
		return g;
	}

	private double gVC(RealVector x, Agent agent, double vplus, double vminus)
	{
		double g = 0;
		int j = 0;
		for (Agent neighbour : agents)
		{
			double vplus_j = x.getEntry(j*2);
			double vminus_j = x.getEntry(j*2+1);
			double y_ij = agent.conductance(neighbour);
			g += (vplus - vplus_j)*y_ij;
			g += (vminus - vminus_j)*y_ij;
			
			++j;
		}
		
		return g;
	}

	private double gCC1(RealVector x, ConstantCurrentAgent agent, double vplus, double vminus)
	{
		double g = 0;
		int j = 0;
		for (Agent neighbour : agents)
		{
			double vplus_j = x.getEntry(j*2);
			double y_ij = agent.conductance(neighbour);
			g += (vplus - vplus_j)*y_ij;
			
			++j;
		}
		
		g -= agent.getCurrent();
		
		return g;
	}

	private double gCC2(RealVector x, ConstantCurrentAgent agent, double vplus, double vminus)
	{
		double g = 0;
		int j = 0;
		for (Agent neighbour : agents)
		{
			double vminus_j = x.getEntry(j*2+1);
			double y_ij = agent.conductance(neighbour);
			g += (vminus - vminus_j)*y_ij;
			
			++j;
		}
		
		g = -g - agent.getCurrent();
		
		return g;
	}

	private double gCP1(RealVector x, ConstantPowerAgent agent, double vplus, double vminus)
	{
		double g = 0;
		int j = 0;
		for (Agent neighbour : agents)
		{
			double vplus_j = x.getEntry(j*2);
			double y_ij = agent.conductance(neighbour);
			g += (vplus - vplus_j)*y_ij;
			
			++j;
		}
		
		g *= vplus-vminus;
		
		g -= agent.getPower();
		
		return g;
	}

	private double gCP2(RealVector x, ConstantPowerAgent agent, double vplus, double vminus)
	{
		double g = 0;
		int j = 0;
		for (Agent neighbour : agents)
		{
			double vminus_j = x.getEntry(j*2+1);
			double y_ij = agent.conductance(neighbour);
			g += (vminus - vminus_j)*y_ij;
			
			++j;
		}
		
		g *= vplus-vminus;
		
		g = -g - agent.getPower();
		
		return g;
	}

	@Override
	public RealVector costGradient(RealVector x)
	{
		RealVector x2 = vector(x.getDimension(), 0.0);
		int i = 0;
		for (Agent agent : agents)
		{
			if(agent instanceof VoltageControlledAgent)
			{
				double vplus = x.getEntry(i*2);
				x2.setEntry(2*i, 2.0*vplus);
				
				if(!agent.isGrounded())
				{
					double vminus = x.getEntry(i*2+1);
					x2.setEntry(2*i+1, 2.0*vminus);
				}
			}
			++i;
		}
		return x2;
	}

	@Override
	public RealMatrix equalityConstraintGradient(RealVector x)
	{
		RealMatrix gradG = new Array2DRowRealMatrix(x.getDimension(), 2*agents.size()); // 1 row per v, and 1 column per g
		int i = 0; // agent index - for constraint functions
		for (Agent agent : agents)
		{
			double vplus = x.getEntry(i*2);
			double vminus = x.getEntry(i*2+1);
			
			int j = 0; // agent index - for wrt
			for (Agent wrt : agents)
			{
				double vplus_wrt = x.getEntry(j*2);
				double vminus_wrt = x.getEntry(j*2+1);
				boolean grounded = wrt.isGrounded();
				
				if(agent instanceof VoltageControlledAgent)
				{
					gradG.setEntry(j*2,   i*2,   gradGVCwrtVplus (x, agent, wrt));
					if(!grounded)
						gradG.setEntry(j*2+1, i*2,   gradGVCwrtVminus(x, agent, wrt));
				}
				else if(agent instanceof ConstantCurrentAgent)
				{
					gradG.setEntry(j*2,   i*2,   gradGCC1wrtVplus (x, (ConstantCurrentAgent)agent, wrt));
					if(!grounded)
						gradG.setEntry(j*2+1, i*2,   gradGCC1wrtVminus(x, (ConstantCurrentAgent)agent, wrt));
					gradG.setEntry(j*2,   i*2+1, gradGCC2wrtVplus (x, (ConstantCurrentAgent)agent, wrt));
					if(!grounded)
						gradG.setEntry(j*2+1, i*2+1, gradGCC2wrtVminus(x, (ConstantCurrentAgent)agent, wrt));
				}
				else if(agent instanceof ConstantPowerAgent)
				{
					gradG.setEntry(j*2,   i*2,   gradGCP1wrtVplus (x, (ConstantPowerAgent)agent, wrt, vplus, vminus, vplus_wrt, vminus_wrt));
					if(!grounded)
						gradG.setEntry(j*2+1, i*2,   gradGCP1wrtVminus(x, (ConstantPowerAgent)agent, wrt, vplus, vminus, vplus_wrt, vminus_wrt));
					gradG.setEntry(j*2,   i*2+1, gradGCP2wrtVplus (x, (ConstantPowerAgent)agent, wrt, vplus, vminus, vplus_wrt, vminus_wrt));
					if(!grounded)
						gradG.setEntry(j*2+1, i*2+1, gradGCP2wrtVminus(x, (ConstantPowerAgent)agent, wrt, vplus, vminus, vplus_wrt, vminus_wrt));
				}
				else
					throw new RuntimeException("Invalid Agent");
				++j;
			}
			++i;
		}
		
		return gradG;
	}

	private double gradGVCwrtVplus(RealVector x, Agent agent, Agent wrt)
	{
		if(agent == wrt)
			return agent.conductanceSum();
		else
			return -agent.conductance(wrt);
	}

	private double gradGVCwrtVminus(RealVector x, Agent agent, Agent wrt)
	{
		if(wrt.isGrounded())
			return 0.0;
		
		if(agent == wrt)
			return agent.conductanceSum();
		else
			return -agent.conductance(wrt);
	}

	private double gradGCC1wrtVplus(RealVector x, ConstantCurrentAgent agent, Agent wrt)
	{
		if(agent == wrt)
			return agent.conductanceSum();
		else
			return -agent.conductance(wrt);
	}

	private double gradGCC1wrtVminus(RealVector x, ConstantCurrentAgent agent, Agent wrt)
	{
		return 0.0;
	}

	private double gradGCC2wrtVplus(RealVector x, ConstantCurrentAgent agent, Agent wrt)
	{
		return 0.0;
	}

	private double gradGCC2wrtVminus(RealVector x, ConstantCurrentAgent agent, Agent wrt)
	{
		if(agent == wrt)
			return -agent.conductanceSum();
		else
			return agent.conductance(wrt);
	}

	private double gradGCP1wrtVplus(RealVector x, ConstantPowerAgent agent, Agent wrt, double vplus, double vminus, double vplus_wrt, double vminus_wrt)
	{
		if(agent == wrt)
			return (2*vplus - vminus)*agent.conductanceSum() - sumVplus_jY_ij(x, agent);
		else
			return -(vplus - vminus)*agent.conductance(wrt);
	}

	public double sumVplus_jY_ij(RealVector x, ConstantPowerAgent agent)
	{
		double sum = 0;
		int j = 0;
		for (Agent agent_j : agents)
		{
			double vplus_j = x.getEntry(j*2);
			double y_ij = agent.conductance(agent_j);
			sum += vplus_j*y_ij;
		}
		return sum;
	}

	public double sumVminus_jY_ij(RealVector x, ConstantPowerAgent agent)
	{
		double sum = 0;
		int j = 0;
		for (Agent agent_j : agents)
		{
			double vminus_j = x.getEntry(j*2+1);
			double y_ij = agent.conductance(agent_j);
			sum += vminus_j*y_ij;
		}
		return sum;
	}

	private double gradGCP1wrtVminus(RealVector x, ConstantPowerAgent agent, Agent wrt, double vplus, double vminus, double vplus_wrt, double vminus_wrt)
	{
		if(agent == wrt)
			return -vplus*agent.conductanceSum() + sumVplus_jY_ij(x, agent);
		else
			return 0.0;
	}

	private double gradGCP2wrtVplus(RealVector x, ConstantPowerAgent agent, Agent wrt, double vplus, double vminus, double vplus_wrt, double vminus_wrt)
	{
		if(agent == wrt)
			return -vminus*agent.conductanceSum() + sumVminus_jY_ij(x, agent);
		else
			return 0.0;
	}

	private double gradGCP2wrtVminus(RealVector x, ConstantPowerAgent agent, Agent wrt, double vplus, double vminus, double vplus_wrt, double vminus_wrt)
	{
		if(agent == wrt)
			return (2*vminus - vplus)*agent.conductanceSum() - sumVminus_jY_ij(x, agent);
		else
			return (vplus - vminus)*agent.conductance(wrt);
	}
	
	Random rand = new Random(0);
	@Override
	public Set<Integer> nextStochasticIndeces(int k)
	{
		// FIXME testing
//if(k >= 901)
//{
//	Set<Integer> is = new HashSet<>();
//	
//	if(k >= 950)
//	{
//		is.add(0);
//		is.add(1);
//		is.add(2);
//		is.add(3);
//	}
////	else
////	if(k >= 950)
//	{
//		is.add(4);
//		is.add(5);
//		is.add(6);
//		is.add(7);
//		is.add(8);
//		is.add(9);
//		is.add(10);
//		is.add(11);
//	}
//	return is;
//}
//else
//{
//	return null;
//}
		Set<Integer> is = new HashSet<>();
		int i = rand.nextInt(agents.size());
		is.add(i*2);
		is.add(i*2+1);
		return is;
	}
	
	
	//// Testing ////
	
	public static void main(String[] args)
	{
		HEMSProblem problem = new HEMSProblem();
		
		Set<Agent> agentSet = Sandbox020SmartBuilding.testCase02();//testCase03_2node();
		problem.agents = new ArrayList<Agent>(agentSet);
		
		augLagSolve(problem, agentSet);
//		solve(problem, agentSet);
	}

//	private static void solve(HEMSProblem problem, Set<Agent> agentSet)
//	{
//		BOBYQAOptimizer opt = new BOBYQAOptimizer(14);
//		RealVector initialX = Sandbox020SmartBuilding.loadState(agentSet);
//		opt.optimize(
//				new MaxEval(100), 
//				new ObjectiveFunction(new OptProblem(problem, agentSet)), 
//				GoalType.MINIMIZE, 
//				new SimpleBounds(new double[]{0,0,0},new double[]{1,1,1}), 
//				new InitialGuess(initialX.toArray()));
//	}
//	
//	public static class OptProblem implements MultivariateFunction
//	{
//		private HEMSProblem problem;
//
//		public OptProblem(HEMSProblem problem, Set<Agent> agentSet)
//		{
//			this.problem = problem;
//		}
//
//		@Override
//		public double value(double[] point)
//		{
//			return problem.cost(vector(point));
//		}
//		
//	}

	public static void augLagSolve(HEMSProblem problem, Set<Agent> agentSet)
	{
		ProblemConfiguration config = new ProblemConfiguration();
		config.initialAlpha = 1e-3;
		config.alphaMultiplier = 1.001; //1.00625;/*stochastic*/ //1.025;/*Centralised*/
		config.initialEpsilon = 1e-2;//1;//10.0
		config.epsilonMultiplier = 1.0;//0.999;//0.99;
		config.iterations = 12000;
		config.initialState = Sandbox020SmartBuilding.loadState(agentSet);
		config.stepMultiplier = 0.5;
		config.stepScale = 0.5;
		config.lambdaBound = 100.0; // FIXME lambda isn't helping
		config.stochastic = true;//true;
		config.equalityConstraintTarget = 1e-12;
		
		Solution s = Sandbox021_AugLagSolver.solve(problem, config);
		
		Sandbox021_AugLagSolver.printCSV(problem, s, 1000);
	}
}