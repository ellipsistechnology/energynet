package ellipsis.energy.mas;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import ellipsis.energy.calculation.AdmittanceMatrix;
import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.calculation.LoadFlowAnalyser;
import ellipsis.energy.grid.Bus;
import ellipsis.energy.grid.Grid;
import ellipsis.energy.test.IEEE13BusGrid;

public class AgentCoordinator
{
	private AdmittanceMatrix Y;
	private Grid grid;
	private Set<Agent> agents;
	
	
	//// Setup ////
	
	public AgentCoordinator(Grid grid, double basePower, double baseVoltage)
	{
		this.grid = grid;
		
		LoadFlowAnalyser analyser = new LoadFlowAnalyser(grid);
		analyser.setBasePower(basePower);
		analyser.setBaseVoltage(baseVoltage);
		AnalysisResults results = analyser.analyse();
		
		Y = results.getAdmittanceMatrix();
		createAgents(results);
	}

	private void createAgents(AnalysisResults results)
	{
		RealMatrix sensitivities = results.sensitivities();
		
		// Create an agent for each bus:
		Map<String, Integer> numbers = results.getBusNumbers();
		agents = new HashSet<>();
		Map<String, Agent> agentMap = new HashMap<>();
		for (String name : numbers.keySet())
		{
			Bus bus = grid.getBus(name);
			int i = numbers.get(name);
			AgentProcess process = new PowerFlowAgentProcess();
			Agent agent = new Agent(bus, Y.get(i, i), process);
			process.setAgent(agent);
			
			agents.add(agent);
			agentMap.put(name, agent);
		}
		
		// Setup agents:
		for (Agent agent_i : agents)
		{
			Bus bus_i = agent_i.getBus();
			String name_i = bus_i.getName();

			// Setup neighbours:
			{
				int index_i = numbers.get(name_i); // includes slack index
				for (String name_j : numbers.keySet())
				{
					int index_j = numbers.get(name_j);
					Complex y_ij = Y.get(index_i, index_j);
					if(!y_ij.equals(Complex.ZERO) && index_i != index_j)
					{
						Agent agent_j = agentMap.get(name_j);
						agent_i.addNeighbour(y_ij, agent_j);
					}
				}
			}
			
			// Setup sensitivities:
			{
				int index_i = results.indexWithoutSlack(name_i);
				Map<Agent, Complex> neighbours_i = agent_i.getNeighbours();
				int size_i = neighbours_i.size()+1;
				if(agent_i.hasSlackNeighbour())
					--size_i;
				RealMatrix sensitivities_i = new Array2DRowRealMatrix(new double[2][size_i*2]);
				int size = sensitivities.getColumnDimension()/2;
				
				// Neighbour sensitivities:
				{
					int j = 0;
					for (Agent agent_j : neighbours_i.keySet())
					{
						if(agent_j.isSlack())
							continue;
						Bus bus_j = agent_j.getBus();
						String name_j = bus_j.getName();
						int index_j = results.indexWithoutSlack(name_j);
						double sensitivities_ij_argP = sensitivities.getEntry(index_i, index_j);
						double sensitivities_ij_argQ = sensitivities.getEntry(index_i, index_j+size);
						double sensitivities_ij_vP = sensitivities.getEntry(index_i+size, index_j);
						double sensitivities_ij_vQ = sensitivities.getEntry(index_i+size, index_j+size);
						sensitivities_i.setEntry(0, j, sensitivities_ij_argP);
						sensitivities_i.setEntry(0, j+size_i, sensitivities_ij_argQ);
						sensitivities_i.setEntry(1, j, sensitivities_ij_vP);
						sensitivities_i.setEntry(1, j+size_i, sensitivities_ij_vQ);
						++j;
					}
				}
				
				// Self sensitivities (append to the end of the matrix):
				double sensitivities_ii_argP = sensitivities.getEntry(index_i, index_i);
				double sensitivities_ii_argQ = sensitivities.getEntry(index_i, index_i+size);
				double sensitivities_ii_vP = sensitivities.getEntry(index_i+size, index_i);
				double sensitivities_ii_vQ = sensitivities.getEntry(index_i+size, index_i+size);
				sensitivities_i.setEntry(0, size_i-1, sensitivities_ii_argP);
				sensitivities_i.setEntry(0, 2*size_i-1, sensitivities_ii_argQ);
				sensitivities_i.setEntry(1, size_i-1, sensitivities_ii_vP);
				sensitivities_i.setEntry(1, 2*size_i-1, sensitivities_ii_vQ);
				
				agent_i.setSensitivities(sensitivities_i);
			}
		}
	}
	
	
	//// Coordination ////

	private static final int K = 100000;
	public boolean distributedPowerFlow()
	{
		boolean converged = false;
		for(int k = 0; k < K && !converged; ++k)
		{
			converged = k != 0; // first iteration must not converge - deltas will start at zero
System.out.print(k);
//System.out.print(',');
			for (Agent agent : agents)
			{
				if(!agent.getAgentProcess().stepInit())
					converged = false;

//System.out.print(',');
//System.out.print(agent.getV().abs());
//System.out.print(',');
//System.out.print(agent.getV().getArgument());
//System.out.print(',');
//System.out.print(agent.getDeltaS().abs());
//System.out.print(',');
//System.out.print(agent.getDeltaS().getArgument());
			}
System.out.println();

			for (Agent agent : agents)
			{
				agent.getAgentProcess().step();
			}
			
//			debug(k);
		}
		
		return converged;
	}
	
	
	//// Accessors ////

	public void debug(int k)
	{
		// Headings:
		System.out.print("k, Delta P");
		for (@SuppressWarnings("unused") Agent agent_i : agents)
			System.out.print(',');
		System.out.print("Delta Q");
		for (@SuppressWarnings("unused") Agent agent_i : agents)
			System.out.print(',');
		System.out.print("P");
		for (@SuppressWarnings("unused") Agent agent_i : agents)
			System.out.print(',');
		System.out.print("Q");
		for (@SuppressWarnings("unused") Agent agent_i : agents)
			System.out.print(',');
		System.out.print("v arg");
		for (@SuppressWarnings("unused") Agent agent_i : agents)
			System.out.print(',');
		System.out.print("|v|");
		
		System.out.println();
		
		System.out.print(',');
		for (Agent agent_i : agents)
		{
			System.out.print(agent_i);
			System.out.print(',');
		}
		for (Agent agent_i : agents)
		{
			System.out.print(agent_i);
			System.out.print(',');
		}
		for (Agent agent_i : agents)
		{
			System.out.print(agent_i);
			System.out.print(',');
		}
		for (Agent agent_i : agents)
		{
			System.out.print(agent_i);
			System.out.print(',');
		}
		for (Agent agent_i : agents)
		{
			System.out.print(agent_i);
			System.out.print(',');
		}
		for (Agent agent_i : agents)
		{
			System.out.print(agent_i);
			System.out.print(',');
		}
		System.out.println();
		
		// Iteration counter:
		System.out.print(k);
		System.out.print(',');
		
		// Delta P:
		for (Agent agent_i : agents)
		{
			System.out.print(((PowerFlowAgentProcess)agent_i.getAgentProcess()).getDeltaS().getReal());
			System.out.print(',');
		}
		
		// Delta Q:
		for (Agent agent_i : agents)
		{
			System.out.print(((PowerFlowAgentProcess)agent_i.getAgentProcess()).getDeltaS().getImaginary());
			System.out.print(',');
		}
		
		// P:
		for (Agent agent_i : agents)
		{
			System.out.print(agent_i.getBus().getNetPower(false).getReal()/IEEE13BusGrid.BASE_POWER);
			System.out.print(',');
		}
		
		// Q:
		for (Agent agent_i : agents)
		{
			System.out.print(agent_i.getBus().getNetPower(false).getImaginary()/IEEE13BusGrid.BASE_POWER);
			System.out.print(',');
		}
		
		// V arg:
		for (Agent agent_i : agents)
		{
			System.out.print(((PowerFlowAgentProcess)agent_i.getAgentProcess()).getV().getArgument());
			System.out.print(',');
		}
		
		// |v|:
		for (Agent agent_i : agents)
		{
			System.out.print(((PowerFlowAgentProcess)agent_i.getAgentProcess()).getV().abs());
			System.out.print(',');
		}
		
		System.out.println();
	}

	public AdmittanceMatrix getY()
	{
		return Y;
	}

	public Grid getGrid()
	{
		return grid;
	}

	public Set<Agent> getAgents()
	{
		return agents;
	}
	
	
	//// Test ////

	/**
	 * Test consensus.
	 * @param args
	 */
	public static void main(String[] args)
	{
		// Analyse grid centrally throguh Newton-Raphson:
		Grid grid = new IEEE13BusGrid();
		LoadFlowAnalyser analyser = new LoadFlowAnalyser(grid);
		analyser.setBasePower(IEEE13BusGrid.BASE_POWER);
		analyser.setBaseVoltage(IEEE13BusGrid.BASE_VOLTAGE);
		AnalysisResults results = analyser.analyse();
		
		// Apply decentralised analysis:
		AgentCoordinator coordinator = new AgentCoordinator(grid, IEEE13BusGrid.BASE_POWER, IEEE13BusGrid.BASE_VOLTAGE);
		coordinator.distributedPowerFlow();
		
		// Compare and log:
		for (Agent agent : coordinator.agents)
		{
			String name = agent.getBus().getName();
			Complex trueV = results.getBusVoltage(name);
			Complex distV = ((PowerFlowAgentProcess)agent.getAgentProcess()).getV();
			System.out.println(name+","+trueV.abs()+","+trueV.getArgument()+","+distV.abs()+","+distV.getArgument());
		}
	}
}