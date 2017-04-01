package ellipsis.energy.mas;

import java.util.LinkedHashMap;
import java.util.Map;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.linear.RealMatrix;

import ellipsis.energy.grid.Bus;

public class Agent
{
	private Bus bus;
	private Map<Agent, Complex> neighbours = new LinkedHashMap<>();
	private Complex selfAdmittance;
	private boolean slackNeighbour = false;
	private RealMatrix sensitivities;
	private AgentProcess process;
	
	public Agent(Bus bus, Complex selfAdmittance, AgentProcess process)
	{
		this.bus = bus;
		this.selfAdmittance = selfAdmittance;
		this.process = process;
	}

	
	//// Accessors ////
	
	@Override
	public String toString()
	{
		return bus.toString();
	}

	public Bus getBus()
	{
		return bus;
	}

	public void addNeighbour(Complex weight, Agent neighbour)
	{
		neighbours.put(neighbour, weight);
		if(neighbour.isSlack())
			slackNeighbour = true;
	}

	public boolean isSlack()
	{
		return !bus.getSlackVoltage().equals(Complex.ZERO);
	}
	
	public Map<Agent, Complex> getNeighbours()
	{
		return neighbours;
	}
	
	public boolean hasSlackNeighbour()
	{
		return slackNeighbour;
	}
	
	public RealMatrix getSensitivities()
	{
		return sensitivities;
	}
	
	public void setSensitivities(RealMatrix sensitivities)
	{
		this.sensitivities = sensitivities;
	}
	
	public AgentProcess getAgentProcess()
	{
		return process;
	}
	
	public Complex getSelfAdmittance()
	{
		return selfAdmittance;
	}
}