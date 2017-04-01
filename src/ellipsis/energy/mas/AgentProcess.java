package ellipsis.energy.mas;

public interface AgentProcess
{
	boolean stepInit();
	void step();
	void setAgent(Agent agent);
}
