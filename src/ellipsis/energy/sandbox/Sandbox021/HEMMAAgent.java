package ellipsis.energy.sandbox.Sandbox021;

import java.util.LinkedList;
import java.util.Queue;
import java.util.Random;

import ellipsis.energy.sandbox.Sandbox021_AugLagSolver.Problem;
import ellipsis.energy.sandbox.Sandbox021_AugLagSolver.ProblemConfiguration;
import ellipsis.energy.sandbox.sandbox020.Agent;
import static ellipsis.energy.sandbox.Sandbox021.HEMMAAgent.HEMMAState.*;

public class HEMMAAgent implements Runnable
{	
	public static enum HEMMAState
	{
		Idle,
		SessionComplete,
		SessionInitialisation,
		SessionExecution
	}
	
	public static enum HEMMAMessageType
	{
		discoverNeighbours,
		identifyNeighbour,
		startSession,
		cancelSession,
		finishSession,
		variableUpdate,
		
		startSession_reject,
		startSession_accepted,
		finishSession_rejected,
		finishSession_accepted
	}
	
	public static class HEMMAMessage
	{
		HEMMAMessageType type;
		int ttl;
		int sessionId;
		
		public HEMMAMessage(HEMMAMessageType type)
		{
			this.type = type;
		}
		
		public HEMMAMessage(HEMMAMessageType type, int ttl)
		{
			this(type);
			this.ttl = ttl;
		}
		
		public HEMMAMessage(HEMMAMessageType type, int ttl, int sessionId)
		{
			this(type, ttl);
			this.sessionId = sessionId;
		}
	}
	
	private HEMMAState state = Idle;
//	private Agent agent;
//	private Problem problem;
//	private ProblemConfiguration config;
	private Queue<HEMMAMessage> messageQueue = new LinkedList<>();
	
	/**
	 * Used to identify an agents neighbours. 
	 * Can only be sent in the Idle state. 
	 * Prior to sending, the agents neighbour 
	 * list is cleared. Upon receiving a 
	 * DN request an agent will immediately respond
	 * with an IN message addressed to the source 
	 * of the DN request and only on the interface 
	 * on which the DN request was received. DN 
	 * requests are never forwarded (sent with TTL=1).
	 * 
	 * Boradcast.
	 */
	public void discoverNeighbours()
	{
		// TODO send IN
	}
	
	/**
	 * Used to inform neighbouring agents of the senders 
	 * presence. Upon receiving an IN message the sender 
	 * will be added to the recipients neighbour list. 
	 * IN requests are never forwarded (sent with TTL=1).
	 * 
	 * Single recipient.
	 */
	public void identifyNeighbour()
	{
		
	}
	
	/**
	 * Request that a new optimisation session begin. 
	 * This can only be sent if the agent is not currently 
	 * in an optimisation session. Each new request will 
	 * be associated with a randomly generated unique identifier. 
	 * Optimisation will begin once all neighbours have 
	 * accepted the SS request.
	 * 
	 * Neighbours.
	 */
	public void startSession(int sessionId)
	{
		switch(state)
		{
		case Idle:
			state = SessionInitialisation;
		case SessionInitialisation:
			// TODO check if current session is superceded and accept or reject accordingly
			// TODO send SS on to neighbours if accepting the request
			// TODO send accept response if accepting the request
		default:
			break;
		}
	}
	
	/**
	 * Cancel the optimisation associated with a particular 
	 * session ID.
	 * 
	 * Broadcast.
	 */
	public void cancelSession(int sessionId)
	{
		
	}
	
	/**
	 * Request that the session be completed. The optimisation 
	 * will be concluded once all neighbours have accepted the 
	 * FS request.
	 * 
	 * Neighbours.
	 */
	public void finishSession()
	{
		
	}
	
	/**
	 * Exchange state variables with a neighbour. The 
	 * senders optimisation state variables are 
	 * attached to the request and the receiver will 
	 * attach its variables to the response.
	 * 
	 * Single recipient.
	 */
	public void variableUpdate()
	{
		
	}
	
	Random rand = new Random(0);
	public void run()
	{
		init();
		
		while(true)
		{
			switch (state)
			{
			case Idle:
				// TODO handle state change to trigger a session start request
				//startSession(rand.nextInt());
				break;
			case SessionExecution:
				// TODO handle completion criteria met
				// state = HEMMAState.SessionComplete;
				// TODO send FS to neighbours
				break;
			default:
				throw new RuntimeException("Unsupported state found during state engine execution: "+state);
			}
			
			HEMMAMessage message = messageQueue.poll();
			if(message != null)
			{
				switch (message.type)
				{
				case discoverNeighbours:
					discoverNeighbours();
					break;
				case cancelSession:
					cancelSession(message.sessionId);
					break;
				case finishSession:
					
					break;
				case finishSession_accepted:
					
					break;
				case finishSession_rejected:
					
					break;
				case identifyNeighbour:
					
					break;
				case startSession:
					
					break;
				case startSession_accepted:
					
					break;
				case startSession_reject:
					
					break;
				case variableUpdate:
					
					break;
				default:
					break;
				}
			}
		}
	}

	private void init()
	{
		// TODO broadcast DN message
	}
}
