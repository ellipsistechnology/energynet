package ellipsis.energy.sandbox;

import static ellipsis.util.ListUtil.setEach;
import static ellipsis.util.Sum.sum;
import static ellipsis.util.VectorHelper.vector;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import com.mls.util.Util;

import ellipsis.energy.sandbox.sandbox020.Agent;
import ellipsis.energy.sandbox.sandbox020.ConstantCurrentAgent;
import ellipsis.energy.sandbox.sandbox020.ConstantPowerAgent;
import ellipsis.energy.sandbox.sandbox020.VoltageControlledAgent;
import ellipsis.util.TeeOutputStream;

/**
 * Future work:
 *  - Allow for stochastic loads that are variable and not CC or CP.
 *  
 * @author bmillar
 */
public class Sandbox020SmartBuilding
{
    private static final PrintStream OUT;
    
    static
    {
        try
        {
            OUT = new PrintStream(new TeeOutputStream(new FileOutputStream("/tmp/Sandbox020.csv"), System.out));
        }
        catch (FileNotFoundException e)
        {
            throw new RuntimeException(e);
        }
    }
    
    private static final double RESISTIVITY = 0.0132; // ohms/m
    private final static double PS_VOLTAGE = 12;
    private static final double PS_VOLTAGE_MIN = PS_VOLTAGE*0.9;
    private static final double PS_VOLTAGE_MAX = PS_VOLTAGE*1.1;
//    private static final double PS_CURRENT_MAX = 10;
    
    private final static double LED_VOLTAGE = 12;
    private final static double LED_VOLTAGE_MIN = 9;
    private final static double LED_VOLTAGE_MAX = 15;
    
    private final static double CHARGER_VOLTAGE = 12;
    private final static double CHARGER_VOLTAGE_MIN = 10;
    private final static double CHARGER_VOLTAGE_MAX = 18;
    
    private final static int ITERATION_COUNT = 3000;
    private final static int LOG_FREQUENCY = ITERATION_COUNT/1000;
    
    private static Set<Agent> agents;
    
    public static void main(String[] args)
    {
//        agents = testCase01_4bus();
        agents = testCase02();
//        agents = testCase03_2node();
        
        logHeaders(agents);
        logOptimisationStep(-1, agents);
        
        // Optimise:
        solveDistributed();
//        solveCentral();
        
        logResults(agents);
    }

    private static int loopStop = 184; // debug point
    private static String agentStop = "PS1";
    
    public static void solveDistributed()
    {
        for (int k = 0; k < ITERATION_COUNT; k++)
        {
            for (Agent agent : agents)
            {
                if(k == loopStop && agent.getName().equals(agentStop))
                    Util.nullop();
                agent.step();
            }
            
            logOptimisationStep(k, agents);
        }
    }
    
    static Random rand = new Random(0);
    static int index;
    static boolean stochastic = true;
    @SuppressWarnings("deprecation")
    public static void solveCentral()
    {
        // Get optimization parameters from first agent in the list (should all be the same):
        Agent firstAgent = agents.iterator().next();
        double alphaInit = firstAgent.getAlpha(); 
        double alphaMultiplier = firstAgent.getAlphaMultiplier(); 
        double alphaMax = firstAgent.getAlphaMax();
        
        // Init:
        double alpha = alphaInit;
        double epsilon = Agent.epsilonInit;
        
        // Iterate:
        for (int k = 0; k < ITERATION_COUNT; k++)
        {
            if(k == loopStop)
                Util.nullop();
            
            // Random selection for stochastic implementation:
            index = (index+1)%agents.size();//rand.nextInt(agents.size());
            
            // Grad dec WRT v:
            RealVector grad = vector(Agent.STATE_DIMENSION*agents.size(), Double.MAX_VALUE);//grad(k);
            while(grad.getNorm() > epsilon)
            {
                grad = grad(k, stochastic);
                
                centralMinimise(grad);
                
                // Project:
                RealVector projected = project();
                grad = grad(k, stochastic);
                grad = grad.ebeMultiply(projected);
            }
            
            // Dual variable update:
            int i = 0;
            for (Agent agent : agents)
            {
                if(!stochastic || i == index)
                {
                    RealVector g = agent.g();
                    agent.setLambda(Agent.constrain(agent.getLambda().add(g.mapMultiply(alpha*agent.getLambdaMultiplier())), agent.getLambdaMax()));
                }
                ++i;
            }
            
            // Increase alpha:
            alpha *= alphaMultiplier;
            if(alpha > alphaMax)
                alpha = alphaMax;
            for (Agent agent : agents)
                agent.setAlpha(alpha);
            
            // Decrease epsilon:
            epsilon *= Agent.epsilonMultiplier;
            if(epsilon < Agent.epsilonMin)
                epsilon = Agent.epsilonMin;
            for (Agent agent : agents)
                agent.setEpsilon(epsilon);
            
            logOptimisationStep(k, agents);
        }
    }

    private static RealVector project()
    {
        RealVector projected = new ArrayRealVector();
        for (Agent agent : agents)
        {
            projected = projected.append(agent.projectVoltages());
        }
        return projected;
    }

    private static RealVector grad(int k, boolean filter)
    {
        RealVector grad = new ArrayRealVector();
        for (Agent agent : agents)
        {
            if(k == loopStop && agent.getName().equals(agentStop))
                Util.nullop();
            grad = grad.append(agent.lagrangeGrad());
//System.out.print(agent.getVplus()+","+agent.getVminus());
        }
//System.out.println();
        
        // Zero out entries for stochastic gradient descent:
        if(filter)
            grad = zeroNonIndex(grad);
        
        return grad;
    }

    private static void centralMinimise(RealVector grad)
    {
        RealVector v = loadState();
        double stepSize = backtrack(grad, v);
        RealVector step = grad.mapMultiply(-stepSize);
        v = v.add(step);
        storeState(v);
//System.out.println(lagrange(agents));
    }
    
    protected static RealVector zeroNonIndex(RealVector vin)
    {
        RealVector vout = new ArrayRealVector(vin.getDimension());
        int i = index*Agent.STATE_DIMENSION;
        vout.setEntry(i,   vin.getEntry(i));
        vout.setEntry(i+1, vin.getEntry(i+1));
        return vout;
    }

    private static double backtrack(RealVector grad, RealVector vOld)
    {
        double lag = Double.MAX_VALUE;
        final double baseLag = lagrange(agents);
        double t = 1;
        final double gradSquared = grad.dotProduct(grad);
        final double delta = 0.5;
        final double beta = 0.5;
        
        while(lag > baseLag - delta*t*gradSquared)
        {
            t = beta*t;
            RealVector v = vOld.add(grad.mapMultiply(-t));
            storeState(v);
            lag = lagrange(agents);
//double slope = grad(0, true).projection(grad.mapMultiply(-1)).getNorm();
//System.out.println(t+","+lag+","+(baseLag - t*gradSquared)+","+slope+","+grad(0,true).getNorm());
        }
        
        double t2 = 0;//t;
        for(; t2 >= -1e-3; t2 -= 0.5e-3)//t/5.0)
        {
            RealVector v = vOld.add(grad.mapMultiply(-t2));
            storeState(v);
            lag = lagrange(agents);
//double slope = grad(0, true).projection(grad.mapMultiply(-1)).getNorm();
//System.out.println(t2+","+lag+","+(baseLag - t2*gradSquared)+","+slope+","+grad(0,true).getNorm());
        }
        
        storeState(vOld);
        
        return t;
    }

    private static void storeState(RealVector v)
    {
        storeState(v, agents);
    }

	public static void storeState(RealVector v, Set<Agent> agents)
	{
		int i = 0;
		for (Agent agent : agents)
        {
            agent.setVplus(v.getEntry(i++));
            agent.setVminus(v.getEntry(i++));
        }
	}

    private static RealVector loadState()
    {
		return loadState(agents);
    }

	public static RealVector loadState(Set<Agent> agents)
	{
		RealVector v = new ArrayRealVector(agents.size()*Agent.STATE_DIMENSION);
        int i = 0;
        for (Agent agent : agents)
        {
            v.setEntry(i++, agent.getVplus());
            v.setEntry(i++, agent.getVminus());
        }
        return v;
	}
    
    
    //// Test Cases ////

    public static Set<Agent> testCase03_2node()
    {
        // Setup agents:
        Set<Agent> agents = new LinkedHashSet<>();
        
        Agent ps1 = makeVCAgent("PS1", 9);
        ps1.setGrounded(true);
        agents.add(ps1);
        
        Agent cc1 = makeCPAgent("CP1", -10);
        agents.add(cc1);
        
        link(ps1, cc1, RESISTIVITY);
        
        return agents;
    }
    
    public static Set<Agent> testCase02()
    {
        // Setup agents:
        Set<Agent> agents = new LinkedHashSet<>();

        Agent ps1 = makeVCAgent("PS1", 100.0); // FIXME reduce current limit
        ps1.setGrounded(true);
        Agent ps2 = makeVCAgent("PS2", 120.0);
//        ps2.setGrounded(true);
        
        Agent cc1 = makeCCAgent("CC1", 0.0);//-6.0);
        Agent cc2 = makeCCAgent("CC2", 0.0);//-6.0);
        
        Agent cp1 = makeCPAgent("CP1", 0.0);//-40.0);
        Agent cp2 = makeCPAgent("CP2", 0.0);//-60.0);
//        Agent cp3 = makeCPAgent("CP3", 0.0);
//        Agent cp4 = makeCPAgent("CP4", 0.0);
        
        agents.add(ps1);
        agents.add(ps2);
        agents.add(cc1);
        agents.add(cc2);
        agents.add(cp1);
        agents.add(cp2);
//        agents.add(cp3);
//        agents.add(cp4);
        
setEach(agents, Agent::setVplus, 10.8);
setEach(agents, Agent::setVminus, 0.0);
        
        // Setup network connections:
        link(ps1, cc1, RESISTIVITY*1.0);
        link(ps1, cc2, RESISTIVITY*1.0);
        link(ps1, cp1, RESISTIVITY*2.0);
        link(cp1, cp2, RESISTIVITY*1.0);
        link(ps2, cc1, RESISTIVITY*1.0);
        link(ps2, cc2, RESISTIVITY*1.0);
        link(ps2, cp2, RESISTIVITY*2.0);
        
        // Set optimization parameters:
        setEach(agents, Agent::setAlpha,            1.0e-6);//FIXME alpha 1e-6); //2e-0; :-| 1e-3 :-(
        setEach(agents, Agent::setAlphaMax,         2e-6);
        setEach(agents, Agent::setAlphaMultiplier,  1.001); //1.005 :-(
        setEach(agents, Agent::setLambdaMultiplier, 1.0);//0.3);
        setEach(agents, Agent::setMuMultiplier,     1.0);//0.1);
        
        // Initialise:
        agents.forEach(Agent::init);
        
        return agents;
    }

    // FIXME make this a better scenario with loops - otherwise what's the benefit of connecting the PSes?
    public static Set<Agent> testCase02_old()
    {
        // Setup agents:
        Set<Agent> agents = new LinkedHashSet<>();

        Agent ps1 = makeVCAgent("PS1", 10);//12.0);
        ps1.setGrounded(true);
        Agent ps2 = makeVCAgent("PS2", 12);//12.0);
//        Agent ps3 = makeVCAgent("PS3", 100);//10.0);

        Agent cc1 = makeCCAgent("CC1", -12.0);
//        Agent cc2 = makeCCAgent("CC2", -10.0);
        
        Agent cp1 = makeCPAgent("CP1", -100.0);
//        Agent cp2 = makeCPAgent("CP2", 0);//-36.0);
//        Agent cp3 = makeCPAgent("CP3", 0.0);
//        Agent cp4 = makeCPAgent("CP4", 0.0);
        
        agents.add(ps1);
        agents.add(ps2);
//        agents.add(ps3);
        agents.add(cc1);
//        agents.add(cc2);
        agents.add(cp1);
//        agents.add(cp2);
//        agents.add(cp3);
//        agents.add(cp4);
        
        // Setup network connections:
        link(ps1, cc1, RESISTIVITY*1.0);
//        link(ps2, ps3, RESISTIVITY*3);
        link(ps1, cp1, RESISTIVITY*2.0);
        link(ps2, cc1, RESISTIVITY*1.0);
        link(ps2, cp1, RESISTIVITY*2.0);
//        link(cp1, cp2, RESISTIVITY*1);
//        link(cp2, cp3, RESISTIVITY*1);
//        link(cp3, cp4, RESISTIVITY*1);
        
        // Set optimization parameters:
        setEach(agents, Agent::setAlpha,            2e-3);
        setEach(agents, Agent::setAlphaMax,         1e2);
        setEach(agents, Agent::setAlphaMultiplier,  1.0005);
        setEach(agents, Agent::setLambdaMultiplier, 0.8);
        setEach(agents, Agent::setMuMultiplier,     0.5);

        // Initialise:
        agents.forEach(Agent::init);
        
        return agents;
    }

    protected static Agent makeCPAgent(String name, double power)
    {
        Agent cp = new ConstantPowerAgent(power); // Watts (load)
        cp.setName(name);
        cp.setVmax(CHARGER_VOLTAGE_MAX);
        cp.setVmin(CHARGER_VOLTAGE_MIN);
        cp.setVplus(CHARGER_VOLTAGE);
        cp.setVminus(0);
        
        return cp;
    }

    protected static Agent makeCCAgent(String name, double current)
    {
        Agent cc = new ConstantCurrentAgent(current); // Amps (load)
        cc.setName(name);
        cc.setVmax(LED_VOLTAGE_MAX);
        cc.setVmin(LED_VOLTAGE_MIN);
        cc.setVplus(LED_VOLTAGE);
        cc.setVminus(0);
        
        return cc;
    }

    protected static Agent makeVCAgent(String name, double maxCurrent)
    {
        VoltageControlledAgent ps = new VoltageControlledAgent();
        ps.setName(name);
        ps.setVmax(PS_VOLTAGE_MAX);
        ps.setVmin(PS_VOLTAGE_MIN);
        ps.setVplus(PS_VOLTAGE);
        ps.setVminus(0);
        ps.setIMax(maxCurrent);
        return ps;
    }

    protected static Set<Agent> testCase01_4bus()
    {
        // Setup agents:
        Set<Agent> agents = new LinkedHashSet<>();
        
        VoltageControlledAgent ps1 = new VoltageControlledAgent();
        ps1.setName("PS1");
        ps1.setGrounded(true);
        ps1.setVmax(PS_VOLTAGE_MAX);
        ps1.setVmin(PS_VOLTAGE_MIN);
        ps1.setVplus(PS_VOLTAGE);
        ps1.setVminus(0);
        ps1.setIMax(10);
        
        VoltageControlledAgent ps2 = new VoltageControlledAgent();
        ps2.setName("PS2");
        ps2.setVmax(PS_VOLTAGE_MAX);
        ps2.setVmin(PS_VOLTAGE_MIN);
        ps2.setVplus(PS_VOLTAGE);
        ps2.setVminus(0);
        ps2.setIMax(15);
        
        Agent cc1 = new ConstantCurrentAgent(-10.0); // Amps (load)
        cc1.setName("CC1");
        cc1.setVmax(LED_VOLTAGE_MAX);
        cc1.setVmin(LED_VOLTAGE_MIN);
        cc1.setVplus(LED_VOLTAGE);
        cc1.setVminus(0);
        
        Agent cp1 = new ConstantPowerAgent(-150.0); // Watts (load)
        cp1.setName("CP1");
        cp1.setVmax(CHARGER_VOLTAGE_MAX);
        cp1.setVmin(CHARGER_VOLTAGE_MIN);
        cp1.setVplus(CHARGER_VOLTAGE);
        cp1.setVminus(0);
        
        agents.add(ps1);
        agents.add(ps2);
        agents.add(cc1);
        agents.add(cp1);
        
        // Setup network connections:
        link(ps1, cc1, RESISTIVITY*2);
        link(cc1, ps2, RESISTIVITY*3);
        link(ps2, cp1, RESISTIVITY*1);
        link(cp1, ps1, RESISTIVITY*5);
        
        // Set optimization parameters:
        setEach(agents, Agent::setAlpha,           0.5e-2);
        setEach(agents, Agent::setAlphaMax,        1e1);
        setEach(agents, Agent::setAlphaMultiplier, 1.0001);
        
        return agents;
    }

    private static void link(Agent a1, Agent a2, double resistance)
    {
        double conductance = 1.0/resistance;
        a1.addNeighbour(a2, conductance);
        a2.addNeighbour(a1, conductance);
    }

    private static void logHeaders(Set<Agent> agents)
    {
        OUT.print("it,");

        logHeader(agents, "[v+]");
        logHeader(agents, "[v-]");
        logHeader(agents, "[v]" );
        logHeader(agents, "[g1]");
        logHeader(agents, "[g2]");
        logHeader(agents, "[h]" );
        logHeader(agents, "[lambda1]");
        logHeader(agents, "[lambda2]");
        logHeader(agents, "[mu]");
        logHeader(agents, "[alpha]");
        agents.forEach(agent -> OUT.print(agent.getName()+"["+agent.primaryVariable()+"],"));
        logHeader(agents, "[resistance]");
        logHeader(agents, "[current+]");
        logHeader(agents, "[current-]");
        logHeader(agents, "[Lagrange]");
        logHeader(agents, "[grad_v+ L]");
        logHeader(agents, "[grad_v- L]");
        OUT.print("||grad||,");
        logHeader(agents, "[epsilon]");
        OUT.print("Lagrange,");
        logHeader(agents, "[cost]");
        OUT.print("cost,");
        OUT.print("G(v).lambda,");
        OUT.print("(alpha/2)||g(v)||^2");
        
        OUT.println();
    }

    protected static void logHeader(Set<Agent> agents, String label)
    {
        agents.forEach(agent -> OUT.print(agent.getName()+label+","));
    }
    
    private static void logOptimisationStep(int k, Set<Agent> agents)
    {
        if(k%LOG_FREQUENCY != 0 && !forceLog(k))
            return;
        
        OUT.print(k+",");
        
        agents.forEach(a -> OUT.print(a.getVplus()                 +","));
        agents.forEach(a -> OUT.print(a.getVminus()                +","));
        agents.forEach(a -> OUT.print(a.voltage()                  +","));
        agents.forEach(a -> OUT.print(a.g().getEntry(0)            +","));
        agents.forEach(a -> OUT.print(a.g().getEntry(1)            +","));
        agents.forEach(a -> OUT.print(a.h()                        +","));
        agents.forEach(a -> OUT.print(a.getLambda().getEntry(0)    +","));
        agents.forEach(a -> OUT.print(a.getLambda().getEntry(1)    +","));
        agents.forEach(a -> OUT.print(a.getMu()                    +","));
        agents.forEach(a -> OUT.print(a.getAlpha()                 +","));
        agents.forEach(a -> OUT.print(a.primaryValue()             +","));
        agents.forEach(a -> OUT.print(a.resistance()               +","));
        agents.forEach(a -> OUT.print(a.currentPlus()              +","));
        agents.forEach(a -> OUT.print(-a.currentMinus()            +","));
        agents.forEach(a -> OUT.print(a.lagrange()                 +","));
        agents.forEach(a -> OUT.print(a.lagrangeGrad().getEntry(0) +","));
        agents.forEach(a -> OUT.print(a.lagrangeGrad().getEntry(1) +","));
                            OUT.print(grad(0, false).getNorm()     +",");
        agents.forEach(a -> OUT.print(a.getEpsilon()               +","));
                            OUT.print(lagrange(agents)             +",");
        agents.forEach(a -> OUT.print(a.costLocal()                +","));
                            OUT.print(cost(agents)                 +",");
                            OUT.print(lambdaG(agents)              +",");
                            OUT.print(augG(agents)                  +",");
                
        OUT.println();
    }

    private static boolean forceLog(int k)
    {
        return k <= 0;//k >= 1080 && k <= 1090;
    }

    private static double lagrange(Set<Agent> agents)
    {
//        zz__show2dvplusPlot.toString();
        Iterator<Agent> it = agents.iterator();
        double x = it.next().lagrange();
        double x2 = it.next().lagrange();
        return 
                cost(agents) + 
                lambdaG(agents) +
                augG(agents) +
                inequalities(agents);
    }

    private static double augG(Set<Agent> agents)
    {
        return sum(a -> a.getAlpha()*a.g().dotProduct(a.g())/2, agents);
    }

    private static double inequalities(Set<Agent> agents)
    {
        return sum(a -> (a.getMu()+a.getAlpha()*a.h() > 0) ? a.getAlpha()*a.h()*a.h()/2 : -a.getMu()*a.getMu()/(a.getAlpha()*2), agents);
    }

    private static double lambdaG(Set<Agent> agents)
    {
        return sum(a -> a.getLambda().dotProduct(a.g()), agents);
    }

    private static double cost(Set<Agent> agents)
    {
        return sum(Agent::costLocal, agents);
    }

    private static void logResults(Set<Agent> agents)
    {
        // Log Lagrange function around solution:
        double deltaV = 1e-6;
        OUT.println("agent,v^+ - d,v^+,v^+ + d,v^- - d,v^-,v^- + d");
        for (Agent agent : agents)
        {
            OUT.print(agent.getName());
            OUT.print(",");
            
            double oldVPlus = agent.getVplus();
            double oldVMinus = agent.getVminus();
            
            // v+:
            agent.setVplus(oldVPlus-deltaV);
            printLagrange(agents);
            agent.setVplus(oldVPlus);
            printLagrange(agents);
            agent.setVplus(oldVPlus+deltaV);
            printLagrange(agents);
            agent.setVplus(oldVPlus);
            
            // v-:
            agent.setVminus(oldVMinus-deltaV);
            printLagrange(agents);
            agent.setVminus(oldVMinus);
            printLagrange(agents);
            agent.setVminus(oldVMinus+deltaV);
            printLagrange(agents);
            agent.setVminus(oldVMinus);
            
            OUT.println();
        }
    }

    private static void printLagrange(Set<Agent> agents)
    {
        OUT.print(lagrange(agents));
        OUT.print(",");
    }
    
    
    //// Debug ////
    
    @SuppressWarnings("unused")
    private static Object zz__showGrad = new Object()
    {
        public String toString() 
        {
            StringBuffer sb = new StringBuffer();

            for (Agent a : agents)
            {
                RealVector est = a.gradEstimate();
                RealVector grad = a.lagrangeGrad();
                sb.append(a);
                sb.append(grad);
                sb.append(est);
                sb.append('\n');
            }
            
            return sb.toString();
        };
    };
    
    @SuppressWarnings("unused")
    private static Object zz__show2dvplusPlot = new Object()
    {
        public String toString() 
        {
            if(agents.size() != 2)
                return "Only 2 agents are supported for this operation.";
            
            double lowx = 12;
            double highx = 12.5;
            double lowy = 12;
            double highy = 12.5;
            double stepx = (highx-lowx)/10;
            double stepy = (highy-lowy)/10;
            Agent agentX=null, agentY=null;
            
            // Store values and assign agents as x or y:
            Map<Agent, Double> oldValues = new HashMap<>();
            for (Agent a : agents)
            {
                if(agentX == null)
                    agentX = a;
                else
                    agentY = a;
                oldValues.put(a, a.getVplus());
            }
            
            // Overwrite lambda and/or alpha:
            HashMap<Agent, Double> oldAlpha = new HashMap<>();
            HashMap<Agent, RealVector> oldLambda = new HashMap<>();
            for(Agent a : agents)
            {
                oldAlpha.put(a, a.getAlpha());
                oldLambda.put(a, a.getLambda());
//                a.setAlpha(0);
//                a.setLambda(vector(0.0, 0.0));
            }
            
            // Header:
            StringBuffer sb = new StringBuffer();
            for(double x = lowx; x <= highx; x += stepx)
                sb.append(","+x);
            sb.append("\n");
            
            // Values:
            for(double y = lowy; y <= highy; y += stepy)
            {
                agentY.setVplus(y);
                sb.append(y+",");
                for(double x = lowx; x <= highx; x += stepx)
                {
                    agentX.setVplus(x);
                    sb.append(lagrange(agents));
                    sb.append(",");
                }
                sb.append("\n");
            }
            
            for(double y = lowy; y <= highy; y += stepy)
            {
                agentY.setVplus(y);
                sb.append(y+",");
                for(double x = lowx; x <= highx; x += stepx)
                {
                    agentX.setVplus(x);
                    sb.append(cost(agents));
                    sb.append(",");
                }
                sb.append("\n");
            }
            
            for(double y = lowy; y <= highy; y += stepy)
            {
                agentY.setVplus(y);
                sb.append(y+",");
                for(double x = lowx; x <= highx; x += stepx)
                {
                    agentX.setVplus(x);
                    sb.append(lambdaG(agents));
                    sb.append(",");
                }
                sb.append("\n");
            }
            
            for(double y = lowy; y <= highy; y += stepy)
            {
                agentY.setVplus(y);
                sb.append(y+",");
                for(double x = lowx; x <= highx; x += stepx)
                {
                    agentX.setVplus(x);
                    sb.append(augG(agents));
                    sb.append(",");
                }
                sb.append("\n");
            }
            
            
            // Restore values:
            for (Agent a : agents)
            {
                double vp = oldValues.get(a);
                a.setVplus(vp);
                sb.append(a.getName());
                sb.append("=");
                sb.append(vp);
                sb.append("\n");
                
                a.setAlpha(oldAlpha.get(a));
                a.setLambda(oldLambda.get(a));
            }
            
            return sb.toString();
        };
    };
}
