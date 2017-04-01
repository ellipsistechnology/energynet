package ellipsis.energy.sandbox;

import static ellipsis.util.Sum.sum;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;

import com.mls.util.Util;

import ellipsis.energy.sandbox.sandbox019.Agent;
import ellipsis.energy.sandbox.sandbox019.ConstantCurrentAgent;
import ellipsis.energy.sandbox.sandbox019.ConstantPowerAgent;
import ellipsis.energy.sandbox.sandbox019.VoltageControlledAgent;
import ellipsis.util.TeeOutputStream;

/**
 * Future work:
 *  - Allow for stochastic loads that are variable and not CC or CP.
 *  
 * @author bmillar
 */
public class Sandbox019SmartBuilding
{
    private static final PrintStream OUT;
    
    static
    {
        try
        {
            OUT = new PrintStream(new TeeOutputStream(new FileOutputStream("/tmp/Sandbox019.csv"), System.out));
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
    private final static double CHARGER_VOLTAGE_MIN = 6;
    private final static double CHARGER_VOLTAGE_MAX = 18;
    
    private final static int ITERATION_COUNT = 60000;
    private final static int LOG_FREQUENCY = ITERATION_COUNT/1000;
    
    // Debug Point:
    private static int loopStop = 49159;
    private static String agentStop = "PS1";
    private static Set<Agent> agents;
    
    public static void main(String[] args)
    {
        //agents = testCase01_4bus();
        //agents = testCase02();
        agents = testCase03_2node();
        
        // Optimise:
        logHeaders(agents);
        logOptimisationStep(-1, agents);
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
        
        logResults(agents);
    }
    
    protected static Set<Agent> testCase03_2node()
    {
        // Setup agents:
        Set<Agent> agents = new LinkedHashSet<>();
        
        Agent ps1 = makeVCAgent("PS1", 100);
        ps1.setGrounded(true);
        agents.add(ps1);
        
        Agent cc1 = makeCCAgent("CC1", -10);
        agents.add(cc1);
        
        link(ps1, cc1, RESISTIVITY);
        
        return agents;
    }

    protected static Set<Agent> testCase02()
    {
        // Setup agents:
        Set<Agent> agents = new LinkedHashSet<>();

        Agent ps1 = makeVCAgent("PS1", 12.0);
        ps1.setGrounded(true);
        Agent ps2 = makeVCAgent("PS2", 12.0);
        Agent ps3 = makeVCAgent("PS3", 10.0);

        Agent cc1 = makeCCAgent("CC1", -10.0);
        Agent cc2 = makeCCAgent("CC2", -10.0);
        
        Agent cp1 = makeCPAgent("CP1", -36.0);
        Agent cp2 = makeCPAgent("CP2", -36.0);
        Agent cp3 = makeCPAgent("CP3", 0.0);
        Agent cp4 = makeCPAgent("CP4", 0.0);
        
        agents.add(ps1);
        agents.add(ps2);
        agents.add(ps3);
        agents.add(cc1);
        agents.add(cc2);
        agents.add(cp1);
        agents.add(cp2);
        agents.add(cp3);
        agents.add(cp4);
        
        // Setup network connections:
        link(ps1, ps2, RESISTIVITY*2.5);
        link(ps2, ps3, RESISTIVITY*3);
        link(ps1, cc1, RESISTIVITY*0.1);
        link(ps2, cc2, RESISTIVITY*0.1);
        link(ps3, cp1, RESISTIVITY*1);
        link(cp1, cp2, RESISTIVITY*1);
        link(cp2, cp3, RESISTIVITY*1);
        link(cp3, cp4, RESISTIVITY*1);
        
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
        ps1.setIMax(15);//PS_CURRENT_MAX);
        
        VoltageControlledAgent ps2 = new VoltageControlledAgent();
        ps2.setName("PS2");
        ps2.setVmax(PS_VOLTAGE_MAX);
        ps2.setVmin(PS_VOLTAGE_MIN);
        ps2.setVplus(PS_VOLTAGE);
        ps2.setVminus(0);
        ps2.setIMax(15);//PS_CURRENT_MAX);
        
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
        logHeader(agents, "[z]" );
        logHeader(agents, "[g1]");
        logHeader(agents, "[g2]");
        logHeader(agents, "[h]" );
        logHeader(agents, "[h+z^2]" );
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
        logHeader(agents, "[grad_z L]");
        logHeader(agents, "[epsilon]");
        OUT.print("Lagrange,");
        logHeader(agents, "[cost]");
        OUT.print("cost,");
        
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
        agents.forEach(a -> OUT.print(a.getZ()                     +","));
        agents.forEach(a -> OUT.print(a.g().getEntry(0)            +","));
        agents.forEach(a -> OUT.print(a.g().getEntry(1)            +","));
        agents.forEach(a -> OUT.print(a.h()                        +","));
        agents.forEach(a -> OUT.print((a.h()+a.getZ()*a.getZ())    +","));
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
        agents.forEach(a -> OUT.print(a.lagrangeGrad().getEntry(2) +","));
        agents.forEach(a -> OUT.print(a.getEpsilon()               +","));
                            OUT.print(lagrange(agents)             +",");
        agents.forEach(a -> OUT.print(a.costLocal()                +","));
                            OUT.print(cost(agents)                 +",");
                
        OUT.println();
    }

    private static boolean forceLog(int k)
    {
        return k >= 49100 && k <= 49200;
    }

    private static double lagrange(Set<Agent> agents)
    {
//        zz__show2dvplusPlot.toString();
        return 
                cost(agents) + 
                sum(a -> a.getLambda().dotProduct(a.g()), agents) +
                sum(a -> a.getAlpha()*a.g().dotProduct(a.g())/2, agents);
    }

    private static double cost(Set<Agent> agents)
    {
        return sum(Agent::costLocal, agents);
    }

    private static void logResults(Set<Agent> agents)
    {
        // Log Lagrange function around solution:
        double deltaV = 1e-6;
        double deltaZ = 1;
        OUT.println("agent,v^+ - d,v^+,v^+ + d,v^- - d,v^-,v^- + d,z - d,z,z + d");
        for (Agent agent : agents)
        {
            OUT.print(agent.getName());
            OUT.print(",");
            
            double oldVPlus = agent.getVplus();
            double oldVMinus = agent.getVminus();
            double oldZ = agent.getZ();
            
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

            // z:
            agent.setZ(oldZ-deltaZ);
            printLagrange(agents);
            agent.setZ(oldZ);
            printLagrange(agents);
            agent.setZ(oldZ+deltaZ);
            printLagrange(agents);
            agent.setZ(oldZ);
            
            OUT.println();
        }
    }

    private static void printLagrange(Set<Agent> agents)
    {
        OUT.print(lagrange(agents));
        OUT.print(",");
    }
    
    
    //// Debug ////
    
    private static Object zz__show2dvplusPlot = new Object()
    {
        public String toString() 
        {
            if(agents.size() != 2)
                return "Only 2 agents are supported for this operation.";
            
            double lowx = 11;
            double highx = 11.5;
            double lowy = 11;
            double highy = 11.5;
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
            
            // Restore values:
            for (Agent a : agents)
            {
                double vp = oldValues.get(a);
                a.setVplus(vp);
                sb.append(a.getName());
                sb.append("=");
                sb.append(vp);
                sb.append("\n");
            }
            
            return sb.toString();
        };
    };
}
