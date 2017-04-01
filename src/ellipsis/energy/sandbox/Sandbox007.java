package ellipsis.energy.sandbox;

import org.apache.commons.math3.complex.Complex;

import ellipsis.energy.calculation.AnalysisResults;
import ellipsis.energy.calculation.LoadFlowAnalyser;
import ellipsis.energy.grid.Grid;

/**
 * Based on {@link Sandbox004} but uses a simple 5 bus network.
 * @author bmillar
 *
 */
public class Sandbox007 extends Sandbox004
{
    protected static final String SLACK_SOURCE = "Slack Source";
    protected static final String SLACK_BUS = "Slack Bus";
    protected static final double BASE_POWER = 100e6;
    protected static final double BASE_VOLTAGE = 1e3;
    protected static final double BASE_IMPEDANCE = BASE_VOLTAGE*BASE_VOLTAGE/BASE_POWER;
    
    protected static final String LINE_2_3 = "2-3";
    protected static final String LINE_2_4 = "2-4";
    protected static final String LINE_3_5 = "3-5";
    protected static final String LINE_1_2 = "1-2";
    protected static final String LINE_1_3 = "1-3";
    
    protected static final String LOAD_2 = "Load 2";
    protected static final String LOAD_5 = "Load 5";
    
    private static final String DG_3 = "DG 3";
    private static final String DG_4 = "DG 4";
    
    protected static final String BUS_2 = "Bus 2";
    protected static final String BUS_3 = "Bus 3";
    protected static final String BUS_4 = "Bus 4";
    protected static final String BUS_5 = "Bus 5";

	public static void main(String[] args)
	{
		// Modified from Example6_8:
		Grid grid = 
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
        grid.getSource(DG_3).setPowerOutput(650e6, 0);
        grid.getSource(DG_4).setPmax(650e6);
        grid.getSource(DG_4).setPowerOutput(650e6, 0);
        
//        test(grid);
        
        new Sandbox007().run(grid);
	}

	public static void test(Grid grid)
	{
		grid.getLoad(LOAD_2).setLoad(new Complex(1.5*grid.getBasePower(), 0*grid.getBasePower()));
		grid.getLoad(LOAD_5).setLoad(new Complex(0.8*grid.getBasePower(), 0*grid.getBasePower()));
		
		grid.getSource(DG_3).setPowerOutput(0*grid.getBasePower(), 0.05*grid.getBasePower());
		grid.getSource(DG_4).setPowerOutput(0*grid.getBasePower(), 0.05*grid.getBasePower());
		
		analyse(grid);
	}
	
	public static LoadFlowAnalyser analyse(Grid grid)
    {
        // Set up analysis:
        LoadFlowAnalyser lfa = new LoadFlowAnalyser(grid);
        lfa.setBasePower(grid.getBasePower());
        lfa.setBaseVoltage(grid.getBaseVoltage());
        
        // Analyse:
        AnalysisResults results = lfa.analyse();
        if(!results.getDidConverge())
            System.out.println("WARNING: Alanysis did not converge!\n");
        
        // 6.7 a)
        System.out.println("6.7 a)");
        System.out.println("Voltage:");
        System.out.println(SLACK_BUS+": "+c2polar(results.getBusVoltage(SLACK_BUS))+" p.u.");
        System.out.println(BUS_2+": "+c2polar(results.getBusVoltage(BUS_2))+" p.u.");
        System.out.println(BUS_3+": "+c2polar(results.getBusVoltage(BUS_3))+" p.u.");
        System.out.println(BUS_4+": "+c2polar(results.getBusVoltage(BUS_4))+" p.u.");
        System.out.println(BUS_5+": "+c2polar(results.getBusVoltage(BUS_5))+" p.u.");
        
        // 6.7 b)
        System.out.println("\n6.7 b)");
        System.out.println("Power:");
        System.out.println(SLACK_BUS+": "+c2s(results.getBusPower(SLACK_BUS))+" p.u.");
        System.out.println(BUS_2+": "+c2s(results.getBusPower(BUS_2))+" p.u.");
        System.out.println(BUS_3+": "+c2s(results.getBusPower(BUS_3))+" p.u.");
        System.out.println(BUS_4+": "+c2s(results.getBusPower(BUS_4))+" p.u.");
        System.out.println(BUS_5+": "+c2s(results.getBusPower(BUS_5))+" p.u.");
        
        // 6.7 c)
        System.out.println("\n6.7 c)");
        
        System.out.println(LINE_1_2+":\t"+
                "I12="+c2s(results.getLineCurrent(LINE_1_2))+",\t"+
                "S12="+c2s(results.getLinePower(LINE_1_2))+",\t"+
                "S21="+c2s(results.getLinePowerReverse(LINE_1_2))+",\t"+
                "SL="+c2s(results.getLineLoss(LINE_1_2)));
        
        System.out.println(LINE_1_3+":\t"+
                "I13="+c2s(results.getLineCurrent(LINE_1_3))+",\t"+
                "S13="+c2s(results.getLinePower(LINE_1_3))+",\t"+
                "S31="+c2s(results.getLinePowerReverse(LINE_1_3))+",\t"+
                "SL="+c2s(results.getLineLoss(LINE_1_3)));
        
        System.out.println(LINE_2_3+":\t"+
                "I23="+c2s(results.getLineCurrent(LINE_2_3))+",\t"+
                "S23="+c2s(results.getLinePower(LINE_2_3))+",\t"+
                "S32="+c2s(results.getLinePowerReverse(LINE_2_3))+",\t"+
                "SL="+c2s(results.getLineLoss(LINE_2_3)));
        
        return lfa;
    }
	
	public static String c2s(Complex c)
    {
        if(c.getImaginary() < 0)
            return String.format("%.4f - j%.4f", c.getReal(), (-c.getImaginary()));
        else
            return String.format("%.4f + j%.4f", c.getReal(), c.getImaginary());
    }
	
	public static String c2polar(Complex c)
    {
        return String.format("%.4f/_%.4f", c.abs(), (c.getArgument()));
    }
}