package ellipsis.energy.sandbox;

import org.apache.commons.math3.complex.Complex;

import ellipsis.energy.grid.Grid;

import static ellipsis.energy.sandbox.GridConstants.*;
public class Sandbox011 extends Sandbox009
{
	protected static final String LINE_5_6 = "5-6";
	protected static final String BUS_6 = "Bus 6";
	protected static final String LOAD_6 = "Load 6";

	protected static final String LINE_6_7 = "6-7";
	protected static final String BUS_7 = "Bus 7";
	protected static final String LOAD_7 = "Load 7";
	
	public static void main(String[] args)
	{
		pStepSize = 0.01;
	    qStepSize = 0.01;
	    vabsStepSize = 0.005; // 0.012 seems best
	    vargStepSize = 0.01;
	    lambdaPStepSize = 0.006;
	    lambdaQStepSize = 0.006;
	    
		new Sandbox011().run();
	}
	
	@Override
	public void initGrid()
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
                    
                    Line(LINE_5_6, 1, 0.01*BASE_IMPEDANCE, 0.03*BASE_IMPEDANCE).
                        Bus(BUS_6).
                        	Load(LOAD_6).
                        	Line(LINE_6_7, 1, 0.01*BASE_IMPEDANCE, 0.03*BASE_IMPEDANCE).
	                            Bus(BUS_7).
	                            	Load(LOAD_7).
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
		grid.getLoad(LOAD_6).setLoad(new Complex(150e6, 50e6));
		grid.getLoad(LOAD_7).setLoad(new Complex(80e6, 25e6));
		
        grid.getSource(DG_3).setPmax(650e6);
        grid.getSource(DG_3).setPowerOutput(100e6, 0);
        grid.getSource(DG_4).setPmax(650e6);
        grid.getSource(DG_4).setPowerOutput(100e6, 0);
	}
}