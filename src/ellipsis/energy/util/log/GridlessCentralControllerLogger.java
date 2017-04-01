package ellipsis.energy.util.log;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;

import com.mls.util.Timer;

import ellipsis.energy.util.GridlessCentralCoordinator;

public class GridlessCentralControllerLogger
{
	private static final String LOOP_TIMER = "loop";
	private static final String UPDATE_POWER_AND_VOLTAGE_TIMER = "updatePowerAndVoltage";
	private static final String CONTROL_SCHEDULE_TIMER = "controlSchedule";
	private static final String OPTIMISE_TIMER = "optimise";
	private static final String MAX_FORECAST_POWER_CHANGE_TIMER = "maxForecastPowerChange";
	private static final String EBE_LAMBDA_DELTA_S_PRODUCT_TIMER = "ebeLambdaDeltaSProduct";
	private static final String SELECT_SUBSETS_TIMER = "selectSubsets";
	private static final String MAKE_GRIDLESS_ADP_TIMER = "makeGridlessADP";
	private static final String DECOMPOSE_TIMER = "decompose";
	
	private static File decompositionFile = new File("/opt/energynet/cs03/decomposition.csv");
	private static File iterationsFile = new File("/opt/energynet/cs03/iterations.csv");
	private static File ctgFile = new File("/opt/energynet/cs03/ctg.csv");
	static
	{
		try
		{
			PrintStream out = new PrintStream(new FileOutputStream(decompositionFile, false)); // New file.
			out.println("Epsilon Decomposition Times\nControls,Max Power,Product,Select Subset,Make Grid,Decomposition");
			out.close();
		} 
		catch (IOException e)
		{
			throw new RuntimeException(e);
		}
	}
	static
	{
		try
		{
			PrintStream out = new PrintStream(new FileOutputStream(iterationsFile, false)); // New file.
			out.println("Central Iterations\nControls,Local Optimisation (total),Control Amalgamation,Voltage Updates,Loop");
			out.close();
		} 
		catch (IOException e)
		{
			throw new RuntimeException(e);
		}
	}
	static
	{
		try
		{
			PrintStream out = new PrintStream(new FileOutputStream(ctgFile, false)); // New file.
			out.println("Costs-to-go\nControls,Central Iterations...");
			out.close();
		} 
		catch (IOException e)
		{
			throw new RuntimeException(e);
		}
	}
	
	private GridlessCentralCoordinator cc;
	
	public GridlessCentralControllerLogger(GridlessCentralCoordinator cc)
	{
		this.cc = cc;
	}
	
	
	//// Decomposition Logging ////
	
	public void log_makeGridlessADP_start()
	{
		Timer.getTimer(MAKE_GRIDLESS_ADP_TIMER).start();
	}
	
	public void log_makeGridlessADP_end()
	{
		Timer.getTimer(MAKE_GRIDLESS_ADP_TIMER).stop();
	}

	public void log_selectSubsets_end()
	{
		Timer.getTimer(SELECT_SUBSETS_TIMER).stop();
	}

	public void log_selectSubsets_start()
	{
		Timer.getTimer(SELECT_SUBSETS_TIMER).start();
	}

	public void log_ebeLambdaDeltaSProduct_start()
	{
		Timer.getTimer(EBE_LAMBDA_DELTA_S_PRODUCT_TIMER).start();
	}

	public void log_ebeLambdaDeltaSProduct_end()
	{
		Timer.getTimer(EBE_LAMBDA_DELTA_S_PRODUCT_TIMER).stop();
	}

	public void log_maxForecastPowerChange_start()
	{
		Timer.getTimer(MAX_FORECAST_POWER_CHANGE_TIMER).start();
	}

	public void log_maxForecastPowerChange_end()
	{
		Timer.getTimer(MAX_FORECAST_POWER_CHANGE_TIMER).stop();
	}

	public void log_decompose_start()
	{
		Timer.getTimer(DECOMPOSE_TIMER).start();
	}

	/**
	 * Stops the timer then adds a line to the output log file.
	 */
	public void log_decompose_end()
	{
		Timer.getTimer(DECOMPOSE_TIMER).stop();
		try
		{
			log_decompose();
		} 
		catch (IOException e)
		{
			throw new RuntimeException(e);
		}
	}

	public void log_decompose() throws IOException
	{
		// Output data:
		PrintStream out = new PrintStream(new FileOutputStream(decompositionFile, true)); // Append.

		out.print(cc.gridlessADP().controlDimension());
		out.print(',');
		out.print(Timer.getTimer(MAX_FORECAST_POWER_CHANGE_TIMER).reset());
		out.print(',');
		out.print(Timer.getTimer(EBE_LAMBDA_DELTA_S_PRODUCT_TIMER).reset());
		out.print(',');
		out.print(Timer.getTimer(SELECT_SUBSETS_TIMER).reset());
		out.print(',');
		out.print(Timer.getTimer(MAKE_GRIDLESS_ADP_TIMER).reset());
		out.print(',');
		out.print(Timer.getTimer(DECOMPOSE_TIMER).reset());
		out.print(',');
		out.println();
		
		out.close();
	}
	
	
	//// Iterations Logging ////
	
	public void log_optimise_start()
	{
		Timer.getTimer(OPTIMISE_TIMER).start();
	}

	public void log_optimise_end()
	{
		Timer.getTimer(OPTIMISE_TIMER).stop();
	}

	public void log_controlSchedule_start()
	{
		Timer.getTimer(CONTROL_SCHEDULE_TIMER).start();
	}

	public void log_controlSchedule_end()
	{
		Timer.getTimer(CONTROL_SCHEDULE_TIMER).stop();
	}

	public void log_updatePowerAndVoltage_start()
	{
		Timer.getTimer(UPDATE_POWER_AND_VOLTAGE_TIMER).start();
	}

	public void log_updatePowerAndVoltage_end()
	{
		Timer.getTimer(UPDATE_POWER_AND_VOLTAGE_TIMER).stop();
	}

	public void log_loop_start()
	{
		Timer.getTimer(LOOP_TIMER).start();
	}

	public void log_loop_end()
	{
		Timer.getTimer(LOOP_TIMER).stop();
		try
		{
			log_loop();
		} catch (FileNotFoundException e)
		{
			throw new RuntimeException(e);
		}
	}


	private void log_loop() throws FileNotFoundException
	{
		// Output data for iterations:
		{
			PrintStream out = new PrintStream(new FileOutputStream(iterationsFile, true)); // Append.
	
			out.print(cc.gridlessADP().controlDimension());
			out.print(',');
			out.print(Timer.getTimer(OPTIMISE_TIMER).reset());
			out.print(',');
			out.print(Timer.getTimer(CONTROL_SCHEDULE_TIMER).reset());
			out.print(',');
			out.print(Timer.getTimer(UPDATE_POWER_AND_VOLTAGE_TIMER).reset());
			out.print(',');
			out.print(Timer.getTimer(LOOP_TIMER).reset());
			out.print(',');
			out.println();
			
			out.close();
		}
		
		// Output data for costs-to-go:
		log_loop_ctg();
	}

	
	//// Costs-to-go ////

	private ArrayList<Double> ctgOverIterations = new ArrayList<>();

	public void log_iteration_end(double[] ctg)
	{
		ctgOverIterations.add(ctg[0]); // CTG at t=0.
	}

	public void log_loop_ctg() throws FileNotFoundException
	{
		PrintStream out = new PrintStream(new FileOutputStream(ctgFile, true)); // Append.
		
		out.print(cc.gridlessADP().controlDimension());
		out.print(',');
		for (Double ctg : ctgOverIterations)
		{
			out.print(ctg);
			out.print(',');
		}
		out.println();
		out.close();
	}
}