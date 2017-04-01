package ellipsis.energy.util.log;

import java.util.List;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import ellipsis.energy.util.GridlessCentralCoordinator;

public class LoggedGridlessCentralController extends GridlessCentralCoordinator
{
	private GridlessCentralControllerLogger logger = new GridlessCentralControllerLogger(this);
	
	
	//// Decomposition ////
	
	@Override
	public RealVector maxForecastPowerChange()
	{
		logger.log_maxForecastPowerChange_start();
		RealVector maxForecastPowerChange = super.maxForecastPowerChange();
		logger.log_maxForecastPowerChange_end();
		return maxForecastPowerChange;
	}
	
	@Override
	public RealMatrix ebeLambdaDeltaSProduct(RealVector DeltaS_t_max)
	{
		logger.log_ebeLambdaDeltaSProduct_start();
		RealMatrix ebeLambdaDeltaSProduct = super.ebeLambdaDeltaSProduct(DeltaS_t_max);
		logger.log_ebeLambdaDeltaSProduct_end();
		return ebeLambdaDeltaSProduct;
	}
	
	@Override
	public List<List<Integer>> selectSubsets(RealMatrix LambdaProduct)
	{
		logger.log_selectSubsets_start();
		List<List<Integer>> Bs = super.selectSubsets(LambdaProduct);
		logger.log_selectSubsets_end();
		return Bs;
	}
	
	@Override
	public void makeGridlessADP(List<List<Integer>> Bs)
	{
		logger.log_makeGridlessADP_start();
		super.makeGridlessADP(Bs);
		logger.log_makeGridlessADP_end();
	}
	
	@Override
	public void decompose()
	{
		logger.log_decompose_start();
		super.decompose();
		logger.log_decompose_end();
	}
	
	
	////Iterations ////
	
	@Override
	public void optimise(double bw)
	{
		logger.log_optimise_start();
		super.optimise(bw);
		logger.log_optimise_end();
	}
	
	@Override
	public RealVector[] controlSchedule()
	{
		logger.log_controlSchedule_start();
		RealVector[] controlSchedule = super.controlSchedule();
		logger.log_controlSchedule_end();
		return controlSchedule;
	}
	
	@Override
	public void update(RealVector[] u)
	{
		logger.log_updatePowerAndVoltage_start();
		super.update(u);
		logger.log_updatePowerAndVoltage_end();
	}
	
	@Override
	public void loop(int iterations, double bandwidth)
	{
		logger.log_loop_start();
		super.loop(iterations, bandwidth);
		logger.log_loop_end();
	}
	
	@Override
	public RealVector[] iteration(double bandwidth)
	{
		RealVector[] iteration = super.iteration(bandwidth);
		
		// Logging:
		RealVector[] u = controlSchedule();
		double[] ctg = gridlessADP().ctgs(u);
		logger.log_iteration_end(ctg);
		
		return iteration;
	}
}