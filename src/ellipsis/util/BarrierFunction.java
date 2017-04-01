package ellipsis.util;

import java.util.List;

public interface BarrierFunction extends TwiceDifferentiableFunction
{
	List<TwiceDifferentiableFunction> getConstraints();
}