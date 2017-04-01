package ellipsis.energy.util.log;

import ellipsis.energy.util.ADP;
import ellipsis.energy.util.DiscreteGridlessADP;

public class LoggedGridlessADP extends DiscreteGridlessADP
{
	@Override
	public ADP makeNew()
	{
		return new LoggedGridlessADP();
	}
}