package ellipsis.energy.smartgrid.global_local;

public class NoSolutionException extends RuntimeException
{
	private static final long serialVersionUID = -955650238896140850L;

	public NoSolutionException()
	{
		super();
	}

	public NoSolutionException(String arg0)
	{
		super(arg0);
	}

	public NoSolutionException(Throwable arg0)
	{
		super(arg0);
	}

	public NoSolutionException(String arg0, Throwable arg1)
	{
		super(arg0, arg1);
	}
}
