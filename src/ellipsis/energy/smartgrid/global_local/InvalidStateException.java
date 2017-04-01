package ellipsis.energy.smartgrid.global_local;

public class InvalidStateException extends RuntimeException
{
	private static final long serialVersionUID = 9208814078298833664L;

	public InvalidStateException()
	{
		super();
	}

	public InvalidStateException(String arg0, Throwable arg1)
	{
		super(arg0, arg1);
	}

	public InvalidStateException(String arg0)
	{
		super(arg0);
	}

	public InvalidStateException(Throwable arg0)
	{
		super(arg0);
	}
}
