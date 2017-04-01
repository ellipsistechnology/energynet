package ellipsis.util;

import java.io.IOException;
import java.io.OutputStream;

public class TeeOutputStream extends OutputStream
{
	OutputStream out1;
	OutputStream out2;
	
	public TeeOutputStream(OutputStream out1, OutputStream out2)
	{
		this.out1 = out1;
		this.out2 = out2;
	}

	@Override
	public void write(int b) throws IOException
	{
		out1.write(b);
		out2.write(b);
	}

	@Override
	public void write(byte[] b) throws IOException
	{
		out1.write(b);
		out2.write(b);
	}

	@Override
	public void write(byte[] b, int off, int len)
			throws IOException
	{
		out1.write(b, off, len);
		out2.write(b, off, len);
	}
}