package ellipsis.util;

import java.io.IOException;
import java.io.OutputStream;

/**
 * Swallows any data and does nothing with it.
 * @author bmillar
 *
 */
public class EmptyOutputStream extends OutputStream
{
    @Override
    public void write(int b) throws IOException
    {
        // Do nothing.
    }
}