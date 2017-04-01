package ellipsis.util;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;

public class ContinuousTimeProperties extends ContinuousProperties
{
    private static final long serialVersionUID = -5249303595350811290L;

    /**
     * Key must be a Date value (or string representation of a date in the format HHmm)
     * and will be converted into milliseconds before calling super.put().
     * Value will be treated as per {@link ContinuousProperties#put(Object, Object)}.
     */
    @Override
    public synchronized Object put(Object key, Object value)
    {
        return super.put(convertToMilliseconds(key), value);
    }
    
    /**
     * Key must be a Date or string representation of a date.
     */
    @Override
    public synchronized Object get(Object key)
    {
        return super.get(convertToMilliseconds(key));
    }

    private Long convertToMilliseconds(Object key)
    {
        Date date;
        if(!(key instanceof Date))
        {
            SimpleDateFormat df = new SimpleDateFormat("HHmm");
            try
            {
                date = df.parse(key.toString());
            } 
            catch (ParseException e)
            {
                throw new RuntimeException(e);
            }
        }
        else
        {
            date = (Date)key;
        }
        
        return date.getTime();
    }
    
    public static void main(String[] args) throws IOException
    {
        FileInputStream in = new FileInputStream(new File("/svn/energynet/casestudies/src/ellipsis/energy/casestudies/ieee/wind24.properties"));
        ContinuousTimeProperties p = new ContinuousTimeProperties();
        p.load(in);
        p.sort();
        for (int i = 0; i < 2400;) // 15 minute intervals
        {
            String key = zeroPad(i);
            System.out.println(key+"="+p.get(key));
            if(i%100 == 45)
                i += 55; // Go to next hour
            else
                i += 15; // Next 15 minutes
        }
    }

    private static String zeroPad(int i)
    {
        String ret = ""+i;
        while(ret.length() < 4)
            ret = "0"+ret;
        return ret;
    }
}
