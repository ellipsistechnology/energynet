package ellipsis.util;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.mls.util.Util;

public class CSVLogger
{
	public static final int LOG_LEVEL_ERROR = 3;
	public static final int LOG_LEVEL_WARNING = 2;
	public static final int LOG_LEVEL_INFO = 1;
	public static final int LOG_LEVEL_DEBUG = 0;
	
	private int logLevel = LOG_LEVEL_WARNING;
	private PrintStream out = System.out;
	private boolean lineStart = true;
	protected Map<String, List<Object>> columnData;
	private Map<String, Method> methodMapping;

	public void setOut(String filePath)
	{
		try
		{
			File file = new File(filePath);
			file.createNewFile();
			setOut(new PrintStream(file));
		}
		catch (IOException e)
		{
			throw new RuntimeException(e);
		}
	}
	
	
	//// Basic output ////
	
	protected boolean invalidLogLevel(int logLevel)
	{
		return logLevel < this.logLevel;
	}
	
	public void printCell(String value, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		if(lineStart)
			lineStart = false;
		else
			out.print(',');
		out.print(value);
	}
	
	public void endRow(int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		out.println();
		lineStart = true;
	}

	public void printHeading(String heading, int logLevel)
	{
		if(invalidLogLevel(logLevel))
			return;
		
		if(!lineStart)
		{
			out.println(); // Finish previous line
			lineStart = true;
		}
		
		out.println('-'); // Empty line
		out.println(heading);
	}
	
	
	//// Data-set printing ////
	
	public void resetData()
	{
		columnData = new LinkedHashMap<String, List<Object>>();
		methodMapping = new HashMap<String, Method>();
	}
	
	public void resetDataWithColumns(String... columnNames)
	{
		resetData();
		for (int i = 0; i < columnNames.length; i++)
		{
			columnData.put(columnNames[i], new ArrayList<Object>());
		}
	}
	
	public void addData(String column, Object value)
	{
		List<Object> values = columnData.get(column);
		if(values == null)
		{
			values = new ArrayList<Object>();
			columnData.put(column, values);
		}
		values.add(value);
	}
	
	public void addData(Object... data)
	{
		Iterator<String> columns = columnData.keySet().iterator();
		for (int i = 0; i < data.length; i++)
		{
			String column = columns.next();
			Object value = data[i];
			addData(column, value);
		}
	}

	public void fillData(Object defaultValue)
	{
		// Find longest column:
		int length = 0;
		for (String column : columnData.keySet())
		{
			List<Object> values = columnData.get(column);
			int size = values.size();
			if(size > length)
				length = size;
		}

		// Fill all columns that are shorter than the longest:
		for (String column : columnData.keySet())
		{
			List<Object> values = columnData.get(column);
			int size = values.size();
			while(size < length)
			{
				values.add(defaultValue);
				size = values.size();
			}
		}
	}

	public void addColumn(String name)
	{
		columnData.put(name, new ArrayList<Object>());
	}
	
	public void printData(int logLevel)
	{
		printData(logLevel, true);
	}


	private void printData(int logLevel, boolean includeHeadings)
	{
		// Setup data for fast access and print header:
		int size = columnData.size();
		List<?>[] values = new List[size];
		int j = 0;
		int rowCount = 0;
		for (String column : columnData.keySet())
		{
			if(includeHeadings)
				printCell(column, logLevel); // header
			
			values[j] = columnData.get(column);
			rowCount = Math.max(values[j].size(), rowCount);
			++j;
		}
		endRow(logLevel); // end header row
		
		// Print data:
		for (int i = 0; i < rowCount; i++)
		{
			for (j = 0; j < values.length; j++)
			{
				printCell(get(values[j], i), logLevel);
			}
			
			endRow(logLevel);
		}
	}
	
	private String get(List<?> list, int i)
	{
		if(i >= list.size())
			return "";
		
		Object value = list.get(i);
		
		if(value == null)
			return "";
		
		return value.toString();
	}

	public void printDataHorizontally(int logLevel)
	{
		for (String column : columnData.keySet())
		{
			printCell(column, logLevel); // header
			List<Object> values = columnData.get(column);
			for (Object value : values)
			{
				if(value == null)
					printCell("", logLevel);
				else
					printCell(value.toString(), logLevel);
			}
			endRow(logLevel);
		}
	}
	
	public void printDataAndReset(int logLevel)
	{
		printData(logLevel);
		resetData();
	}
	
	public void printDataHorizontallyAndReset(int logLevel)
	{
		printDataHorizontally(logLevel);
		resetData();
	}

	public void printDataAndContinue(int logLevel, boolean includeHeadings)
	{
		printData(logLevel, includeHeadings);
		for (String key : columnData.keySet())
		{
			columnData.get(key).clear();
		}
	}
	
	public void mapColumnsToClass(Class<?> c, String... fieldName)
	{
		Iterator<String> columns = columnData.keySet().iterator();
		for (int i = 0; i < fieldName.length; i++)
		{
			String column = columns.next();
			Method getter = Util.getGetter(c, fieldName[i]);
			methodMapping.put(column, getter);
		}
	}
	
	public void addObject(Object o)
	{
		for (String column : columnData.keySet())
		{
			try
			{
				Method getter = methodMapping.get(column);
				Object value = getter.invoke(o, (Object[])null);
				addData(column, value);
			} 
			catch (IllegalArgumentException e)
			{
				logError(e);
			} 
			catch (IllegalAccessException e)
			{
				logError(e);
			} 
			catch (InvocationTargetException e)
			{
				logError(e);
			}
		}
	}

	public void logError(Exception e)
	{
		throw new RuntimeException(e);
	}

	public int getLogLevel()
	{
		return logLevel;
	}

	public void setLogLevel(int logLevel)
	{
		this.logLevel = logLevel;
	}

	public PrintStream getOut()
	{
		return out;
	}

	public void setOut(PrintStream out)
	{
		this.out = out;
	}
}
