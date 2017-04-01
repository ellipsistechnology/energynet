package ellipsis.energy.grid;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.complex.ComplexUtils;

import com.mls.io.CSVReader;

import ellipsis.energy.grid.res.SolarPanel;
import ellipsis.energy.smartgrid.ControllableDemand;

public class GridGenerator
{
	public double slackVoltage = 4160; // V
	
	public int maxLineCountFromSlack = 8;
	public int minLineCountFromSlack = 4;
	public int maxLineCount = 3;
	public int maxBusDepth = 6;
	public int minBusDepth = 1;
	
	public double minLineResistance = 0.06;
	public double maxLineResistance = 0.06;
	public double minLineInductance = 0.03;
	public double maxLineInductance = 0.03;
	public double minLineLength = 1; // km
	public double maxLineLength = 5; // km
	
	public double loadDensity = 0.9; // Probability of having a load on a bus
	public int maxLoadsPerBus = 1; // Loads per bus
	public int minLoadsPerBus = 1; // Loads per bus
	public double maxInitialLoad = 1000; // VA
	public double minInitialLoad = 100; // VA
	public double maxInitialLoadArg = 0.32; // radians
	public double minInitialLoadArg = 0; // radians
	
	public double pvDensity = 0.25; // Probability of having a PV on a bus
	public int maxPvPerBus = 1; // PVs per bus
	public int minPvPerBus = 1; // PVs per bus
	public double maxPvCapacity = 30e3; // VA
	public double minPvCapacity = 1500; // VA
	public double minPvCurrentToIrradianceRatio = 0.0075; // A/(W/m^2)
	public double maxPvCurrentToIrradianceRatio = 0.0075; // A/(W/m^2)
	public double minPvVoltage = 24; // V
	public double maxPvVoltage = 24; // V
	public double minPvUnitCount = 1000;
	public double maxPvUnitCount = 12000;
	public double initialIrradiance = 1000; // W/m^2
	
	public double storageDensity = 0.25; // Probability of having storage on a bus
	public int maxStoragePerBus = 1; // Storage units per bus
	public int minStoragePerBus = 1; // Storage units per bus
	public int maxStorageTotal = Integer.MAX_VALUE;
	public double maxStorageCapacity = 50000; // Wh
	public double minStorageCapacity = 2000; // Wh
	public double initialStorageCapacityPercent = 0.5; // [0,1]
	public double maxStorageChargeRate = 2000; // W
	public double minStorageChargeRate = 1000; // W
	
	private int busCounter = 1;
	private int storageCounter = 0;
	
	public Grid createGrid(String name)
	{
		return createGrid(name, 0);
	}
	
	public Grid createGrid(String name, int seed)
	{
		busCounter = 1;
		Random rand = new Random(seed);
		
		// Grid:
		Grid grid = new Grid();
		grid.setName(name);
		
		// Start with the slack bus and bus 1:
		Bus slackBus = new Bus(grid);
		slackBus.setName("Slack Bus");
		SlackSource slack = new SlackSource();
		slack.setVoltage(new Complex(slackVoltage, 0));
		slack.setName("Slack");
		slackBus.addChild(slack);
		grid.add(slackBus);
		
		Bus bus1 = createBus(grid, rand);
		Line line = new Line();
		line.setName(slackBus.getName()+"-"+bus1.getName());
		line.setFromBus(slackBus);
		line.setToBus(bus1);
		line.setResistancePerMetre(rand(rand, minLineResistance, maxLineResistance));
		line.setInductancePerMetre(rand(rand, minLineInductance, maxLineInductance));
		line.setLength(rand(rand, minLineLength, maxLineLength));
		grid.add(bus1);
		
		// Create the tree of busses and lines:
		createBusTree(grid, bus1, rand, 1);
		
		return grid;
	}

	private void createBusTree(Grid grid, Bus fromBus, Random rand, int depth)
	{
		if(depth > maxBusDepth)
			return;
		
		int minLineCount = 
				depth == 1 ? minLineCountFromSlack :
				depth < minBusDepth ? 1 : 
				0;
		int max = depth == 1 ? maxLineCountFromSlack : maxLineCount;
		int childCount = (int)(minLineCount + rand.nextDouble()*(max-minLineCount));
		for(int i = 0; i < childCount; ++i)
		{
			// Next bus:
			Bus toBus = createBus(grid, rand);
			
			// Line:
			Line line = new Line();
			line.setName(fromBus.getName()+"-"+toBus.getName());
			line.setFromBus(fromBus);
			line.setToBus(toBus);
			line.setResistancePerMetre(rand(rand, minLineResistance, maxLineResistance));
			line.setInductancePerMetre(rand(rand, minLineInductance, maxLineInductance));
			line.setLength(rand(rand, minLineLength, maxLineLength));
			
			fromBus.addChild(line);
			
			createBusTree(grid, toBus, rand, depth+1);
		}
	}

	protected Bus createBus(Grid grid, Random rand)
	{
		Bus bus = new Bus(grid);
		bus.setName("Bus"+busCounter);
		++busCounter;
		
		// Load:
		if(rand.nextDouble() < loadDensity)
		{
			int count = (int)rand(rand, minLoadsPerBus, maxLoadsPerBus);
			for(int i = 0; i < count; ++i)
			{
				Load load = new Load();
				load.setName("Load-"+bus.getName()+"-"+(i+1));
				
				double magnitude = rand(rand, minInitialLoad, maxInitialLoad);
				double arg = rand(rand, minInitialLoadArg, maxInitialLoadArg);
				load.setLoad(ComplexUtils.polar2Complex(magnitude, arg));
				bus.addChild(load);
			}
		}
		
		// PV:
		if(rand.nextDouble() < pvDensity)
		{
			int count = (int)(minPvPerBus + rand.nextDouble()*(maxPvPerBus-minPvPerBus));
			for(int i = 0; i < count; ++i)
			{
				SolarPanel pv = new SolarPanel();
				pv.setName("PV-"+bus.getName()+"-"+(i+1));
				
				pv.setCurrentToIrradianceRatio(rand(rand, minPvCurrentToIrradianceRatio, maxPvCurrentToIrradianceRatio));
				pv.setIrradiance(initialIrradiance);
				pv.setPmax(rand(rand, minPvCapacity, maxPvCapacity));
				pv.setVoltage(rand(rand, minPvVoltage, maxPvVoltage));
				pv.setPvUnitCount(rand(rand, minPvUnitCount, maxPvUnitCount));
				
				bus.addChild(pv);
			}
		}
		
		// Storage:
		if(rand.nextDouble() < storageDensity)
		{
			int count = (int)(rand(rand, minStoragePerBus, maxStoragePerBus));
			for(int i = 0; i < count; ++i)
			{
				if(storageCounter >= maxStorageTotal)
					break;
				
				ControllableDemand bat = new ControllableDemand();
				bat.setName("Storage-"+bus.getName()+"-"+(i+1));
				
				double capacity = rand(rand, minStorageCapacity, maxStorageCapacity);
				bat.setMaxCapacity(capacity);
				bat.setCapacity(capacity*initialStorageCapacityPercent);
				
				double chargeDischarge = rand(rand, minStorageChargeRate, maxStorageChargeRate);
				bat.setMaxChargeRate(chargeDischarge);
				bat.setMaxDischargeRate(chargeDischarge);
				bat.setChargeRate(new Complex(rand.nextBoolean() ? chargeDischarge : -chargeDischarge));
				
				bus.addChild(bat);
				++storageCounter;
			}
		}
		
		return bus;
	}

	protected double rand(Random rand, double min, double max)
	{
		return min + rand.nextDouble()*(max-min);
	}
	
	
	//// Load Data To Generate Grid ////
	
	public static class LineConfiguration
	{
		public int id;
		public Complex impedance;
	}
	
	public static final String BUS_PREFIX = "";
	public static final String LOAD_PREFIX = "L";
	public static final String LINE_JOIN = "-";
	public static final String SWITCH_JOIN = "/";
	public static final String LINE_CONFIGURATION_PREFIX = "Configuration";
	public static final String LINE_CONFIGURATION_HEADER = "Z (R +jX) in ohms per mile";

	/**
	 * Reads grid data from the given streams and creates a Grid. Note that this will not
	 * create a slack source.
	 * @param name
	 * @param lineConfiguration
	 * @param lineData
	 * @param switchData
	 * @param loadData
	 * @return
	 */
    public static Grid loadGrid(
            String name, 
            InputStream lineConfiguration, 
            InputStream lineData,
            InputStream switchData,
            InputStream loadData)
    {
        Map<Integer, LineConfiguration> config = readLineConfiguration(lineConfiguration);
        Grid grid = new Grid();
        grid.setName(name);
        readLineData(lineData, grid, config);
        readSwitchData(switchData, grid);
        readLoadData(loadData, grid);
        return grid;
    }
	
	/**
	 * Reads line data and creates lines and busses as needed, and adds them to the grid.
	 * @param in
	 * @param config
	 * @param grid
	 * @return
	 */
	public static void readLineData(InputStream in, Grid grid, Map<Integer, LineConfiguration> config)
	{
		CSVReader csv = new CSVReader();
		csv.openStream(in);
		
		// Parse header and find column order:
		String[] header = csv.readRow();
		int nodeAIndex = -1;
		int nodeBIndex = -1;
		int lengthIndex = -1;
		int configIndex = -1;
		for (int i = 0; i < header.length; i++)
		{
			String heading = header[i].trim();
            if("Node A".equals(heading))
				nodeAIndex = i;
			else if("Node B".equals(heading))
				nodeBIndex = i;
			else if("Length (ft.)".equals(heading))
				lengthIndex = i;
			else if("Config.".equals(heading))
				configIndex = i;
		}
		if(nodeAIndex == -1 || nodeBIndex == -1 || lengthIndex == -1 || configIndex == -1)
			throw new RuntimeException("Incorrect heading in line configuration data.");
		
		String[] row = csv.readRow();
		while(row != null)
		{
			try
			{
			    // Parse row:
    			int nodeA = Integer.parseInt(row[nodeAIndex]);
    			int nodeB = Integer.parseInt(row[nodeBIndex]);
    			double length = Double.parseDouble(row[lengthIndex]);
    			int configId = Integer.parseInt(row[configIndex]);
    			
    			// Get busses:
    			String busAName = BUS_PREFIX+nodeA;
    			String busBName = BUS_PREFIX+nodeB;
    			Bus busA = grid.getOrCreateBus(busAName);
    			Bus busB = grid.getOrCreateBus(busBName);
    			
    			// Create and add line:
    			Line line = new Line();
    			line.setName(busAName+LINE_JOIN+busBName);
    			line.setImpedencePerMetre(config.get(configId).impedance);
    			line.setLength(feetToMetres(length));
    			line.setFromBus(busA);
    			line.setToBus(busB);
			}
			catch(NumberFormatException e){} // Just ignore poorly formed lines.

            // Next row:
            row = csv.readRow();
		}
	}
	
	private static double feetToMetres(double feet)
	{
		return feet/3.2808;
	}

	/**
	 * e.g.
	 * 
     *       Configuration 1:
     *       
     *                  Z (R +jX) in ohms per mile
     *        0.4576  1.0780   0.1560  0.5017   0.1535  0.3849
     *                         0.4666  1.0482   0.1580  0.4236
     *                                          0.4615  1.0651
     *
	 * @param lineConfiguration
	 * @return
	 */
    private static Map<Integer, LineConfiguration> readLineConfiguration(InputStream in)
	{
		BufferedReader reader = new BufferedReader(new InputStreamReader(in));
		Map<Integer, LineConfiguration> config = new HashMap<>();
		while(true) // continues until the end of the stream
		{
            // Find next configuration and get its ID:
    		String line;
    		do
    		{
    		    line = readNextNonEmptyLine(reader);
    		} while(line != null && !line.trim().startsWith(LINE_CONFIGURATION_PREFIX));
    		if(line == null)
    		    break;
    		
    		line = line.trim();
    		int idStartIndex = line.indexOf(LINE_CONFIGURATION_PREFIX) + LINE_CONFIGURATION_PREFIX.length()+1;
    		int colonIndex = line.lastIndexOf(':');
    		int id = Integer.parseInt(line.substring(idStartIndex, colonIndex));
    		
    		// Skip over header:
    		do
            {
                line = readNextNonEmptyLine(reader);
            } while(line != null && !line.trim().startsWith(LINE_CONFIGURATION_HEADER));
            if(line == null)
                break;
    		
    		// Line 1:
            line = readNextNonEmptyLine(reader);
            if(line == null)
                break;
    		double[] ds = lineToNumbers(line);
    		Complex Z_AA = new Complex(ds[0], ds[1]);
//    		Complex Z_AB = new Complex(ds[2], ds[3]);
//    		Complex Z_AC = new Complex(ds[4], ds[5]);
            
            line = readNextNonEmptyLine(reader);
            if(line == null)
                break;
            ds = lineToNumbers(line);
            Complex Z_BB = new Complex(ds[0], ds[1]);
//            Complex Z_BC = new Complex(ds[2], ds[3]);
            
            line = readNextNonEmptyLine(reader);
            if(line == null)
                break;
            ds = lineToNumbers(line);
            Complex Z_CC = new Complex(ds[0], ds[1]);
            
            // Create config:
            LineConfiguration c = new LineConfiguration();
            c.impedance = perMileToPerMetre(Z_AA.add(Z_BB).add(Z_CC).divide(3)); // FIXME averaging for now since I don't really understand these values
            c.impedance = c.impedance.divide(30);//200); //(30);// FIXME just trying to get sensible results here
            config.put(id, c);
		}
		
		return config;
	}

    protected static String readNextNonEmptyLine(BufferedReader reader)
    {
        String line;
        do
        { 
            line = readLoud(reader); 
        } while(line != null && line.trim().isEmpty());
        
        return line;
    }

    private static Complex perMileToPerMetre(Complex c)
    {
        return c.divide(1609.344);
    }

    private static double[] lineToNumbers(String line)
    {
        String[] toks = line.trim().split("\\s+");
        double[] ds = new double[toks.length];
        for (int i = 0; i < ds.length; i++)
        {
            ds[i] = Double.parseDouble(toks[i]);
        }
        return ds;
    }

    protected static String readLoud(BufferedReader reader)
    {
        try
        {
            return reader.readLine();
        }
        catch (IOException e)
        {
            throw new RuntimeException(e);
        }
    }

	private static void readSwitchData(InputStream in, Grid grid)
	{
		CSVReader csv = new CSVReader();
		csv.openStream(in);
		
		// Parse header and find column order 
		// (Node A	Node B	Normal):
		String[] header = csv.readRow();
		int nodeAIndex = -1;
		int nodeBIndex = -1;
		int normalIndex = -1;
		for (int i = 0; i < header.length; i++)
		{
			String heading = header[i].trim();
            if("Node A".equals(heading))
				nodeAIndex = i;
			else if("Node B".equals(heading))
				nodeBIndex = i;
			else if("Normal".equals(heading))
				normalIndex = i;
		}
		if(    nodeAIndex    == -1
			|| nodeBIndex    == -1
			|| normalIndex   == -1)
			throw new RuntimeException("Incorrect heading in switch data.");
		
		// Read and parse load data:
		String[] row = csv.readRow();
		while(row != null)
		{
			if(row.length == 3)
			{
    			// Parse row:
    			int nodeA = Integer.parseInt(row[nodeAIndex]);
    			int nodeB = Integer.parseInt(row[nodeBIndex]);
    			boolean closed = row[normalIndex].trim().equals("closed");
    			
    			// Get busses:
    			String busAName = BUS_PREFIX+nodeA;
    			String busBName = BUS_PREFIX+nodeB;
    			Bus busA = grid.getOrCreateBus(busAName);
    			Bus busB = grid.getOrCreateBus(busBName);
    			
    			// Create switch:
    			Switch sw = new Switch();
    			sw.setName(nodeA+SWITCH_JOIN+nodeB);
    			sw.setFromBus(busA);
    			sw.setToBus(busB);
    			sw.setOn(closed);
			}

            // Read next row:
            row = csv.readRow();
		}
	}

	private static void readLoadData(InputStream in, Grid grid)
	{
		CSVReader csv = new CSVReader();
		csv.openStream(in);
		
		// Parse header and find column order 
		// (Node	Load	Ph-1	Ph-1	Ph-2	Ph-2	Ph-3	Ph-4):
		String[] header = csv.readRow();
		int nodeIndex = -1;
		int loadIndex = -1;
		int Ph1ActiveIndex = -1;
		int Ph1ReactiveIndex = -1;
		int Ph2ActiveIndex = -1;
		int Ph2ReactiveIndex = -1;
		int Ph3ActiveIndex = -1;
		int Ph3ReactiveIndex = -1;
		for (int i = 0; i < header.length; i++)
		{
			String heading = header[i].trim();
            if("Node".equals(heading))
				nodeIndex = i;
			else if("Load".equals(heading))
				loadIndex = i;
			else if("Ph-1".equals(heading))
			{
				if(Ph1ActiveIndex == -1)
					Ph1ActiveIndex = i;
				else
					Ph1ReactiveIndex = i;
			}
			else if("Ph-2".equals(heading))
			{
				if(Ph2ActiveIndex == -1)
					Ph2ActiveIndex = i;
				else
					Ph2ReactiveIndex = i;
			}
			else if("Ph-3".equals(heading))
			{
				if(Ph3ActiveIndex == -1)
					Ph3ActiveIndex = i;
				else
					Ph3ReactiveIndex = i;
			}
		}
		if(    nodeIndex        == -1
			|| loadIndex        == -1
			|| Ph1ActiveIndex   == -1
			|| Ph1ReactiveIndex == -1
			|| Ph2ActiveIndex   == -1
			|| Ph2ReactiveIndex == -1
			|| Ph3ActiveIndex   == -1
			|| Ph3ReactiveIndex == -1)
			throw new RuntimeException("Incorrect heading in load data.");
		
		// Read second header (WARNING: This assumes kW comes before kVAr):
		csv.readRow();
		
		// Read and parse load data:
		String[] row = csv.readRow();
		while(row != null)
		{
		    try
		    {
    			// Parse row:
    			int node = Integer.parseInt(row[nodeIndex]);
    			double ph1Active = Double.parseDouble(row[Ph1ActiveIndex]);
    			double ph1Reactive = Double.parseDouble(row[Ph1ReactiveIndex]);
    			double ph2Active = Double.parseDouble(row[Ph2ActiveIndex]);
    			double ph2Reactive = Double.parseDouble(row[Ph2ReactiveIndex]);
    			double ph3Active = Double.parseDouble(row[Ph3ActiveIndex]);
    			double ph3Reactive = Double.parseDouble(row[Ph3ReactiveIndex]);
    			
    			// Get load power:
    			Complex p1 = new Complex(ph1Active, ph1Reactive);
    			Complex p2 = new Complex(ph2Active, ph2Reactive);
    			Complex p3 = new Complex(ph3Active, ph3Reactive);
    			
    			// Create load:
    			Load load = new Load();
    			load.setName(LOAD_PREFIX+node);
    			load.setLoad(p1.add(p2).add(p3).multiply(1e3)); // data is given in kW
    			grid.getOrCreateBus(BUS_PREFIX+node).addChild(load);
		    }
		    catch(NumberFormatException e) {} // Ignore poorly formed rows.
			
			// Next row:
			row = csv.readRow();
		}
	}
}
