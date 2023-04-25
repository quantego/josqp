package com.quantego.josqp;


import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.FileHandler;
import java.util.logging.Formatter;
import java.util.logging.LogRecord;
import java.util.logging.SimpleFormatter;


public class Parser {

	int nRows;
	int nCols;
	int Anz;
	int Pnz;
	public int[] Ap;
	public int[] Ai;
	public double[] Ax;
	public int[] Pp;
	public int[] Pi;
	public double[] Px;
	public double[] q;
	public double[] l;
	public double[] u;

	/**
	 * Create a new Problem instance based on vector inputs (CSC formatted).
	 * @param q linear objective coefficients
	 * @param P_x quadratic objective (coefficients)
	 * @param P_p quadratic objective (next column starts)
	 * @param P_i quadratic objective (row indices)
	 * @param A_x constraints (coefficients)
	 * @param A_p constraints (next column starts)
	 * @param A_i constraints (row indices)
	 * @param l constraint lower bounds
	 * @param u constraint upper bounds
	 */
	public Parser(double[] q, double[] P_x, int[] P_p, int[] P_i, double[] A_x, int[] A_p, int[] A_i, double[] l, double[] u) {
		this.q = q;
		this.Px = P_x;
		this.Pp = P_p;
		this.Pi = P_i;
		this.Ax = A_x;
		this.Ap = A_p;
		this.Ai = A_i;
		this.l = l;
		this.u = u;
		nCols = q.length;
		nRows = l.length;
		Anz = Ax.length;
		Pnz = Px.length;
	}

	/**
	 * Read problem as a collection of files containing the individual vectors as CSV (without header but a row index).
	 * Requires files 'prefix q', 'prefix P_x', 'prefix P_p', 'prefix P_i', 'prefix l', 'prefix u', , 'prefix A_x', 'prefix A_p', 'prefix A_i' without (file ending).
	 * @param prefix
	 * @return
	 */
	public static Parser readVectors(String prefix, int skipRows, int skipCols) {
		Parser p = new Parser(
				Utils.readDoubleColumn(prefix+ " q", ",", skipRows, skipCols),
				Utils.readDoubleColumn(prefix+ " P_x", ",", skipRows, skipCols),
				Utils.readIntColumn(prefix+ " P_p", ",", skipRows, skipCols),
				Utils.readIntColumn(prefix+ " P_i", ",", skipRows, skipCols),
				Utils.readDoubleColumn(prefix+ " A_x", ",", skipRows, skipCols),
				Utils.readIntColumn(prefix+ " A_p", ",", skipRows, skipCols),
				Utils.readIntColumn(prefix+ " A_i", ",", skipRows, skipCols),
				Utils.readDoubleColumn(prefix+ " l", ",", skipRows, skipCols),
				Utils.readDoubleColumn(prefix+ " u", ",", skipRows, skipCols)
		);
		return p;
	}

	/**
	 * Prints this problem as vectors that can be directly dumped into C code (probably easiest to test small problems before writing a parser in C).
	 */
	public void printToC() {
		System.out.println(String.format("c_int A_p[%d] = %s;",Ap.length,Arrays.toString(Ap).replace("[", "{").replace("]", "}")));
		System.out.println(String.format("c_int A_i[%d] = %s;",Ai.length,Arrays.toString(Ai).replace("[", "{").replace("]", "}")));
		System.out.println(String.format("c_float A_x[%d] = %s;",Ax.length,Arrays.toString(Ax).replace("[", "{").replace("]", "}")));
		System.out.println(String.format("c_int P_p[%d] = %s;",Pp.length,Arrays.toString(Pp).replace("[", "{").replace("]", "}")));
		System.out.println(String.format("c_int P_i[%d] = %s;",Pi.length,Arrays.toString(Pi).replace("[", "{").replace("]", "}")));
		System.out.println(String.format("c_float P_x[%d] = %s;",Px.length,Arrays.toString(Px).replace("[", "{").replace("]", "}")));
		System.out.println(String.format("c_float l[%d] = %s;",l.length,Arrays.toString(l).replace("[", "{").replace("]", "}").replace("Infinity","OSQP_INFTY")));
		System.out.println(String.format("c_float u[%d] = %s;",u.length,Arrays.toString(u).replace("[", "{").replace("]", "}").replace("Infinity","OSQP_INFTY")));
		System.out.println(String.format("c_float q[%d] = %s;",q.length,Arrays.toString(q).replace("[", "{").replace("]", "}")));
		System.out.println(String.format("c_int m=%d;",nRows));
		System.out.println(String.format("c_int n=%d;",nCols));
		System.out.println(String.format("c_int A_nnz=%d;",Anz));
		System.out.println(String.format("c_int P_nnz=%d;",Pnz));
	}

	public CSCMatrix getA( ) {
		return new CSCMatrix(nRows, nCols, Anz, Ap, Ai, Ax);
	}

	public CSCMatrix getP() {
		return new CSCMatrix(nCols, nCols, Pnz, Pp, Pi, Px);
	}

	public double[] getl() {
		return l;
	}

	public double[] getu() {
		return u;
	}

	public double[] getq() {
		return q;
	}

	public OSQP.Data getData() {
		return new OSQP.Data(nCols,nRows,getP(),getA(),q,l,u);
	}

	private enum Section {
		ROWS, COLUMNS, RHS, BOUNDS, HEAD, SENSE, OBJ, RANGES;
	}

	private static void parseRow(int[] shape, Map<String, Integer> rows, List<Double> l, List<Double> u, String[] tokens) {
		String rowName = tokens[2];
		switch(tokens[1]) {
			case "L": u.add(0.); l.add(-OSQP.OSQP_INFTY); rows.put(rowName, shape[0]++); break;
			case "G": u.add(OSQP.OSQP_INFTY); l.add(0.); rows.put(rowName, shape[0]++); break;
			case "E": u.add(0.); l.add(0.); rows.put(rowName, shape[0]++); break;
			default: break;
		}
	}

	private static void parseBnd(Map<String, Integer> cols, List<Integer> Ai, List<Integer> Ap, List<Double> l, List<Double> u, String[] tokens) {
		String col = tokens[3];
		int rowIndex = Ai.get(Ap.get(cols.get(col)));
		double bnd;
		switch(tokens[1]) {
			case "LO":
				bnd = Double.parseDouble(tokens[4]);
				l.set(rowIndex, bnd);
				break;
			case "UP":
				bnd = Double.parseDouble(tokens[4]);
				u.set(rowIndex, bnd);
				break;
			case "FX":
				bnd = Double.parseDouble(tokens[4]);
				l.set(rowIndex, bnd); u.set(rowIndex, bnd);
				break;
			//TODO: Complete the following bound keys.
			case "MI":
				break;
			case "PI":
				break;
			case "FR":
				l.set(rowIndex, -OSQP.OSQP_INFTY); u.set(rowIndex, OSQP.OSQP_INFTY);
				break;
			default: throw new IllegalStateException("Unkown bound key: "+tokens[1]);
		}
	}

	private static void parseCol(int[] shape, Map<String, Integer> rows, Map<String, Integer> cols,
								 List<Integer> Ai, List<Integer> Ap, List<Double> Ax,
								 List<Double> l, List<Double> u, List<Double> q,
								 String objName, String[] tokens, double sign) {
		String colName = tokens[1];
		if (!cols.containsKey(colName)) {
			q.add(0.);
			Ap.set(shape[1],Ai.size());
			Ai.add(shape[0]++);
			Ax.add(1.);
			Ap.add(Ai.size());
			l.add(0.);
			u.add(OSQP.OSQP_INFTY);
			cols.put(colName, shape[1]++);
		}
		int colIndex = cols.get(colName); //TODO: this won't work if columns are not ordered
		for (int i=2; i<tokens.length; i+=2) {
			String rowName = tokens[i];
			if (!rowName.matches(objName)) {
				Ai.add(rows.get(rowName));
				Ax.add(Double.parseDouble(tokens[i+1]));
				Ap.set(colIndex+1, Ai.size());
			} else {
				q.set(colIndex, sign*Double.parseDouble(tokens[i+1]));
			}
		}
	}

	private static void parseRhs(Map<String, Integer> rows, List<Double> l, List<Double> u, String[] tokens) {
		for (int i=2; i<tokens.length; i+=2) {
			int row = rows.get(tokens[i]);
			if (u.get(row)==0.)
				u.set(row, Double.parseDouble(tokens[i+1]));
			if (l.get(row)==0.)
				l.set(row, Double.parseDouble(tokens[i+1]));
		}
	}


	private static void parseRanges(Map<String, Integer> rows, List<Double> l, List<Double>  u, String[] tokens) {
		if (!rows.containsKey(tokens[2]))
			throw new IllegalStateException(String.format("Error in reading ranges."));
		Integer rowNum = rows.get(tokens[2]);
		Double range = Double.parseDouble(tokens[3]);
		if (u.get(rowNum) == OSQP.OSQP_INFTY)
			u.set(rowNum, l.get(rowNum) + range);
		else if (l.get(rowNum) == -OSQP.OSQP_INFTY)
			l.set(rowNum, u.get(rowNum) - range);
		else
			throw new IllegalStateException(String.format("Error in reading ranges."));
	}


	public static Parser readMps(String filename) {
		OSQP.LOG.info(String.format("Begin parsing file %s.",filename));
		double tme = System.currentTimeMillis();
		String objname = null;
		boolean maximize = false;
		int[] shape = new int[2];
		Map<String, Integer> rows = new HashMap<>();
		Map<String, Integer> cols = new HashMap<>();
		List<Integer> Ai = new ArrayList<>();
		List<Integer> Ap = new ArrayList<>();
		Ap.add(0);
		List<Double> Ax = new ArrayList<>();
		List<Double> u = new ArrayList<>();
		List<Double> l = new ArrayList<>();
		List<Double> q = new ArrayList<>();
		try {
			FileInputStream in = new FileInputStream(filename);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			Section section = Section.HEAD;
			String line;
			while((line = br.readLine()) != null) {
				String[] tokens = line.split("\\s+");
				if (tokens.length > 1 || !tokens[0].isEmpty())
					try {
						if (tokens[0].isEmpty()) {
							if (tokens[1].charAt(0) == '*') // skipping commented lines
								continue;
							switch(section) {
								case ROWS:
									if (objname==null && tokens[1].matches("N"))
										objname = tokens[2];
									else
										parseRow(shape, rows, l, u, tokens);
									break;
								case COLUMNS:
									parseCol(shape, rows, cols, Ai, Ap, Ax, l, u, q, objname, tokens, maximize?-1.:1.);
									break;
								case RHS:
									parseRhs(rows, l, u, tokens);
									break;
								case RANGES:
									parseRanges(rows, l, u, tokens);
									break;
								case BOUNDS:
									parseBnd(cols, Ai, Ap, l, u, tokens);
									break;
								case OBJ:
									objname = tokens[1];
									break;
								case SENSE:
									maximize = tokens[1].matches("MAX") || tokens[1].matches("MAXIMIZE");
									break;
								default:
									throw new IllegalStateException(
											String.format("Line %s started with an empty character but contains no data.",line)
									);
							}
						} else {
							// skipping the commented lines
							if (tokens[0].charAt(0) == '*')
								continue;
							// At this stage, it is for sure that the line is a header line.
							switch(tokens[0]) {
								case "NAME":     section = Section.HEAD;    break;
								case "OBJSENSE": section = Section.SENSE;   break;
								case "OBJNAME":  section = Section.OBJ;     break;
								case "ROWS":     section = Section.ROWS;    break;
								case "COLUMNS":  section = Section.COLUMNS; break;
								case "RHS":      section = Section.RHS;     break;
								case "RANGES":   section = Section.RANGES;  break;
								case "BOUNDS":   section = Section.BOUNDS;  break;
								case "ENDATA":                              break; //end of file
								default: throw new IllegalStateException(
										String.format("Unknown section name: %s",tokens[0])
								);
							}
						}
					} catch (ArrayIndexOutOfBoundsException e) {
						e.printStackTrace();
						throw new IllegalStateException(String.format("Error in line: '%s'", line));
					}
			}
			br.close();
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		Utils.toDoubleArray(q);
		double[] P_x = new double[shape[1]];
		int[] P_p = new int[shape[1]+1];
		int[] P_i = new int[shape[1]];
		for (int i=0; i<shape[1]; i++) {
			P_i[i] = i;
			P_p[i+1] = i+1;
//			P_x[i] = 1.e-4;
		}
		OSQP.LOG.info(String.format("Read MPS in %.2fsec.",(System.currentTimeMillis()-tme)/1000.));
		return new Parser(Utils.toDoubleArray(q),
				P_x, P_p, P_i,
				Utils.toDoubleArray(Ax), Utils.toIntArray(Ap), Utils.toIntArray(Ai),
				Utils.toDoubleArray(l), Utils.toDoubleArray(u));
	}

	public static void main(String... args) {

		try {
			FileHandler fh = new FileHandler("josqp.log");
			OSQP.LOG.addHandler(fh);
			Formatter formatter = new SimpleFormatter();
			fh.setFormatter(formatter);
		} catch (IOException e) {
			e.printStackTrace();
		}

		//Parser p = Parser.readMps("src/test/resources/neos-3025225_lp.mps"); //Optimal objective  5.257405203e-02
			//Parser p = Parser.readMps("src/test/resources/s82_lp.mps"); //Optimal objective -3.457971082e+01
			//Parser p = Parser.readMps("src/test/resources/qap15.mps"); //Optimal objective  1.040994041e+03
			//Parser p = Parser.readMps("src/test/resources/irish-electricity.mps"); //Optimal objective  2.546254563e+06
//		Parser p = Parser.readMps("src/test/resources/supportcase10.mps"); //Optimal objective  3.383923666e+00
			//Parser p = Parser.readMps("src/test/resources/ex10.mps"); //Optimal objective 100
//		Parser p = Parser.readMps("src/test/resources/savsched1.mps"); //Optimal objective 2.1740357143e+02
			//Parser p = Parser.readMps("src/test/resources/sample1.mps");
			String mpsFileDir = "src/test/resources/sample1.mps";
			if (args.length >= 1)
				mpsFileDir = args[0];
			else {
				OSQP.LOG.info("Usage: java -jar josqp.jar <mps_file_name>");
				return;
			}
			//String mpsFileDir = "src/test/resources/sample1.mps";
			Parser p = Parser.readMps(mpsFileDir);

			OSQP.Data data = p.getData();
			OSQP.Settings settings = new OSQP.Settings();
			//settings.eps_rel = 1.e-6;
			//settings.alpha = 1.66666667;
			//settings.sigma = 1.e-4;
			settings.polish = true;
			//settings.adaptive_rho = false;
			//settings.rho = 0.000001;
			//settings.eps_rel = 1.e-6;
			settings.verbose = true;
			OSQP opt = new OSQP(data,settings);
			opt.solve();

		}

}
