package com.quantego.josqp;


import java.io.*;
import java.util.*;
import java.util.logging.FileHandler;
import java.util.logging.Formatter;
import java.util.logging.LogRecord;
import java.util.logging.SimpleFormatter;
import java.util.regex.Pattern;

/**
 * A class to store a matrix in the Row-Column-Value (RCV) format.
 */
class RcvMat {
	protected int cap, nnz;
	protected int[] rows, cols;
	protected double[] vals;

	/**
	 * Constructor of the class RcvMat.
	 *
	 * @param  cap  initial capacity of the class
	 */
	protected RcvMat(int cap) {
		this.cap = cap;
		nnz = 0;
		rows = new int[cap];
		cols = new int[cap];
		vals = new double[cap];
	}

	/**
	 * Adds a new entry to its corresponding RCV-formatted matrix
	 *
	 * @param  r  row number of the new element
	 * @param  c  column number of the new element
	 * @param  v  value of the new element
	 */
	protected void addEntry(int r, int c, double v) {
		// Check capacity and add space if needed
		if (nnz >= cap) {
			rows = Arrays.copyOf(rows, 2 * cap);
			cols = Arrays.copyOf(cols, 2 * cap);
			vals = Arrays.copyOf(vals, 2 * cap);
			cap *= 2;
		}
		// add the new entry
		rows[nnz] = r;
		cols[nnz] = c;
		vals[nnz] = v;
		nnz++;
	}
}

/**
 * This class is used to parse MPS- and QPS-formatted files and
 * solve the corresponding QP problems with jOSQP solver
 *
 * @author Farid Alavi
 */
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
	public double offset;

	/**
	 * Create a new Problem instance based on vector inputs (CSC formatted).
	 * @param q       linear objective coefficients
	 * @param P_x     quadratic objective (coefficients)
	 * @param P_p     quadratic objective (next column starts)
	 * @param P_i     quadratic objective (row indices)
	 * @param A_x     constraints (coefficients)
	 * @param A_p     constraints (next column starts)
	 * @param A_i     constraints (row indices)
	 * @param l       constraint lower bounds
	 * @param u       constraint upper bounds
	 * @param offset  offset of the objective function
	 */
	public Parser(double[] q, double[] P_x, int[] P_p, int[] P_i, double[] A_x, int[] A_p, int[] A_i, double[] l, double[] u,
								double offset) {
		this.q = q;
		this.Px = P_x;
		this.Pp = P_p;
		this.Pi = P_i;
		this.Ax = A_x;
		this.Ap = A_p;
		this.Ai = A_i;
		this.l = l;
		this.u = u;
		this.offset = offset;
		nCols = q.length;
		nRows = l.length;
		Anz = Ax.length;
		Pnz = Px.length;
	}
	/**
	 * Create a new Problem instance based on vector inputs (CSC formatted).
	 * @param q       linear objective coefficients
	 * @param P_x     quadratic objective (coefficients)
	 * @param P_p     quadratic objective (next column starts)
	 * @param P_i     quadratic objective (row indices)
	 * @param A_x     constraints (coefficients)
	 * @param A_p     constraints (next column starts)
	 * @param A_i     constraints (row indices)
	 * @param l       constraint lower bounds
	 * @param u       constraint upper bounds
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
		this.offset = offset;
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
	private static Parser readVectors(String prefix, int skipRows, int skipCols) {
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

	/**
	 * Get the constraint matrix (A)
	 *
	 * @return    The constraint matrix (A) in the CSC format
	 */
	public CSCMatrix getA( ) {
		return new CSCMatrix(nRows, nCols, Anz, Ap, Ai, Ax);
	}

	/**
	 * Get the quadratic objective matrix (P)
	 *
	 * @return    The quadratic objective matrix (P) in the CSC format
	 */
	public CSCMatrix getP() {
		return new CSCMatrix(nCols, nCols, Pnz, Pp, Pi, Px);
	}

	/**
	 * Get the lower-bound vector (l)
	 */
	public double[] getl() {
		return l;
	}

	/**
	 * Get the upper-bound vector (u)
	 */
	public double[] getu() {
		return u;
	}

	/**
	 * Get the linear objective coefficients vector (q)
	 */
	public double[] getq() {
		return q;
	}

	/**
	 * Get a copy of the problem data.
	 */
	public OSQP.Data getData() {
		return new OSQP.Data(nCols, nRows, getP(), getA(), q, l, u, offset);
	}

	private enum Section {
		ROWS, COLUMNS, RHS, BOUNDS, HEAD, SENSE, OBJ, RANGES, QUADOBJ, ENDDATA;
	}

	private static void parseRow(int[] shape, Map<String, Integer> rows, List<Double> l, List<Double> u,
								 List<Integer> rowSense, String[] tokens) {
		String rowName = removeSpecialChars(tokens[2]);
		int rowN;
		switch(tokens[1]) {
			case "L":
				rowN = shape[0]++;
				rows.put(rowName, rowN);
				rowSense.add(rowN, -1);
				u.add(rowN, 0.0);
				l.add(rowN, -OSQP.OSQP_INFTY);
				break;
			case "G":
				rowN = shape[0]++;
				rows.put(rowName, rowN);
				rowSense.add(rowN, 1);
				u.add(rowN, OSQP.OSQP_INFTY);
				l.add(rowN, 0.0);
				break;
			case "E":
				rowN = shape[0]++;
				rows.put(rowName, rowN);
				rowSense.add(rowN, 0);
				u.add(rowN, 0.0);
				l.add(rowN, 0.0);
				break;
			default:
				break;
		}
	}

	private static void parseBnd(Map<String, Integer> cols, Map<String, Integer> rows, List<Double> l, List<Double> u, String[] tokens) {
		int placeColName = 3;
		if (tokens.length % 2 == 0) {
			placeColName = 2;
		}
		String colName = removeSpecialChars(tokens[placeColName]);
		int rowIndex = rows.get(colName + "_bnd");
		double bnd;
		switch(tokens[1]) {
			case "LO":
				bnd = Double.parseDouble(tokens[placeColName + 1]);
				l.set(rowIndex, bnd);
				break;
			case "UP":
				bnd = Double.parseDouble(tokens[placeColName + 1]);
				u.set(rowIndex, bnd);
				break;
			case "FX":
				bnd = Double.parseDouble(tokens[placeColName + 1]);
				l.set(rowIndex, bnd);
				u.set(rowIndex, bnd);
				break;
			case "MI": // inf lower bound
				l.set(rowIndex, -OSQP.OSQP_INFTY);
				break;
			case "PL": // inf upper bound
				u.set(rowIndex, OSQP.OSQP_INFTY);
				break;
			case "FR": // free variable, no upper or lower bound
				l.set(rowIndex, -OSQP.OSQP_INFTY);
				u.set(rowIndex,  OSQP.OSQP_INFTY);
				break;
			default:
				throw new IllegalStateException("Unkown bound key: " + tokens[1]);
		}
	}

	private static void parseCol(int[] shape, Map<String, Integer> rows, Map<String, Integer> cols,
								 RcvMat Arcv, List<Double> l, List<Double> u, List<Double> q,
								 String objName, String[] tokens, double sign) {
		String colName = removeSpecialChars(tokens[1]);
		int colN, rowN;
		if (!cols.containsKey(colName)) {
			colN = shape[1]++;
			rowN = shape[0]++;
			String rowName = colName + "_bnd";
			cols.put(colName, colN);
			q.add(0.0);
			rows.put(rowName, rowN);
			Arcv.addEntry(rowN, colN, 1.0);
			l.add(rowN, 0.0);
			u.add(rowN, OSQP.OSQP_INFTY);
		}
		int colIndex = cols.get(colName);
		for (int i=2; i<tokens.length; i+=2) {
			String rowName = removeSpecialChars(tokens[i]);
			if (!Pattern.matches(rowName, objName)) { // initially there was rowName.matches method, but that was showing a strange behavior in some examples (e.g. qbandm.qps)! It could be a bug in the Java language!
				Arcv.addEntry(rows.get(rowName), colIndex, Double.parseDouble(tokens[i+1]));
			} else {
				q.set(colIndex, sign*Double.parseDouble(tokens[i+1]));
			}
		}
	}

	private static String removeSpecialChars(String input) {
		String output = input.replace(".", "!_");
		output = output.replace("+", "!!");
		return output;
	}

	private static void parseRhs(Map<String, Integer> rows, List<Integer> rowSense, List<Double> l, List<Double> u,
															 String objName, List<Double> offset, String[] tokens) {
		int rowN;
		String rowName;
		for (int i=2; i<tokens.length; i+=2) {
			rowName = removeSpecialChars(tokens[i]);
			if (Pattern.matches(rowName, objName)) {
				offset.add(0, -Double.parseDouble(tokens[i + 1]));
				continue;
			}
			rowN = rows.get(rowName);
			assert u.size() == l.size();
			assert rowN <= u.size() : "Error in parsing RHS!";
			switch (rowSense.get(rowN)) {
				case -1:
					u.set(rowN, Double.parseDouble(tokens[i+1]));
					break;
				case 1:
					l.set(rowN, Double.parseDouble(tokens[i+1]));
					break;
				case 0:
					u.set(rowN, Double.parseDouble(tokens[i+1]));
					l.set(rowN, Double.parseDouble(tokens[i+1]));
					break;
				default:
					throw new IllegalStateException(String.format("Error in reading RHS."));
			}
		}
	}


	private static void parseRanges(Map<String, Integer> rows, List<Double> l, List<Double>  u, String[] tokens) {
		if (!rows.containsKey(removeSpecialChars(tokens[2])))
			throw new IllegalStateException(String.format("Error in reading ranges."));
		Integer rowNum = rows.get(removeSpecialChars(tokens[2]));
		Double range = Double.parseDouble(tokens[3]);
		if (u.get(rowNum) == OSQP.OSQP_INFTY)
			u.set(rowNum, l.get(rowNum) + range);
		else if (l.get(rowNum) == -OSQP.OSQP_INFTY)
			l.set(rowNum, u.get(rowNum) - range);
		else
			throw new IllegalStateException(String.format("Error in reading ranges."));
	}

	private static void parseQuadObj(Map<String, Integer> cols, RcvMat Prcv, String[] tokens) {
		// First constructing a matrix in RCV format and then translating it to the CSC format.
		ArrayList<Integer> rows, columns;
		ArrayList<Double> values;
		int col1, col2;
		double val;
		try {
			col1 = cols.get(removeSpecialChars(tokens[1]));
			col2 = cols.get(removeSpecialChars(tokens[2]));
			val  = Double.parseDouble(tokens[3]);
			if (col1 < col2) {
				Prcv.addEntry(col1, col2, val);
			} else if (col1 > col2) {
				Prcv.addEntry(col2, col1, val);
			} else {
				Prcv.addEntry(col1, col2, val);
			}
		} catch (Exception e) {
		System.out.println("Error in reading the QUADOBJ section!");
		}
	}


	/**
	 * This function reads the content of a QPS file and constructs
	 * the problem that is described therein.
	 *
	 * @param  filename  The absolute path to a QPS-formatted file
	 * @return      An instance of the Parser class representing the
	 *              problem described in the input file
	 */
	public static Parser readQmps(String filename) {
		OSQP.LOG.info(String.format("Begin parsing file %s.",filename));
		double tme = System.currentTimeMillis();
		String objname = null;
		boolean maximize = false;
		int[] shape = new int[2];
		Map<String, Integer> rows = new HashMap<>();
		Map<String, Integer> cols = new HashMap<>();
		RcvMat Arcv = new RcvMat(1024);
		CSCMatrix Acsc;
		RcvMat Prcv = new RcvMat(1024);
		CSCMatrix Pcsc;
		List<Double>  u          = new ArrayList<>();
		List<Double>  l          = new ArrayList<>();
		List<Double>  q          = new ArrayList<>();
		List<Integer> rowSense   = new ArrayList<>();
		List<Double>  offset     = new ArrayList<>();
		try {
			FileInputStream in = new FileInputStream(filename);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			Section currentSection = Section.HEAD;
			String line;
			while((line = br.readLine()) != null) {
				String[] tokens = line.split("\\s+");
				if (tokens.length > 1 || !tokens[0].isEmpty())
					try {
						if (tokens[0].isEmpty()) {
							if (tokens[1].charAt(0) == '*') // skipping commented lines
								continue;
							switch(currentSection) {
								case ROWS:
									if (objname==null && tokens[1].matches("N"))
										objname = removeSpecialChars(tokens[2]);
									else
										parseRow(shape, rows, l, u, rowSense, tokens);
									break;
								case COLUMNS:
									parseCol(shape, rows, cols, Arcv, l, u, q, objname, tokens, maximize?-1.:1.);
									break;
								case RHS:
									parseRhs(rows, rowSense, l, u, objname, offset, tokens);
									break;
								case RANGES:
									parseRanges(rows, l, u, tokens);
									break;
								case BOUNDS:
									parseBnd(cols, rows, l, u, tokens);
									break;
								case OBJ:
									objname = tokens[1];
									break;
								case SENSE:
									maximize = tokens[1].matches("MAX") || tokens[1].matches("MAXIMIZE");
									break;
								case QUADOBJ:
									parseQuadObj(cols, Prcv, tokens);
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
								case "NAME":
									currentSection = Section.HEAD;
									break;
								case "OBJSENSE":
									currentSection = Section.SENSE;
									break;
								case "OBJNAME":
									currentSection = Section.OBJ;
									break;
								case "ROWS":
									currentSection = Section.ROWS;
									break;
								case "COLUMNS":
									currentSection = Section.COLUMNS;
									break;
								case "RHS":
									currentSection = Section.RHS;
									break;
								case "RANGES":
									currentSection = Section.RANGES;
									break;
								case "BOUNDS":
									currentSection = Section.BOUNDS;
									break;
								case "QUADOBJ":
									currentSection = Section.QUADOBJ;
									break;
								case "ENDATA":
									currentSection = Section.ENDDATA;
									break; //end of file
								default:
									throw new IllegalStateException(
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
		OSQP.LOG.info(String.format("Read MPS in %.2fsec.",(System.currentTimeMillis()-tme)/1000.));
		Acsc = CSCMatrix.triplet_to_csc(rows.size(), cols.size(), Arcv.nnz, Arcv.rows, Arcv.cols, Arcv.vals, null);
		Pcsc = CSCMatrix.triplet_to_csc(cols.size(), cols.size(), Prcv.nnz, Prcv.rows, Prcv.cols, Prcv.vals, null);
		if (offset.size() == 0)
			offset.add(0.0);
		return new Parser(Utils.toDoubleArray(q),
				Pcsc.Ax, Pcsc.Ap, Pcsc.Ai,
				Acsc.Ax, Acsc.Ap, Acsc.Ai,
				Utils.toDoubleArray(l), Utils.toDoubleArray(u), offset.get(0));
	}

	/**
	 * This serves as the program's entry point. An address for a QPS
	 * file should be provided as an argument. This function extracts
	 * the problem described in the given file and solves it using
	 * the specific algorithm of the JOSQP solver.
	 *
	 * @param  args  A String representing the address of the qps file
	 *               For example: "src/test/resources/sample1.qps"
	 */
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
		//Parser p = Parser.readMps("src/test/resources/supportcase10.mps"); //Optimal objective  3.383923666e+00
		//Parser p = Parser.readMps("src/test/resources/ex10.mps"); //Optimal objective 100
		//Parser p = Parser.readMps("src/test/resources/savsched1.mps"); //Optimal objective 2.1740357143e+02
		//Parser p = Parser.readMps("src/test/resources/sample1.mps");
		//	String mpsFileDir = "src/test/resources/sample1.mps";
		//String qpsFileDir = "src/test/resources/sample1.qps";
		String qpsFileDir;
		if (args.length >= 1)
			qpsFileDir = args[0];
		else {
			OSQP.LOG.info("Usage: java -jar josqp.jar <qps_file_name>");
			return;
		}
		//Parser p = Parser.readMps(mpsFileDir);
		Parser p = Parser.readQmps(qpsFileDir);

		OSQP.Data data = p.getData();
		OSQP.Settings settings = new OSQP.Settings();
		settings.max_iter = 100000;
		//settings.eps_rel = 1.e-6;
		//settings.alpha = 1.66666667;
		//settings.sigma = 1.e-4;
		//settings.polish = true;
		//settings.adaptive_rho = false;
		//settings.rho = 0.000001;
		//settings.eps_rel = 1.e-6;
		settings.verbose = true;
		OSQP opt = new OSQP(data,settings);
		opt.solve();

	}

}
