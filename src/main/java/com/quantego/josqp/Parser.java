package com.quantego.josqp;


import java.io.*;
import java.util.*;
import java.util.logging.FileHandler;
import java.util.logging.Formatter;
import java.util.logging.LogRecord;
import java.util.logging.SimpleFormatter;
import java.util.regex.Pattern;

class RcvMat {
	public int cap, nnz;
	public int[] rows, cols;
	public double[] vals;

	public RcvMat(int cap) {
		this.cap = cap;
		nnz = 0;
		rows = new int[cap];
		cols = new int[cap];
		vals = new double[cap];
	}

	public void addEntry(int r, int c, double v) {
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
		ROWS, COLUMNS, RHS, BOUNDS, HEAD, SENSE, OBJ, RANGES, QUADOBJ, ENDDATA;
	}

	private static void parseRow(int[] shape, Map<String, Integer> rows, List<Double> l, List<Double> u,
								 List<Integer> rowSense, String[] tokens) {
		String rowName = tokens[2].replace(".", "\\.");
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
		String colName = tokens[3].replace(".", "\\.");
		int rowIndex = rows.get(colName + "_bnd");
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
				l.set(rowIndex, bnd);
				u.set(rowIndex, bnd);
				break;
			case "MI": // inf lower bound
				l.set(rowIndex, -OSQP.OSQP_INFTY);
				break;
			case "PI": // inf upper bound
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
		String colName = tokens[1].replace(".", "\\.");
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
			String rowName = tokens[i].replace(".", "\\.");
			if (!Pattern.matches(rowName, objName)) { // initially there was rowName.matches method, but that was showing a strange behavior in some examples (e.g. qbandm.qps)! It could be a bug in the Java language!
				Arcv.addEntry(rows.get(rowName), colIndex, Double.parseDouble(tokens[i+1]));
			} else {
				q.set(colIndex, sign*Double.parseDouble(tokens[i+1]));
			}
		}
	}

	private static void parseRhs(Map<String, Integer> rows, List<Integer> rowSense, List<Double> l, List<Double> u, String[] tokens) {
		int rowN;
		String rowName;
		for (int i=2; i<tokens.length; i+=2) {
			rowName = tokens[i].replace(".", "\\.");
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

	private static void parseQuadObj(Map<String, Integer> cols, RcvMat Prcv, String[] tokens) {
		// First constructing a matrix in RCV format and then translating it to the CSC format.
		ArrayList<Integer> rows, columns;
		ArrayList<Double> values;
		int col1, col2;
		double val;
		try {
			col1 = cols.get(tokens[1]);
			col2 = cols.get(tokens[2]);
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

	private static void extractUpperTriangle(ArrayList<Integer> rcvR, ArrayList<Integer> rcvC, ArrayList<Double> rcvV) {
		for (int i = 0; i < rcvR.size(); i++) {
			if (rcvR.get(i) > rcvC.get(i)) {
				rcvR.remove(i);
				rcvC.remove(i);
				rcvV.remove(i);
			}
		}
	}

	private static void convertRcvToCsc(ArrayList<Integer> rcvR, ArrayList<Integer> rcvC, ArrayList<Double> rcvV, int NCol,
										List<Integer> cscI, List<Integer> cscP, List<Double> cscX) {
		// convert RCV format to CSC format.

		// sorting the RCV format (first based on the column, then based on the row).
		int r, c, col;
		double v;
		// step 1: sorting based on the column
		for (int i = 0; i < rcvC.size(); i++) {
			for (int j = 0; j < i; j++) {
				if (rcvC.get(i) < rcvC.get(j)) {
					Collections.swap(rcvR, i, j);
					Collections.swap(rcvC, i, j);
					Collections.swap(rcvV, i, j);
					//r = rcvR.get(i);   c = rcvC.get(i);   v = rcvV.get(i);
					//rcvR.remove(i);    rcvC.remove(i);    rcvV.remove(i);
					//rcvR.add(j, r);    rcvC.add(j, c);    rcvV.add(j, v);
					break;
				}
			}
		}
		// step 2: sorting each column based on the rows.
		int idxStart = 0;
		int idxEnd = 0;
		for (int colN = 0; colN < Collections.max(rcvC); colN++) {
			if (rcvC.get(idxStart) > colN)
				continue;
			while (rcvC.get(idxEnd + 1) == colN)
				idxEnd++;
			for (int i = idxStart; i <= idxEnd; i++) {
				for (int j = idxStart; j < i; j++) {
					if (rcvR.get(i) < rcvR.get(j)) {
						Collections.swap(rcvR, i, j);
						Collections.swap(rcvC, i, j);
						Collections.swap(rcvV, i, j);
						//r = rcvR.get(i);    c = rcvC.get(i);    v = rcvV.get(i);
						//rcvR.remove(i);     rcvC.remove(i);     v = rcvV.remove(i);
						//rcvR.add(j, r);     rcvC.add(j, c);     rcvV.add(j, v);
						break;
					}
				}
			}
			idxStart = ++idxEnd;
		}

		// constructing the CSC vectors based on the RCV vectors.
		for (int i = 0; i < rcvV.size(); i++) {
			cscX.add(rcvV.get(i));
			cscI.add(rcvR.get(i));
		}
		cscP.add(0);
		int totNCol = 0, currColParsed = 0;
		for (c = 0; c < rcvC.size(); c++) {
			col = rcvC.get(c);
			while (col > currColParsed) { // This while loop is to handle the empty columns at the beginning.
				cscP.add(totNCol);
				currColParsed++;
			}
			totNCol++;
			if (c != rcvC.size() - 1) {
				if (rcvC.get(c + 1) != col) {
					cscP.add(totNCol);
					currColParsed++;
				}
			} else {
				cscP.add(totNCol);
				while (currColParsed < NCol - 1) { // This while loop is to handle the empty columns at the end.
					cscP.add(totNCol);
					currColParsed++;
				}
			}
		}
	}


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
		try {
			FileInputStream in = new FileInputStream(filename);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			Section currentSection = Section.HEAD;
			Section prevSection = currentSection;
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
										objname = tokens[2];
									else
										parseRow(shape, rows, l, u, rowSense, tokens);
									break;
								case COLUMNS:
									parseCol(shape, rows, cols, Arcv, l, u, q, objname, tokens, maximize?-1.:1.);
									break;
								case RHS:
									parseRhs(rows, rowSense, l, u, tokens);
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
							if (prevSection != currentSection) { // postprocessing phase of parsing each section
								switch (prevSection) {
									case COLUMNS:
										//convertRcvToCsc(ArcvR, ArcvC, ArcvV, cols.size(), Ai, Ap, Ax);
										break;
									case QUADOBJ:
										//extractUpperTriangle(PrcvR, PrcvC, PrcvV);
										//convertRcvToCsc(PrcvR, PrcvC, PrcvV, cols.size(), Pi, Pp, Px);
										break;
									default:
										break;
								}
								prevSection = currentSection;
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
		return new Parser(Utils.toDoubleArray(q),
				Pcsc.Ax, Pcsc.Ap, Pcsc.Ai,
				Acsc.Ax, Acsc.Ap, Acsc.Ai,
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
        //Parser p = Parser.readMps("src/test/resources/supportcase10.mps"); //Optimal objective  3.383923666e+00
		//Parser p = Parser.readMps("src/test/resources/ex10.mps"); //Optimal objective 100
		//Parser p = Parser.readMps("src/test/resources/savsched1.mps"); //Optimal objective 2.1740357143e+02
		//Parser p = Parser.readMps("src/test/resources/sample1.mps");
		//	String mpsFileDir = "src/test/resources/sample1.mps";
		//String qpsFileDir = "src/test/resources/sample1.qps";
		//String qpsFileDir = "src/test/resources/qafiro.qps";
		//String qpsFileDir = "src/test/resources/qbandm.qps";
		//String qpsFileDir = "src/test/resources/qsctap1.qps";
		//String qpsFileDir = "src/test/resources/cvxqp1_l.qps";
		//String qpsFileDir = "src/test/resources/boyd1.qps";
		String qpsFileDir = "src/test/resources/boyd2.qps";
		//String qpsFileDir;
		if (args.length >= 1)
			qpsFileDir = args[0];
		else {
			OSQP.LOG.info("Usage: java -jar josqp.jar <qps_file_name>");
			//return;
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
