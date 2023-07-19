package com.quantego.josqp;


import java.io.*;
import java.util.*;
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
		ROWS, COLUMNS, RHS, BOUNDS, HEAD, SENSE, OBJ, RANGES, QUADOBJ;
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

	private static void parseQuadObj(Map<String, Integer> cols, List<Integer> PrcvR, List<Integer> PrcvC, List<Double> PrcvV, String[] tokens) {

		// First constructing a matrix in RCV format and then translating it to the CSC format.

		ArrayList<Integer> rows, columns;
		ArrayList<Double> values;
		int col1, col2;
		double val;
		try {
			col1 = cols.get(tokens[1]);
			col2 = cols.get(tokens[2]);
			val  = Double.parseDouble(tokens[3]);
			if (col1 != col2) {
				PrcvR.add(col1);
				PrcvC.add(col2);
				PrcvV.add(val);
				PrcvR.add(col2);
				PrcvC.add(col1);
				PrcvV.add(val);
			} else {
				PrcvR.add(col1);
				PrcvC.add(col2);
				PrcvV.add(2 * val);
			}
		}
		catch (Exception e) {
			System.out.println("Error in reading the QUADOBJ section!");
		}
	}

	private static void extractUpperTriangle(List<Integer> rcvR, List<Integer> rcvC, List<Double> rcvV) {
		for (int i = 0; i < rcvR.size(); i++) {
			if (rcvR.get(i) > rcvC.get(i)) {
				rcvR.remove(i);
				rcvC.remove(i);
				rcvV.remove(i);
			}
		}
	}

	private static void convertRcvToCsc(List<Integer> rcvR, List<Integer> rcvC, List<Double> rcvV, int NCol,
										List<Integer> cscI, List<Integer> cscP, List<Double> cscX) {
		// convert RCV format to CSC format.

		// sorting the RCV format (first based on the column, then based on the row).
		int r, c, col;
		double v;
			// step 1: sorting based on the column
		for (int i = 0; i < rcvC.size(); i++) {
			for (int j = 0; j < i; j++) {
				if (rcvC.get(i) < rcvC.get(j)) {
					r = rcvR.get(i);   c = rcvC.get(i);   v = rcvV.get(i);
					rcvR.remove(i);    rcvC.remove(i);    rcvV.remove(i);
					rcvR.add(j, r);    rcvC.add(j, c);    rcvV.add(j, v);
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
						r = rcvR.get(i);    c = rcvC.get(i);    v = rcvV.get(i);
						rcvR.remove(i);     rcvC.remove(i);     v = rcvV.remove(i);
						rcvR.add(j, r);     rcvC.add(j, c);     rcvV.add(j, v);
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
		for (int i = 0; i < cscP.size(); i++) { // Because the indexing in OSQP starts from 1.
			cscP.set(i, cscP.get(i) + 1);
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
		List<Integer> Ai = new ArrayList<>();
		List<Integer> Ap = new ArrayList<>();
		Ap.add(0);
		List<Double> Ax = new ArrayList<>();
		List<Integer> Pi = new ArrayList<>(); // row indices
		List<Integer> Pp = new ArrayList<>(); // column pointers
		List<Double>  Px = new ArrayList<>();
		List<Integer> PrcvR = new ArrayList<>();
		List<Integer> PrcvC = new ArrayList<>();
		List<Double> PrcvV = new ArrayList<>();
		List<Double> u = new ArrayList<>();
		List<Double> l = new ArrayList<>();
		List<Double> q = new ArrayList<>();
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
										parseRow(shape, rows, l, u, tokens);
									break; case COLUMNS: parseCol(shape, rows, cols, Ai, Ap, Ax, l, u, q, objname, tokens, maximize?-1.:1.);
									break; case RHS: parseRhs(rows, l, u, tokens);
									break; case RANGES: parseRanges(rows, l, u, tokens);
									break; case BOUNDS: parseBnd(cols, Ai, Ap, l, u, tokens);
									break; case OBJ: objname = tokens[1];
									break; case SENSE: maximize = tokens[1].matches("MAX") || tokens[1].matches("MAXIMIZE");
									break; case QUADOBJ: parseQuadObj(cols, PrcvR, PrcvC, PrcvV, tokens);
									break; default: throw new IllegalStateException(
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
									prevSection = currentSection;
									break;
								case "OBJSENSE":
									currentSection = Section.SENSE;
									prevSection = currentSection;
									break;
								case "OBJNAME":
									currentSection = Section.OBJ;
									break;
								case "ROWS":
									currentSection = Section.ROWS;
									prevSection = currentSection;
									break;
								case "COLUMNS":
									currentSection = Section.COLUMNS;
									prevSection = currentSection;
									break;
								case "RHS":
									currentSection = Section.RHS;
									prevSection = currentSection;
									break;
								case "RANGES":
									currentSection = Section.RANGES;
									prevSection = currentSection;
									break;
								case "BOUNDS":
									currentSection = Section.BOUNDS;
									prevSection = currentSection;
									break;
								case "QUADOBJ":
									currentSection = Section.QUADOBJ;
									prevSection = currentSection;
									break;
								case "ENDATA":
									if (prevSection == Section.QUADOBJ) // postprocessing of the quadratic objective section
										extractUpperTriangle(PrcvR, PrcvC, PrcvV);
										convertRcvToCsc(PrcvR, PrcvC, PrcvV, cols.size(), Pi, Pp, Px);
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
		double[] P_x = new double[shape[1]];
		int[] P_p = new int[shape[1]+1];
		int[] P_i = new int[shape[1]];
		for (int i = 0; i < Pi.size(); i++) {
			P_i[i] = Pi.get(i);
			P_x[i] = Px.get(i);
		}
		for (int i = 0; i < Pp.size(); i++) {
			P_p[i] = Pp.get(i);
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
        //Parser p = Parser.readMps("src/test/resources/supportcase10.mps"); //Optimal objective  3.383923666e+00
		//Parser p = Parser.readMps("src/test/resources/ex10.mps"); //Optimal objective 100
		//Parser p = Parser.readMps("src/test/resources/savsched1.mps"); //Optimal objective 2.1740357143e+02
		//Parser p = Parser.readMps("src/test/resources/sample1.mps");
		//	String mpsFileDir = "src/test/resources/sample1.mps";
		String qpsFileDir = "src/test/resources/sample1.qps";
		if (args.length >= 1)
			qpsFileDir = args[0];
		else {
			OSQP.LOG.info("Usage: java -jar josqp.jar <mps_file_name>");
			//return;
		}
		//Parser p = Parser.readMps(mpsFileDir);
		Parser p = Parser.readQmps(qpsFileDir);

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
