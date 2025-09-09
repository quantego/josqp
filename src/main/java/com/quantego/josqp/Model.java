package com.quantego.josqp;

import java.util.*;
import java.util.logging.Handler;

/**
 * The class Model holds the formulation of a quadratic programming model. A model can be built algebraically by either
 * calling the static function getBuilder, or by calling the constructor directly. The builder provides an algebraic
 * way to build a model by adding variables, constraints, and an objective function. Use this object to call the
 * solver, retrieve its solution, and change model parameters before a resolve.
 *
 * @author Nils Loehndorf
 *
 */
public class Model {

    /**
     * Algebraic model builder
     */
    public static class Builder {
        Model _m;
        Builder(Model m) {
            _m = m;
        }

        public Variable addVariable() {
            return new Variable(_m, _m.addVar());
        }

        /**
         * Add n variables to the model.
         * @param n number of variables
         * @return  array of variables
         */
        public Variable[] addVariables(int n) {
            Variable[] X = new Variable[n];
            for (int i=0; i<n; i++)
                X[i] = addVariable();
            return X;
        }

        /**
         * To add the constraint to the model, the functions leq, eq, or geq must be called.
         * @return constraint builder
         */
        public Constraint addConstraint() {
            return new Constraint(_m);
        }

        /**
         * To set the objective function, the functions maximize or minimize must be called.
         * @return objective builder
         */
        public Objective setObjective() {
            return new Objective(_m);
        }

        /**
         * Returns the constructed model.
         * @return model
         */
        public Model build() {
            return _m;
        }

    }

    /**
     * Decision variables
     */
    public static class Variable {
       Model _m;
       int _x;
        Variable(Model m, int i) {
            _m = m;
            _x = i;
        }

        /**
         * Variable name
         * @param name
         * @return builder
         */
        Variable name(String name) {
            _m.setName(_x, name);
            return this;
        }

        /**
         * Objective coefficient of this variable.
         * @param value
         * @return builder
         */
        public Variable obj(double value) {
            _m.setObj(_x,value);
            return this;
        }

        /**
         * Set lower bound. By default variables are 0 <= x <= inf.
         * @param value
         * @return builder
         */
        public Variable lb(double value) {
            _m.setLb(_x,value);
            return this;
        }

        /**
         * Set upper bound. By default variables are 0 <= x <= inf.
         * @param value
         * @return builder
         */
        public Variable ub(double value) {
            _m.setUb(_x, value);
            return this;
        }

        /**
         * Defines this variable as free, -inf <= x <= inf.
         * @return builder
         */
        public Variable free() {
            _m.setLb(_x, -OSQP.OSQP_INFTY);
            _m.setUb(_x, OSQP.OSQP_INFTY);
            return this;
        }

        /**
         * Returns the index of this variable in the model.
         * @return  index
         */
        public int getIndex() {
            return _x;
        }

        /**
         * Returns the solution value of this variable.
         * @return solution value
         */
        public double getSolution() {
            return _m.getSolution(_x);
        }
    }

    /**
     * Model objective function.
     */
    public static class Objective {

        static class Pair {
            final Variable x1;
            final Variable x2;
            Pair(Variable x1, Variable x2) {
                this.x1=x1;
                this.x2=x2;
            }

            /**
             * Pairs are equal if they have the same variables.
             * @param o
             * @return
             */
            @Override
            public boolean equals(Object o) {
                Pair pair = (Pair) o;
                return x1==pair.x1 && x2==pair.x2;
            }

            @Override
            public int hashCode() {
                return Objects.hash(x1, x2);
            }
        }

        Model _m;
        Map<Variable, Double> _a = new HashMap<>();
        Map<Pair,Double> _p = new HashMap<>();
        Double _c = 0.;

        Objective(Model m) {
            _m = m;
        }

        /**
         * Add a vector of variables
         * @param A vector of coefficients
         * @param X vector of variables
         * @return builder
         */
        public Objective add(double[] A, Variable[] X) {
            for (int i=0; i<X.length; i++)
                this.add(A[i], X[i]);
            return this;
        }

        /**
         * Add a vector of variables
         * @param a identical coefficient
         * @param X vector of variables
         * @return builder
         */
        public Objective add(double a, Variable[] X) {
            for (Variable x : X)
                this.add(a, x);
            return this;
        }

        /**
         * Add a constant offset to the objective.
         * @param c constant
         * @return builder
         */
        public Objective add(double c) {
            if (_a==null) throw new IllegalStateException("Objective has already been built. Use set instead.");
            _c += c;
            return this;
        }

        /**
         * Add a variable with a coefficient.
         * @param a coefficient
         * @param x variable
         * @return
         */
        public Objective add(double a, Variable x) {
            if (_a==null) throw new IllegalStateException("Objective has already been built. Use set instead.");
            _a.compute(x,(key,val)->val==null?a:val+a);
            return this;
        }

        /**
         * Add quadratic term a*x1*x2
         * @param a coefficient
         * @param x1 variable
         * @param x2 variable
         * @return builder
         */
        public Objective add(double a, Variable x1, Variable x2) {
            if (_a==null) throw new IllegalStateException("Objective has already been built. Use set instead.");
            _p.compute(new Pair(x1,x2),(key,val)->val==null?a:val+a);
            return this;
        }

        /**
         * Set the objective term a*x1*x2
         * @param p coefficient
         * @param x1 variable
         * @param x2 variable
         */
        public void set(double p, Variable x1, Variable x2) {
            if (_a!=null) throw new IllegalStateException("Objective function must first be built by calling maximize or minimize.");
            _m.setQbj(x1._x,x2._x,p);
        }

        /**
         * Set the objective term q*x
         * @param q coefficient
         * @param x variable
         */
        public void set(double q, Variable x) {
            if (_a!=null) throw new IllegalStateException("Objective function must first be built by calling maximize or minimize.");
            _m.setObj(x._x,q);
        }

        /**
         * Set the objective term c
         * @param c constant
         */
        public void set(double c) {
            _m.setOffset(c);
        }

        /**
         * Set the objective sense to maximize
         * @param maximize true if problem is a maximization problem
         */
        public void setSense(boolean maximize) {
            if (_a!=null) throw new IllegalStateException("Objective function must first be built by calling maximize or minimize.");
            _m.setSense(maximize);
        }

        /**
         * Build the objective function and add it to the model as a maximization problem.
         */
        public void maximize() {
            if (_a==null) throw new IllegalStateException("Objective has already been built.");
            addToModel();
            _m.setSense(true);
            _m.setOffset(_c);
            _a = null; _p=null; _c=null;
        }

        /**
         * Build the objective function and add it to the model as a minimization problem.
         */
        public void minimize() {
            if (_a==null) throw new IllegalStateException("Objective has already been built.");
            addToModel();
            _m.setSense(false);
            _m.setOffset(_c);
            _a = null; _p=null; _c=null;
        }

        private void addToModel() {
            for (Map.Entry<Variable,Double> e : _a.entrySet()) {
                _m.setObj(e.getKey()._x,e.getValue());
            }
            for (Map.Entry<Pair,Double> e : _p.entrySet()) {
                _m.setQbj(e.getKey().x1._x,e.getKey().x2._x,e.getValue());
            }
        }

        /**
         * Returns the objective value of the solution.
         * @return objective value
         */
        public double getValue() {
            return _m.getObj();
        }
    }

    /**
     * Model constraint. Once functions leq, eq, or geq are called, the constraint is linked with the model, and parameters
     * can be changed.
     */
    public static class Constraint {

        /**
         * Constraint types
         * EQ: Ax == b
         * GEQ: Ax >= b
         * LEQ: Ax <= b
         * NEQ: Ax != b
         */
        public static enum Type {
            EQ, GEQ, LEQ, NEQ
        }

        Map<Variable, Double> _a = new HashMap<>();
        Model _m;
        Integer _y = null;

        Constraint(Model m) {
            _m = m;
        }

        /**
         * Add a vector of variables
         * @param A
         * @param X
         * @return
         */
        public Constraint add(double[] A, Variable[] X) {
            for (int i=0; i<X.length; i++)
                this.add(A[i], X[i]);
            return this;
        }

        /**
         * Add a vector of variables
         * @param a identical coefficient
         * @param X vector of variables
         * @return
         */
        public Constraint add(double a, Variable[] X) {
            for (Variable x : X)
                this.add(a, x);
            return this;
        }

        /**
         * Add a constant offset to the constraint.
         * @param a constant
         * @param x variable
         * @return builder
         */
        public Constraint add(double a, Variable x) {
            if (_y!=null) throw new IllegalStateException("Constraint has already been added to the model by function geq.");
            _a.compute(x,(key,val)->val==null?a:val+a);
            return this;
        }

        /**
         * Add the constraint to the model as a <= b.
         * @param b right hand side
         */
        public void leq(double b) {
            if (_y!=null) throw new IllegalStateException("Constraint has already been added to the model by function leq.");
            _y = _m.addCtr(Constraint.Type.LEQ,b);
            for (Map.Entry<Variable,Double> e : _a.entrySet()) {
                _m.setLhs(_y,e.getKey()._x,e.getValue());
            }
            _a = null;
        }

        /**
         * Add the constraint to the model as a == b.
         * @param b right hand side
         */
        public void eq(double b) {
            if (_y!=null) throw new IllegalStateException("Constraint has already been added to the model by function eq.");
            _y = _m.addCtr(Constraint.Type.EQ,b);
            for (Map.Entry<Variable,Double> e : _a.entrySet()) {
                _m.setLhs(_y,e.getKey()._x,e.getValue());
            }
            _a = null;
        }

        /**
         * Add the constraint to the model as a >= b.
         * @param b right hand side
         */
        public void geq(double b) {
            if (_y!=null) throw new IllegalStateException("Constraint has already been added to the model by function geq.");
            _y = _m.addCtr(Constraint.Type.GEQ,b);
            for (Map.Entry<Variable,Double> e : _a.entrySet()) {
                _m.setLhs(_y,e.getKey()._x,e.getValue());
            }
            _a = null;
        }

        /**
         * Sets the left-hand side of the constraint for the given variable.
         * @param x variable
         * @param a coefficient
         */
        public void setLhs(Variable x, double a) {
            if (_y==null) throw new IllegalStateException("Constraint has not yet been added to model. Call leq, eq, geq, before calling lhs or rhs.");
            _m.setLhs(_y,x._x,a);
        }

        /**
         * Sets the right-hand side of the constraint.
         * @param b right-hand side
         */
        public void setRhs(double b) {
            if (_y==null) throw new IllegalStateException("Constraint has not yet been added to model. Call leq, eq, geq, before calling lhs or rhs.");
            _m.setRhs(_y,b);
        }

        /**
         * Sets the type of the constraint.
         * @param type constraint type
         */
        public void setType(Constraint.Type type) {
            if (_y==null) throw new IllegalStateException("Constraint has not yet been added to model. Call leq, eq, geq, before calling lhs or rhs.");
            _m.setType(_y,type);
        }

        /**
         * Returns the index of the constraint in the model.
         * @return index
         */
        public int getIndex() {
            return _y;
        }

        /**
         * Returns the dual solution value associated with this constraint.
         * @return dual solution value
         */
        public double getSolution() {
            return _m.getDual(_y);
        }

    }

    OSQP _solver;
    OSQP.Settings _param;
    OSQP.Data _data;
    CSCMatrixBuilder _A;
    CSCMatrixBuilder _P;
    List<Integer> _cr;
    List<Integer> _rr;
    List<Double> _ub;
    List<Double> _lb;
    List<Double> _q;
    List<Constraint.Type> _type;
    List<Double> _rhs;
    Map<Integer, String> _varNames;
    int _nCols;
    int _nRows;
    boolean _maximize;
    boolean _newQ;
    boolean _newP;
    boolean _newL;
    boolean _newU;
    boolean _newA;
    boolean _newCols;
    boolean _newRows;
    double[] _x;
    double[] _y;
    double _offset;
    public Model() {
        _A = new CSCMatrixBuilder();
        _P = new CSCMatrixBuilder();
        _ub = new ArrayList<>();
        _lb = new ArrayList<>();
        _q = new ArrayList<>();
        _cr = new ArrayList<>();
        _rr = new ArrayList<>();
        _type = new ArrayList<>();
        _rhs = new ArrayList<>();
        _param = new OSQP.Settings();
    }

    /**
     * Returns a builder object that can be used to build a model algebraically by adding variables, constraints, and
     * an objective function.
     * @return builder
     */
    public static Model.Builder getBuilder() {
        return new Model.Builder(new Model());
    }

    /**
     * Add a variable to the model and return its index.
     * @return index of the variable
     */
    public int addVar() {
        _q.add(0.);
        int i = addCtr(Constraint.Type.GEQ, 0.0);
        _A.set(i,_nCols,1.0);
        _P.set(_nCols,_nCols,0.0);
        _cr.add(i);
        _newCols = true;
        return _nCols++;
    }

    /**
     * Set the name of the variable at the given index
     * @param var index of the variable
     * @param name name of the variable
     */
    public void setName(int var, String name) {
        if (_varNames==null)
            _varNames = new HashMap<>();
        _varNames.put(var,name);
    }

    /**
     * Add a constraint to the model and return its index.
     * @param type constraint type
     * @param rhs right-hand side
     * @return index of the constraint
     */
    public int addCtr(Constraint.Type type, double rhs) {
        switch (type) {
            case EQ: _lb.add(rhs); _ub.add(rhs); break;
            case LEQ:  _lb.add(-OSQP.OSQP_INFTY); _ub.add(rhs); break;
            case GEQ:  _lb.add(rhs); _ub.add(OSQP.OSQP_INFTY);break;
            case NEQ:  _lb.add(-OSQP.OSQP_INFTY); _ub.add(OSQP.OSQP_INFTY); break;
        }
        _rr.add(_nRows);
        _rhs.add(rhs);
        _type.add(type);
        _newRows = true;
        return _nRows++;
    }

    /**
     * Set the type of the constraint at the given index.
     * @param ctr index of the constraint
     * @param type constraint type
     */
    public void setType(int ctr, Constraint.Type type) {
        if (ctr>=_nRows) throw new IllegalArgumentException(String.format("Constraint %d is unknown.",ctr));
        if (type==_type.get(ctr)) return;
        double rhs = _rhs.get(ctr);
        switch (type) {
            case EQ: _lb.set(ctr,rhs); _ub.set(ctr,rhs); break;
            case LEQ:  _lb.set(ctr,-OSQP.OSQP_INFTY); _ub.set(ctr,rhs); break;
            case GEQ:  _lb.set(ctr,rhs); _ub.set(ctr,OSQP.OSQP_INFTY);break;
            case NEQ:  _lb.set(ctr,-OSQP.OSQP_INFTY); _ub.set(ctr,OSQP.OSQP_INFTY); break;
        }
        _type.set(ctr, type);
        _newL = true; _newU = true;
    }

    /**
     * Set the right-hand side of the constraint at the given index.
     * @param ctr index of the constraint
     * @param value right-hand side
     */
    public void setRhs(int ctr, double value) {
        if (ctr>=_nRows) throw new IllegalArgumentException(String.format("Constraint %d is unknown.",ctr));
        if (value==_rhs.get(ctr)) return;
        switch (_type.get(ctr)) {
            case EQ: _lb.set(ctr,value); _ub.set(ctr,value); _newL = true; _newU = true; break;
            case LEQ:  _lb.set(ctr,-OSQP.OSQP_INFTY); _ub.set(ctr,value); _newU=true; break;
            case GEQ:  _lb.set(ctr,value); _ub.set(ctr,OSQP.OSQP_INFTY); _newL=true; break;
            case NEQ:  _lb.set(ctr,-OSQP.OSQP_INFTY); _ub.set(ctr,OSQP.OSQP_INFTY); break;
        }
        _newL = true; _newU = true;
    }

    /**
     * Set the left-hand side of the constraint at the given index.
     * @param ctr index of the constraint
     * @param var index of the variable
     * @param value coefficient
     */
    public void setLhs(int ctr, int var, double value) {
        if (ctr>=_nRows) throw new IllegalArgumentException(String.format("Constraint %d is unknown.",ctr));
        if (var>=_nCols) throw new IllegalArgumentException(String.format("Variable %d is unknown.",var));
        _A.set(ctr,var,value);
        _newA = true;
    }

    /**
     * Set the quadratic coefficient of the objective function.
     * @param var index of the variable
     * @param value coefficient
     */
    public void setObj(int var, double value) {
        if (var >= _nCols) throw new IllegalArgumentException(String.format("Variable %d is unknown.",var));
        _q.set(var,value);
        _newQ = true;
    }

    /**
     * Set the quadratic coefficient of the objective function.
     * @param var1 index of the first variable
     * @param var2 index of the second variable
     * @param value coefficient
     */
    public void setQbj(int var1, int var2, double value) {
        if (var1 >= _nCols) throw new IllegalArgumentException(String.format("Variable %d is unknown.",var1));
        if (var2 >= _nCols) throw new IllegalArgumentException(String.format("Variable %d is unknown.",var2));
        _P.set(var1,var2,2*value);
        _newP = true;
    }

    /**
     * Set the upper bound of the variable at the given index.
     * @param var index of the variable
     * @param value upper bound
     */
    public void setUb(int var, double value) {
        if (var >= _nCols) throw new IllegalArgumentException(String.format("Variable %d is unknown.",var));
        if (value>OSQP.OSQP_INFTY) value=OSQP.OSQP_INFTY;
        _ub.set(_cr.get(var), value);
        _newU = true;
    }

    /**
     * Set the lower bound of the variable at the given index.
     * @param var index of the variable
     * @param value lower bound
     */
    public void setLb(int var, double value) {
        if (var >= _nCols) throw new IllegalArgumentException(String.format("Variable %d is unknown.",var));
        if (value<-OSQP.OSQP_INFTY) value=-OSQP.OSQP_INFTY;
        _lb.set(_cr.get(var), value);
        _newL = true;
    }

    /**
     * Set the offset of the objective function.
     * @param value offset
     */
    public void setOffset(double value) {
        _offset = value;
    }

    /**
     * Set the sense of the objective function.
     * @param maximize true if problem is a maximization problem
     */
    public void setSense(boolean maximize) {
        _maximize = maximize;
    }

    /**
     * Returns the settings of the solver.
     * @return settings
     */
    public OSQP.Settings getParam() {
        return _param;
    }

    void init(OSQP.Settings settings) {
        if (settings != null)
            _param = settings;
        else
            _param = new OSQP.Settings();

        _data = new OSQP.Data(_nCols,
                _nRows,
                _P.build(_maximize),
                _A.build(false),
                _maximize ? Utils.toDoubleArrayNeg(_q) : Utils.toDoubleArray(_q),
                Utils.toDoubleArray(_lb),
                Utils.toDoubleArray(_ub),
                0
        );
        _solver = new OSQP(_data,_param);
    }

    /**
     * Solve the model using the default parameter settings.
     * @return status of the solver
     */
    public OSQP.Status solve() {
        return solve(new OSQP.Settings());
    }

    /**
     * Solve the model using customized parameter settings.
     * @param settings the customized solver settings (or {@code null} for the default settings).
     * @return status of the solver
     */
    public OSQP.Status solve(OSQP.Settings settings) {
        if (_solver==null) {
            init(settings);
        }
        else {
            if (_newCols || _newRows) {
                init(settings);
                if (_newCols)
                    _x = Arrays.copyOf(_x,_nCols);
                if (_newRows)
                    _y = Arrays.copyOf(_y,_nRows);
                _solver.warm_start(_x,_y);
            }
            else {
                if (_newQ) {
                    _solver.update_lin_cost(_maximize ? Utils.toDoubleArrayNeg(_q) : Utils.toDoubleArray(_q));
                }
                if (_newL || _newU)
                    _solver.update_bounds(Utils.toDoubleArray(_lb), Utils.toDoubleArray(_ub));
                if (_newA) {
                    _A.update(false);
                    _solver.update_A(_A.getElements(), _A.getIndices());
                }
                if (_newP) {
                    _P.update(_maximize);
                    _solver.update_P(_P.getElements(), _P.getIndices());
                }
            }
        }
        _newA = false; _newP=false; _newQ = false; _newL = false; _newU = false; _newCols = false; _newRows=false; 
        _solver.solve();
        _x = _solver.getPrimalSolution();
        _y = _solver.getDualSolution();
        return _solver.getInfo().status;
    }

    /**
     * Add another log handler to the solver.
     * @param h log handler
     */
    public void addLogHandler(Handler h) {
        OSQP.LOG.addHandler(h);
    }

    /**
     * Return the solution of the variable at the given index.
     * @param var index of the variable
     * @return solution value
     */
    public double getSolution(int var) {
        if (_x==null) throw new IllegalStateException("Solver does not have a primal solution. Call solve first.");
        return _x[var];
    }

    /**
     * Return the dual solution of the constraint at the given index.
     * @param ctr index of the constraint
     * @return dual solution
     */
    public double getDual(int ctr) {
        if (_x==null) throw new IllegalStateException("Solver does not have a dual solution. Call solve first.");
        return (_maximize?-1.0:1.0)*_y[ctr];
    }

    /**
     * Return the objective value of the solution.
     * @return objective value
     */
    public double getObj() {
        return (_maximize?-1:1)*_solver.getObjectiveValue()+_offset;
    }

    static String termToString(double d, String s) {
		if (d==0.0) return "";
		if (d==1.0)
			return " + "+s;
		if (d==-1.0)
			return " - "+s;
		if (d>0.0)
			return String.format(" + %.5f %s",d,s);
		return String.format(" - %.5f %s",-d,s);
	}

    /**
     * Returns a string of the full model in .lp format
     * @return model string
     */
    @Override
	public String toString() {
		StringBuilder modelString = new StringBuilder();
		modelString.append(_maximize ? "Maximize" : "Minimize").append("\nobj:");
		if (_offset != 0.0) modelString.append(" "+_offset);
		//objective function
		for (int col=0; col<_nCols; col++) {
			double c = _q.get(col);
			modelString.append(termToString(c,getName(col)));
		}
		if (_P != null) {
			modelString.append(" + [");
            _P.update(_maximize);
            int[] starts = _P.getStarts();
            double[] elements = _P.getElements();
            int[] index = _P.getIndices();
            for (int col=0; col<_nCols; col++) {
                int begin = starts[col];
                int end = starts[col+1];
                for (int j=begin; j<end; j++) {
                    int col2 = index[j];
                    double element = elements[j];
                    if (element != 0.0) {
                        if (col!=col2)
                            modelString.append(termToString(element,getName(col)+" "+getName(col2)));
                        else
                            modelString.append(termToString(element,getName(col)+ "^2"));

                    }

                }
            }
			modelString.append(" ] / 2");
		}
		modelString.append("\nSubject To\n");
		//constraints
		List<StringBuilder> constraintStrings = new ArrayList<>(_nRows);
		for (int row=0; row<_nRows; row++)
			constraintStrings.add(new StringBuilder());
        _A.update(false);
		int[] starts = _A.getStarts();
        double[] elements = _A.getElements();
		int[] index = _A.getIndices();
		for (int col=0; col<_nCols; col++) {
			int begin = starts[col];
			int end = starts[col+1];
			for (int j=begin; j<end; j++) {
				int row = index[j];
				double element = elements[j];
				constraintStrings.get(row).append(termToString(element,getName(col)));
			}
		}
		for (Integer row : _rr) {
			double lb = _lb.get(row);
			double ub = _ub.get(row);
			if (Double.compare(lb,ub)==0)
				modelString.append(constraintStrings.get(row).append(" = ").append(lb).append("\n"));
			else if (lb==-OSQP.OSQP_INFTY && ub<OSQP.OSQP_INFTY)
				modelString.append(constraintStrings.get(row).append(" <= ").append(ub).append("\n"));
			else if (ub==OSQP.OSQP_INFTY && lb>-OSQP.OSQP_INFTY)
				modelString.append(constraintStrings.get(row).append(" >= ").append(lb).append("\n"));
		}
//		//bounds
		modelString.append("Bounds\n");
		for (int col=0; col<_nCols; col++) {
			double lb = _lb.get(_cr.get(col));
			double ub = _ub.get(_cr.get(col));
			if (lb==0 && ub<OSQP.OSQP_INFTY)
				modelString.append(getName(col)).append(" <= ").append(ub).append("\n");
			else if (lb<=-OSQP.OSQP_INFTY && ub>=OSQP.OSQP_INFTY)
				modelString.append("-inf <= ").append(getName(col)).append(" <= inf\n");
			else if (lb!=0 && lb>-OSQP.OSQP_INFTY && ub>=OSQP.OSQP_INFTY)
				modelString.append(lb).append(" <= ").append(getName(col)).append(" <= inf").append("\n");
			else if (lb<=-OSQP.OSQP_INFTY && ub<OSQP.OSQP_INFTY)
				modelString.append("-inf <= ").append(getName(col)).append(" <= ").append(ub).append("\n");
			else if (lb>-OSQP.OSQP_INFTY && ub<OSQP.OSQP_INFTY)
				modelString.append(lb).append(" <= ").append(getName(col)).append(" <= ").append(ub).append("\n");
		}
		modelString.append("End");
		return modelString.toString();
	}

    /**
     * Returns the name of the variable at the given index.
     * @param var index of the variable
     * @return variable name
     */
    public String getName(int var) {
        if (_varNames==null || !_varNames.containsKey(var))
            return "x_"+var;
        return _varNames.get(var);
    }

}
