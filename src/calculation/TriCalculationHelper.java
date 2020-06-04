package calculation;

import Jama.Matrix;

public class TriCalculationHelper{
	public static final int STRAIN = 0;
	public static final int STRESS = 1;
	
	private static TriCalculationHelper calc;
	
	private TriCalculationHelper() {
	}
	
	public static TriCalculationHelper getInstance() {
		if (calc==null) {
			calc = new TriCalculationHelper();
		}
		return calc;
	}
	
	public Matrix computeConstitutiveMatrix(double youngModulus,
			double shearConstant, int kindOfProblem) {
		double c1;
		Matrix res=null;
		
		switch (kindOfProblem) {
		case STRESS:
			c1 = youngModulus/(1f-shearConstant*shearConstant);
			res = new Matrix(new double[][] {
					new double[] {1d, shearConstant, 0d},
					new double[] {shearConstant, 1d, 0d},
					new double[] {0d, 0d, (1d-shearConstant)/2d}
			});
			res=res.times(c1);
			break;
		case STRAIN:
			c1 = youngModulus/((1+shearConstant)*(1-2*shearConstant));
			res = new Matrix(new double[][] {
					new double[] {1d-shearConstant, shearConstant, 0d},
					new double[] {shearConstant, 1d-shearConstant, 0d},
					new double[] {0d, 0d, (1d-2d*shearConstant)/2d}
			});
			res=res.times(c1);
			break;

		}
		return res;
	}
	
	public Matrix constraints(int[] bcdof, Matrix stiffness) {
		Matrix res = new Matrix(stiffness.getArrayCopy());
		
		for (int i = 0; i < bcdof.length; i++) {
			int c = bcdof[i]-1;
			for (int j = 0; j < stiffness.getRowDimension(); j++) {
				res.set(c, j, 0);
			}
			res.set(c, c, 1);
		}
		
		return res;
	}
	
	public Matrix computeDisplacementMatrix(int[] bcdof, Matrix stiffness,
			Matrix forces) {
		Matrix U = new Matrix(forces.getRowDimension(), 1);
		Matrix ssm = constraints(bcdof, stiffness);
		U = ssm.inverse().times(forces);
		return U;
	}
	
	
	public double computeAreaElements(double x[], double y[]) {
		return ( (x[0]-x[2])*(y[1]-y[2]) - (x[1]-x[2])*(y[0]-y[2]) )/2f;
	}
	
	public Matrix computeDifferentOperator(double x[], double y[]) {
		Matrix B = new Matrix(3, 6);
		B.set(0, 0, y[1]-y[2]);
		B.set(2, 1, y[1]-y[2]);
		
		B.set(0, 2, y[2]-y[0]);
		B.set(2, 3, y[2]-y[0]);
		
		B.set(0, 4, y[0]-y[1]);
		B.set(2, 5, y[0]-y[1]);
		
		B.set(1, 1, x[2]-x[1]);
		B.set(2, 0, x[2]-x[1]);
		
		B.set(1, 3, x[0]-x[2]);
		B.set(2, 2, x[0]-x[2]);
		
		B.set(1, 5, x[1]-x[0]);
		B.set(2, 4, x[1]-x[0]);
		
		return B.times( 1f/(2f*computeAreaElements(x, y)) );
	}
	
	public Matrix computeElementStiffness(Matrix C, Matrix vertex, double thickness) {
		Matrix K = new Matrix(6, 6);
		double x[] = new double[]{
				vertex.get(0, 0), vertex.get(1, 0), vertex.get(2, 0) 
		};
		double y[] = new double[]{
				vertex.get(0, 1), vertex.get(1, 1), vertex.get(2, 1) 
		};
		Matrix B = computeDifferentOperator(x, y);
		K = B.transpose().times(C).times(B).times(computeAreaElements(x, y)).times(thickness);
		
		return K;
	}
	
	public Matrix getIndexMatrix(double[] nd, int nnel, int ndof) {
		Matrix indexMatrix = new Matrix(nnel*ndof,1);
		int k=0;
		for (int i = 0; i < nnel; i++) {
			int start = ((int)nd[i]-1)*ndof;
			for (int j = 0; j < ndof; j++) {
				indexMatrix.set(k,0, start+j);
				k++;
			}
		}
		return indexMatrix;
	}
	
	public void assembleStiffnessMatrix (int nn, Matrix stiffness, Matrix K, Matrix index) {
		if (stiffness==null) {
			return;
		}
		
		int edof = index.getRowDimension();
		for (int i = 0; i < edof; i++) {
			int ii = (int)index.get(i,0);
			for (int j = 0; j < edof; j++) {
				int jj = (int)index.get(j,0);
				stiffness.set(ii, jj, stiffness.get(ii, jj)+K.get(i, j));
			}
		}
	}
	
	public Matrix getVertexOfElement(Matrix coordinates, Matrix elements, int ioe) {
		int nrow = elements.getColumnDimension();
		int ncol = coordinates.getColumnDimension();
		Matrix vertex = new Matrix(nrow, ncol);
		
		double[] nd = elements.getArray()[ioe];
		for (int i = 0; i < nd.length; i++) {
			vertex.set(i, 0, coordinates.get((int)nd[i]-1, 0));
			vertex.set(i, 1, coordinates.get((int)nd[i]-1, 1));
		}
		
		return vertex;
	}
	
	public void computeStressStrainElements(Matrix C, Matrix U, 
			int ioe, Matrix elements, Matrix coordinates,
			Matrix stress, Matrix strain) {	 // output
		Matrix ue = new Matrix(6, 1);
		double[] element = elements.getArray()[ioe];
		int uei = 0;
		for (int i = 0; i < element.length; i++) {
			int Ui = (int)element[i]*2-2;
			ue.set(uei++, 0, U.get(Ui, 0));
			ue.set(uei++, 0, U.get(Ui+1, 0));
		}
		Matrix vertex = getVertexOfElement(coordinates, elements, ioe);
		double x[] = new double[]{
				vertex.get(0, 0), vertex.get(1, 0), vertex.get(2, 0) 
		};
		double y[] = new double[]{
				vertex.get(0, 1), vertex.get(1, 1), vertex.get(2, 1) 
		};
		Matrix B = computeDifferentOperator(x, y);
		Matrix temp = new Matrix(3,1);
		strain.setMatrix(0, 2, 0, 0, B.times(ue));
		stress.setMatrix(0, 2, 0, 0, C.times(strain));
	}
	
	public static void main(String[] args) {
		TriCalculationHelper calc = TriCalculationHelper.getInstance();
		double youngModulus = 2e11;
		double shearConstant = 0.3;
		double thickness = 0.1;
		
//		Matrix Cev = calc.computeConstitutiveMatrix(youngModulus, shearConstant, CalculationHelper.STRAIN);
		Matrix Ces = calc.computeConstitutiveMatrix(youngModulus, shearConstant, TriCalculationHelper.STRESS);
		
		int nnel = 3;
		int nn = 4;
		int ne = 2;
		int ndof = 2;
		
		Matrix coordinates = new Matrix(new double[][] {
				new double[] {0d,0d},
				new double[] {4d,0d},
				new double[] {4d,2d},
				new double[] {0d,2d},
		});
		Matrix elements = new Matrix(new double[][] {
				new double[] {1,2,4},
				new double[] {2,3,4},
		});
		int[] bcdof = new int[]{1,2,3,4,7};
		Matrix forces = new Matrix(nn*2,1);
		forces.set(5, 0, 2000);
		forces.set(7, 0, 2000);
		
		Matrix stiffness = new Matrix(2*nn, 2*nn);
		for (int i = 0; i < ne; i++) {
			Matrix vertex = calc.getVertexOfElement(coordinates, elements, i);
			Matrix k = calc.computeElementStiffness(Ces, vertex, thickness);
			
			double[] nd = elements.getArray()[i];
			
			Matrix index = calc.getIndexMatrix(nd, nnel, ndof);
			calc.assembleStiffnessMatrix(nn, stiffness, k, index);
		}
//		stiffness = stiffness.times(1e6); 
//		stiffness.print(2*nn, 2*nn);
//		
		Matrix U = calc.computeDisplacementMatrix(bcdof, stiffness, forces);
//		U = U.times(1e11);
//		U.print(U.getRowDimension(), 1);
		for (int i = 0; i < elements.getRowDimension(); i++) {
			Matrix strain = new Matrix(3,1);
			Matrix stress = new Matrix(3,1);
			calc.computeStressStrainElements(Ces, U, i, elements, coordinates, stress, strain);
//			strain = strain.times(1e11);
//			stress = stress.times(1e11);
			strain.print(3, 1);
			stress.print(3, 1);
		}
	}
}
