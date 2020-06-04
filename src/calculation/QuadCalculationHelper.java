package calculation;

import Jama.Matrix;

public class QuadCalculationHelper {
	public static final int STRAIN = 0;
	public static final int STRESS = 1;

	private static QuadCalculationHelper calc;

	private QuadCalculationHelper() {
	}

	public static QuadCalculationHelper getInstance() {
		if (calc==null) {
			calc = new QuadCalculationHelper();
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

	public void computeShapeFunction(double xi, double yi,
			Matrix shape, Matrix dhdu, Matrix dhdv) {
		/*shape = new Matrix(4,1);
		dhdu = new Matrix(4,1);
		dhdv = new Matrix(4,1);*/

		shape.set(0, 0, 0.25d*(1d-xi)*(1d-yi) );
		shape.set(1, 0, 0.25*(1d+xi)*(1d-yi) );
		shape.set(2, 0, 0.25d*(1d+xi)*(1d+yi) );
		shape.set(3, 0, 0.25d*(1d-xi)*(1d+yi) );

		dhdu.set(0, 0, -0.25d*(1d-yi) );
		dhdu.set(1, 0,  0.25d*(1d-yi) );
		dhdu.set(2, 0,  0.25d*(1d+yi) );
		dhdu.set(3, 0, -0.25d*(1d+yi) );

		dhdv.set(0, 0, -0.25d*(1d-xi) );
		dhdv.set(1, 0, -0.25d*(1d+xi) );
		dhdv.set(2, 0,  0.25d*(1d+xi) );
		dhdv.set(3, 0,  0.25d*(1d-xi) );

	}

	public Matrix computeJacobian(int nnel, Matrix dhdu, Matrix dhdv,
			Matrix vertex) {
		Matrix jacobian = new Matrix(2,2);

		for (int i = 0; i < nnel; i++) {
			jacobian.set(0, 0, jacobian.get(0, 0)+dhdu.get(i, 0)*vertex.get(i, 0));
			jacobian.set(0, 1, jacobian.get(0, 1)+dhdu.get(i, 0)*vertex.get(i, 1));
			jacobian.set(1, 0, jacobian.get(1, 0)+dhdv.get(i, 0)*vertex.get(i, 0));
			jacobian.set(1, 1, jacobian.get(1, 1)+dhdv.get(i, 0)*vertex.get(i, 1));
		}

		return jacobian;
	}

	public void computeShapeDerivatives(int nnel, Matrix dhdu, Matrix dhdv, Matrix invJacobian,
			Matrix dhdx, Matrix dhdy) {
		if (nnel==4) {
			for (int i = 0; i < nnel; i++) {
				dhdx.set(i, 0, invJacobian.get(0, 0)*dhdu.get(i, 0) + invJacobian.get(0, 1)*dhdv.get(i, 0));
				dhdy.set(i, 0, invJacobian.get(1, 0)*dhdu.get(i, 0) + invJacobian.get(1, 1)*dhdv.get(i, 0));
			}
		} else {
			System.out.println("nnel == "+nnel);
		}
	}

	public Matrix computeDifferentOperator(int nnel, Matrix dhdx, Matrix dhdy) {
		Matrix B = new Matrix(3, 2*nnel);
		for (int i = 0; i < nnel; i++) {
			int i1 = i*2;
			int i2 = i1+1;

			B.set(0, i1, dhdx.get(i, 0));
			B.set(1, i2, dhdy.get(i, 0));
			B.set(2, i1, dhdy.get(i, 0));
			B.set(2, i2, dhdx.get(i, 0));
		}
		return B;
	}

	public Matrix computeElementStiffness(Matrix C, Matrix vertex, double[] gausspoint, double thickness) {
		int nnel = 4;
		Matrix K = new Matrix(2*nnel, 2*nnel);
		Matrix dhdu = new Matrix(nnel, 1);
		Matrix dhdv = new Matrix(nnel, 1);
		Matrix shape = new Matrix(nnel, 1);

		for (int i=0; i<gausspoint.length; i++) {
			double xi = gausspoint[i];
			for (int j=0; j<gausspoint.length; j++) {
				double yi = gausspoint[j];
				computeShapeFunction(xi, yi, shape, dhdu, dhdv);
				Matrix jacobian = computeJacobian(nnel, dhdu, dhdv, vertex);

				Matrix dhdx = new Matrix(nnel, 1);
				Matrix dhdy = new Matrix(nnel, 1);
				computeShapeDerivatives(nnel, dhdu, dhdv, jacobian.inverse(), dhdx, dhdy);

				Matrix B = computeDifferentOperator(nnel, dhdx, dhdy);
//				Matrix test = new Matrix(B.getArrayCopy());
//				test.print(B.getRowDimension(), B.getColumnDimension());

				K = K.plus(B.transpose().times(C).times(B).times(jacobian.det())).times(thickness);
			}

		}
//		System.out.println("-------------------------------");

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

	public Matrix filterByNode(Matrix elements, int nodeIndex) {
		int erow = elements.getRowDimension();
		double[][] e = elements.getArrayCopy();
		double[][] ne = new double[erow][4];

		int nerow = 0;
		for (int i = 0; i < erow; i++) {
			for (int j = 0; j < 4; j++) {
				if (e[i][j]==nodeIndex) {
					ne[nerow] = e[i];
					nerow++;
					break;
				}
			}
		}
		Matrix eles = new Matrix(nerow, 4);
		for (int i = 0; i < nerow; i++) {
			for (int j = 0; j < 4; j++) {
				eles.set(i, j, ne[i][j]);
			}
		}
		return eles;
	}

	public Matrix computeNodeDiffOperator(double[] gausspoint,
			Matrix elements, Matrix coordinates,
			int ioe, int ni) throws IllegalArgumentException {
		Matrix B = new Matrix(3, 8);
		Matrix dhdu = new Matrix(4, 1);
		Matrix dhdv = new Matrix(4, 1);
		Matrix shape = new Matrix(4, 1);
		double[][] gpValue = new double[][] {
				new double[] {gausspoint[0],gausspoint[0]},
				new double[] {gausspoint[0],gausspoint[1]},
				new double[] {gausspoint[1],gausspoint[0]},
				new double[] {gausspoint[1],gausspoint[1]}
		};
		double[] e = elements.getArray()[ioe];
		int noe = 0;
		while(noe<4 && e[noe]!=ni) noe++;

		if (noe==4) throw new IllegalArgumentException("node not found - ni="+ni+
				" e="+String.format("%.0f,%.0f,%.0f,%.0f",
						e[0],e[1],e[2],e[3]));

		computeShapeFunction(gpValue[noe][0], gpValue[noe][1], shape, dhdu, dhdv);
		Matrix vertex = calc.getVertexOfElement(coordinates, elements, ioe);
		Matrix jacobian = computeJacobian(4, dhdu, dhdv, vertex);
		Matrix dhdx = new Matrix(4, 1);
		Matrix dhdy = new Matrix(4, 1);
		computeShapeDerivatives(4, dhdu, dhdv, jacobian.inverse(), dhdx, dhdy);

		B = computeDifferentOperator(4, dhdx, dhdy);
		return B;
	}

	public void computeStressStrainNodeInElement(
			Matrix U, Matrix C,	double[] gausspoint,
			Matrix coordinates, Matrix elements, int ioe, int ni,
			Matrix stress, Matrix strain) {	 // output
		Matrix ue = new Matrix(8, 1);
		double[] element = elements.getArray()[ioe];
		int uei = 0;
		for (int i = 0; i < element.length; i++) {
			int Ui = (int)element[i]*2-2;
			ue.set(uei++, 0, U.get(Ui, 0));
			ue.set(uei++, 0, U.get(Ui+1, 0));
		}

		Matrix B = computeNodeDiffOperator(gausspoint, elements, coordinates, ioe, ni);
		strain.setMatrix(0, 2, 0, 0, B.times(ue));
		stress.setMatrix(0, 2, 0, 0, C.times(strain));

	}

	public void computeStressStrainNode(
			Matrix U, Matrix C, double[] gausspoint,
			Matrix coordinates, Matrix elements, int ni,
			Matrix stress, Matrix strain) {
		Matrix es = filterByNode(elements, ni);
		Matrix tstress = new Matrix(3,1);
		Matrix tstrain = new Matrix(3,1);
		for (int i = 0; i < es.getRowDimension(); i++) {
			Matrix nodeStre = new Matrix(3,1);
			Matrix nodeStra = new Matrix(3,1);
			computeStressStrainNodeInElement(U, C, gausspoint, coordinates, es, i, ni, nodeStre, nodeStra);
			tstress = stress.plus(nodeStre);
			tstrain = strain.plus(nodeStra);
		}
		tstress = tstress.times(1f/es.getRowDimension());
		tstrain = tstrain.times(1f/es.getRowDimension());

		for (int i = 0; i < 3; i++) {
			stress.set(i, 0, tstress.get(i, 0));
			strain.set(i, 0, tstrain.get(i, 0));
		}
		
	}

	public static void main(String[] args) {
		QuadCalculationHelper calc = QuadCalculationHelper.getInstance();
		double youngModulus = 2e11;
		double shearConstant = 0.3;
		double thickness = 0.1;

//		Matrix Cev = calc.computeConstitutiveMatrix(youngModulus, shearConstant, CalculationHelper.STRAIN);
		Matrix Ces = calc.computeConstitutiveMatrix(youngModulus, shearConstant, QuadCalculationHelper.STRESS);

		int nnel = 4;
		int nn = 9;
		int ne = 4;
		int ndof = 2;

		Matrix coordinates = new Matrix(new double[][] {
				new double[] {0d,0d},
				new double[] {10d,0d},
				new double[] {20d,0d},
				new double[] {0d,10d},
				new double[] {10d,10d},
				new double[] {20d,10d},
				new double[] {0d,20d},
				new double[] {10d,20d},
				new double[] {20d,20d}

		});
		Matrix elements = new Matrix(new double[][] {
				new double[] {1,2,5,4},
				new double[] {2,3,6,5},
				new double[] {4,5,8,7},
				new double[] {5,6,9,8}
		});

		double[] gausspoint = new double[] {
				-1d/Math.sqrt(3),
				1d/Math.sqrt(3),
		};
		int[] bcdof = new int[]{1,2,7,8,13,14};
		Matrix forces = new Matrix(nn*2,1);
		forces.set( 5, 0, 10000);
		forces.set( 11, 0, 10000);
		forces.set(17, 0, 10000);

		Matrix stiffness = new Matrix(2*nn, 2*nn);
		for (int i = 0; i < ne; i++) {
			Matrix vertex = calc.getVertexOfElement(coordinates, elements, i);
			Matrix k = calc.computeElementStiffness(Ces, vertex, gausspoint, thickness);

			double[] nd = elements.getArray()[i];

			Matrix index = calc.getIndexMatrix(nd, nnel, ndof);
			calc.assembleStiffnessMatrix(nn, stiffness, k, index);
		}
//		stiffness = stiffness.times(1d/1e11); 
//		stiffness.print(2*nn, 2*nn);
	
		Matrix U = calc.computeDisplacementMatrix(bcdof, stiffness, forces);
		U = U.times(1e7);
		U.print(U.getRowDimension(), 1);
		
		Matrix stress = new Matrix(3,1);
		Matrix strain = new Matrix(3,1);
		calc.computeStressStrainNode(U, Ces, gausspoint, coordinates, elements, 9, stress, strain);
		
		Matrix temp = stress.times(1e0);
		temp.print(3, 1);
		temp = strain.times(1e11);
		temp.print(3, 1);
		
	}
}
