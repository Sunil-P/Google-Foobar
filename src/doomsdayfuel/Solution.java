/*Please be sure to comment out all print statements including unnecessary .print() calls*/

package doomsdayfuel;

import java.util.ArrayList;
import java.util.List;

public class Solution {
		
	public static int transitionCount = 0;
	public static int absorptionCount = 0;
	
	public static class Fraction {
		public int num;
		public int den;
		
		public Fraction (int num, int den) {
			
			if (den == 0) {
				this.num = 0;
				this.den = 1;
			}
			else {
				this.num = num;
				this.den = den;
			}
		}
		
		public float absoluteValue() {
			return (float)num/(float)den;
		}
		
		public boolean isWhole() {
			
			if (((float)num/(float)den)%1 == 0) {
				return true;
			}
			
			return false;
		}
		
		public int gcd(int a, int b){
			   if (b==0) return a;
			   return gcd(b,a%b);
		}
		
		public void lowest() {
			Fraction f = lowest(num, den);
			this.num = f.num;
			this.den = f.den;
			
			this.fixSign();
		}
		
		public Fraction lowest(int num, int den)
		{
		    int common_factor = gcd(num, den);
		 
		    den = den / common_factor;
		    num = num / common_factor;
		    
		    Fraction result = new Fraction(num, den);
		    return result;
		}
		
		public Fraction add(Fraction b) {
			
			if (this.num == 0) {
				return b;
			} else if (b.num == 0) {
				return this;
			}
			
			this.fixSign();
			b.fixSign();
			
			int den3 = gcd(this.den, b.den);
			den3 = (this.den * b.den)/den3;
			int num3 = (this.num) * (den3/this.den) + (b.num) * (den3/b.den);
			
			Fraction result = new Fraction(num3, den3);
			result.lowest();
			
			return result;
		}
		
		public Fraction substract(Fraction b) {
			
			if (this.num == 0) {
				return new Fraction(-b.num, b.den);
			} else if (b.num == 0) {
				return this;
			}
			
			b.fixSign();
			this.fixSign();
			
			int den3 = gcd(this.den, b.den);
			den3 = (this.den * b.den)/den3;
			int num3 = (this.num) * (den3/this.den) - (b.num) * (den3/b.den);
			
			if (den3 < 0 && num3 < 0) {
				den3 = -den3;
				num3 = -num3;
			}
			
			Fraction result = new Fraction(num3, den3);
			result.lowest();
			
			return result;
		}
		
		public Fraction multiply(Fraction b) {
			
			if (this.num == 0 || b.num == 0) {
				return new Fraction(0,1);
			}
			
			this.fixSign();
			b.fixSign();
			
			int num3 = this.num * b.num;
			int den3 = this.den * b.den;
			
			Fraction result = new Fraction(num3, den3);
			result.lowest();
			
			return result;
		}
		
		public Fraction multiply(int b) {
			if (this.num == 0 || b == 0) {
				return new Fraction(0,1);
			}
			
			this.fixSign();
			
			int num3 = (int) (this.num * b);
			int den3 = this.den;
			
			Fraction result = new Fraction(num3, den3);
			result.lowest();
			
			return result;
		}
		
		public Fraction divide(Fraction b) {
			if (this.num == 0 || b.num == 0) {
				return new Fraction(0,1);
			}
			
			this.fixSign();
			b.fixSign();
			
			int num3 = (int) this.num * b.den;
			int den3 = this.den * b.num;
			
			Fraction result = new Fraction(num3, den3);
			result.lowest();
			
			return result;
		}
		
		public Fraction divide(int b) {
			if (this.num == 0 || b == 0) {
				return new Fraction(0,1);
			}
			
			this.fixSign();
			
			int num3 = (int) this.num;
			int den3 = this.den * b;
			
			Fraction result = new Fraction(num3, den3);
			result.lowest();
			
			return result;
		}
		
		public Fraction divideTo(int b) {
			if (this.num == 0 || b == 0) {
				return new Fraction(0,1);
			}
			
			this.fixSign();
			
			int num3 = b * this.den;
			int den3 = this.num;
			
			Fraction result = new Fraction(num3, den3);
			result.lowest();
			
			return result;
		}
		
		public void print() {
			if (this.isWhole()) {
				System.out.print(num/den);
			} else {
			System.out.print(num+"/"+den);
			}
		}
		
		public void fixSign() {
			if (this.num < 0 && this.den < 0) {
				this.num = -this.num;
				this.den = -this.den;
			} else if (this.den < 0 && this.num >0) {
				this.num = -this.num;
				this.den = -this.den;
			}
		}
	}
	
	public static class Matrix {
		public Fraction[][] matrix;
		
		public Matrix(int rows, int cols) {
			matrix = new Fraction[rows][cols];
		}
		
		public Matrix(int rows, int cols, Fraction defaultValue) {
			matrix = new Fraction[rows][cols];
			for(int i=0;i<rows;i++) {
				for(int j=0;j<cols;j++) {
					matrix[i][j] = defaultValue;
				}
			}			
		}
		
		public static Matrix formIdentityMatrix(int size) {
			
			Matrix m = new Matrix(size,size);
							
				for(int i=0;i<size;i++) {
					for (int j=0;j<size;j++) {
						if (i == j) {
							m.matrix[i][j] = new Fraction(1,1);
						} else {
							m.matrix[i][j] = new Fraction(0,1);
						}
					}
				}
		
			return m;
		}
		
		public static Matrix add (Matrix mat1, Matrix mat2) {
			
			if (mat1.matrix.length == 0 || mat2.matrix.length == 0) {
				//System.out.println("Argument matrix has zero size");
				return null;
			}
			
			if (mat1.matrix.length != mat2.matrix.length || mat1.matrix[0].length != mat2.matrix[0].length) {
				//System.out.println("Invalid size matrix multiplication operation");
				return null;
			}
			
			Matrix result = new Matrix(mat1.matrix.length,mat1.matrix[0].length);
			
			for (int i=0;i<mat1.matrix.length; i++) {
				for (int j=0; j<mat1.matrix[0].length; j++) {
					result.matrix[i][j] = mat1.matrix[i][j].add(mat2.matrix[i][j]);
				}
			}

			return result;
		}
		
		public void add(Matrix mat2) {
			this.matrix = add(this, mat2).matrix;
		}
		
		public static Matrix substract (Matrix mat1, Matrix mat2) {
			
			if (mat1.matrix.length == 0 || mat2.matrix.length == 0) {
				//System.out.println("Argument matrix has zero size");
				return null;
			}
			
			if (mat1.matrix.length != mat2.matrix.length || mat1.matrix[0].length != mat2.matrix[0].length) {
				//System.out.println("Invalid size matrix multiplication operation");
				return null;
			}
			
			Matrix result = new Matrix(mat1.matrix.length,mat1.matrix[0].length);
			
			for (int i=0;i<mat1.matrix.length; i++) {
				for (int j=0; j<mat1.matrix[0].length; j++) {
					result.matrix[i][j] = mat1.matrix[i][j].substract(mat2.matrix[i][j]);
				}
			}
	
			return result;
		}
		
		public void substract(Matrix mat2) {
			this.matrix = substract(this, mat2).matrix;
		}
		
		public static Matrix dotProduct (Matrix mat1, Matrix mat2) {
			
			int row1 = mat1.matrix.length;
			int col1 = mat1.matrix[0].length;
			
			int row2 = mat2.matrix.length;
			int col2 = mat2.matrix[0].length;
			
			if(col1 != row2){    
	            //System.out.println("Matrices cannot be multiplied"); 
	            return null;
	        }    
			
			Matrix result = new Matrix(row1,col2);
			
			for(int i = 0; i < row1; i++){    
                for(int j = 0; j < col2; j++){
                	result.matrix[i][j] = new Fraction(0,1);
                    for(int k = 0; k < row2; k++){    
                    	result.matrix[i][j] = result.matrix[i][j].add(mat1.matrix[i][k].multiply(mat2.matrix[k][j]));     
                    }    
                }    
            }    

			return result;
		}
		
		public void dotProduct(Matrix mat2) {
			matrix = dotProduct(this, mat2).matrix;
		}
		
		public static Matrix createSubMatrix(Matrix m, int excluding_row, int excluding_col) {
		    Matrix mat = new Matrix(m.matrix.length-1,m.matrix[0].length-1);
		    int r = -1;
		    for (int i=0;i<m.matrix.length;i++) {
		        if (i==excluding_row)
		            continue;
		            r++;
		            int c = -1;
		        for (int j=0;j<m.matrix[0].length;j++) {
		            if (j==excluding_col)
		                continue;
		            mat.matrix[r][++c] = m.matrix[i][j];
		        }
		    }
		    return mat;
		} 
		
		public static Fraction determinant(Matrix m) {
		    if (m.matrix.length != m.matrix[0].length)
		    {
		    	return null;
		    }
		    
		    
		    if (m.matrix.length == 1) {
			return m.matrix[0][0];
		    }
		    if (m.matrix.length==2) {
		        return m.matrix[0][0].multiply(m.matrix[1][1]).substract(m.matrix[0][1].multiply(m.matrix[1][0]));
		    }
		    Fraction sum = new Fraction(0,1);
		    for (int i=0; i<m.matrix[0].length; i++) {
		    	
		    	Fraction d = determinant(createSubMatrix(m, 0, i));
		    	
		        sum =  sum.add(m.matrix[0][i].multiply(d.multiply(Math.changeSign(i))));
		    }
		    return sum;
		} 
		
		
		public static Matrix cofactor(Matrix m) {
		    Matrix result = new Matrix(m.matrix.length,m.matrix[0].length);
		    for (int i=0;i<m.matrix.length;i++) {
		        for (int j=0; j<m.matrix[0].length;j++) {
		        	
		        	int finalSign = Math.changeSign(i) * Math.changeSign(j);
		            result.matrix[i][j] =  
	                        determinant(createSubMatrix(m, i, j)).multiply(finalSign);
		        }
		    }
		    
		    return result;
		}
		
		public static Matrix transpose(Matrix m) {
		    Matrix result =  new Matrix(m.matrix.length,m.matrix[0].length);
		    for (int i=0;i<m.matrix.length;i++) {
		        for (int j=0;j<m.matrix[0].length;j++) {
		            result.matrix[j][i]= m.matrix[i][j];
		        }
		    }
		    return result;
		} 
		
		public static Matrix multiplyByConstant(Matrix mat, Fraction constant){
			for(int i=0;i<mat.matrix.length;i++) {
				for(int j=0;j<mat.matrix[0].length;j++) {
					mat.matrix[i][j] = mat.matrix[i][j].multiply(constant);
				}
			}
			
			return mat;
		}
		
		public void multiplyByConstant(Fraction constant) {
			this.matrix = multiplyByConstant(this, constant).matrix;
		}
		
		public static Matrix inverse(Matrix m) {
		    return (multiplyByConstant(transpose(cofactor(m)), determinant(m).divideTo(1)));
		}
		
		public void invert() {
			this.matrix = inverse(this).matrix;
		}
		
		public void print() {
			
			for (int i=0;i<matrix.length;i++) {
				for (int j=0;j<matrix[0].length;j++) {
					matrix[i][j].print();
					System.out.print("\t");
				}
				
				System.out.print("\n");
			}
		}
	}
	
	public static class Math {
		
		public static int gcd(int a, int b)
		{
		    if (b == 0)
		        return a;
		    return gcd(b, a % b);
		}
		
		public static int findlcm(int arr[], int n)
		{
		    int ans = arr[0];
		 
		    for (int i = 1; i < n; i++)
		        ans = (((arr[i] * ans)) /
		                (gcd(arr[i], ans)));
		 
		    return ans;
		}
		
		public static int changeSign(int i) {
			if (i%2 == 0) {
				return 1;
			}
			return -1;
		}
	}	
	
	public static Matrix convertToProbability(int[][] matrix){
		
		if (matrix.length == 0) {
			return null;
		}
		
		Matrix result = new Matrix(matrix.length,matrix[0].length);
		
		for (int i=0;i<matrix.length;i++) {
			int total = 0;
			for (int j=0;j<matrix[i].length;j++) {
				total = total + matrix[i][j];
			}
			
			System.out.println("Probability Row : "+i + " Total = "+total);
			
			for (int j=0;j<matrix[i].length;j++) {
				result.matrix[i][j] = new Fraction(matrix[i][j], total);
			}
		}
		
		return result;
	}
	
	public static int[][] getOrderedMatrix(int[][] matrix){
		
		if (matrix.length == 0) {
			return null;
		}
		
		List<Integer> transitionRows = new ArrayList<Integer>();
		List<Integer> absorptionRows = new ArrayList<Integer>();
		
		for (int i=0;i<matrix.length;i++) {
			int total = 0;
			for (int j=0;j<matrix[i].length;j++) {
				total = total + matrix[i][j];
			}
			
			System.out.println("Row : "+i + " total = "+total);
			if (total == 0) {
				System.out.println("Row : "+i+" is absorbing state");
				absorptionRows.add(i);
			} else {
				System.out.println("Row : "+i+" is transition state");
				transitionRows.add(i);
			}
		}
		
		transitionCount = transitionRows.size();
		absorptionCount = absorptionRows.size();
		
		int[][] resultmatrices = new int[transitionRows.size()][matrix[0].length];
		
		for (int i=0;i<transitionRows.size();i++) {
			int transitionRow = transitionRows.get(i);
			int j = 0;
			for (j=0;j<absorptionRows.size();j++) {
				int absorptionRow = absorptionRows.get(j);
				resultmatrices[i][j] = matrix[transitionRow][absorptionRow];
				System.out.println("matrix filled with absorption values ["+i+","+j+"] = "+resultmatrices[i][j]);
			}
			
			for (int k=0;k<transitionRows.size();k++,j++) {
				int l = transitionRows.get(k);
				resultmatrices[i][j] = matrix[transitionRow][l];
				System.out.println("matrix filled with transition values ["+i+","+j+"] = "+resultmatrices[i][j]);
			}
			
		}
				
		return resultmatrices;
	}
	
	
	
	public static Matrix[] getRQSplitMatrices(Matrix m){
		
		Matrix R = new Matrix(transitionCount,absorptionCount);
		Matrix Q = new Matrix(transitionCount,transitionCount);

		System.out.println("Transition Count :"+transitionCount);
		System.out.println("Absorption Count :"+transitionCount);

		for(int i=0;i<transitionCount;i++) {
			
			int j =0;
			for (j=0;j<absorptionCount;j++) {
				R.matrix[i][j] = m.matrix[i][j];
				System.out.println("R matrix i="+i+", j="+j+" = "+m.matrix[i][j].num + "/"+ m.matrix[i][j].den);
			}
			
			for(int k=0;j<m.matrix[0].length;j++,k++) {
				Q.matrix[i][k] = m.matrix[i][j]; 
				System.out.println("Q matrix i="+i+", k="+k+" = "+Q.matrix[i][k].num + "/"+ Q.matrix[i][k].den);

			}
		}
		
		System.out.println("Transition To Absorption Matrix");
		R.print();
		
		System.out.println("Transition To Transition Matrix");
		Q.print();
		
		Matrix[] result = {R,Q};
		
		return result;
	}
	


	public static int[] solution(int[][] matrix) {

		int sum = 0;
		for(int i=0;i<matrix[0].length;i++) {
			sum += matrix[0][i];
		}
		
		int size = matrix.length;
		
		if (size == 1) {
			size++;
		}
		
		
		if (sum == 0) {
			int[] result = new int[size];
			result[0] = 1;
			result[size-1] = 1;
			return result;
		}
	
		System.out.println("Printing Ordered matrix");
		int[][] transitionMatrix = getOrderedMatrix(matrix);
		for(int i=0;i<transitionMatrix.length;i++) {
			for (int j=0;j<transitionMatrix[0].length;j++) {
				System.out.print(transitionMatrix[i][j] + "\t");
			}
			System.out.print("\n");
		}
		
		
		Matrix probMatrix = convertToProbability(transitionMatrix);
		Matrix[] RQmatrix = getRQSplitMatrices(probMatrix);
		Matrix R = RQmatrix[0];
		Matrix Q = RQmatrix[1];
		
		Matrix I = Matrix.formIdentityMatrix(Q.matrix.length);
		System.out.println("Print Identity matrix");
		I.print();
		
		Matrix IQ = Matrix.substract(I, Q);
		System.out.println("I-Q Matrix :");
		IQ.print();

		Matrix F = Matrix.inverse(IQ);
		System.out.println("I-Q Inverse or F Matrix :");
		F.print();
		
		Matrix FR = Matrix.dotProduct(F, R);
		System.out.println("FR Matrix :");
		FR.print();
		
		
		List<Fraction> ansList = new ArrayList<Fraction>();
		int[] denoms = new int[FR.matrix[0].length];
		for(int i=0;i<FR.matrix[0].length;i++) {
			Fraction val =FR.matrix[0][i];
			ansList.add(val);
			denoms[i] = val.den;
			System.out.println(val.num+"/"+val.den);
		}
		
		int lcm = Math.findlcm(denoms, denoms.length);
		int[] result = new int[denoms.length+1];
		for(int i=0;i<ansList.size();i++) {
			Fraction f = ansList.get(i);
			int factor = lcm/f.den;
			result[i] = factor*f.num;
		}
		
		result[ansList.size()] = lcm;
		return result;
		
	}
	
	public static void main(String[] args) {

		int[][] matrix = {{0, 1, 0, 0, 0, 1}, {4, 0, 0, 3, 2, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};
		int[][] matrix2 = {{1, 1, 1, 1, 1},  {0, 0, 0, 0, 0}, {1, 1, 1, 1, 1}, {0, 0, 0, 0, 0}, {1, 1, 1, 1, 1}};
		
		
		int []sol = solution(matrix);
		
		for(int i=0;i<sol.length;i++) {
			System.out.print(sol[i]+" ");
		}
	}

}
