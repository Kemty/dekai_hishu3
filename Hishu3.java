package jishu;

import java.io.FileReader;
import java.io.StreamTokenizer;

public class Hishu3 {
	public static int Tdata=414;
	public static double data[][] = new double[Tdata+7][3];
	
	private static void load_data() {
		int i, j, n = 0;
		try {
			FileReader fr = new FileReader("CO2month.txt");
			StreamTokenizer st = new  StreamTokenizer(fr);
			while (st.nextToken() != StreamTokenizer.TT_EOF) {
				i = n/3;
				j = n%3;
				data[i][j] = st.nval;
				n++; 
			}
			//System.out.println("データ数=\t"+n/2);
			fr.close();
		} catch (Exception e) {
			System.out.println(e);
		}
	}
		
	public static void main(String args[]) {
		load_data();
		int t;
		double y[] = new double[Tdata+7];
		y[0] = 350;
		for (t = 1; t <= Tdata; t++) y[t] = data[t-1][2];
		
		int N = 13;
		double xp[][][] = new double[Tdata+7][N][1];
		double Vp[][][] = new double[Tdata+7][N][N];
		double x[][][] = new double[Tdata+7][N][1];
		double V[][][] = new double[Tdata+7][N][N];
		double K[][] = new double[N][1];
		double F[][] = {
				{ 2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
				{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
				{ 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
				{ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
				{ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
				{ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
				{ 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
				{ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
				{ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
				{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0}
		};
		double tF[][] = turned(F);
		double d;
		double s = 1, am = 1, as = 1;
		double Q[][] = {
				{am*s, 0},
				{0, as*s}
		};
		double tG[][] = {
				{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
				{ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
		};
		double G[][] = turned(tG);
		
		for(int i = 0; i < N; i++) {
			x[0][i][0] = y[0];
			for(int j = 0; j < N; j++) V[0][i][j] = i==j ? 100: 0;
		}
		
		for (t = 1; t <= 396; t++) {
			// 一期先予測
			double[][] tmp = mxm(F, x[t-1]);
			for(int i = 0; i < N; i++) xp[t][i][0] = tmp[i][0];
			tmp = mpm(mxm(F, mxm(V[t-1], tF)), mxm(G, mxm(Q, tG)));
			for (int i = 0; i < N; i++) System.arraycopy(tmp[i], 0, Vp[t][i], 0, N);
			
			if(y[t] == -1) {
				t++;
				tmp = mxm(F, xp[t-1]);
				for(int i = 0; i < N; i++) xp[t][i][0] = tmp[i][0];
				tmp = mpm(mxm(F, mxm(Vp[t-1], tF)), mxm(G, mxm(Q, tG)));
				for(int i = 0; i < N; i++) System.arraycopy(tmp[i], 0, Vp[t][i], 0, N);
				d = Vp[t][0][0]+s;
			}
			// フィルタ
			d = Vp[t][0][0]+Vp[t][0][2]+Vp[t][2][0]+Vp[t][2][2]+s;
			for(int i = 0; i < N; i++) K[i][0] = (Vp[t][i][0]+Vp[t][i][2])/d;
			tmp = mpm(xp[t], sxm(y[t]-xp[t][0][0]-xp[t][2][0], K));
			for(int i = 0; i < N; i++) x[t][i][0] = tmp[i][0];
			tmp = mpm(Vp[t], sxm(-d, mxm(K, turned(K))));
			for(int i = 0; i < N; i++) System.arraycopy(tmp[i], 0, V[t][i], 0, N);
			//System.out.println(t+"\t"+data[t-1][0]+"\t"+data[t-1][1]+"\t"+y[t]+"\t"+(x[t][0][0]+x[t][2][0])+"\t"+Math.sqrt(d));
		}
		
		for (t = 397; t <= 420; t++) {
			// 一期先予測
			double[][] tmp = mxm(F, xp[t-1]);
			for(int i = 0; i < N; i++) xp[t][i][0] = tmp[i][0];
			tmp = mpm(mxm(F, mxm(Vp[t-1], tF)), mxm(G, mxm(Q, tG)));
			for(int i = 0; i < N; i++) System.arraycopy(tmp[i], 0, Vp[t][i], 0, N);
			d = Vp[t][0][0]+s;
			System.out.println((t-397)+"\t"+data[t-1][0]+"\t"+data[t-1][1]+"\t"+y[t]+"\t"+(xp[t][0][0]+xp[t][2][0])+"\t"+Math.sqrt(d));
			//System.out.println((t-397)+"\t"+(xp[t][0][0]+xp[t][2][0]));
			//if(t-397 < 18) System.out.println((t-397)+"\t"+y[t]);
		}		

	}
	
	public static double[][] mxm(double[][] matrix1, double[][] matrix2) {
		int n1 = matrix1.length;
		int m1 = matrix1[0].length;
		int m2 = matrix2[0].length;
		double newmatrix[][] = new double[n1][m2];
		
		for(int i = 0;i < n1; i++) {
			for(int j = 0; j < m2; j++) {
				double matIJ = 0; //i,j成分を求める
				for(int k = 0; k < m1; k++) {
					matIJ += matrix1[i][k]*matrix2[k][j];
				}
				newmatrix[i][j] = matIJ;
			}
		}
		return newmatrix;
	}
	
	public static double[][] mpm(double[][] mat1, double[][] mat2) {
		int sizeh = mat1.length, sizew = mat1[0].length;
		double[][] ret = new double[sizeh][sizew];
		for(int i = 0; i < sizeh; i++) {
			for(int j = 0; j < sizew; j++) ret[i][j] = mat1[i][j]+mat2[i][j];
		}
		return ret;
	}
	
	public static double[][] sxm(double sch, double[][] mat) {
		int sizeh = mat.length, sizew = mat[0].length;
		double[][] ret = new double[sizeh][sizew];
		for(int i = 0; i < sizeh; i++) {
			for(int j = 0; j < sizew; j++) ret[i][j] = sch*mat[i][j];
		}
		return ret;
	}
	
	public static double tr(double[][] matrix) {
		int n = matrix.length;
		double s = 0;
		for(int i = 0; i < n; i++) {
			s += matrix[i][i];
		}
		return s;
	}
	
	public static double[][] turned(double[][] matrix) {
		int sizeh = matrix.length, sizew = matrix[0].length;
		double[][] ret = new double[sizew][sizeh];
		for(int i = 0; i < sizeh; i++) {
			for(int j = 0; j < sizew; j++) ret[j][i] = matrix[i][j];
		}
		return ret;
	}
	
	public static void printm(double[][] matrix) {
		int n = matrix.length;
		int m = matrix[0].length;
		for(int i = 0; i < n; i++) {
			for(int j = 0; j < m; j++) {
				System.out.print(matrix[i][j]);
				if(j < m-1) System.out.print("\t");
			}
			System.out.println();
		}
	}
}

