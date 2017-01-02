import java.io.*;
import java.util.*;

public class heirarchical {
	public static int[] minimum_similarity;
	public static void main(String[] args)
	{
		// set no_of_clusters = 1 in order to run hierarchical clustering until completion. Or else set it depending upon number of clusters you need in the final result.
		int no_of_clusters = 1;
		String file_name = "iyer.txt";
		String[] s = file_name.split("\\.");
		String pca = "PCA_Heirarchical_";
		String st = pca + s[0] + ".csv";
		String[] splitter;
		int rows=0,cols=0;
		double[][] gene_cluster;
		try{
			BufferedReader buff_rdr = new BufferedReader(new FileReader(new File(file_name)));
			String line = buff_rdr.readLine();
			while(line != null)
			{
				splitter = line.split("\t");
				cols = splitter.length;
				rows++;
				line = buff_rdr.readLine();
			}
			gene_cluster = new double[rows][cols];
			buff_rdr = new BufferedReader(new FileReader(new File(file_name)));
				for(int i=0;i<rows;i++)
				{
					line = buff_rdr.readLine();
					splitter = line.split("\t");
					for(int j=0;j<cols;j++)
					{
						gene_cluster[i][j] = Double.parseDouble(splitter[j]);
					}
				}
				int length = rows - (no_of_clusters-1);
				heirarchical_clustering(length, rows, cols, gene_cluster);
				
				try{
				    BufferedWriter bfw = new BufferedWriter(new FileWriter(pca+file_name));
					for(int i=0;i<rows;i++)
					{
						for(int j=0;j<cols;j++)
						{
							if(j == 1)
								bfw.write(Double.toString(minimum_similarity[i]));
							else
							    bfw.write(Double.toString(gene_cluster[i][j]));
							if(j!=cols-1)
							  bfw.write("\t");
						}
						bfw.newLine();
					}
				    bfw.close();
				}catch(IOException i)
				{
					i.printStackTrace();
				}
		}catch(IOException i)
		{
			i.printStackTrace();
		}
	}
	
	public static void heirarchical_clustering(int iterator, int rows, int cols, double[][] gene_cluster)
	{
	double[][] similarity_matrix = new double[rows][rows];
	int i=0;
	double min = Double.MAX_VALUE;
	minimum_similarity = new int[rows];
	System.out.print("Iter:" + iterator);
	while(i<rows)
	{
		int j=0;
		while(j<rows)
		{
			int k=2;
			double distance = 0.0;
			while(k<cols)
			{	
			distance = distance + compute_euclidean_distance(gene_cluster[j][k],gene_cluster[i][k]);
			k++;
		    }
			distance = Math.sqrt(distance);
			if(i == j)
				similarity_matrix[i][j] = Double.MAX_VALUE;
			else
				similarity_matrix[i][j] = distance;
			//System.out.println("[" + i + "]" + "[" + j + "]" + " : " + distance);
			j++;
		}
		i++;
	}
	
	i = 0;
	while(i<rows)
	{
		int j = 0;
		min = Double.MAX_VALUE;
		while(j<rows)
		{
			if(similarity_matrix[i][j] < min)
			{
				minimum_similarity[i] = j;
				min = similarity_matrix[i][j];
			}
			j++;
		}
		i++;
	}
	
	double minimum_value = Double.MAX_VALUE;
	int number1 = 0,number2 = minimum_similarity[0];
	int t;
	i = 0;
	while(i<iterator)
	{
	t = 0;
	for(int m=0;m<rows;m++)
	{
		int value1 = minimum_similarity[i];
		int v = minimum_similarity[t];
		if(similarity_matrix[i][value1] < similarity_matrix[t][v])
		{
			t = i;
			number1 = t;
			number2 = minimum_similarity[t];
		}
	}
	
	int p = 0;
	while(p<rows)
	{
		if(similarity_matrix[p][number2] < similarity_matrix[p][number1])
		{
			similarity_matrix[p][number1] = similarity_matrix[p][number2];
			similarity_matrix[number1][p] = similarity_matrix[p][number2];
		}
		p++;
	}
	
	int q = 0;
	while(q<rows)
	{
		similarity_matrix[q][number2] = Double.MAX_VALUE;
		similarity_matrix[number2][q] = Double.MAX_VALUE;
		q++;
	}
	
	int r = 0;
	while(r<rows)
	{
		int val1 = minimum_similarity[r];
 		if(val1 == number2)
 		{
 			minimum_similarity[r] = number1;
 		}
 		if(similarity_matrix[number1][r] < similarity_matrix[number1][number2])
		{
			minimum_similarity[number1] = r;
		}
		r++;
	}
	
	i++;
	}
	
	// Construct ground truth and clustering incidence matrix inorder to calculate Jaccard Coeeficient
	int[][] ground_truth_incidence_matrix = new int[rows][rows];
	int[][] clustering_incidence_matrix = new int[rows][rows];
	i = 0;
	while(i<rows)
	{
		int j = 0;
		while(j<rows)
		{
			if(gene_cluster[i][1] == gene_cluster[j][1])
			{
				ground_truth_incidence_matrix[i][j] = 1;
			}
			else
			{
				ground_truth_incidence_matrix[i][j] = 0;
			}
			j++;
		}
		i++;
	}
	
	i = 0;
	while(i<rows)
	{
		int j = 0;
		while(j<rows)
		{
			if(minimum_similarity[i] == minimum_similarity[j])
			{
				clustering_incidence_matrix[i][j] = 1;
			}
			else
			{
				clustering_incidence_matrix[i][j] = 0;
			}
			j++;
		}
		i++;
	}
	
	// Computing Jaccard Coefficient based on ground truth and clustering incidence matrix
	int m11=0,m10=0,m01=0,m00=0;
	for(int p=0;p<rows;p++)
	{
		for(int q=0;q<rows;q++)
		{
			if((clustering_incidence_matrix[p][q] == 1) && (ground_truth_incidence_matrix[p][q] == 1))
				m11 = m11 + 1;
			if((clustering_incidence_matrix[p][q] == 1) && (ground_truth_incidence_matrix[p][q] == 0))
				m10 = m10 + 1;
			if((clustering_incidence_matrix[p][q] == 0) && (ground_truth_incidence_matrix[p][q] == 1))
				m01 = m01 + 1;
			if((clustering_incidence_matrix[p][q] == 0) && (ground_truth_incidence_matrix[p][q] == 0))
				m00 = m00 + 1;
			
		}
	}
	System.out.println("M11:" + m11);
	System.out.println("M10:" + m10);
	System.out.println("M01:" + m01);
	System.out.println("M00:" + m00);
	int denominator = m11 + m10 + m01;
	int numerator = m11 + m00;
	int deno1 = m11 + m10 + m00 + m01;
	double jaccard_coefficient = ((double)m11 / (double) denominator);
	double rand_index = ((double)numerator / (double) deno1);
	System.out.println("K-Means Jaccard Coefficient : " + jaccard_coefficient);
	System.out.println("K-Means Rand Index : " + rand_index);
}
	
	// Computes Euclidean distance between the data points.
	public static double compute_euclidean_distance(double x, double y)
	{
		double result, temp;
		result = Math.pow(x - y, 2);
		return result;
	}

}
