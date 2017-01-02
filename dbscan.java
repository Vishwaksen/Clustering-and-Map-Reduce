import java.io.*;
import java.util.*;

public class dbscan {
	public static ArrayList<ArrayList<Integer>> cluster = new ArrayList<ArrayList<Integer>>();
	public static int[][] resultant_matrix;
	public static void main(String[] args)
	{
		double eps = 2.0;
		int min_pts = 1;
		String file_name = "iyer.txt";
		String pca = "PCA_Dbscan_";
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
	
				dbscan_clustering(eps,min_pts, rows, cols, gene_cluster);
				
				try{
				    BufferedWriter bfw = new BufferedWriter(new FileWriter(pca+file_name));
					for(int i=0;i<rows;i++)
					{
						for(int j=0;j<cols;j++)
						{
							if(j == 1)
								bfw.write(Double.toString(resultant_matrix[i][1]));
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
	
	public static void dbscan_clustering(double eps, int min_pts, int rows, int cols, double[][] gene_cluster)
	{
	int[] visited = new int[rows];
	int[] noise = new int[rows];
	for(int i=0;i<rows;i++)
	{
		visited[i] = 0;
		noise[i] = 0;
	}
	
	for(int i=0;i<rows;i++)
	{
		//ArrayList<Integer> neighbouring_pts = new ArrayList<Integer>();
		if(visited[i] == 0)
		{
			visited[i] = 1;
			ArrayList<Integer> neighbouring_pts = find_regional_points(i,eps,gene_cluster,rows,cols);
			if(neighbouring_pts.size() < min_pts)
				noise[i] = 1;
			else
			{
				ArrayList<Integer> temp = new ArrayList<Integer>();
				temp = expand_cluster(i,neighbouring_pts,temp,eps,min_pts,gene_cluster,rows,cols,visited);
				cluster.add(temp);
			}
		}
	}
	
    // Construct ground truth and clustering incidence matrix inorder to calculate Jaccard Coeeficient
	int[][] ground_truth_incidence_matrix = new int[rows][rows];
	int[][] clustering_incidence_matrix = new int[rows][rows];
    resultant_matrix  = new int[rows][2];
	int i = 0;
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
	HashMap<Integer,Integer> hmap = new HashMap<Integer,Integer>();
	int counter = 1;
	for(ArrayList<Integer> c : cluster)
	{
		for(Integer n : c)
		{
			hmap.put(n,counter);
		}
		counter++;
	}
	int m=0;
	while(m<rows)
	{
		if(hmap.containsKey(m+1))
		{
			int temp = hmap.get(m+1);
			resultant_matrix[m][1] = temp;
		}
		else
			resultant_matrix[m][1] = 0;
		m++;
	}
	
	i = 0;
	while(i<rows)
	{
		int j = 0;
		while(j<rows)
		{
			if(resultant_matrix[i][1] == resultant_matrix[j][1])
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
		double result;
		result = Math.pow(x - y, 2);
		return result;
	}

	public static ArrayList<Integer> find_regional_points(int pt, double eps, double[][] gene_cluster, int rows, int cols)
	{
		double distance;
		ArrayList<Integer> neighbouring_pts = new ArrayList<Integer>();
			for(int j=0;j<rows;j++)
			{
				distance = 0.0;
				for(int k=2;k<cols;k++)
				{
					distance = distance + compute_euclidean_distance(gene_cluster[j][k],gene_cluster[pt][k]);
				}
				distance = Math.sqrt(distance);
				if(distance <= eps)
				{
					neighbouring_pts.add(j);
				}
			}
		return neighbouring_pts;	
	}
	
	public static ArrayList<Integer> expand_cluster(int pt,ArrayList<Integer> neighbouring_pts,ArrayList<Integer>temp, double eps, int min_pts, double[][] gene_cluster, int rows, int cols, int[] visited)
	{
		int flag1, flag2;
		temp.add(pt);
		for(int k=0;k<neighbouring_pts.size();k++)
		{
			int i = neighbouring_pts.get(k);
			if(visited[i] == 0)
			{
				visited[i] = 1;
				System.out.println("I:" + i);
				ArrayList<Integer> temp_neighbouring_pts = find_regional_points(i,eps,gene_cluster,rows,cols);
				if(temp_neighbouring_pts.size() >= min_pts)
					neighbouring_pts = merge_neighbouring_pts(neighbouring_pts, temp_neighbouring_pts);
			}
			flag1 = 0;
			flag2 = 0;
			if(!temp.contains(i))
				flag1 = 1;
			for(ArrayList<Integer> j : cluster)
			{
				if(!j.contains(i))
					flag2 = 1;
			}
			if(flag1 == 1 && flag2 == 1)
				temp.add(i);
			
		}
		return temp;
	}
	
	public static ArrayList<Integer> merge_neighbouring_pts(ArrayList<Integer> neighbouring_pts, ArrayList<Integer> temp_neighbouring_pts)
	{	
		for(int i: temp_neighbouring_pts)
		{
			if (!neighbouring_pts.contains(i))
			{
			   neighbouring_pts.add(i);
			}
		}
	  return neighbouring_pts;
	}
	
}
