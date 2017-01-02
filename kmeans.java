import java.io.*;
import java.util.*;
import java.util.concurrent.ThreadLocalRandom;

public class kmeans {
	public static int[] assigned_clusters;
	public static void main(String[] args)
	{
		// Specify the number of clusters
		int no_of_clusters = 5;
		//Specify the number of iterations; if you set itr=0 then kmeans algo will run until completion without any restrictions on number of iterations.
		int itr = 0;
		String file_name = "iyer.txt";
		String pca = "PCA_Kmeans_";
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

			kmeans_clustering(no_of_clusters, rows, cols, gene_cluster,itr);
			
			try{
			    BufferedWriter bfw = new BufferedWriter(new FileWriter(pca+file_name));
				for(int i=0;i<rows;i++)
				{
					for(int j=0;j<cols;j++)
					{
						if(j == 1)
							bfw.write(Double.toString(assigned_clusters[i]));
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
	
	public static void kmeans_clustering(int no_of_clusters, int rows, int cols, double[][] gene_cluster, int itr)
	{
	double clusters[][] = new double[no_of_clusters][cols-2];
	Random r = new Random();
	int counter = 0;
	int flag = 0;
	int it=0;
	double[] prev_total = new double[no_of_clusters];
	int[] prev_no = new int[no_of_clusters];
	assigned_clusters = new int[rows];
	double min = Double.MAX_VALUE;
	double distance = 0.0;	
	// Random Initialization of Centroids.
	/*for(int i=0; i<no_of_clusters;i++)
	{
		int random_cluster_value = ThreadLocalRandom.current().nextInt(0, rows-1);
		System.out.println("Random:" + random_cluster_value);
		
		for(int j=2;j<cols;j++)
		{
			clusters[i][j-2] = gene_cluster[random_cluster_value][j];
		}
	}
	*/
	
	//Hardcoded Initial Clusters
		int[] random_cluster_value = {0,96,192,288,385};
		for(int i=0; i<no_of_clusters;i++)
		{
			//int random_cluster_value = ThreadLocalRandom.current().nextInt(0, rows-1);
			System.out.println("Random:" + random_cluster_value[i]);
			System.out.println("GENER:" + gene_cluster[random_cluster_value[i]][0]);
			
			for(int j=2;j<cols;j++)
			{
				clusters[i][j-2] = gene_cluster[random_cluster_value[i]][j];
			}
		}
    
		
	while((flag != no_of_clusters) && ((it < itr)|| (itr == 0)))
	{
    int m=0;
	while(m < rows)
	{
		if(counter <= 0)
     		min = Double.MAX_VALUE;
		else 
		{
		distance = 0.0;	
		int n = 2;
		while(n<cols)
		{
			int index = assigned_clusters[m]; 
			index = index - 1;
			distance = distance + compute_euclidean_distance(gene_cluster[m][n],clusters[index][n-2]);
			n++;
		}
		distance = Math.sqrt(distance);
		min = distance;
		}
		
		int p = 0;
		while(p < no_of_clusters)
		{
			int q = 2;
			distance = 0;
			while(q < cols)
			{
				distance = distance + compute_euclidean_distance(gene_cluster[m][q], clusters[p][q-2]);
				q++;
			}
			distance = Math.sqrt(distance);
			if(distance < min)
			{
				assigned_clusters[m] = p + 1;
				min = distance;
			}
			p++;
		}
		m++;
	}
	
	int i = 0,no=0;
	flag = 0;
	double total=0.0;
	double temp_total=0.0;
	int temp_no=0;
	System.out.println("Iteration " + (counter+1));
	while(i < no_of_clusters)
	{
		ArrayList<Double> gene_classification = new ArrayList<Double>();
		int j = 2;
		while(j < cols)
		{
			int k = 0;
			total = 0.0;
			no = 0;	
			while(k < rows)
			{
				if((assigned_clusters[k]-1) == i)
				{
					no++;
					total = total + gene_cluster[k][j];
					gene_classification.add(gene_cluster[k][0]);
				}
				k++;
			}
			if(no!=0){
				temp_total = total;
			    temp_no = no;
			    clusters[i][j-2] = total/no;
			    
			}
			j++;
		}
		
		if((prev_total[i] == total) && (prev_no[i] == no))
		{
			flag = flag + 1;
		}
		prev_total[i] = temp_total;
		prev_no[i] = temp_no;
		System.out.println("Genes Belonging to Cluster " + (i+1));
		for(double d : gene_classification)
		{
		System.out.print(d + " ");
		}
		System.out.println("\nTotal Number of Genes in Cluster " + (i+1) + " : " + no);
		System.out.println();
		i++;
	}
	counter = counter + 1;
	it++;
	System.out.println();
	}
	
	// Construct ground truth and clustering incidence matrix inorder to calculate Jaccard Coeeficient
	int[][] ground_truth_incidence_matrix = new int[rows][rows];
	int[][] clustering_incidence_matrix = new int[rows][rows];
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
	
	i = 0;
	while(i<rows)
	{
		int j = 0;
		while(j<rows)
		{
			if(assigned_clusters[i] == assigned_clusters[j])
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
	//System.out.println("M11:" + m11);
	//System.out.println("M10:" + m10);
	//System.out.println("M01:" + m01);
	//System.out.println("M00:" + m00);
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
