import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.conf.Configured;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;
import org.apache.hadoop.util.Tool;
import org.apache.hadoop.util.ToolRunner;

/**
 * @author Vishwaksen *
 */
class KMeansMapper extends Mapper<Object, Text, Text, Text> {

	public void map(Object key, Text eachGenevalue, Context context) throws IOException, InterruptedException {
		// find the distance and assign it to a particular centroid
		// The key is basically the cluster Number which ranges from 1 - 5
		// The value is the gene it self
		context.write(getCluster(eachGenevalue, context), eachGenevalue);
	}

	private Text getCluster(Text value, Context context) {
		String[] currentGeneExpValues = value.toString().split("\\s");
		ArrayList<String[]> centroids = new ArrayList<>();
		int centroidCounter = 1;
		while (context.getConfiguration().get(String.valueOf(centroidCounter)) != null) {
			centroids.add(context.getConfiguration().get(String.valueOf(centroidCounter)).split("\\s"));
			centroidCounter++;
		}
		centroidCounter--; // since we have exceeded by 1
		double[] centrDistances = new double[centroidCounter];
		for (int i = 2; i < currentGeneExpValues.length; i++) {
			for (int j = 0; j < centroidCounter; j++) {
				centrDistances[j] += Math
						.pow(Float.valueOf(currentGeneExpValues[i]) - Float.valueOf(centroids.get(j)[i]), 2);
			}
		}
		// finding min index -- min index + 1 will be the cluster to which the
		// current gene belongs
		int minIndex = 0;
		for (int i = 0; i < centroidCounter; i++) {
			centrDistances[i] = Math.sqrt(centrDistances[i]);
			if (centrDistances[minIndex] > centrDistances[i]) {
				minIndex = i;
			}
		}
		return new Text(String.valueOf(minIndex + 1));// return the cluster
														// number (1 - 5)
	}
}

class KMeansReducer extends Reducer<Text, Text, Text, Text> {

	public void reduce(Text key, Iterable<Text> allGenesForThisCentroid, Context context)
			throws IOException, InterruptedException {
		ArrayList<String> allGenes = new ArrayList<>();
		StringBuilder newCentroid = new StringBuilder();
		double[] newCentroidValues = new double[20]; // hard coded 20, should
														// have
		// been 18
		int totalGenesForThisCentroid = 0;
		int geneExpLength = 0;
		for (Text eachGene : allGenesForThisCentroid) {
			totalGenesForThisCentroid++;
			String[] currentGeneExpValues = eachGene.toString().split("\\s");
			allGenes.add(currentGeneExpValues[0]);
			geneExpLength = currentGeneExpValues.length;
			for (int i = 2; i < currentGeneExpValues.length; i++) {
				newCentroidValues[i] += Float.valueOf(currentGeneExpValues[i]);
			}
		}
		// finding the new mean for each expression
		for (int i = 2; i < geneExpLength; i++) {
			newCentroidValues[i] /= totalGenesForThisCentroid;
			newCentroid.append(" " + newCentroidValues[i]);
		}
		//// key corresponds to the cluster number which ranges from 1 - 5
		// The below line does not work as the job driver is not able to get
		//// this new set data
		// context.getConfiguration().set(key.toString(),
		// newCentroid.toString());
		context.write(key,
				new Text(newCentroid.toString().trim() + "\n" + Arrays.toString(allGenes.toArray(new String[0]))));
		System.out.println("Total Genes in Cluster " + key.toString() + " : " + totalGenesForThisCentroid);
	}
}

public class KMeansMR extends Configured implements Tool {

	private static final int maxIteration = 0;
	// 2,15 or 3,2 or try 4,5
	private static final int noOfClusters = 5;
	// "new_dataset_2.txt" //cho.txt //iyer.txt
	private static final String inputFile = "iyer.txt";
	private static final Map<String, String> masterGeneDataMap = new HashMap<>();
	private static final Map<Integer, Integer> groundTruthGeneMapping = new HashMap<>();
	private static final String[] geneIdForInitialCentroids = new String[] { "1", "25", "50", "75", "100" };
	private static int totalIterations = 0;

	public int run(String[] args) throws Exception {

		// Randomly selecting the first five genes as the centroids from the
		// input file;
		BufferedReader br = new BufferedReader(new FileReader(new File(args[0])));

		String eachLine = "";
		while ((eachLine = br.readLine()) != null) {
			masterGeneDataMap.put(eachLine.split("\\s")[0], eachLine);
			groundTruthGeneMapping.put(Integer.valueOf(eachLine.split("\\s")[0]),
					Integer.valueOf(eachLine.split("\\s")[1]));
		}

		Configuration customClusterConfig = new Configuration();
		String[] previousClusterInfo = new String[noOfClusters];
		int clusterIndex = 1;
		for (String geneId : geneIdForInitialCentroids) {
			customClusterConfig.set(String.valueOf(clusterIndex), masterGeneDataMap.get(geneId));
			previousClusterInfo[clusterIndex - 1] = masterGeneDataMap.get(geneId);
			clusterIndex++;
		}

		boolean shouldIterate = false;
		boolean success;
		int iteration = 1;
		String[] newClusterInfo = new String[noOfClusters];
		do {
			shouldIterate = false;
			// updating configuration values from the generated output
			if (iteration > 1) {
				// update the previous cluster info and empty the
				processMROutput(newClusterInfo, previousClusterInfo);

			}

			System.out.println("Iteration: " + iteration);
			Job job = Job.getInstance(customClusterConfig, "K MEans");
			job.setJarByClass(KMeansMR.class);

			FileInputFormat.addInputPath(job, new Path(args[0]));
			FileOutputFormat.setOutputPath(job, new Path(args[1] + iteration));

			job.setMapperClass(KMeansMapper.class);
			job.setReducerClass(KMeansReducer.class);

			job.setOutputKeyClass(Text.class);
			job.setOutputValueClass(Text.class);

			success = job.waitForCompletion(true);

			// this will be used by MR as starting clusters for each
			// iteration except for the first one
			updateConfiguration(newClusterInfo, customClusterConfig, iteration);

			// check if previous and new values match
			for (int m = 0; m < noOfClusters; m++) {
				if (previousClusterInfo[m] == null && newClusterInfo[m] == null) {
					continue;
				}
				if (previousClusterInfo[m] == null || newClusterInfo[m] == null
						|| !previousClusterInfo[m].equals(newClusterInfo[m])) {
					shouldIterate = true;
					break;
				}
			}
			totalIterations = iteration++;
		} while ((maxIteration == 0 && shouldIterate) || iteration < maxIteration);
		br.close();
		return success ? 0 : 1;
	}

	private void updateConfiguration(String[] newClusterInfo, Configuration customClusterConfig, int iteration)
			throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File("output" + iteration + "/part-r-00000"))); // hard
		// coded
		String eachGeneratedCluster = "";
		int counter = 0;
		int lineCounter = 0;
		while ((eachGeneratedCluster = br.readLine()) != null) {
			if (lineCounter % 2 == 0) {
				// we divide by 2 as every alternate line has the cluster info
				StringBuilder newCentroid = new StringBuilder();
				String[] splits = eachGeneratedCluster.split("\\s");
				// the first value is the key
				for (int i = 1; i < splits.length; i++) {
					newCentroid.append(splits[i] + " ");
				}
				newClusterInfo[counter++] = newCentroid.toString().trim();
				customClusterConfig.set(String.valueOf(counter), "blank blank " + newCentroid.toString().trim());
			}
			lineCounter++;
		}
		br.close();
	}

	private void processMROutput(String[] newClusterInfo, String[] previousClusterInfo) {
		// copying newClusterInfo to previous
		for (int i = 0; i < newClusterInfo.length; i++) {
			previousClusterInfo[i] = newClusterInfo[i];
		}

	}

	/**
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		long startTime = System.currentTimeMillis();
		if (args == null || args.length == 0) {
			args = new String[] { inputFile, "output" };
		}
		// deleting the output folder if already present
		try {
			int outputCounter = 1;
			while (outputCounter <= 100) {
				File outputFolder = new File("output" + outputCounter);
				// emptying the outputFolder
				File[] files = outputFolder.listFiles();
				for (File eachOutputFile : files) {
					// directly deleting the file as we know it is a file and
					// not a directory
					eachOutputFile.delete();
				}
				outputFolder.delete();
				outputCounter++;
			}
		} catch (Exception e) {
			//e.printStackTrace();
		}
		// call the MapReduce method
		KMeansMR kmeansmr = new KMeansMR();
		ToolRunner.run(kmeansmr, args);
		System.out.println("MR took: " + (System.currentTimeMillis() - startTime) / (double) 1000 + " secs");

		// calculating the Rand & Jaccard Coefficient
		Map<Integer, Integer> clusterGeneMapping = new HashMap<>();
		populateGeneMapping(clusterGeneMapping);
		int[][] groundTruthIncidenceMatrix = new int[masterGeneDataMap.size()][masterGeneDataMap.size()];
		int[][] clusterIncidenceMatrix = new int[masterGeneDataMap.size()][masterGeneDataMap.size()];
		for (int m = 0; m < groundTruthIncidenceMatrix.length; m++) {
			for (int n = 0; n < groundTruthIncidenceMatrix.length; n++) {
				if (groundTruthGeneMapping.get(m) == groundTruthGeneMapping.get(n)) {
					groundTruthIncidenceMatrix[m][n] = 1;
				} else {
					groundTruthIncidenceMatrix[m][n] = 0;
				}
				if (clusterGeneMapping.get(m) == clusterGeneMapping.get(n)) {
					clusterIncidenceMatrix[m][n] = 1;
				} else {
					clusterIncidenceMatrix[m][n] = 0;
				}
			}
		}

		// Computing Jaccard Coefficient based on ground truth and clustering
		// incidence matrix
		int m11 = 0, m10 = 0, m01 = 0, m00 = 0;
		for (int p = 0; p < masterGeneDataMap.size(); p++) {
			for (int q = 0; q < masterGeneDataMap.size(); q++) {
				if ((clusterIncidenceMatrix[p][q] == 1) && (groundTruthIncidenceMatrix[p][q] == 1))
					m11 = m11 + 1;
				if ((clusterIncidenceMatrix[p][q] == 1) && (groundTruthIncidenceMatrix[p][q] == 0))
					m10 = m10 + 1;
				if ((clusterIncidenceMatrix[p][q] == 0) && (groundTruthIncidenceMatrix[p][q] == 1))
					m01 = m01 + 1;
				if ((clusterIncidenceMatrix[p][q] == 0) && (groundTruthIncidenceMatrix[p][q] == 0))
					m00 = m00 + 1;
			}
		}
		int denominator = m11 + m10 + m01;
		double jaccard_coefficient = ((double) m11 / (double) denominator);
		double rand_index = (((double) m11 + (double) m00) / ((double) denominator + m00));

		System.out.println("Map Reduce K-Means Jaccard Coefficient : " + jaccard_coefficient);
		System.out.println("Map Reduce K-Means Rand Index : " + rand_index);
	}

	private static void populateGeneMapping(Map<Integer, Integer> geneMapping) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File("output" + totalIterations + "/part-r-00000")));
		int lineCounter = 1;
		String clusterGenes = "";
		int clusterId = 1;
		while ((clusterGenes = br.readLine()) != null) {
			if (lineCounter % 2 == 0) {
				clusterGenes = clusterGenes.substring(1, clusterGenes.length() - 1);
				for (String eachGeneId : clusterGenes.split(",")) {
					geneMapping.put(Integer.valueOf(eachGeneId.trim()), clusterId);
				}
				clusterId++;
			}
			lineCounter++;
		}
		br.close();
	}
}
