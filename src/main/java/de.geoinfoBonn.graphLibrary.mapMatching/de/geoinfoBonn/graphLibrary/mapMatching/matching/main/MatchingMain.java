package de.geoinfoBonn.graphLibrary.mapMatching.matching.main;

import java.awt.geom.Point2D;
import java.io.*;
import java.util.*;

import org.apache.commons.lang3.StringUtils;
import org.geotools.api.feature.simple.SimpleFeature;
import org.tinylog.Logger;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicLong;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import com.google.common.collect.Iterables;
import com.google.common.math.LongMath;
import de.geoinfoBonn.graphLibrary.mapMatching.core.generic.DiGraph;
import de.geoinfoBonn.graphLibrary.mapMatching.core.generic.DiGraph.DiGraphArc;
import de.geoinfoBonn.graphLibrary.mapMatching.core.generic.DiGraph.DiGraphNode;
import de.geoinfoBonn.graphLibrary.mapMatching.core.generic.DoubleWeightDataWithInfo;
import de.geoinfoBonn.graphLibrary.mapMatching.io.Road;
import de.geoinfoBonn.graphLibrary.mapMatching.io.Road.RoadInfo;
import de.geoinfoBonn.graphLibrary.mapMatching.io.RoadReader;
import de.geoinfoBonn.graphLibrary.mapMatching.matching.Matching;
import de.geoinfoBonn.graphLibrary.mapMatching.matching.Track;
import org.geotools.api.referencing.crs.CoordinateReferenceSystem;
import org.geotools.geopkg.FeatureEntry;
import org.geotools.geopkg.GeoPackage;

import static de.geoinfoBonn.graphLibrary.mapMatching.matching.Track.getIdAsLong;
import static de.geoinfoBonn.graphLibrary.mapMatching.matching.main.AbstractMain.getOptionalArg;
import static de.geoinfoBonn.graphLibrary.mapMatching.matching.main.AbstractMain.containsOptionalArg;

//@formatter:off
/**
 * Executable for map matching.
 *
 * @author Jan-Henrik Haunert
 * @author Axel Forsch (forsch@igg.uni-bonn.de)
 * @version %I%, %G%
 *
 */
//@formatter:on
public class MatchingMain {

	private static final int DEFAULT_PARTITION_SIZE = 1000;
	private static final int DEFAULT_NUMBER_OF_THREADS = 6;

	public static void main(String[] args) {

		if (containsOptionalArg(args, "-h")) {
			printHelpText();
			System.exit(0);
		}

		// PROCESS ARGUMENTS
		// Core matching parameters (radius, max candidates, candidate cost weight)
		Matching.VERBOSE = containsOptionalArg(args, "-v"); // todo: verbose/debug doesn't work
		if (getOptionalArg(args, "-r") != null)
			Matching.RADIUS = Double.parseDouble(getOptionalArg(args, "-r"));
		if (getOptionalArg(args, "-k") != null)
			Matching.MAX_CAND_N = Integer.parseInt(getOptionalArg(args, "-k"));
		if (getOptionalArg(args, "-c") != null)
			Matching.CANDIDATE_COST_WEIGHT = Double.parseDouble(getOptionalArg(args, "-c"));

		// Regarding offroad candidates
		if (getOptionalArg(args, "-off") != null)
			Matching.ADD_OFFROAD_CANDIDATE = Boolean.parseBoolean(getOptionalArg(args, "-off")); // was -o
		if (getOptionalArg(args, "-offwt") != null)
			Matching.OFF_ROAD_WEIGHT = Double.parseDouble(getOptionalArg(args, "-offwt")); // was -w

		// Regarding which outputs to write
		boolean writeMatches = containsOptionalArg(args, "-ma");
		boolean writeChunks = containsOptionalArg(args, "-ch");
		boolean writeChunkPaths = containsOptionalArg(args, "-chp");
		boolean writeGlobalPaths = containsOptionalArg(args, "-gp");
		boolean writeSegments = containsOptionalArg(args, "-s");

		// The input network attribute with the road's
		// weight (e.g., geometric length, travel time)
		// If left unspecified, the geometric length computed from coordinates is used
		// as the road's weight
		String linkDistId = getOptionalArg(args, "-nwt"); // was "-l"

		// The input network attribute unique link ID
		String linkIdName = getOptionalArg(args, "-nid"); // was "-t"

		// Regarding partitioning and multithreading
		String threadsInput = getOptionalArg(args, "-t");
		String partitionsInput = getOptionalArg(args, "-n");
		int numberOfThreads = threadsInput == null ? DEFAULT_NUMBER_OF_THREADS : Integer.parseInt(threadsInput);
		int partitionSize = partitionsInput == null ? DEFAULT_PARTITION_SIZE : Integer.parseInt(partitionsInput);

		// Deviation and distance penalty factors
		if (getOptionalArg(args, "-pdev") != null)
			Matching.DEVIATION_PENALTY_FACTOR = Double.parseDouble(getOptionalArg(args, "-pdev"));
		if (getOptionalArg(args, "-pdist") != null)
			Matching.DISTANCE_PENALTY_FACTOR = Double.parseDouble(getOptionalArg(args, "-pdist"));

		// Link weight adjustments
		String weightAdjustmentsFile = getOptionalArg(args, "-plink");
		final Map<String,Double> weightAdjustments = (weightAdjustmentsFile == null) ? null : readWeightAdjustments(weightAdjustmentsFile);

		// InfoGenerator
		Function<SimpleFeature, RoadInfo> roadInfoGenerator = feature -> {

			// RoadId
			long roadId = getIdAsLong(linkIdName, feature);

			// Weight
			double weight = 1.;
			if(weightAdjustments != null) {
				for(Map.Entry<String,Double> e : weightAdjustments.entrySet()) {

					String criteria = e.getKey();

					int eqlIdx = criteria.indexOf('=');
					String attr = criteria.substring(0, eqlIdx);
					String testVal = criteria.substring(eqlIdx+1);

					Object val = feature.getAttribute(attr);
					if(val != null) {
						if(testCriteria(val, testVal)) {
							weight *= e.getValue();
						}
					}
				}
			}

			// Results
			return new RoadInfo(roadId,weight);
		};

		// Read roads and CRS
		LinkedList<Road<RoadInfo>> roads = RoadReader.importFromGpkg(args[0],roadInfoGenerator,linkDistId);
		CoordinateReferenceSystem crs = RoadReader.readCRS(args[0]);

		// Read trajectories
		String trajectoryIdName = getOptionalArg(args, "-tid");
		List<Track> trajectories = Track.importTrajectories(args[1],trajectoryIdName);

		// Filter trajectories
		String selectedIds = getOptionalArg(args, "-select");
		String minInput = getOptionalArg(args, "-min");
		String maxInput = getOptionalArg(args, "-max");

		if(selectedIds != null) {
			Set<Long> includedIds = Arrays.stream(selectedIds.split(",")).map(Long::parseLong).collect(Collectors.toSet());
			trajectories = trajectories.stream().filter(f -> includedIds.contains(f.getId())).collect(Collectors.toList());
			if (minInput != null || maxInput != null) {
				throw new RuntimeException("For filtering IDs, shouldn't combine -select and -min/max");
			}
		} else if(minInput != null || maxInput != null) {
			long min = minInput == null ? Long.MIN_VALUE : Long.parseLong(minInput);
			long max = maxInput == null ? Long.MAX_VALUE : Long.parseLong(maxInput);
			trajectories = trajectories.stream().filter(f -> f.getId() >= min && f.getId() <= max).collect(Collectors.toList());
		}

		// Prepare trajectory data
		Iterable<List<Track>> partitions = Iterables.partition(trajectories, partitionSize);
		long numberOfPartitions = StreamSupport.stream(partitions.spliterator(), false).count();


		// STATUS UPDATE BEFORE MATCHING
		Logger.info("Number of trajectories in input file: " + trajectories.size());
		Logger.info("Program arguments:");
		// Core parameters
		Logger.info(" - search radius:           " + Matching.RADIUS);
		Logger.info(" - max. num. of candidates: " + Matching.MAX_CAND_N);
		Logger.info(" - candidate cost weight:   " + Matching.CANDIDATE_COST_WEIGHT);
		Logger.info(" - use offroad candidates?  " + Matching.ADD_OFFROAD_CANDIDATE);
		if (Matching.ADD_OFFROAD_CANDIDATE)
			Logger.info(" - offroad weight:      " + Matching.OFF_ROAD_WEIGHT);
		// Other penalty factors
		Logger.info(" - distance penalty factor: " + Matching.DISTANCE_PENALTY_FACTOR);
		Logger.info(" - deviation penalty factor: " + Matching.DEVIATION_PENALTY_FACTOR);
		// Computational parameters
		Logger.info(" - number of threads:       " + numberOfThreads);
		Logger.info(" - batch size:              " + partitionSize);
		// Inputs
		Logger.info(" - road input data:         " + args[0]);
		Logger.info(" - trajectory input data:   " + args[1]);
		// Input column names
		Logger.info(" - network ID column:       " + linkIdName);
		Logger.info(" - network weight column:   " + linkDistId);
		Logger.info(" - trajectory ID column:    " + trajectoryIdName);
		// Output details
		Logger.info(" - output Matches?      " + writeMatches);
		Logger.info(" - output Chunks?       " + writeChunks);
		Logger.info(" - output ChunkPaths?   " + writeChunkPaths);
		Logger.info(" - output GlobalPaths?  " + writeGlobalPaths);
		Logger.info(" - output Segments?     " + writeSegments);
		// Debugging / other
		if (selectedIds != null) {
			Logger.info(" - selected " + trajectories.size() + " trajectories with IDs: " + trajectories.stream().map(t -> Long.toString(t.getId())).collect(Collectors.joining(",")));
		}
		if (minInput != null || maxInput != null) {
			Logger.info(" - selected " + trajectories.size() + " trajectories between ID range " + minInput + " and " + maxInput);
		}
		Logger.info("- split trajectories into " + numberOfPartitions + " partitions containing " + partitionSize + " trajectories each.");

		// Initialise output lists
		ArrayList<Track> paths = new ArrayList<>(partitionSize);
		ArrayList<Track> matches = writeMatches ? new ArrayList<>(partitionSize) : null;
		ArrayList<Track> chunks = writeChunks ? new ArrayList<>(partitionSize) : null;
		ArrayList<Track> chunkPaths = writeChunkPaths ? new ArrayList<>(partitionSize) : null;
		ArrayList<Track> globalPaths = writeGlobalPaths ? new ArrayList<>(partitionSize) : null;
		ArrayList<Track> segments = writeSegments ? new ArrayList<>(partitionSize) : null;

		// Initialise output arrays
		TrajectoryWorker[] workers = new TrajectoryWorker[numberOfThreads];
		Thread[] threads = new Thread[numberOfThreads];

		// Specify output file (and delete if already exists)
		File outputFile = new File(args[2]);
		if(outputFile.delete()) {
			Logger.warn("File " + outputFile.getAbsolutePath() + " already exists. Deleting...");
		}

		try {

			// Crete new output file
			Logger.info("Creating and initialising new geopackage file: " + outputFile.getAbsolutePath());
			GeoPackage out = new GeoPackage(outputFile);
			out.init();

			// Create shutdown hook (to run at exit)
			Thread closeGpkgHook = new Thread(() -> {
				Logger.info("Initiating shutdown...");

				// Stop running threads
				for(Thread thread : threads) {
					thread.interrupt();
				}

				// Create spatial indices
				try {
					for(FeatureEntry entry : out.features()) {
						Logger.info("Creating spatial index for entry: " + entry.getDescription() + "...");
						long starttime = System.currentTimeMillis();
						out.createSpatialIndex(entry);
						Logger.info("Spatial index created in " + (System.currentTimeMillis() - starttime) / 1000.0 + "s.");
					}
				} catch (IOException e) {
					System.out.println(e.getMessage());
					Logger.warn("Could not create spatial indexes!");
				}

				// Close file
				Logger.info("Closing geopackage...");
				out.close();

				// Log message
				Logger.info("Wrote results to " + outputFile.getAbsolutePath());
			});
			Runtime.getRuntime().addShutdownHook(closeGpkgHook);

			int currPartition = 0;

			for(final List<Track> partition : partitions) {

				currPartition++;
				Logger.info("COMPUTING PARTITION " + currPartition + " OF " + numberOfPartitions + "...");

				// Initiate matching threads
				ConcurrentLinkedQueue<Track> trajectoriesQueue = new ConcurrentLinkedQueue<>(partition);
				AtomicLong counter = new AtomicLong();
				for(int i = 0; i < numberOfThreads; i++) {
					workers[i] = new TrajectoryWorker(trajectoriesQueue,counter,linkIdName,roads,
							writeMatches,writeChunks,writeChunkPaths,writeGlobalPaths,writeSegments);
					threads[i] = new Thread(workers[i],"Worker # " + i);
					threads[i].start();
				}

				// wait until al threads have finished
				for(Thread thread : threads) {
					try {
						thread.join();
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
				}

				// Paths
				Arrays.stream(workers).forEach(w -> paths.addAll(w.getPaths()));
				paths.sort(Comparator.comparingLong(Track::getId).thenComparingInt(Track::getSubtrack));
				Track.writeToGpkg(out, "paths", paths,false, false, crs);
				paths.clear();

				// Matches
				if(writeMatches) {
					Arrays.stream(workers).forEach(w -> matches.addAll(w.getMatches()));
					matches.sort(Comparator.comparingLong(Track::getId).thenComparingInt(Track::getSubtrack));
					Track.writeToGpkg(out,"matches", matches, true, false, crs);
					matches.clear();
				}

				// Chunks
				if(writeChunks) {
					Arrays.stream(workers).forEach(w -> chunks.addAll(w.getChunks()));
					chunks.sort(Comparator.comparingLong(Track::getId).thenComparingInt(Track::getSubtrack).thenComparingInt(Track::getSection));
					Track.writeToGpkg(out,"chunks", chunks, true, false, crs);
					chunks.clear();
				}

				// Chunk paths
				if(writeChunkPaths) {
					Arrays.stream(workers).forEach(w -> chunkPaths.addAll(w.getChunkPaths()));
					chunkPaths.sort(Comparator.comparingLong(Track::getId).thenComparingInt(Track::getSubtrack));
					Track.writeToGpkg(out,"chunkPaths", chunkPaths, false, false, crs);
					chunkPaths.clear();
				}

				// Global paths
				if(writeGlobalPaths) {
					Arrays.stream(workers).forEach(w -> globalPaths.addAll(w.getGlobalPaths()));
					globalPaths.sort(Comparator.comparingLong(Track::getId).thenComparingInt(Track::getSubtrack));
					Track.writeToGpkg(out, "globalPaths",globalPaths, false, false, crs);
					globalPaths.clear();
				}


				// Segments
				if(writeSegments) {
					Arrays.stream(workers).forEach(w -> segments.addAll(w.getSegments()));
					segments.sort(Comparator.comparingLong(Track::getId).thenComparingInt(Track::getSubtrack).thenComparingInt(Track::getSection));
					Track.writeToGpkg(out,"segments", segments, true, true, crs);
					segments.clear();
				}

				System.gc();
			}

			// Normal exit
			System.exit(0);

		} catch (IOException e) {
			throw new RuntimeException(e);
		}

	}

	private static void printHelpText() {
		System.out.println("Executable for map matching.");
		System.out.println();
		System.out.println("Mandatory arguments:");
		System.out.println("1) path to road shapefile (must be program first argument) [string] ");
		System.out.println("2) path to trajectory shapefile (must be second program argument) [string] ");
		System.out.println();
		System.out.println("Optional arguments:");
		// Core matching parameters
		System.out.println("-r  , search radius (default: 100.0) [double]");
		System.out.println("-k  , maximum number of candidates per trajectory point (default: INTEGER MAX VALUE) [int]");
		System.out.println("-c  , candidate cost weight (default: 0.01) [double]");
		System.out.println("-off  , use offroad candidates? (default: true) [boolean]");
		System.out.println("-offwt  , parameter for offroad weight (default: 15) [double]");
		// Computational parameters
		System.out.println("-t  , number of threads for multithreaded matching (default: 6)");
		System.out.println("-n  , number of trajectories in each batch (default: 1000)");
		// Additional penalties
		System.out.println("-pdev  , Penalty for deviation between subsequent matching lines (default: 1.4)");
		System.out.println("-pdist  , Penalty for deviation between euclidean and network distance  (default: 0.6) [double < 1.0]");
		System.out.println("-plink , Path to file specifying link-type penalties (default: none) [string]" );
		// Debugging / other
		System.out.println("-h  , flag to (just) print help");
		System.out.println("-v  , flag to print verbose output");
		System.out.println("-select, select trajectories by ID (comma-separated list)");
		System.out.println("-min, minimum trajectory ID to compute");
		System.out.println("-max, maximum trajectory ID to compute");
		// Input data details
		System.out.println("-nid, name of link ID column in road data [string]");
		System.out.println("-nwt, name of weight column in road data [string]");
		System.out.println("-tid, name of trajectory ID column in trajectory data [string]");
		// Output details
		System.out.println("-ma , include matches in output");
		System.out.println("-ch , include chunks in output");
		System.out.println("-chp, include chunk paths in output");
		System.out.println("-gp , include global paths in output");
		System.out.println("-s  , include segments in output");
	}

	public static <I> ArrayList<Point2D> extractPointsFromPath(
			List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>>> path,
			boolean unique) {
		ArrayList<Point2D> pointList = new ArrayList<>();
		for (DiGraphNode<Point2D, DoubleWeightDataWithInfo<I>> node : path) {
			if(unique & !pointList.isEmpty()) {
				if(node.getNodeData().equals(pointList.getLast())) {
					continue;
				}
			}
			pointList.add(node.getNodeData());
		}
		return pointList;
	}

	public static Map<String, Double> readWeightAdjustments(String weightAdjustmentsFile) {

		Map<String,Double> weightAdjustments = new LinkedHashMap<>();

		try {
			BufferedReader in = new BufferedReader(new FileReader(weightAdjustmentsFile));
			String recString;
			while((recString = in.readLine()) != null) {

				int sepIdx = recString.lastIndexOf(' ');
				double weight = Double.parseDouble(recString.substring(sepIdx+1));
				String criteria = recString.substring(0, sepIdx).replaceAll("(\"[\\w\\h]+\")|\\s*","$1").replace("\"", ""); // This regex removes white space unless within quotes

				if(StringUtils.countMatches(criteria,'=') != 1) {
					throw new RuntimeException("Equals symbol \"=\" should appear exactly once per line in weight adjustments file!");
				}

				weightAdjustments.put(criteria, weight);
				Logger.info("Incorporating weight factor of " + weight + " when " + criteria);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		return weightAdjustments;
	}

	private static boolean testCriteria(Object val, String testVal) {
		Class<?> c = val.getClass();

		if(c.equals(String.class)) {
			return val.equals(testVal);
		} else if (c.equals(Integer.class)) {
			return val.equals(Integer.parseInt(testVal));
		} else if (c.equals(Double.class)) {
			return val.equals(Double.parseDouble(testVal));
		} else if (c.equals(Long.class)) {
			return val.equals(Long.parseLong(testVal));
		} else if (c.equals(Boolean.class)) {
			return val.equals(Boolean.parseBoolean(testVal));
		} else {
			return false;
		}
	}

	static class TrajectoryWorker implements Runnable {

		private final ConcurrentLinkedQueue<Track> tracksQueue;
		private final DiGraph<Point2D, DoubleWeightDataWithInfo<RoadInfo>> g;
		private final AtomicLong counter;
		private final String typeColName;

		private final boolean saveMatches;
		private final boolean saveChunks;
		private final boolean saveChunkPaths;
		private final boolean saveGlobalPaths;
		private final boolean saveSegments;

		private final ArrayList<Track> paths;
		private final ArrayList<Track> matches;
		private final ArrayList<Track> chunks;
		private final ArrayList<Track> chunkPaths;
		private final ArrayList<Track> globalPaths;
		private final ArrayList<Track> segments;

		TrajectoryWorker(ConcurrentLinkedQueue<Track> tracksQueue, AtomicLong counter, String typeColName, LinkedList<Road<RoadInfo>> roads,
						 boolean saveMatches, boolean saveChunks, boolean saveChunkPaths, boolean saveGlobalPaths, boolean saveSegments) {
			this.tracksQueue = tracksQueue;
			this.counter = counter;
			this.g = Road.buildGraph(roads);
			this.typeColName = typeColName;
			this.saveMatches = saveMatches;
			this.saveChunks = saveChunks;
			this.saveChunkPaths = saveChunkPaths;
			this.saveGlobalPaths = saveGlobalPaths;
			this.saveSegments = saveSegments;
			this.paths = new ArrayList<>();
			this.matches = saveMatches ? new ArrayList<>() : null;
			this.chunks = saveChunks ? new ArrayList<>() : null;
			this.chunkPaths = saveChunkPaths ? new ArrayList<>() : null;
			this.globalPaths = saveGlobalPaths ? new ArrayList<>() : null;
			this.segments = saveSegments ? new ArrayList<>() : null;
		}

		@Override
		public void run() {
			while(true) {
				Track track = tracksQueue.poll();
				if(track == null) {
					return;
				}

				long id = counter.incrementAndGet();
				if(LongMath.isPowerOfTwo(id)) {
					Logger.info("Computing trajectory " + id);
				}

				// PERFORM MATCHING
				Matching<RoadInfo> m = new Matching<>(g, track, new RoadInfo(-1,Matching.OFF_ROAD_WEIGHT));

				// SAVE RESULTS
				// Paths
				paths.add(new Track(track.getId(), track.getSubtrack(), extractPointsFromPath(m.getPath(),true)));

				// Segments
				if(saveSegments) {
					int segmentCounter = 1;
					for (DiGraphArc<Point2D, DoubleWeightDataWithInfo<RoadInfo>> arc : m.getPathArcs()) {

						Point2D sourcePoint = arc.getSource().getNodeData();
						Point2D targetPoint = arc.getTarget().getNodeData();

						if(sourcePoint.equals(targetPoint)) {
							throw new RuntimeException("Source and target points are the same!");
						}

						ArrayList<Point2D> segment = new ArrayList<>();
						segment.add(sourcePoint);
						segment.add(targetPoint);
						Long type = typeColName != null ? arc.getArcData().getInfo().getId() : null;
						segments.add(new Track(track.getId(), track.getSubtrack(),segmentCounter++,type, segment));
					}
				}

				// Matches
				if(saveMatches) {
					ArrayList<Point2D> track_matchedPoints = extractPointsFromPath(m.getMatches(),false);
					LinkedList<Integer> segmentIdx = m.getMatchesArcCount();
					for (int i = 0; i < track.getTrackPoints().size(); i++) {
						ArrayList<Point2D> match = new ArrayList<>();
						match.add(track.getTrackPoints().get(i));
						match.add(track_matchedPoints.get(i));
						matches.add(new Track(track.getId(), track.getSubtrack(),segmentIdx.get(i), match));
					}
				}

				// Chunks
				if(saveChunks) {
					int chunkCounter = 1;
					for (ArrayList<DiGraphNode<Point2D, DoubleWeightDataWithInfo<RoadInfo>>> chunk : m.getChunks()) {
						chunks.add(new Track(track.getId(), track.getSubtrack(), chunkCounter++, extractPointsFromPath(chunk,false)));
					}
				}

				// Chunk paths
				if(saveChunkPaths) {
					for (List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<RoadInfo>>> chunkPath : m
							.getShortestPathsForChunks()) {
						chunkPaths.add(new Track(track.getId(), track.getSubtrack(), extractPointsFromPath(chunkPath,false)));
					}
				}

				// Global paths
				if(saveGlobalPaths) {
					List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<RoadInfo>>> p = m.getShortestPathForWholeTrajectory();
					if (p != null && p.size() > 1) {
						ArrayList<Point2D> p_points = extractPointsFromPath(p,false);
						globalPaths.add(new Track(track.getId(), track.getSubtrack(), p_points));
					}
				}

				// RESTORE GRAPH
				m.restoreGraph();
			}
		}

		public ArrayList<Track> getPaths() {
			return paths;
		}

		public ArrayList<Track> getMatches() {
			return matches;
		}

		public ArrayList<Track> getChunks() {
			return chunks;
		}

		public ArrayList<Track> getChunkPaths() {
			return chunkPaths;
		}

		public ArrayList<Track> getGlobalPaths() {
			return globalPaths;
		}

		public ArrayList<Track> getSegments() {
			return segments;
		}

	}
}
