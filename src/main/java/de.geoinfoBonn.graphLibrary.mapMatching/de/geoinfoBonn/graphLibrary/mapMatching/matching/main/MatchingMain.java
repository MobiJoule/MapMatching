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

import static de.geoinfoBonn.graphLibrary.mapMatching.matching.main.AbstractMain.getOptionalArg;
import static de.geoinfoBonn.graphLibrary.mapMatching.matching.main.AbstractMain.containsOptionalArg;

//@formatter:off
/**
 * Executable for map matching.
 *
 * Mandatory arguments:
 * 1) path to road shapefile (must be program first argument) [string]
 * 2) path to trajectory shapefile (must be second program argument) [string]
 *
 * Optional arguments:
 *  -t  , number of threads for multithreaded matching (default: 6)
 *  -n  , number of trajectories in each batch (default: 1000)
 *  -h  , flag to (just) print help
 *  -r  , parameter for radius (default: 25.0) [double]
 *  -k  , parameter for maximum number of candidates (default: INF) [int]
 *  -c  , parameter for candidate cost weight (default: 0.01) [double]
 *  -o  , use offroad candidates? (default: true) [boolean]
 *  -w  , parameter for offroad weight (default: 1.5) [double]
 *  -v  , flag to print verbose output
 *  -tid, name of trajectory ID column in trajectory data [string]
 *  -select, for processing only selected trajectory IDs (comma-separated)
 *  -min, minimum trajectory ID to compute
 *  -max, maximum trajectory ID to compute
 *  -nid, name of link ID column in road data [string]
 *  -nwt, name of weight column in road data [string]
 *  -ma , include matches in output
 *  -ch , include chunks in output
 *  -chp, include chunk paths in output
 *  -gp , include global paths in output
 *  -s  , include segments in output
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
		Matching.VERBOSE = containsOptionalArg(args, "-v");
		if (getOptionalArg(args, "-r") != null)
			Matching.RADIUS = Double.parseDouble(getOptionalArg(args, "-r"));
		if (getOptionalArg(args, "-k") != null)
			Matching.MAX_CAND_N = Integer.parseInt(getOptionalArg(args, "-k"));
		if (getOptionalArg(args, "-w") != null)
			Matching.OFF_ROAD_WEIGHT = Double.parseDouble(getOptionalArg(args, "-w"));
		if (getOptionalArg(args, "-c") != null)
			Matching.CANDIDATE_COST_WEIGHT = Double.parseDouble(getOptionalArg(args, "-c"));
		if (getOptionalArg(args, "-o") != null)
			Matching.ADD_OFFROAD_CANDIDATE = Boolean.parseBoolean(getOptionalArg(args, "-o"));

		// Which outputs to write
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

		// Partitions and threads
		String threadsInput = getOptionalArg(args, "-t");
		String partitionsInput = getOptionalArg(args, "-n");
		int numberOfThreads = threadsInput == null ? DEFAULT_NUMBER_OF_THREADS : Integer.parseInt(threadsInput);
		int partitionSize = partitionsInput == null ? DEFAULT_PARTITION_SIZE : Integer.parseInt(partitionsInput);

		// Weight adjustments
		String weightAdjustmentsFile = getOptionalArg(args, "-adj");
		final Map<String,Double> weightAdjustments = (weightAdjustmentsFile == null) ? null : readWeightAdjustments(weightAdjustmentsFile);

		// InfoGenerator
		Function<SimpleFeature, RoadInfo> roadInfoGenerator = feature -> {

			// RoadId
			long roadId = (long) feature.getAttribute(linkIdName);

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
		Logger.info("Number of trajectories in input file: " + trajectories.size());

		// Filter trajectories
		String selectedIds = getOptionalArg(args, "-select");
		String minInput = getOptionalArg(args, "-min");
		String maxInput = getOptionalArg(args, "-max");

		if(selectedIds != null) {
			Set<Long> includedIds = Arrays.stream(selectedIds.split(",")).map(Long::parseLong).collect(Collectors.toSet());
			trajectories = trajectories.stream().filter(f -> includedIds.contains(f.getId())).collect(Collectors.toList());
			Logger.info("Selected " + trajectories.size() + " trajectories with IDs: " + trajectories.size());
			if (minInput != null || maxInput != null) {
				throw new RuntimeException("For filtering IDs, shouldn't combine -select and -min/max");
			}
		} else if(minInput != null || maxInput != null) {
			long min = minInput == null ? Long.MIN_VALUE : Long.parseLong(minInput);
			long max = maxInput == null ? Long.MAX_VALUE : Long.parseLong(maxInput);
			trajectories = trajectories.stream().filter(f -> f.getId() >= min && f.getId() <= max).collect(Collectors.toList());
			Logger.info("Trajectories after filtering trajectory IDs: " + trajectories.size());
		}

		// Print arguments
		printArguments(args, partitionSize, numberOfThreads, linkDistId, linkIdName, trajectoryIdName, minInput, maxInput,
				writeMatches,writeChunks,writeChunkPaths,writeGlobalPaths,writeSegments);

		// Prepare trajectory data
		Iterable<List<Track>> partitions = Iterables.partition(trajectories, partitionSize);
		long numberOfPartitions = StreamSupport.stream(partitions.spliterator(), false).count();
		Logger.info("Split trajectories file into " + numberOfPartitions + " partitions containing " + partitionSize + " trajectories each.");

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
					Logger.info("Could not create spatial indexes!");
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

	private static void printArguments(String[] args, int partitionSize, int numberOfThreads, String netWtColName,
									   String netIdColName, String trajIdColName, String minInput, String maxInput,
									   boolean writeMatches, boolean writeChunks, boolean writeChunkPaths, boolean writeGlobalPaths, boolean writeSegments) {
		System.out.println("Program arguments:");
		System.out.println(" - batch size:              " + partitionSize);
		System.out.println(" - number of threads:       " + numberOfThreads);
		System.out.println(" - road data:               " + args[0]);
		System.out.println(" - trajectory data:         " + args[1]);
		System.out.println(" - radius:                  " + Matching.RADIUS);
		System.out.println(" - max. num. of candidates: " + Matching.MAX_CAND_N);
		System.out.println(" - candidate cost weight:   " + Matching.CANDIDATE_COST_WEIGHT);
		System.out.println(" - use offroad candidates?  " + Matching.ADD_OFFROAD_CANDIDATE);
		if (Matching.ADD_OFFROAD_CANDIDATE)
			System.out.println(" - offroad weight:         " + Matching.OFF_ROAD_WEIGHT);
		System.out.println(" - network weight column:   " + netWtColName);
		System.out.println(" - network ID column:       " + netIdColName);
		System.out.println(" - trajectory ID column:    " + trajIdColName);
		System.out.println(" - min trajectory ID to process  :   " + minInput);
		System.out.println(" - max trajectory ID to process :   " + maxInput);
		System.out.println(" - output Matches?      " + writeMatches);
		System.out.println(" - output Chunks?       " + writeChunks);
		System.out.println(" - output ChunkPaths?   " + writeChunkPaths);
		System.out.println(" - output GlobalPaths?  " + writeGlobalPaths);
		System.out.println(" - output Segments?     " + writeSegments);
	}

	private static void printHelpText() {
		System.out.println("Executable for map matching.");
		System.out.println();
		System.out.println("Mandatory arguments:");
		System.out.println("1) path to road shapefile (must be program first argument) [string] ");
		System.out.println("2) path to trajectory shapefile (must be second program argument) [string] ");
		System.out.println();
		System.out.println("Optional arguments:");
		System.out.println("-t  , number of threads for multithreaded matching (default: 6)");
		System.out.println("-n  , number of trajectories in each batch (default: 1000)");
		System.out.println("-h  , flag to (just) print help");
		System.out.println("-r  , parameter for radius (default: 25.0) [double]");
		System.out.println("-k  , parameter for maximum number of candidates (default: INF) [int]");
		System.out.println("-c  , parameter for candidate cost weight (default: 0.01) [double]");
		System.out.println("-o  , use offroad candidates? (default: true) [boolean]");
		System.out.println("-w  , parameter for offroad weight (default: 1.5) [double]");
		System.out.println("-v  , flag to print verbose output");
		System.out.println("-tid, name of trajectory ID column in trajectory data [string]");
		System.out.println("-min, minimum trajectory ID to compute");
		System.out.println("-max, maximum trajectory ID to compute");
		System.out.println("-nid, name of link ID column in road data [string]");
		System.out.println("-nwt, name of weight column in road data [string]");
		System.out.println("-ma , include matches in output");
		System.out.println("-ch , include chunks in output");
		System.out.println("-chp, include chunk paths in output");
		System.out.println("-gp , include global paths in output");
		System.out.println("-s  , include segments in output");
		System.out.println();
		System.out.println("Developed by Prof. Dr.-Ing. Jan-Henrik Haunert and Axel Forsch");
		System.out.println("Contact: forsch@igg.uni-bonn.de");
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
				String criteria = recString.substring(0, sepIdx).replaceAll("(\"[\\w\\h]+\")|\\s*","$1").replaceAll("\"", ""); // This regex removes white space unless within quotes

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
