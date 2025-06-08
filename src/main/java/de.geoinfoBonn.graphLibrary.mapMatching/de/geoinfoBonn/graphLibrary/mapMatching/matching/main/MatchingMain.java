package de.geoinfoBonn.graphLibrary.mapMatching.matching.main;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import org.opengis.referencing.crs.CoordinateReferenceSystem;

import de.geoinfoBonn.graphLibrary.core.generic.DiGraph;
import de.geoinfoBonn.graphLibrary.core.generic.DiGraph.DiGraphArc;
import de.geoinfoBonn.graphLibrary.core.generic.DiGraph.DiGraphNode;
import de.geoinfoBonn.graphLibrary.core.generic.DoubleWeightDataWithInfo;
import de.geoinfoBonn.graphLibrary.mapMatching.io.Road;
import de.geoinfoBonn.graphLibrary.mapMatching.io.RoadReader;
import de.geoinfoBonn.graphLibrary.mapMatching.matching.Matching;
import de.geoinfoBonn.graphLibrary.mapMatching.matching.Track;

//@formatter:off
/**
 * Executable for map matching.
 * 
 * Mandatory arguments:
 * 1) path to road shapefile (must be program first argument) [string] 
 * 2) path to trajectory shapefile (must be second program argument) [string] 
 * 
 * Optional arguments:
 * -h  , flag to (just) print help
 * -r  , parameter for radius (default: 25.0) [double]
 * -k  , parameter for maximum number of candidates (default: INF) [int]
 * -c  , parameter for candidate cost weight (default: 0.01) [double]
 * -o  , add offroad candidates? (default: true) [boolean]
 * -w  , parameter for offroad weight (default: 1.5) [double]
 * -v  , flag to print verbose output
 * -l  , name of weight column in road data [string]
 * -t  , name of type column in road data [string]
 * -ma , output path for matches [string]
 * -ch , output path for chunks [string]
 * -chp, output path for chunk paths [string]
 * -gp , output path for global paths [string]
 * -s  , output path for segments [string]
 * 
 * @author Jan-Henrik Haunert
 * @author Axel Forsch (forsch@igg.uni-bonn.de)
 * @version %I%, %G%
 *
 */
//@formatter:on
public class MatchingMain {
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

		// The shapefile's column that contains the double attribute with the road's
		// weight (e.g., geometric length, travel time)
		// If left unspecified, the geometric length computed from coordinates is used
		// as the road's weight
		String weightColName = containsOptionalArg(args, "-l") ? getOptionalArg(args, "-l") : null;
		String typeColName = containsOptionalArg(args, "-t") ? getOptionalArg(args, "-t") : null;

		if (Matching.VERBOSE)
			printArguments(args, weightColName, typeColName);

		LinkedList<Road<String>> roads = RoadReader.importFromShapefile(args[0],
				f -> (String) f.getAttribute(typeColName), weightColName);

		CoordinateReferenceSystem crs = null;
		try {
			crs = RoadReader.readCRS(args[1].split("\\.")[0] + ".prj");
		} catch (FileNotFoundException e) {
			System.err.println(
					"Unable to read coordinate reference system of input data. No .prj file found. Defaulting to WGS84.");
		}

		DiGraph<Point2D, DoubleWeightDataWithInfo<String>> g = Road.<String>buildGraph(roads);
		ArrayList<Track> trajectories = Track.importFromShapefile(args[1]);
		System.out.println("Number of trajectories in input file: " + trajectories.size());
		ArrayList<Track> paths = new ArrayList<Track>();
		ArrayList<Track> matches = new ArrayList<Track>();
		ArrayList<Track> chunks = new ArrayList<Track>();
		ArrayList<Integer> chunkIds = new ArrayList<Integer>();
		ArrayList<Track> chunkPaths = new ArrayList<Track>();
		ArrayList<Track> globalPaths = new ArrayList<Track>();
		ArrayList<Track> segments = new ArrayList<Track>();
		ArrayList<Integer> segmentIds = new ArrayList<Integer>();
		ArrayList<String> types = null;
		if (typeColName != null) {
			types = new ArrayList<String>();
		}

		int counter = 1;
		for (Track track : trajectories) {
			System.out.print("Starting track #" + counter++ + " (id=" + track.getId() + ", subtrack="
					+ track.getSubtrack() + ")");
			if (Matching.VERBOSE)
				System.out.println();
			else
				System.out.print(" ... ");
			Matching<String> m = new Matching<String>(g, track.getTrackPoints());
			paths.add(new Track(track.getId(), track.getSubtrack(), extractPointsFromPath(m.getPath())));

			int segmentCounter = 1;
			for (DiGraphArc<Point2D, DoubleWeightDataWithInfo<String>> arc : m.getPathArcs()) {
				ArrayList<Point2D> segment = new ArrayList<Point2D>();
				segment.add(arc.getSource().getNodeData());
				segment.add(arc.getTarget().getNodeData());
				segments.add(new Track(track.getId(), track.getSubtrack(), segment));
				segmentIds.add(segmentCounter++);
				if (typeColName != null) {
					types.add(arc.getArcData().getInfo());
				}
			}

			ArrayList<Point2D> track_matchedPoints = extractPointsFromPath(m.getMatches());
			for (int i = 0; i < track.getTrackPoints().size(); i++) {
				ArrayList<Point2D> match = new ArrayList<Point2D>();
				match.add(track.getTrackPoints().get(i));
				match.add(track_matchedPoints.get(i));
				matches.add(new Track(track.getId(), track.getSubtrack(), match));
			}

			int chunkCounter = 1;
			for (ArrayList<DiGraphNode<Point2D, DoubleWeightDataWithInfo<String>>> chunk : m.getChunks()) {
				chunks.add(new Track(track.getId(), track.getSubtrack(), extractPointsFromPath(chunk)));
				chunkIds.add(chunkCounter++);
			}
			for (List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<String>>> chunkPath : m
					.getShortestPathsForChunks()) {
				chunkPaths.add(new Track(track.getId(), track.getSubtrack(), extractPointsFromPath(chunkPath)));
			}

			List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<String>>> p = m.getShortestPathForWholeTrajectory();
			if (p != null && p.size() > 1) {
				ArrayList<Point2D> p_points = extractPointsFromPath(p);
				globalPaths.add(new Track(track.getId(), track.getSubtrack(), p_points));
			}
			m.restoreGraph();
			if (Matching.VERBOSE)
				System.out.print("Track ");
			System.out.println("done!");
		}

		File pathsFile = new File(args[2]);
		Track.writeToShapeFile(pathsFile, paths, null, null, crs);

		// Optional outputs
		File matchesFile = null;
		if (getOptionalArg(args, "-ma") != null) {
			matchesFile = new File(getOptionalArg(args, "-ma"));
			Track.writeToShapeFile(matchesFile, matches, null, null, crs);
		}

		File chunksFile = null;
		if (getOptionalArg(args, "-ch") != null) {
			chunksFile = new File(getOptionalArg(args, "-ch"));
			Track.writeToShapeFile(chunksFile, chunks, chunkIds, null, crs);
		}

		File chunksPathFile = null;
		if (getOptionalArg(args, "-chp") != null) {
			chunksPathFile = new File(getOptionalArg(args, "-chp"));
			Track.writeToShapeFile(chunksPathFile, chunkPaths, chunkIds, null, crs);
		}

		File globalPathsFile = null;
		if (getOptionalArg(args, "-gp") != null) {
			globalPathsFile = new File(getOptionalArg(args, "-gp"));
			Track.writeToShapeFile(globalPathsFile, globalPaths, null, null, crs);
		}

		File segmentsFile = null;
		if (getOptionalArg(args, "-s") != null) {
			segmentsFile = new File(getOptionalArg(args, "-s"));
			Track.writeToShapeFile(segmentsFile, segments, segmentIds, types, crs);
		}
	}

	private static void printArguments(String[] args, String weightColName, String typeColName) {
		System.out.println("Program arguments:");
		System.out.println(" - road data:               " + args[0]);
		System.out.println(" - trajectory data:         " + args[1]);
		System.out.println(" - radius:                  " + Matching.RADIUS);
		System.out.println(" - max. num. of candidates: " + Matching.MAX_CAND_N);
		System.out.println(" - candidate cost weight:   " + Matching.CANDIDATE_COST_WEIGHT);
		System.out.println(" - use offroad candidates?  " + Matching.ADD_OFFROAD_CANDIDATE);
		if (Matching.ADD_OFFROAD_CANDIDATE)
			System.out.println(" - offroad weight:          " + Matching.OFF_ROAD_WEIGHT);
		System.out.println(" - weight column:           " + weightColName);
		System.out.println(" - type column:             " + typeColName);
		System.out.println(" - output matches:          " + getOptionalArg(args, "-ma"));
		System.out.println(" - output chunks:           " + getOptionalArg(args, "-ch"));
		System.out.println(" - output chunk paths:      " + getOptionalArg(args, "-chp"));
		System.out.println(" - output global paths:     " + getOptionalArg(args, "-gp"));
		System.out.println(" - output segments:         " + getOptionalArg(args, "-s"));
		System.out.println();
	}

	private static void printHelpText() {
		System.out.println("Executable for map matching.");
		System.out.println();
		System.out.println("Mandatory arguments:");
		System.out.println("1) path to road shapefile (must be program first argument) [string] ");
		System.out.println("2) path to trajectory shapefile (must be second program argument) [string] ");
		System.out.println();
		System.out.println("Optional arguments:");
		System.out.println("-h  , flag to (just) print help");
		System.out.println("-r  , parameter for radius (default: 25.0) [double]");
		System.out.println("-k  , parameter for maximum number of candidates (default: INF) [int]");
		System.out.println("-c  , parameter for candidate cost weight (default: 0.01) [double]");
		System.out.println("-o  , use offroad candidates? (default: true) [boolean]");
		System.out.println("-w  , parameter for offroad weight (default: 1.5) [double]");
		System.out.println("-v  , flag to print verbose output");
		System.out.println("-l  , name of weight column in road data [string]");
		System.out.println("-t  , name of type column in road data [string]");
		System.out.println("-ma , output path for matches [string]");
		System.out.println("-ch , output path for chunks [string]");
		System.out.println("-chp, output path for chunk paths [string]");
		System.out.println("-gp , output path for global paths [string]");
		System.out.println("-s  , output path for segments [string]");
		System.out.println();
		System.out.println("Developed by Prof. Dr.-Ing. Jan-Henrik Haunert and Axel Forsch");
		System.out.println("Contact: forsch@igg.uni-bonn.de");
	}

	public static String getOptionalArg(String[] args, String identifier) {
		for (int i = 0; i < args.length - 1; i++) {
			if (args[i].equals(identifier))
				return args[i + 1];
		}
		return null;
	}

	public static boolean containsOptionalArg(String[] args, String identifier) {
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals(identifier))
				return true;
		}
		return false;
	}

	public static ArrayList<Point2D> extractPointsFromPath(
			List<DiGraphNode<Point2D, DoubleWeightDataWithInfo<String>>> path) {
		ArrayList<Point2D> pointList = new ArrayList<Point2D>();
		for (DiGraphNode<Point2D, DoubleWeightDataWithInfo<String>> node : path) {
			pointList.add(node.getNodeData());
		}
		return pointList;
	}
}
