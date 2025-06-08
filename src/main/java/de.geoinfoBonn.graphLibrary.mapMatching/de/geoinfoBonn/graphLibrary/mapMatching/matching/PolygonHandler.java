package de.geoinfoBonn.graphLibrary.mapMatching.matching;

import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.geotools.geometry.jts.JTSFactoryFinder;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.index.strtree.STRtree;
import org.locationtech.jts.io.ParseException;
import org.locationtech.jts.io.WKTReader;

import de.geoinfoBonn.graphLibrary.core.structures.Feature;

/**
 * Provides indexed access to the polygons.
 * 
 * @author Axel Forsch
 *
 */
public class PolygonHandler {

	private STRtree index;

	/**
	 * Creates a new, empty PolygonHandler
	 */
	public PolygonHandler() {
		this.index = new STRtree();
	}

	/**
	 * Creates a new PolygonHandler and initializes polygons with the given csv-file
	 * 
	 * @param polygonCsv path to the polygon file
	 */
	public PolygonHandler(String polygonCsv) {
		this();
		this.readPolygons(polygonCsv);
	}

	/**
	 * Returns a list of all polygons, whose bounding box point
	 * 
	 * @param p point p to check for intersection with polygons bboxes
	 * @return all bounding boxes, whose bounding box cover <code>p</code>
	 */
	@SuppressWarnings("unchecked")
	public List<Feature> queryPolygons(Point2D p) {
		return index.query(new Envelope(p.getX(), p.getX(), p.getY(), p.getY()));
	}

	/**
	 * Returns a list of all polygons, whose bounding box point
	 * 
	 * @param p point p to check for intersection with polygons bboxes
	 * @return all bounding boxes, whose bounding box cover <code>p</code>
	 */
	@SuppressWarnings("unchecked")
	public List<Feature> queryPolygons(Point p) {
		return index.query(p.getEnvelopeInternal());
	}

	/**
	 * Returns a list of all polygons, whose bounding box intersect the given
	 * bounding box
	 * 
	 * @param env bounding box to check for intersection with polygons bboxes
	 * @return all bounding boxes, whose bounding box intersect with
	 *         <code>env</code>
	 */
	@SuppressWarnings("unchecked")
	public List<Feature> queryPolygons(Envelope env) {
		return index.query(env);
	}

	/**
	 * Returns a list of all polygons, whose bounding box intersect the given
	 * bounding box
	 * 
	 * @param env bounding box to check for intersection with polygons bboxes
	 * @return all bounding boxes, whose bounding box intersect with
	 *         <code>env</code>
	 */
	@SuppressWarnings("unchecked")
	public List<Feature> queryPolygons(de.geoinfoBonn.graphLibrary.core.geometry.Envelope env) {
		return index.query(new Envelope(env.getxMin(), env.getxMax(), env.getyMin(), env.getyMax()));
	}

	/**
	 * Returns the number of polygons handled by this object.
	 * 
	 * @return number of polygons
	 */
	public int numPolygons() {
		return index.size();
	}

	/**
	 * Returns the depth of the STRtree used to index the polygons
	 * 
	 * @return depth of the STRtree
	 */
	public int indexDepth() {
		return index.depth();
	}

	/**
	 * Loades polygons from csv file.
	 * 
	 * @param polygonCsv
	 */
	public void readPolygons(String polygonCsv) {
		GeometryFactory gf = JTSFactoryFinder.getGeometryFactory();
		WKTReader wktReader = new WKTReader(gf);

		HashMap<String, Object> attributes;
		String[] tokens = null;

		try (BufferedReader in = new BufferedReader(new FileReader(polygonCsv))) {
			String line = in.readLine(); // skip header line
			while ((line = in.readLine()) != null) {
				tokens = line.split(";");

				attributes = new HashMap<>();
				attributes.put("id", Integer.parseInt(tokens[0]));
				attributes.put("osm_id", Integer.parseInt(tokens[1]));

//				tokens[2] = tokens[2].replace(", ()", "");

				Polygon p = (Polygon) wktReader.read(tokens[2]);
				Feature f = new Feature(p, attributes);
				index.insert(p.getEnvelopeInternal(), f);
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e1) {
			e1.printStackTrace();
		} catch (ParseException e) {
			System.err.println(tokens[2]);
			e.printStackTrace();
		} catch (IllegalArgumentException e) {
			System.err.println(tokens[2]);
			e.printStackTrace();
		}
	}

	public ArrayList<Feature> getAllPolygons() {
		return getAllPolygonsRecursive(index.itemsTree());
	}

	@SuppressWarnings({ "rawtypes", "unchecked" })
	public ArrayList<Feature> getAllPolygonsRecursive(List treeItem) {
		if (treeItem.get(0) instanceof Feature) {
			return (ArrayList<Feature>) treeItem;
		}
		ArrayList collected = new ArrayList();
		for (var item : treeItem) {
			collected.addAll(getAllPolygonsRecursive((List) item));
		}
		return collected;
	}
}
