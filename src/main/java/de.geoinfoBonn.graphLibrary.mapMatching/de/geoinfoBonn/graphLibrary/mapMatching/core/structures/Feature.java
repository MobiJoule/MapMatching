package de.geoinfoBonn.graphLibrary.mapMatching.core.structures;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;

/***
 * 
 * Class representing a feature with geometry and attributes
 *
 */
public class Feature {

	private Geometry geometry;
	private HashMap<String, Object> attributes;

	/**
	 * Initializes this feature without attributes.
	 * 
	 * @param geom the geometry of this feature
	 */
	public Feature(Geometry geom) {
		this(geom, null);
	}

	/**
	 * Initializes this feature with geometry and attributes.
	 * 
	 * @param geom       the geometry of this feature
	 * @param attributes the map of attributes for this feature
	 */
	public Feature(Geometry geom, HashMap<String, Object> attributes) {
		this.geometry = geom;
		if (attributes == null)
			this.attributes = new HashMap<String, Object>();
		else
			this.attributes = attributes;
	}

	/**
	 * Returns the (JTS-) geometry of this feature
	 * 
	 * @return a JTS geometry
	 */
	public Geometry getGeometry() {
		return geometry;
	}

	public void setGeometry(Geometry geometry) {
		this.geometry = geometry;
	}

	/**
	 * Sets the attribute for this geometry. Returns true if former attribute value
	 * got replaced, false if attribute is a new attribute.
	 * 
	 * @param key   name of the attribute
	 * @param value value of the attribute
	 * @return true if feature previously had a value for this attribute
	 */
	public boolean setAttribute(String key, Object value) {
		boolean replace = hasAttribute(key);
		attributes.put(key, value);
		return replace;
	}

	/**
	 * Returns if the feature currently has a value for the given attribute.
	 * 
	 * @param name name of the attribute
	 * @return boolean value for the result
	 */
	public boolean hasAttribute(String name) {
		return attributes.containsKey(name);
	}

	/**
	 * Returns the value of the attribute, or null if feature does not have this
	 * attribute.
	 * 
	 * @param name name of the attribute
	 * @return attribute value or null if absent
	 */
	public Object getAttribute(String name) {
		return attributes.get(name);
	}

	/**
	 * Returns a String representation of this features geometry type.
	 * 
	 * @return possible values: GeometryCollection, LinearRing, LineString,
	 *         MultiLineString, MultiPoint, MultiPolygon, Point, Polygon
	 */
	public String getGeometryType() {
		return geometry.getGeometryType();
	}

	/**
	 * Returns a list of all the attribute names this feature has.
	 * 
	 * @return ArrayList of the names
	 */
	public List<String> getAttributeNames() {
		return new ArrayList<>(attributes.keySet());
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Feature [geometryType=");
		sb.append(getGeometryType());
		for (Entry<String, Object> attribute : attributes.entrySet()) {
			sb.append(",");
			sb.append(attribute.getKey());
			sb.append("=");
			sb.append(attribute.getValue());
		}
		sb.append("]");

		return sb.toString();
	}

	public static final class Builder {

		private GeometryFactory gf;

		private LinkedList<Coordinate> _coordinates;
		private HashMap<String, Object> _attributes;

		public Builder() {
			gf = new GeometryFactory();
			_coordinates = new LinkedList<>();
			_attributes = new HashMap<>();
		}

		public Builder coordinate(double x, double y) {
			_coordinates.add(new Coordinate(x, y));
			return this;
		}

		public Builder coordinate(Coordinate c) {
			_coordinates.add(c);
			return this;
		}

		public Builder coordinates(Coordinate[] c) {
			_coordinates.addAll(Arrays.asList(c));
			return this;
		}

		public Builder attribute(String key, Object value) {
			_attributes.put(key, value);
			return this;
		}

		public Feature buildPointFeature() {
			Geometry geom;
			if (_coordinates.size() == 1) {
				geom = gf.createPoint(_coordinates.getFirst());
			} else if (_coordinates.size() > 1) {
				geom = gf.createMultiPointFromCoords(_coordinates.toArray(new Coordinate[_coordinates.size()]));
			} else {
				throw new IllegalArgumentException("Number of points is zero");
			}
			return new Feature(geom, _attributes);
		}

		public Feature buildPointFeature(Point point) {
			if (_coordinates.size() > 0)
				System.err.println("Coordinates ignored.");
			return new Feature(point, _attributes);
		}

		public Feature buildLineFeature() {
			Geometry geom = gf.createLineString(_coordinates.toArray(new Coordinate[_coordinates.size()]));
			return new Feature(geom, _attributes);
		}

		public Feature buildLineFeature(LineString line) {
			if (_coordinates.size() > 0)
				System.err.println("Coordinates ignored.");
			return new Feature(line, _attributes);
		}

		public Feature buildPolygonFeature() {
			Geometry geom = gf.createPolygon(_coordinates.toArray(new Coordinate[_coordinates.size()]));
			return new Feature(geom, _attributes);
		}

		public Feature buildPolygonFeature(Polygon polygon) {
			if (_coordinates.size() > 0)
				System.err.println("Coordinates ignored.");
			return new Feature(polygon, _attributes);
		}

		public Feature buildCopy(Feature feature) {
			Geometry geom = feature.geometry.copy();
			for (Entry<String, Object> attribute : feature.attributes.entrySet()) {
				_attributes.put(attribute.getKey(), attribute.getValue());
			}
			return new Feature(geom, _attributes);
		}
	}

	public static Builder builder() {
		return new Builder();
	}
}
