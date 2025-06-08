package de.geoinfoBonn.graphLibrary.mapMatching.io;

import java.awt.geom.Point2D;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Function;

import org.tinylog.Logger;

import org.geotools.api.data.SimpleFeatureReader;
import org.geotools.api.referencing.crs.CoordinateReferenceSystem;
import org.geotools.geopkg.GeoPackage;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.geotools.api.feature.simple.SimpleFeature;

public class RoadReader<I> {

	/**
	 * method for reading a shapefile with roads
	 *
	 * @param filename:      the name of the shapefile
	 * @param infoGenerator:   function for retrieving link ID
	 * @param weightColName: the name of the column containing the weight data
	 * @return a list of roads
	 * @throws IOException
	 */
	public static <I> LinkedList<Road<I>> importFromGpkg(String filename, Function<SimpleFeature, I> infoGenerator,
														 String weightColName, Map<String,Double> weightAdjustments) {
		LinkedList<Road<I>> roads = new LinkedList<>();
		int featureCount = 0;
		int importedFeatureCount = 0;
		Map<String,Integer> adjustedFeatureCount = new HashMap<>();
		if(weightAdjustments != null) {
			weightAdjustments.forEach((k,v) -> adjustedFeatureCount.put(k,0));
		}

		try {
			if (filename.endsWith(".gpkg")) {

				GeoPackage geopkg = new GeoPackage(new File(filename));
				SimpleFeatureReader reader = geopkg.reader(geopkg.features().get(0), null, null);

				while(reader.hasNext()) {
					SimpleFeature feature = reader.next();

					I type = null;
					if (infoGenerator != null) {
						type = infoGenerator.apply(feature);
					}
					Geometry myGeometry = (Geometry) feature.getDefaultGeometry();
					double totalLength = myGeometry.getLength();

					// Compute weight (or use length as default)
					double weight = weightColName != null ? (double) feature.getAttribute(weightColName) : totalLength;

					// Adjust weights
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
									adjustedFeatureCount.merge(e.getKey(), 1, Integer::sum);
								}
							}
						}
					}

					if (myGeometry.getGeometryType().equals("LineString")) {
						ArrayList<Point2D> points = new ArrayList<>();
						Coordinate[] xyz = myGeometry.getCoordinates();
                        for (Coordinate coordinate : xyz) {
                            Point2D p = new Point2D.Double(coordinate.x, coordinate.y);
                            points.add(p);
                        }
						Road<I> r = new Road<>(points, type, weight);
						roads.add(r);
						importedFeatureCount++;
					} else if (myGeometry.getGeometryType().equals("MultiLineString")) {
						for (int geoIndex = 0; geoIndex < myGeometry.getNumGeometries(); geoIndex++) {
							Geometry myGeometryPart = myGeometry.getGeometryN(geoIndex);
							double partLength = myGeometryPart.getLength();
							ArrayList<Point2D> points = new ArrayList<>();
							Coordinate[] xyz = myGeometryPart.getCoordinates();
                            for (Coordinate coordinate : xyz) {
                                Point2D p = new Point2D.Double(coordinate.x, coordinate.y);
                                points.add(p);
                            }
							Road<I> r = new Road<>(points, type, weight * partLength / totalLength);
							roads.add(r);
							importedFeatureCount++;
						}
					} else {
						Logger.warn("Feature is not a polyline but a " + myGeometry.getGeometryType());
					}
					featureCount++;
				}
				reader.close();
				geopkg.close();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		// Print info to logger
		Logger.info(featureCount + " features read.");
		Logger.info(importedFeatureCount + " features imported to the road graph.");
		for(Map.Entry<String,Integer> entry : adjustedFeatureCount.entrySet()) {
			Logger.info(entry.getValue() + " features with attribute \"" + entry.getKey() + "\" adjusted with factor of " + weightAdjustments.get(entry.getKey()));
		}




		return roads;
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

	public static CoordinateReferenceSystem readCRS(String filename) {

		GeoPackage geopkg;

		try {
			geopkg = new GeoPackage(new File(filename));
			SimpleFeatureReader reader = geopkg.reader(geopkg.features().get(0), null, null);
			CoordinateReferenceSystem crs = reader.getFeatureType().getCoordinateReferenceSystem();
			reader.close();
			geopkg.close();
			return crs;
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

	}
}
