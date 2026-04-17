package de.geoinfoBonn.graphLibrary.mapMatching.core.geometry;

import java.awt.geom.Line2D;
import java.awt.geom.Point2D;

public class Calculations {

	public static Point2D nearestPointOnSegment(Point2D p, Line2D seg) {
		double sx = seg.getX1();
		double sy = seg.getY1();
		double ex = seg.getX2();
		double ey = seg.getY2();

		// Laenge von Vektor a
		double aLength = Math.hypot(ex - sx, ey - sy);
		// Einheitsvektor zu a
		double[] aNormalized = { (ex - sx) / aLength, (ey - sy) / aLength };
		// Vektor b
		double[] b = { p.getX() - sx, p.getY() - sy };

		// Skalarprodukt aus Einheitsvektor und b
		double product = b[0] * aNormalized[0] + b[1] * aNormalized[1];

		// Berechnung von Vektor c
		double[] c = { (product * aNormalized[0]), (product * aNormalized[1]) };

		// Berechnung des Fusspunkts
		Point2D.Double basePoint = new Point2D.Double(sx + c[0], sy + c[1]);

		// Pruefen, ob der Fusspunkt auf seg liegt
		double distanceToStart = Math.hypot(basePoint.getX() - sx, basePoint.getY() - sy);
		double distanceToEnd = Math.hypot(basePoint.getX() - ex, basePoint.getY() - ey);
		if (distanceToEnd <= aLength && distanceToStart <= aLength) {
			return basePoint;
		} else if (distanceToEnd > distanceToStart) {
			return new Point2D.Double(sx, sy);
		} else {
			return new Point2D.Double(ex, ey);
		}
	}
}
