package de.geoinfoBonn.graphLibrary.mapMatching.core.generic.data;

import java.io.Serializable;

public class DoubleWeightData implements WeightedArcData, Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1306575331185980551L;
	private double value;

	@Override
	public double getValue() {
		return value;
	}

	public void setValue(double value) {
		this.value = value;
	}

	public DoubleWeightData(double weight) {
		super();
		this.value = weight;
	}

	public DoubleWeightData() {

	}

	@Override
	public String toString() {
		return "DoubleWeightData [value=" + String.format("%.2f", value) + "]";
	}
}