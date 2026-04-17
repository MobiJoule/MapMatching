package de.geoinfoBonn.graphLibrary.mapMatching.core.generic;

import de.geoinfoBonn.graphLibrary.mapMatching.core.generic.data.DoubleWeightData;

public class DoubleWeightDataWithInfo<I> extends DoubleWeightData {

	private static final long serialVersionUID = -4952717897982831931L;

	protected I s;

	public DoubleWeightDataWithInfo(double w, I s) {
		super(w);
		this.s = s;
	}

	public I getInfo() {
		return s;
	}

	@Override
	public String toString() {
		return super.toString() + " " + s;
	}
}
