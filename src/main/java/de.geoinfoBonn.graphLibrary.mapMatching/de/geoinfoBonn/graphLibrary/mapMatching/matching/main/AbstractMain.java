package de.geoinfoBonn.graphLibrary.mapMatching.matching.main;

public class AbstractMain {

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
}
