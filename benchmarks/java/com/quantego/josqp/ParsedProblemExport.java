package com.quantego.josqp;

/** Export exactly the matrices read by JOSQP, without scaling or solving. */
public final class ParsedProblemExport {
    public static void main(String[] args) {
        StabilityBenchmark.exportData(Parser.readQmps(args[0]).getData(),args[1]);
    }
}
