package irreducibleshifts;

import java.io.*;
import java.util.*;
import udslibz.*;

/**
 * Checks the properties of irreducible sequence classes.
 * @author UndeadScythes
 */
public final class IrreducibleShifts3 {
    private IrreducibleShifts3() {}

    public static void main(final String[] args) throws IOException {
        final int lowN = ArgUtils.getInt(args, "-l", 0);
        final int highN = ArgUtils.getInt(args, "-h", 0);
        for(int n = lowN; n <= highN; n++) { // current test degree n
            final BufferedWriter out = new BufferedWriter(new FileWriter("IrreducibleShifts2-" + n + ".csv"));
            out.write("n,f,d,q,g,dec seqs,class reps,trace seqs,dec classes,dec phases\n");
            final int N = (1 << n) - 1; // m-sequence length N
            Polynomial f = PolynomialUtils.getPrimitive(n, 0); // primitive polynomial f
            while(f != null) {
                final List<Integer> factors = IntegerUtils.getFactors(N); // factors of N
                if(factors.isEmpty()) {
                    break;
                }
                final GaloisLFSR lfsr = new GaloisLFSR(n, f, 1);
                for(int d : factors) {
                    Polynomial g = PolynomialUtils.getStrictIrreducible(n, d, 0); // irreducible polynomial g
                    if(g == null) {
                        continue;
                    }
                    while(g != null) {
                        lfsr.reset(1);
                        final int q = N / d; // length of sub sequences q
                        final Sequence[] seqs = new Sequence[q]; // sub sequences
                        for(int i = 0; i < q; i++) {
                            seqs[i] = new Sequence(2 * d);
                        }
                        for(int i = 0; i < 2 * d; i++) {
                            for(int j = 0; j < q; j++) {
                                seqs[j].setElement(i, lfsr.getBit(0));
                                lfsr.clock();
                            }
                        }
                        final Polynomial alpha = g.getPrimitiveRoot(); // alpha a primitive element in terms of beta
                        final List<Polynomial> classReps = g.getClassReps(alpha);
                        final List<Sequence> classSeqs = new ArrayList<Sequence>();
                        for(int i = 0; i < q; i++) {
                            classSeqs.add(Sequence.fromTrace(classReps.get(i), g, n));
                        }
                        boolean allMatch = true;
                        for(int i = 0; i < q; i++) {
                            boolean flag = false;
                            for(int j = 0; j < q; j++) {
                                if(seqs[i].contains(classSeqs.get(j))) {
                                    flag = true;
                                    break;
                                }
                            }
                            if(flag) {
                                continue;
                            } else {
                                allMatch = false;
                            }
                        }
                        if(allMatch) {
                            final Sequence[] seqSnaps = new Sequence[n]; // snapshots of lfsr sequences
                            final Polynomial[] seqClasses = new Polynomial[n]; // sequence class of lfsr
                            final int[] phases = new int[n]; // phases of lfsr sequences
                            final GaloisLFSR lfsr2 = new GaloisLFSR(n, g, 1); // Galois LFSR with feedback polynomial g
                            for(int i = 0; i < n; i++) {
                                for(int j = 0; j < n; j++) {
                                    if(i == 0) {
                                        seqSnaps[j] = new Sequence(n);
                                    }
                                    seqSnaps[j].setElement(i, lfsr2.getBit(j));
                                }
                                lfsr2.clock();
                            }
                            int seqsFound = 0;
                            int phase = 0;
                            while(seqsFound != n) {
                                for(int i = 0; i < n; i++) {
                                    if(seqClasses[i] == null) {
                                        for(int j = 0; j < q; j++) {
                                            if(seqSnaps[i].equalTo(classSeqs.get(j))) {
                                                seqClasses[i] = classReps.get(j);
                                                phases[i] = phase;
                                                seqsFound++;
                                            }
                                        }
                                    }
                                }
                                for(int i = 0; i < n; i++) {
                                    seqSnaps[i].leftShift(1);
                                    seqSnaps[i].setElement(n - 1, lfsr.getBit(i));
                                }
                                lfsr.clock();
                                phase++;
                            }
                            System.out.println("n:" + n + "\tf:" + f.toBinary() + "\td:" + d + "\tq:" + q + "\tg:" + g.toBinary());
                            out.write(n + "," + f.toBinary() + "," + d + "," + q + "," + g.toBinary() + ",");
                            Polynomial[] decSeqClasses = new Polynomial[q];
                            int[] decSeqPhases = new int[q];
                            for(int i = 0; i < q; i++) {
                                boolean flag = false;
                                for(int j = 0; j < q; j++) {
                                    if(seqs[i].contains(classSeqs.get(j))) {
                                        decSeqClasses[i] = classReps.get(j);
                                        decSeqPhases[i] = seqs[i].find(classSeqs.get(j));
                                        System.out.println(" s" + i + "=" + seqs[i].getSubSequence(0, d).toString());
                                        out.write(seqs[i].getSubSequence(0, d).toString() + " ");
                                        flag = true;
                                        break;
                                    }
                                }
                                if(flag) {
                                    continue;
                                }
                            }
                            out.write("," + ArrayUtils.joinBinary(classReps, " ") + "," + ArrayUtils.join(classSeqs, " ") + "," + ArrayUtils.joinBinary(decSeqClasses, " ") + "," + ArrayUtils.join(decSeqPhases, " ") + "\n");
                            System.out.println(" " + "class reps:{" + ArrayUtils.joinBinary(classReps, ", ") + "}\ttrace seqs:{" + ArrayUtils.join(classSeqs, ", ") + "}");
                            System.out.println(" " + "decimation classes:{" + ArrayUtils.joinBinary(decSeqClasses, ", ") + "}\tdecimation phases:{" + ArrayUtils.join(decSeqPhases, ", ") + "}");
                            System.out.println();
                        }
                        g = PolynomialUtils.getStrictIrreducible(n, d, g.nextPoly().toInt());
                    }
                }
                f = PolynomialUtils.getPrimitive(n, f.nextPoly().toInt());
            }
        out.close();
        }
    }
}
