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

        for(int n = lowN; n <= highN; n++) { // degree of the current test polynomial
            final BufferedWriter out = new BufferedWriter(new FileWriter("IrreducibleShifts2-" + n + ".csv"));
            out.write("degree,dec poly,bin poly,order,q,dec alpha,bin alpha,dec class reps,bin class reps,sample seq,dec sample output,bin sample output,class phases\n");
            Polynomial g = PolynomialUtils.getStrictIrreducible(n, 1 << n); // current test polynomial g
            while(g != null) {
                final int d = g.getOrder(); // order of g and length of sequences d
                final int N = (1 << n) - 1; // length of m-sequence N
                final int q = N / d; // number of classes q
                final Polynomial alpha = g.getPrimitiveRoot(); // primitive root alpha represented in beta
                final List<Polynomial> classReps = g.getClassReps(alpha); // class representitives in beta
                final List<Sequence> classSeqs = new ArrayList<Sequence>(); // example class sequences in statndard phase
                for(int i = 0; i < q; i++) {
                    classSeqs.add(Sequence.fromTrace(classReps.get(i), g, n));
                }
                final Sequence[] seqSnaps = new Sequence[n]; // snapshots of lfsr sequences
                final Polynomial[] seqClasses = new Polynomial[n]; // sequence class of lfsr
                final int[] phases = new int[n]; // phases of lfsr sequences
                final GaloisLFSR lfsr = new GaloisLFSR(n, g, 1); // Galois LFSR with feedback polynomial g
                for(int i = 0; i < n; i++) {
                    for(int j = 0; j < n; j++) {
                        if(i == 0) {
                            seqSnaps[j] = new Sequence(n);
                        }
                        seqSnaps[j].setElement(i, lfsr.getBit(j));
                    }
                    lfsr.clock();
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
                System.out.println("n:" + n + "\tg:" + g.toString() + "\td:" + d + "\tq:" + q + "\tclass reps:{" + ArrayUtils.joinBinary(classReps, " ") + "}\tseqs:{" + ArrayUtils.join(classSeqs, " ") + "}");
                out.write(n + "," + g.toInt() + "," + g.toBinary() + "," + d + "," + q + "," + alpha.toInt() + "," + alpha.toBinary() + "," + ArrayUtils.joinDecimal(classReps, " ") + "," + ArrayUtils.joinBinary(classReps, " ") + "," + ArrayUtils.join(classSeqs, " ") + "," + ArrayUtils.joinDecimal(seqClasses, " ") + "," + ArrayUtils.joinBinary(seqClasses, " ") + "," + ArrayUtils.join(phases, " ") + "\n");
                g = PolynomialUtils.getStrictIrreducible(n, g.toInt() + 1);
            }
            out.close();
        }

        for(int n = lowN; n <= highN; n++) {
            final int N = (1 << n) - 1;
            final Polynomial prim = PolynomialUtils.getPrimitive(n, 0);
            final List<Integer> factors = IntegerUtils.getFactors(N);
            if(factors.isEmpty()) {
                continue;
            }
            final GaloisLFSR lfsr = new GaloisLFSR(n, prim, 1);
            for(int d : factors) {
                Polynomial g = PolynomialUtils.getStrictIrreducible(n, d, 0);
                if(g == null) {
                    continue;
                }
                while(g != null) {
                    lfsr.reset(1);
                    final int q = N / d;
                    final Sequence[] seqs = new Sequence[q];
                    for(int i = 0; i < q; i++) {
                        seqs[i] = new Sequence(2 * d);
                    }
                    for(int i = 0; i < 2 * d; i++) {
                        for(int j = 0; j < q; j++) {
                            seqs[j].setElement(i, lfsr.getBit(0));
                            lfsr.clock();
                        }
                    }
                    final Polynomial alpha = g.getPrimitiveRoot();
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
                        System.out.println("n:" + n + "\tf:" + prim.toBinary() + "\td:" + d + "\tq:" + q + "\tg:" + g.toBinary());
                        for(int i = 0; i < q; i++) {
                            boolean flag = false;
                            for(int j = 0; j < q; j++) {
                                if(seqs[i].contains(classSeqs.get(j))) {
                                    System.out.println(" s" + i + "=" + seqs[i].getSubSequence(0, d).toString() + ":" + classReps.get(j).toString() + " >> " + seqs[i].find(classSeqs.get(j)));
                                    flag = true;
                                    break;
                                }
                            }
                            if(flag) {
                                continue;
                            }
                        }
                    }
                    g = PolynomialUtils.getStrictIrreducible(n, d, g.nextPoly().toInt());
                }
            }
        }
    }
}
