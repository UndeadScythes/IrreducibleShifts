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
        final int low = ArgUtils.getInt(args, "-l", 0);
        final int high = ArgUtils.getInt(args, "-h", 0);
        final boolean optimise = ArgUtils.getBool(args, "-o");
        //test1(low, high, optimise);
        test2(low, high);
    }

    private static void test1(final int low, final int high, final boolean optimise) throws IOException {
            for(int n = low; n <= high; n++) {
            final File file = new File("IrreducibleShifts2-" + n + ".csv");
            final BufferedWriter out = new BufferedWriter(new FileWriter(file));
            out.write("n,f,d,q,g,dec seqs,class reps,trace seqs,dec classes,dec phases,seeds,seed count\n");
            final int N = (1 << n) - 1;
            Polynomial f = PolynomialUtils.getPrimitive(n, 0);
            while(f != null) {
                final List<Integer> factorsN = IntegerUtils.getFactors(N);
                if(factorsN.isEmpty()) {
                    out.close();
                    file.delete();
                    break;
                }
                final GaloisLFSR lfsrF = new GaloisLFSR(n, f, 1);
                for(int d : factorsN) {
                    Polynomial g = PolynomialUtils.getStrictIrreducible(n, d, 0);
                    if(g == null) {
                        continue;
                    }
                    while(g != null) {
                        lfsrF.reset(1);
                        final int q = N / d; // q
                        if(optimise && q > n) {
                            break;
                        }
                        final Sequence[] mSeqDecimations = new Sequence[q];
                        for(int i = 0; i < q; i++) {
                            mSeqDecimations[i] = new Sequence(2 * d);
                        }
                        for(int i = 0; i < 2 * d; i++) {
                            for(int j = 0; j < q; j++) {
                                mSeqDecimations[j].setElement(i, lfsrF.getBit(0));
                                lfsrF.clock();
                            }
                        }
                        final Polynomial alpha = g.getPrimitiveRoot();
                        final List<Polynomial> classReps = g.getClassReps(alpha);
                        final List<Sequence> classRepTraces = new ArrayList<Sequence>();
                        for(int i = 0; i < q; i++) {
                            classRepTraces.add(Sequence.fromTrace(classReps.get(i), g, n));
                        }
                        boolean allMatch = true;
                        for(int i = 0; i < q; i++) {
                            boolean flag = false;
                            for(int j = 0; j < q; j++) {
                                if(mSeqDecimations[i].contains(classRepTraces.get(j))) {
                                    flag = true;
                                    break;
                                }
                            }
                            if(flag) {
                                continue;
                            } else {
                                allMatch = false;
                                break;
                            }
                        }
                        Polynomial[] decSeqClasses = new Polynomial[q];
                        int[] decSeqPhases = new int[q];
                        for(int i = 0; i < q; i++) {
                            boolean flag = false;
                            for(int j = 0; j < q; j++) {
                                if(mSeqDecimations[i].contains(classRepTraces.get(j))) {
                                    decSeqClasses[i] = classReps.get(j);
                                    decSeqPhases[i] = mSeqDecimations[i].find(classRepTraces.get(j));
                                    //System.out.println(" s" + i + "=" + mSeqDecimations[i].getSubSequence(0, d).toString());
                                    //out.write(mSeqDecimations[i].getSubSequence(0, d).toString() + " ");
                                    flag = true;
                                    break;
                                }
                            }
                            if(flag) {
                                continue;
                            }
                        }
                        if(allMatch) {
                            System.out.println("n:" + n + "\tf:" + f.toBinary() + "\td:" + d + "\tq:" + q + "\tg:" + g.toBinary() + "\ta:" + alpha.toBinary());
                            System.out.println(" class reps:{" + ArrayUtils.joinBinary(classReps, ", ") + "}\ttrace seqs:{" + ArrayUtils.join(classRepTraces, ", ") + "}");
                            System.out.println(" decimation classes:{" + ArrayUtils.joinBinary(decSeqClasses, ", ") + "}\tdecimation phases:{" + ArrayUtils.join(decSeqPhases, ", ") + "}");
                            out.write(n + "," + f.toBinary() + "," + d + "," + q + "," + g.toBinary() + ",");
                            out.write("," + ArrayUtils.joinBinary(classReps, " ") + "," + ArrayUtils.join(classRepTraces, " ") + "," + ArrayUtils.joinBinary(decSeqClasses, " ") + "," + ArrayUtils.join(decSeqPhases, " ") + ",");
                            Sequence[] lfsrSeqSnap;
                            Polynomial[] lfsrSeqClasses = new Polynomial[n];
                            int[] phases = new int[n];
                            int seed = 1;
                            final int[] requiredSeeds = new int[q];
                            int seedsFound = 0;
                            while(seed <= N && seedsFound != q) {
                                lfsrSeqSnap = new Sequence[n];
                                lfsrSeqClasses = new Polynomial[n];
                                phases = new int[n];
                                final GaloisLFSR lfsrG = new GaloisLFSR(n, g, seed);
                                for(int i = 0; i < n; i++) {
                                    for(int j = 0; j < n; j++) {
                                        if(i == 0) {
                                            lfsrSeqSnap[j] = new Sequence(n);
                                        }
                                        lfsrSeqSnap[j].setElement(i, lfsrG.getBit(j));
                                    }
                                    lfsrG.clock();
                                }
                                int seqsFound = 0;
                                int phase = 0;
                                while(seqsFound != n) {
                                    for(int i = 0; i < n; i++) {
                                        if(lfsrSeqClasses[i] == null) {
                                            for(int j = 0; j < q; j++) {
                                                if(lfsrSeqSnap[i].equalTo(classRepTraces.get(j))) {
                                                    lfsrSeqClasses[i] = classReps.get(j);
                                                    phases[i] = phase;
                                                    seqsFound++;
                                                    for(int k = 0; k < q; k++) {
                                                        if(requiredSeeds[k] == 0 && decSeqClasses[k].equalTo(lfsrSeqClasses[i]) && decSeqPhases[k] == phases[i]) {
                                                            requiredSeeds[k] = seed;
                                                            seedsFound++;
                                                        }
                                                    }
                                                    continue;
                                                }
                                            }
                                        }
                                    }
                                    for(int i = 0; i < n; i++) {
                                        lfsrSeqSnap[i].leftShift(1);
                                        lfsrSeqSnap[i].setElement(n - 1, lfsrG.getBit(i));
                                    }
                                    lfsrG.clock();
                                    phase++;
                                }
                                System.out.println(" seed:" + seed + "\tlfsr classes:{" + ArrayUtils.joinBinary(lfsrSeqClasses, ", ") + "}\tlfsr phases:{" + ArrayUtils.join(phases, ", ") + "}");
                                seed++;
                            }
                            System.out.println(" required seeds:{" + ArrayUtils.join(requiredSeeds, ", ") + "}" + "\tunique:" + IntegerUtils.uniques(requiredSeeds));
                            System.out.println();
                            out.write(ArrayUtils.join(requiredSeeds, " ") + "," + IntegerUtils.uniques(requiredSeeds));
                            out.newLine();
                        }
                        g = PolynomialUtils.getStrictIrreducible(n, d, g.nextPoly().toInt());
                    }
                }
                f = PolynomialUtils.getPrimitive(n, f.nextPoly().toInt());
            }
        out.close();
        }
    }

    private static void test2(final int low, final int high) throws IOException {
        for(int n = low; n <= high; n++) {
            final File file = new File("IrredTest-" + n + ".csv");
            final BufferedWriter out = new BufferedWriter(new FileWriter(file));
            out.write("degree,N,g,d,q,sample,alpha,seqs,seeds,classes,shifts\n");
            final int N = (1 << n) - 1;
            System.out.println("n:" + n + "\tN:" + N);
            Polynomial g = PolynomialUtils.getStrictIrreducible(n, 0);
            if(g == null) {
                out.close();
                file.delete();
            }
            while(g != null) {
                final int d = g.getOrder();
                final int q = N / d;
                System.out.println(" g:" + g.toString("x") + "\td:" + d + "\tq:" + q);
                Polynomial alpha = g.getPrimitiveRoot();
                GaloisLFSR lfsr = new GaloisLFSR(n, g, 1);
                final Sequence[] sampleSeqs = new Sequence[n];
                long start = System.currentTimeMillis();
                for(int i = 0; i < d * 2; i++) {
                    for(int j = 0; j < n; j++) {
                        if(i == 0) {
                            sampleSeqs[j] = new Sequence(d * 2);
                        }
                        sampleSeqs[j].setElement(i, lfsr.getBit(j));
                    }
                    lfsr.clock();
                }
                System.out.println(" sample:{" + ArrayUtils.join(sampleSeqs, d, ", ") + "}\ttime:" + (System.currentTimeMillis() - start) + "ms");
                while(alpha != null) {
                    System.out.println("  a:" + alpha.toString("b"));
                    int seed = 1;
                    final int[] seeds = new int[q];
                    final Sequence[] seqs = new Sequence[q];
                    start = System.currentTimeMillis();
                    for(int seedNo = 0; seedNo < q; seedNo++) {
                        seeds[seedNo] = seed;
                        lfsr = new GaloisLFSR(n, g, seed);
                        seed = 0;
                        for(int i = 0; i < d; i++) {
                            if(i < n) {
                                seed ^= alpha.getCoeff(i) * lfsr.getState();
                            }
                            if(i == 0) {
                                seqs[seedNo] = new Sequence(d);
                            }
                            seqs[seedNo].setElement(i, lfsr.getBit(n - 1));
                            lfsr.clock();
                        }
                    }

                    final int[] shiftSeqs = new int[n];
                    final int[] shifts = new int[n];
                    for(int i = 0; i < n; i++) {
                        for(int j = 0; j < q; j++) {
                            final int temp = sampleSeqs[i].find(seqs[j]);
                            if(temp > -1) {
                                shiftSeqs[i] = j;
                                shifts[i] = (d - temp) % d;
                            }
                        }
                    }
                    System.out.println("   seqs:{" + ArrayUtils.join(seqs, d, ", ") + "}" + "\tseeds:{" + ArrayUtils.joinBinary(seeds, n, ", ") + "}");
                    System.out.println("   classes:{" + ArrayUtils.join(shiftSeqs, ", ") + "}\tshifts(<<):{" + ArrayUtils.join(shifts, ", ") + "}\ttime:" + (System.currentTimeMillis() - start) + "ms");

                    final Sequence mSeq = new Sequence(N);
                    for(int i = 0; i < d; i++) {
                        for(int j = 0; j < q; j++) {
                            mSeq.setElement(i * q + j, seqs[j].getElement(i));
                        }
                    }
                    System.out.println("   interleaved:" + mSeq.toString() + "\tf:" + mSeq.getMinimal().toString("x"));

                    out.write(n + "," + N + "," + g.toBinary() + "," + d + "," + q + "," + ArrayUtils.join(sampleSeqs, n, " ") + "," + alpha.toBinary() + "," + ArrayUtils.join(seqs, n, " ") + "," + ArrayUtils.joinBinary(seeds, n, " ") + "," + ArrayUtils.join(shiftSeqs, " ") + "," + ArrayUtils.join(shifts, " ") + "\n");

                    alpha = g.getPrimitiveRoot(alpha);
                }
                g = PolynomialUtils.getStrictIrreducible(n, g.nextPoly().toInt());
            }
            out.close();
        }
    }
}
