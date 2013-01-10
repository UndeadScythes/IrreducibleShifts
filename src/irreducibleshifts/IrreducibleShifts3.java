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
        for(int n = low; n <= high; n++) {
            final BufferedWriter out = new BufferedWriter(new FileWriter("IrreducibleShifts2-" + n + ".csv"));
            out.write("n,f,d,q,g,dec seqs,class reps,trace seqs,dec classes,dec phases,lfsr seqs, lfsr phases\n");
            final int N = (1 << n) - 1;
            Polynomial f = PolynomialUtils.getPrimitive(n, 0);
            while(f != null) {
                final List<Integer> factorsN = IntegerUtils.getFactors(N);
                if(factorsN.isEmpty()) {
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
                        if(allMatch) {
                            final Sequence[] lfsrSeqSnap = new Sequence[n];
                            final Polynomial[] lfsrSeqClasses = new Polynomial[n];
                            final int[] phases = new int[n];
                            final GaloisLFSR lfsrG = new GaloisLFSR(n, g, 1);
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
                            System.out.println("n:" + n + "\tf:" + f.toBinary() + "\td:" + d + "\tq:" + q + "\tg:" + g.toBinary() + "\ta:" + alpha.toBinary());
                            out.write(n + "," + f.toBinary() + "," + d + "," + q + "," + g.toBinary() + ",");
                            Polynomial[] decSeqClasses = new Polynomial[q];
                            int[] decSeqPhases = new int[q];
                            for(int i = 0; i < q; i++) {
                                boolean flag = false;
                                for(int j = 0; j < q; j++) {
                                    if(mSeqDecimations[i].contains(classRepTraces.get(j))) {
                                        decSeqClasses[i] = classReps.get(j);
                                        decSeqPhases[i] = mSeqDecimations[i].find(classRepTraces.get(j));
                                        System.out.println(" s" + i + "=" + mSeqDecimations[i].getSubSequence(0, d).toString());
                                        //out.write(mSeqDecimations[i].getSubSequence(0, d).toString() + " ");
                                        flag = true;
                                        break;
                                    }
                                }
                                if(flag) {
                                    continue;
                                }
                            }
                            out.write("," + ArrayUtils.joinBinary(classReps, " ") + "," + ArrayUtils.join(classRepTraces, " ") + "," + ArrayUtils.joinBinary(decSeqClasses, " ") + "," + ArrayUtils.join(decSeqPhases, " ") + "," + ArrayUtils.joinBinary(lfsrSeqClasses, " ") + "," + ArrayUtils.join(phases, " ") + "\n");
                            System.out.println(" " + "class reps:{" + ArrayUtils.joinBinary(classReps, ", ") + "}\ttrace seqs:{" + ArrayUtils.join(classRepTraces, ", ") + "}");
                            System.out.println(" " + "sample classes:{" + ArrayUtils.joinBinary(lfsrSeqClasses, ", ") + "}\tsample phases:{" + ArrayUtils.join(phases, ", ") + "}");
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
