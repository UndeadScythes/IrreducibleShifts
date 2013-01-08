package irreducibleshifts;

import java.io.*;
import java.util.*;
import udslibz.algebra.*;
import udslibz.utilities.*;

/**
 * Checks the properties of irreducible sequence classes.
 * @author UndeadScythes
 */
public final class IrreducibleShifts3 {
    private IrreducibleShifts3() {}

    public static void main(final String[] args) throws IOException {
        final int lowN = ArgUtils.getInt(args, "-l", 0);
        final int highN = ArgUtils.getInt(args, "-h", 0);
        for(int n = lowN; n <= highN; n++) {
            final BufferedWriter out = new BufferedWriter(new FileWriter("IrreducibleShifts2-" + n + ".csv"));
            out.write("degree,decimal poly,binary poly,order,q,decimal alpha,binary alpha,decimal class reps,binary class reps\n");
            Polynomial g = Polynomial.getStrictIrreducible(n, 1 << n);
            while(g != null) {
                final int d = g.getOrder();
                final int q = ((1 << n) - 1) / d;
                final Polynomial alpha = g.getPrimitiveRoot();
                final List<Polynomial> classReps = g.getClassReps(alpha);
                String decimalReps = "";
                String reps = "";
                for(int i = 0; i < classReps.size(); i++) {
                    reps = reps.concat(" " + BinaryUtils.toString(classReps.get(i).toInt(), n));
                    decimalReps = decimalReps.concat(" " + classReps.get(i).toInt());
                }
                reps = reps.substring(1);
                decimalReps = decimalReps.substring(1);
                System.out.println("n:" + n + "\tg:" + g.toInt() + "\td:" + d + "\tq:" + q + "\talpha:" + alpha.toInt() + "\t{" + decimalReps.replace(" ", ", ") + "}");
                out.write(n + "," + g.toInt() + "," + g.toBinary() + "," + d + "," + q + "," + alpha.toInt() + "," + alpha.toBinary() + "," + decimalReps + "," + reps + "\n");
                g = Polynomial.getStrictIrreducible(n, g.toInt() + 1);
            }
            out.close();
        }
    }
}
