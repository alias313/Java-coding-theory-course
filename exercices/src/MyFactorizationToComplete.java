
import java.util.Vector;
import jmoreira.pfc.blockcode.reedsolomon.koettervardy.KoetterVardyException;
import jmoreira.pfc.galois.GFBiPolynomial;
import jmoreira.pfc.galois.GFException;
import jmoreira.pfc.galois.GFPolynomial;
import jmoreira.pfc.galois.GaloisField;

import java.util.logging.Level;
import java.util.logging.Logger;
import jmoreira.pfc.IntegerUtils;
import jmoreira.pfc.blockcode.BlockCodeException;
import jmoreira.pfc.blockcode.reedsolomon.RSInterpolationPoint;
import jmoreira.pfc.blockcode.reedsolomon.RSInterpolationPointSet;
import jmoreira.pfc.blockcode.reedsolomon.koettervardy.KoetterVardy;
import jmoreira.pfc.galois.ExtendedGaloisField;
import jmoreira.pfc.galois.GFBiPolynomial;
import jmoreira.pfc.galois.GFException;
import jmoreira.pfc.galois.GFPolynomial;
import jmoreira.pfc.galois.GFVector;
import jmoreira.pfc.galois.GaloisField;
import jmoreira.pfc.galois.RootGaloisField;


public class MyFactorizationToComplete {
    public static Vector<GFPolynomial> factorization(GFBiPolynomial Q, int k)
                     throws KoetterVardyException {

        GaloisField.Element[] coefs = new GaloisField.Element[k];
        Vector<GFPolynomial> factors = new Vector<GFPolynomial>();

        _factor(Q, k, 0, coefs, factors);
        
        System.out.print("COEFS: ");
        for (GaloisField.Element e: coefs) System.out.print(e + ", ");
        System.out.println();

        return factors;
    }
    
    // To complete _factor
    // check additional notes about Roth-Ruckenstein factoring
    private static void _factor(GFBiPolynomial Q, int k, int i,
            GaloisField.Element[] coefs, Vector<GFPolynomial> factors)
            throws KoetterVardyException {
        try {
            GFBiPolynomial M, M_shift, M_change_xy, xBiPoly, xyBiPoly;
            int numRoots;
            GaloisField field;
            GaloisField.Element gamma;
            GFPolynomial zeroPoly;
            // check additional notes about Roth-Ruckenstein factoring
            

            field = Q.field();

            zeroPoly = new GFPolynomial(new GaloisField.Element[]{field.zeroElement()}, field);

            xBiPoly = new GFBiPolynomial(new GaloisField.Element[][]{
                new GaloisField.Element[]{field.zeroElement()},
                new GaloisField.Element[]{field.oneElement()}}, field);

            xyBiPoly = new GFBiPolynomial(new GaloisField.Element[][]{
                new GaloisField.Element[]{field.zeroElement(),
                    field.zeroElement()},
                new GaloisField.Element[]{field.zeroElement(),
                    field.oneElement()}}, field);

            // M(x,y) = <<Q(x,y)>> <- Q(x,y)/x^r
            // use brackets function below.
            // that also needs to be implemented
            M = brackets(Q,field);
            System.out.println("M polynomial iteration " + i + ": " + M + "\n");
	    // find all the roots in F of the univariate polynomial
            // M(0,y);
            Vector<GaloisField.Element> roots = (M.xEval(field.zeroElement())).roots();

            //  for each of the distinct roots gamma of M(0,y) do {
            numRoots = roots.size();
            for (int j = 0; j < numRoots; j++) {
                // coefs[i] = gamma;
                // Do it in 2 steps
                gamma = roots.elementAt(j);
                coefs[i] = gamma;

                // if i = k-1 then output coefs[0],...,coefs[k-1];
                if (i == (k - 1)) {

                    //McElice Corollary 12 pag 35. Exit condition
                    // If Q_k(x,0) = 0 then
                    // then f(x) = coefs[0]+coefs[1]x + ...+ coefs[k-1]x^{k-1}
                    // is an y-root of Q(x,y)
                    //
                    // Note that at this stage
                    // Q_k(x,y) will be the following polynomial:
                    //
                    // Q_k(x,y) = M(x,y) <- <<M(x, xy+gamma)>>
                    //
                    // To compute this M(x,y)  
                    // Use polynomials M, M_shift, M_change_xy 
                    // Do it in 3 steps
                    // First use method shift to
                    // M_shift(x,y) <- M(x, y+gamma);
                    M_shift = M.shift(field.zeroElement(), gamma);                
                    
                    // Second: 
                    // M_change_xy(x,y) <- M_shift(x,xy);
                    // You can use method eval and xBiPoly, xyBiPoly
                    M_change_xy = M_shift.eval(xBiPoly, xyBiPoly);

                    // Third
                    // To compute
                    //M(x,y) <- <<M(x, xy+gamma)>>
                    // Use M and brackets method
                    M = brackets(M_change_xy, field);

                    GFPolynomial eval = M.yEval(field.zeroElement());
                    if (eval.equals(zeroPoly)) {
                        factors.addElement(new GFPolynomial(coefs, field));
                    }
                } else {
                    
                // M_{i+1}(x, y) <- M(x, xy+gamma)
                // Use polynomials M, M_change_xy, M_shift
                // Do it in 2 steps
                // First use method shift to
                // M_shift(x,y) <- M(x, y+gamma);
                M_shift = M.shift(field.zeroElement(), gamma);                
                
                // Second
                // M_change_xy(x,y) <- M_shift(x,xy);
                // You can use method eval and xBiPoly, xyBiPoly
                M_change_xy = M_shift.eval(xBiPoly, xyBiPoly);
                
                // _factor(M_{i+1}(x,y), k, i+1 coefs, factors);
                // Make recursive call to _factor
                _factor(M_change_xy, k, i+1, coefs, factors);

                }
            }            
            
        } catch (Exception e) {
            throw new KoetterVardyException("Factorization fails.",
                    KoetterVardyException.FACTORIZATION, e);
        }

    }
    
    static GFBiPolynomial brackets ( GFBiPolynomial Q, GaloisField GF) throws GFException{
        GFBiPolynomial Qn = new GFBiPolynomial(new GaloisField.Element[][]{
            new GaloisField.Element[]{GF.zeroElement(),
                GF.zeroElement()},
            new GaloisField.Element[]{GF.zeroElement(),
                GF.oneElement()}}, GF);
        //System.out.println("Input poly: " + Q);
        // First check if degree of x is more than 0
        if (Q.isZeroPoly()) return Q;
        // Evaluate x part of bivariate polynomial, if 0 divide by x, else return Qn;
        Qn = Q.copy();

        while (Qn.xEval(GF.zeroElement()).isZeroPoly()) {
            Qn = Qn.mDiv(GF.oneElement(), 1, 0);
        }
        //System.out.println("Normalized poly: " + Qn);
        //System.out.println("Zero evaluation: " + Qn.xEval(GF.zeroElement()));
        System.out.println("Roots of Zero eval: " + Qn.xEval(GF.zeroElement()).roots());
        return Qn;
    }
}
