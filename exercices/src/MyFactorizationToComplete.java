
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
    public static void main(String[] args) {
        try {
            // Reed-Solomon code of length 8 and dimension 2 over field gf8
            // See also the end of the file
            // k is the dimension of the code
            int k = 2;
            // Build the field GF8
            GaloisField gf2 = new RootGaloisField(2);
            ExtendedGaloisField gf8 = new ExtendedGaloisField(gf2, 3);
            // Acces the elements of the field GaloisField.Element using the attribute element
            // element is an array of GaloisField.Element
            GaloisField.Element[] elementsGF8 = gf8.element;

            // Vector to store the code
            Vector<GFVector> code = new Vector();

            // Coefficients of the ``message polynomial''
            GaloisField.Element[] coefs = new GaloisField.Element[k];
            GaloisField.Element[] codewordVec = new GaloisField.Element[elementsGF8.length];
            int first_degree_coef = -1;
            for (int i=0; i<64; i++) {
               // iterate over all 64 different input polynomials (including 0)
               coefs[0] = gf8.element[i%8];
               if (i%8 == 0) first_degree_coef++;
               coefs[1] = gf8.element[first_degree_coef];
               GFPolynomial polM = new GFPolynomial(coefs, gf8);
               for (int j=0; j<8; j++) {
                   codewordVec[j] = elementsGF8[j].mul(coefs[1]).add(coefs[0]); // For k=2 each output letter should be coefs[1]*elementsGF8[i]+coefs[0]
               }
               GFVector codeword = new GFVector(codewordVec,gf8);
               code.add(i, codeword);
               //System.out.println(String.format("%15s %5s", polM + " ->", codeword));
           }
            // Error correction
            GFVector c1 = code.elementAt(10);
            //System.out.println("codeword : " + c1);
            // Introduce 3 errors
            GFVector c2 = c1.copy();
            c2.setElementAt(c1.elementAt(2).add(gf8.oneElement()), 2);
            c2.setElementAt(c1.elementAt(4).add(gf8.oneElement()), 4);
            c2.setElementAt(c1.elementAt(6).add(gf8.oneElement()), 6);

            //System.out.println("corrupted codeword (3 errors) : " + c2);

            //  Check how many codewords at distance 3
            for (GFVector c : code) {
                if (c.dist(c2) <= 3) {
                    //System.out.println("decoded codeword (dist. 3): " + c);
                }
            }
            // Introduce two more errors
            c2.setElementAt(c1.elementAt(0).add(gf8.oneElement()), 4);
            c2.setElementAt(c1.elementAt(7).add(gf8.oneElement()), 6);

            //System.out.println("corrupted codeword (5 errors) : " + c2);
            //  Check how many codewords at distance 5
            for (GFVector c : code) {
                if (c.dist(c2) <= 5) {
                    //System.out.println("decoded codeword (dist. 5): " + c);
                }
            }
            // Decoding 

            int m = 7;
            RSInterpolationPoint[] interpolPoints = new RSInterpolationPoint[elementsGF8.length];
            RSInterpolationPointSet interpolSet = new RSInterpolationPointSet(gf8);
            for (int t = 0; t < elementsGF8.length; t++) {
                interpolPoints[t] = new RSInterpolationPoint(elementsGF8[t],
                        c2.elementAt(t), m);
                try {
                    interpolSet.addPoint(interpolPoints[t]);
                } catch (BlockCodeException ex) {
                    Logger.getLogger(MainReedSolomonToComplete.class.getName()).log(Level.SEVERE, null, ex);
                }
            }

            int cost = gf8.cardinality() * IntegerUtils.comb(m + 1, 2);

            int dy = KoetterVardy.dy_revlex(cost, k);
            System.out.println("dy : " + dy);

            GFBiPolynomial p = KoetterVardy.interpolate(interpolSet, dy, k, 0);

            System.out.println("pol interpolated =\n" + p + "\n");

            //Vector<GFPolynomial> v = KoetterVardy.reconstruct(p, 2);
            //Vector<GFPolynomial> v = MyFactorizationToComplete.factorization(p, 2);

            //System.out.println("Factors=\n" + v + "\n");

            brackets(p, gf8);

        } catch (GFException ex) {
        } catch (KoetterVardyException ex) {
        }

    }
    public static Vector<GFPolynomial> factorization(GFBiPolynomial Q, int k)
                     throws KoetterVardyException {

             GaloisField.Element[] coefs = new GaloisField.Element[k];
             Vector<GFPolynomial> factors = new Vector<GFPolynomial>();

             _factor(Q, k, 0, coefs, factors);

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
                if (i == k - 1) {

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
                    
                    
                    
                    // Second: 
                    // M_change_xy(x,y) <- M_shift(x,xy);
                    // You can use method eval and xBiPoly, xyBiPoly


                    // Third
                    // To compute
                    //M(x,y) <- <<M(x, xy+gamma)>>
                    // Use M and brackets method
                    

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
                
                
                
                // Second
                // M_change_xy(x,y) <- M_shift(x,xy);
                // You can use method eval and xBiPoly, xyBiPoly
                
                
                
                
                
                // _factor(M_{i+1}(x,y), k, i+1 coefs, factors);
                // Make recursive call to _factor
                

                }
            }
            

            
            
            // Remove this line after completing
            throw new RuntimeException("To complete");
            
            
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
        // Use methods xEval and mDiv
        System.out.println(Q.mDiv(GF.oneElement(), Q.multiDegree()));
        return Qn;
    }
}
