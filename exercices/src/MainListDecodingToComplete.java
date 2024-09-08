
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import jmoreira.pfc.IntegerUtils;
import jmoreira.pfc.blockcode.BlockCodeException;
import jmoreira.pfc.blockcode.reedsolomon.RSInterpolationPoint;
import jmoreira.pfc.blockcode.reedsolomon.RSInterpolationPointSet;
import jmoreira.pfc.blockcode.reedsolomon.koettervardy.KoetterVardy;
import jmoreira.pfc.blockcode.reedsolomon.koettervardy.KoetterVardyException;
import jmoreira.pfc.galois.ExtendedGaloisField;
import jmoreira.pfc.galois.GFBiPolynomial;
import jmoreira.pfc.galois.GFException;
import jmoreira.pfc.galois.GFPolynomial;
import jmoreira.pfc.galois.GFVector;
import jmoreira.pfc.galois.GaloisField;
import jmoreira.pfc.galois.RootGaloisField;



public class MainListDecodingToComplete {
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

            /**
             * ******************************************************************
             * To create the code iterate over all polynomials of degree less
             * than k and then use:
             *
             * Generate the message polynomials
             * 
             * GFPolynomial polM = new GFPolynomial(coefs, gf8);
             * 
             * Then
             * 
             * GaloisField.Element[] codewordVec = new GaloisField.Element[elements.length]; 
             * 
             * The evaluate the polynomial for all elements 
             * 
             * GFVector codeword = new GFVector(codewordVec,gf8);
             */

             GaloisField.Element[] codewordVec = new GaloisField.Element[elementsGF8.length];
             int numberOfCodewords = (int)Math.pow(elementsGF8.length, k);
             GFPolynomial[] messagePolynomialsVec = new GFPolynomial[numberOfCodewords];
             int first_degree_coef = -1;
             for (int i=0; i<numberOfCodewords; i++) {
                // iterate over all numberOfCodewords different input polynomials (including 0)
                coefs[0] = gf8.element[i%8];
                if (i%8 == 0) first_degree_coef++;
                coefs[1] = gf8.element[first_degree_coef];
                for (int j=0; j<8; j++) {
                    codewordVec[j] = elementsGF8[j].mul(coefs[1]).add(coefs[0]); // For k=2 each output letter should be coefs[1]*elementsGF8[i]+coefs[0]
                }
                GFVector codeword = new GFVector(codewordVec,gf8);
                code.add(i, codeword);
                GFPolynomial polM = new GFPolynomial(coefs, gf8);
                messagePolynomialsVec[i] = polM;
                System.out.println(String.format("%15s %5s", polM + " ->", codeword));
            }

            // Compute the distance of the code and check with the theoretical value
            System.out.println("Theoretical min distance: " + (elementsGF8.length-k+1));
            int minDist = 8;
            int dist;
            for (GFVector vec1: code) {
                for (GFVector vec2: code) {
                    if (!vec1.sub(vec2).isZeroVector()) {
                        dist = vec1.dist(vec2);
                        if (dist < minDist) minDist = dist;
                    }
                }
            }
            System.out.println("Computed min distance: " + minDist);

            // Compute the minimum weight and check with the theoretical value
            // and with the distance value
            System.out.println("Theoretical min weight: " + (elementsGF8.length-k+1));
            int minWeight = 8;
            int weight;
            for (GFVector vec: code) {
                weight = vec.weight();
                if (weight < minWeight && !vec.isZeroVector()) minWeight = weight;
            }
            System.out.println("Computed min weight: " + minWeight);
            
            // Error correction
            GFVector c1 = code.elementAt(10);
            System.out.println("codeword : " + messagePolynomialsVec[10] + "-> " + c1);
            // Introduce 3 errors
            GFVector c2 = c1.copy();
            c2.setElementAt(c1.elementAt(2).add(gf8.oneElement()), 2);
            c2.setElementAt(c1.elementAt(4).add(gf8.oneElement()), 4);
            c2.setElementAt(c1.elementAt(6).add(gf8.oneElement()), 6);

            System.out.println("corrupted codeword : " + c2);



            // Introduce two more errors
            c2.setElementAt(c1.elementAt(0).add(gf8.oneElement()), 4);
            c2.setElementAt(c1.elementAt(7).add(gf8.oneElement()), 6);

            System.out.println("corrupted codeword (5 errors) : " + c2);

            //  Check how many codewords at distance 5
            for (GFVector c : code) {
                if (c.dist(c2) <= 5) {
                    System.out.println("decoded codeword : " + c);
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
            Vector<GFPolynomial> v = MyFactorizationToComplete.factorization(p, 2);

            System.out.println("Factors: " + v);

            GFPolynomial closestPolynomial = new GFPolynomial(elementsGF8, gf8);
            int minDistanceToC2 = 8;
            for (GFPolynomial gfp : v) {
                GaloisField.Element[] ve = new GaloisField.Element[elementsGF8.length];
                for (int s = 0; s < elementsGF8.length; s++) {
                    ve[s] = gfp.eval(elementsGF8[s]);
                }
                GFVector cv = new GFVector(ve, gf8);
                if (minDistanceToC2 > c2.dist(cv)) {
                    closestPolynomial = gfp;
                    minDistanceToC2 = c2.dist(cv);
                }
                System.out.println(String.format("%15s %15s", gfp + " ->", "codeword: " + cv + " <-> dist to c2: " + c2.dist(cv)));
            }

            System.out.println("Decoded polynomial: " + closestPolynomial);

        } catch (GFException ex) {
            Logger.getLogger(MainReedSolomonToComplete.class.getName()).log(Level.SEVERE, null, ex);
        } catch (KoetterVardyException ex) {
            Logger.getLogger(MainReedSolomonToComplete.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}            
            