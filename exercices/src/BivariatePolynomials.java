import java.util.logging.Level;
import java.util.logging.Logger;

import jmoreira.pfc.galois.ExtendedGaloisField;
import jmoreira.pfc.galois.GFBiPolynomial;
import jmoreira.pfc.galois.GFException;
import jmoreira.pfc.galois.GaloisField;
import jmoreira.pfc.galois.RootGaloisField;

public class BivariatePolynomials {
    public static void main(String[] args) {
        try {
            GaloisField gf2 = new RootGaloisField(2);
            ExtendedGaloisField gf8 = new ExtendedGaloisField(gf2, 8);
            // Acces the elements of the field GaloisField.Element using the attribute element
            // element is an array of GaloisField.Element
            GaloisField.Element[] elementsGF8 = gf8.element;
            GFBiPolynomial oneBiPoly = new GFBiPolynomial(new GaloisField.Element[][]{{gf8.oneElement()}}, gf8);
            System.out.println("oneBiPoly: " + oneBiPoly);
            GFBiPolynomial xBiPoly = new GFBiPolynomial(new GaloisField.Element[][]{{gf8.zeroElement()}, {gf8.oneElement()}}, gf8);
            System.out.println("xBiPoly: " + xBiPoly);
            /*GFBiPolynomial[] bip = new GFBiPolynomial[4];
            for (int l = 0; l <= 5; ++l) {
                bip[l] = oneBiPoly.mMul(gf8.oneElement(), 5-l, l);
                System.out.println("bip: " + bip[l]);
            }
            bip[0] = oneBiPoly.mMul(elementsGF8[43], 5, 3);
            System.out.println("bip: " + bip[0]);
            bip[1] = oneBiPoly.mMul(elementsGF8[22], 3, 6);
            System.out.println("bip: " + bip[1]);
            bip[2] = oneBiPoly.mMul(elementsGF8[14], 3, 3);
            System.out.println("bip: " + bip[2]);
            bip[3] = oneBiPoly.mMul(elementsGF8[16], 0, 0);
            System.out.println("bip: " + bip[3]);

            GFBiPolynomial[] bip2 = new GFBiPolynomial[1];
            bip2[0] = oneBiPoly.mMul(elementsGF8[48], 3, 6);
            System.out.println("bip2: " + bip2[0]);

            GFBiPolynomial res = bip[0].add(bip[1]);
            res = res.add(bip[2]);
            res = res.add(bip[3]);*/
            GaloisField.Element[][] bip_coef = new GaloisField.Element[6][7];
            for (int i = 0; i <= 5; i++) {
                for (int j = 0; j <= 6; j++) {
                    bip_coef[i][j] = gf8.zeroElement();
                }
            }
            bip_coef[0][0] = elementsGF8[15];
            bip_coef[3][3] = elementsGF8[13];
            bip_coef[3][6] = elementsGF8[21];
            bip_coef[5][3] = elementsGF8[42];
            GFBiPolynomial bip = new GFBiPolynomial(bip_coef, gf8);
            GaloisField.Element[][] bip2_coef = new GaloisField.Element[6][7];
            for (int i = 0; i <= 5; i++) {
                for (int j = 0; j <= 6; j++) {
                    bip2_coef[i][j] = gf8.zeroElement();
                }
            }
            bip2_coef[3][6] = elementsGF8[47];
            GFBiPolynomial bip2 = new GFBiPolynomial(bip2_coef, gf8);
            GFBiPolynomial res = bip.add(bip2);
            System.out.println(res.toString());
        } catch (GFException ex) {
            Logger.getLogger(MainReedSolomonToComplete.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
