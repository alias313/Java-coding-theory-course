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
            ExtendedGaloisField gf8 = new ExtendedGaloisField(gf2, 3);
            // Acces the elements of the field GaloisField.Element using the attribute element
            // element is an array of GaloisField.Element
            GaloisField.Element[] elementsGF8 = gf8.element;
            GFBiPolynomial oneBiPoly = new GFBiPolynomial(new GaloisField.Element[][]{{gf8.oneElement()}}, gf8);
            System.out.println("oneBiPoly: " + oneBiPoly);
            GFBiPolynomial xBiPoly = new GFBiPolynomial(new GaloisField.Element[][]{{gf8.zeroElement()}, {gf8.oneElement()}}, gf8);
            System.out.println("xBiPoly: " + xBiPoly);
            GFBiPolynomial[] bip = new GFBiPolynomial[13];
            for (int l = 0; l <= 12; ++l) {
                bip[l] = oneBiPoly.mMul(gf8.oneElement(), 0, l);
                System.out.println("bip: " + bip[l]);
            }
        } catch (GFException ex) {
            Logger.getLogger(MainReedSolomonToComplete.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
