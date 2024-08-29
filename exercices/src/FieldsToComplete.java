import java.util.Vector;
import jmoreira.pfc.galois.ExtendedGaloisField;
import jmoreira.pfc.galois.GFPolynomial;
import jmoreira.pfc.galois.GFPolynomialFactorSet;
import jmoreira.pfc.galois.GaloisField;
import jmoreira.pfc.galois.RootGaloisField;

public class FieldsToComplete {

    public static java.io.PrintStream output;

    /**
     * M�todo principal.
     *
     * @param args argumentos pasados al programa.
     * @throws Exception Si se produce alguna excepci�n.
     */
    public static void main(String[] args) throws Exception {

        /* Reference to pages and exercises are from:
            Finite Fields for Computer Scientists and Engineers
            Authors: McEliece, Robert J.
         */
        

        // output.print and output.println show content! 
        output = System.out;

        output.println("Implement the addition and multiplication tables of the field of 5 elements. "
                + "Compute the inverses and check them using the previous multiplication table.");
        GaloisField gf5 = new RootGaloisField(5);

        GaloisField.Element elemP5 = gf5.primitiveElement();
        output.println(elemP5);

        GaloisField.Element[] elems5 = gf5.element;

        // add table                
        output.print("\t");
        for (GaloisField.Element e : elems5) {
            output.print(e + "\t");
        }
        output.println();

        // complete
        output.println("addition table");
        output.print("\t");
        for (GaloisField.Element e1 : elems5) {
            
            output.print(e1 + "\t");
        }
        output.println();
        for(GaloisField.Element e2 : elems5){
            output.print(e2 + "\t");
            for (GaloisField.Element e1 : elems5) {
                output.print(e2.add(e1) + "\t");
            }
            output.println();

        }
                
        
        
        // mult table 
        output.println("mult table");
        output.print("\t");
        for (GaloisField.Element e1 : elems5) {
            
            output.print(e1 + "\t");
        }
        output.println();
        for(GaloisField.Element e2 : elems5){
            output.print(e2 + "\t");
            for (GaloisField.Element e1 : elems5) {
                output.print(e2.mul(e1) + "\t");
            }
            output.println();

        }        
        // complete
        
        
        // Compute inverses and check them using the mult table
        // complete

        output.println("inverses table");
        output.println("\t");
        for (GaloisField.Element e1 : elems5) {
            if(!e1.isZeroElement()){           
                output.println(e1 + "\t" + e1.inv());
            }
        }
        output.println();        


        output.println("Implement the addition and multiplication tables of the field of 4 elements. "
                + "Compute the inverses and check them using the previous multiplication table.");

//		GaloisField gf3 = new RootGaloisField(4);
        GaloisField gf2 = new RootGaloisField(2);
        GaloisField gf4 = new ExtendedGaloisField(gf2, 2);
        GaloisField.Element elemP = gf4.primitiveElement();

        output.println(elemP);

        GaloisField.Element[] elems4 = gf4.element;

        // add table                
        output.print("\t");
        for (GaloisField.Element e : elems4) {
            output.print(e + "\t");
        }
        output.println();

        // complete
        output.println("addition table");
        output.print("\t");
        for (GaloisField.Element e1 : elems4) {
            
            output.print(e1 + "\t");
        }
        output.println();
        for(GaloisField.Element e2 : elems4){
            output.print(e2 + "\t");
            for (GaloisField.Element e1 : elems4) {
                output.print(e2.add(e1) + "\t");
            }
            output.println();

        }
                
        
        
        // mult table 
        output.println("mult table");
        output.print("\t");
        for (GaloisField.Element e1 : elems4) {
            
            output.print(e1 + "\t");
        }
        output.println();
        for(GaloisField.Element e2 : elems4){
            output.print(e2 + "\t");
            for (GaloisField.Element e1 : elems4) {
                output.print(e2.mul(e1) + "\t");
            }
            output.println();

        }        
        // complete
        
        
        // Compute inverses and check them using the mult table
        // complete

        output.println("inverses table");
        output.println("\t");
        for (GaloisField.Element e1 : elems4) {
            if(!e1.isZeroElement()){           
                output.println(e1 + "\t" + e1.inv());
            }
        }
        output.println();        



        // gcd Corollary 2.4 pag 7
        // Exercise warm up
        output.println("");
        output.println("gcd Corollary 2.4 pag 7");
        output.println("Check that $gcd(x^{15}-1,x^{20}-1)=x^{5}-1$. ");
        output.println("");


        // Complete 
GaloisField.Element[] coefsP1 = {gf2.oneElement(),gf2.zeroElement(),gf2.oneElement()};
                GFPolynomial p1 = new GFPolynomial(coefsP1, gf2);
                output.println("p1(x) = " + p1);
                
                GaloisField.Element[] coefsP2 = {gf2.oneElement(),gf2.zeroElement(),gf2.zeroElement(),gf2.oneElement()};
                GFPolynomial p2 = new GFPolynomial(coefsP2, gf2);
                output.println("p2(x) = " + p2); 
                
                GFPolynomial pGCD = p1.gcd(p2);
                output.println("gcd(p1(x),p2(x)) = " + pGCD);

                output.println("");
                output.println("Example pag 7");
           
                output.println("");                
                // Example pag 7
                GaloisField.Element[] coefsP3 = {gf2.oneElement(),gf2.zeroElement(),gf2.zeroElement(),gf2.zeroElement(),gf2.zeroElement(),gf2.zeroElement(),
                                                                    gf2.zeroElement(),gf2.zeroElement(),gf2.zeroElement(),gf2.zeroElement(),gf2.zeroElement(),
                                                                    gf2.zeroElement(),gf2.zeroElement(),gf2.zeroElement(),gf2.zeroElement() ,gf2.oneElement()};
                GFPolynomial p3 = new GFPolynomial(coefsP3, gf2);
                output.println("p3(x) = " + p3);
                
                GaloisField.Element[] coefsP4 = {gf2.oneElement(),gf2.zeroElement(),gf2.zeroElement(),gf2.zeroElement(),gf2.zeroElement(),gf2.zeroElement(),
                                                                    gf2.zeroElement(),gf2.zeroElement(),gf2.zeroElement(),gf2.zeroElement(),gf2.zeroElement(),
                                                                    gf2.zeroElement(),gf2.zeroElement(),gf2.zeroElement(),gf2.zeroElement() ,
                                                                    gf2.zeroElement(),gf2.zeroElement(),gf2.zeroElement(),gf2.zeroElement(),gf2.zeroElement(),gf2.oneElement()};
                GFPolynomial p4 = new GFPolynomial(coefsP4, gf2);
                output.println("p4(x) = " + p4); 
                
                GFPolynomial pGCD1 = p3.gcd(p4);
                output.println("gcd(p3(x),p4(x)) = " + pGCD1);            
        
        
        
        // Factorization Exercise 9 page 18
        output.println("");
        output.println("Factorization Exercise 9 page 18");
        output.println("In $F_2[x]$, factor the following polynomials $x^2+1, x^2+x+1,x^2+x$.");
        output.println("");
        GaloisField.Element[] coefsP5 = {gf2.oneElement(), gf2.oneElement(), gf2.oneElement()};
        GFPolynomial p5 = new GFPolynomial(coefsP5, gf2);
        output.println("p5(x) = " + p5);
        Vector<GaloisField.Element> rootsP5 = p5.roots();
        for (GaloisField.Element e : rootsP5) {
            output.println("root p5(x) = " + e);
        }

        GFPolynomialFactorSet factorsP5 = p5.squareFree();
        output.println("square free p5(x) = " + factorsP5);

        // Complete
        
        
        
        
        
        // division algorithm in F_p[x] page 24 example 4.7
        output.println("");
        output.println("division algorithm in F_p[x] page 24 example 4.7");
        output.println("In $F_{13}[x]$, let $a(x)=x^8+x^6+10x^4+10x^3+8x^2+2x+8$ "
                + "and $b(x)=3x^6+5x^4+9x^2+4x+8$ Find two polynomials $q(x)$ "
                + "and $r(x)$ such that\n"
                + "  $$a(x)=q(x)b(x)+r(x).$$");
        output.println("");
        GaloisField gf13 = new RootGaloisField(13);
        GaloisField.Element[] elems13 = gf13.element;
        GaloisField.Element[] coefsPa = {elems13[8], elems13[2], elems13[8], elems13[10], elems13[10], elems13[0], elems13[1], elems13[0], elems13[1]};
        GFPolynomial pa = new GFPolynomial(coefsPa, gf13);
        output.println("pa(x) = " + pa);

        GaloisField.Element[] coefsPb = {elems13[8], elems13[4], elems13[9], elems13[0], elems13[5], elems13[0], elems13[3]};
        GFPolynomial pb = new GFPolynomial(coefsPb, gf13);
        output.println("pb(x) = " + pb);

        // Complete
        
        
        
        
        
        // Exercise 11 page 28
        output.println("");
        output.println("Exercise 11 page 28");
        output.println("Given  $p(x)=x^2-x-1$  in $F_3[x]$, construct the field $F_3[x] \\mod p(x)$, "
                + "by giving the addition and multiplication tables. "
                + "Do the same with $p(x)=x^2 + 2$. "
                + "What do you observe?");
        output.println("");
        GaloisField gf3 = new RootGaloisField(3);
        GaloisField.Element[] elems3 = gf3.element;
        GaloisField.Element[] coefsPF3 = {elems3[2], elems3[2], gf3.oneElement()};
        GFPolynomial pf3 = new GFPolynomial(coefsPF3, gf3);
        ExtendedGaloisField gf9 = new ExtendedGaloisField(pf3);

        output.println("pf3(x) = " + pf3);

        GaloisField.Element[] coefsX3 = {gf3.zeroElement(), gf3.oneElement()};
        GFPolynomial px3 = new GFPolynomial(coefsX3, gf3);

        // add table   
        
        
        // mult table 
        
                GaloisField.Element[] coefsPF3b = {elems3[2],gf3.zeroElement(),gf3.oneElement()};  
                GFPolynomial pf3b = new GFPolynomial(coefsPF3b, gf3); 
//                ExtendedGaloisField gf9b = new ExtendedGaloisField(pf3b);
                
                output.println("pf3(x) = " + pf3b);  
                
                GaloisField.Element[] coefsX3b = {gf3.zeroElement(),gf3.oneElement()};  
                GFPolynomial px3b = new GFPolynomial(coefsX3b, gf3);  

                
                GFPolynomial[] elemsgf9Pb = new GFPolynomial[8];
                for(int i = 0; i < 8 ; i++){
                    elemsgf9Pb[i] = px3b.pow(i).mod(pf3b);
                    output.println("x^"+i+ " = " + elemsgf9Pb[i]);
                } 
                
                // add table                
                output.print("\t");
                for(GFPolynomial e : elemsgf9Pb){
                    output.print(e + "\t");
                }
                output.println();
                

                for(GFPolynomial e1 : elemsgf9Pb){
                    output.print(e1.toString() + "\t");
                    for(GFPolynomial e2 : elemsgf9Pb){
                        output.print(e1.add(e2) + "\t");
                    }
                    output.println();
                }
                output.println();
                
                // mult table 
                output.print("\t");
                for(GFPolynomial e : elemsgf9Pb){
                    output.print(e + "\t" );
                }
                output.println();                
                for(GFPolynomial e1 : elemsgf9Pb){
                    output.print(e1.toString() + "\t");
                    for(GFPolynomial e2 : elemsgf9Pb){
                        output.print((e1.mul(e2)).mod(pf3b) + "\t");
                    }
                    output.println();
                }
                output.println();  
                        
        
        output.println();

        // abstract propereties of finite fields
        // Theorem 5.2 pag 31 order
        output.println("");
        output.println("Theorem 5.2 pag 31 order");
        output.println("Using $F_{13}$ show that if $t$ is the order or $\\alpha$, "
                + "then $t$ divides $q-1=12$.");
        output.println("");

        // complete
               output.println("order F13");
                for(GaloisField.Element e : elems13){
                    if(!e.isZeroElement()){
                        GaloisField.Element tmp = e;
                        int i = 1;
                        output.print(e + " ");
                        while(!tmp.isOneElement()){
                            i++;
                            tmp = tmp.mul(e);
                            output.print(e + "^" + i + " = " + tmp + "  ");


                        }
                        output.println();
                    }
                }
        
        
        
        
        // Teorem 5.8 example 5.3 
        output.println("");
        output.println("Teorem 5.8 example 5.3 ");
        output.println("Using  $F_{16}$ show that there are $\\phi(t)$ "
                + "elements of order $t$, for appropriate $t$.");
        output.println("");

        ExtendedGaloisField gf16 = new ExtendedGaloisField(gf2, 4);
        GaloisField.Element[] elems16 = gf16.element;

        // Complete
        
        
        
        
        
        //  example 5.5 pag 37
        output.println("");
        output.println(" example 5.5 pag 37 ");
        output.println("In $F_7$ find a primitive root and compute its powers.");
        output.println("");
        GaloisField gf7 = new RootGaloisField(7);
        GaloisField.Element[] elems7 = gf7.element;

        // Complete     
        
        
        
        
        // Example 5.7 pag 42 minimal polynomial
        output.println("");
        output.println(" Example 5.7 pag 42 minimal polynomial ");
        output.println("in $F_{25}$ find a primitive root and its primitive polynomial.");
        output.println("");

        GaloisField.Element[] coefsPol5 = {elems5[3], gf5.zeroElement(), gf5.oneElement()};
        GFPolynomial pol5 = new GFPolynomial(coefsPol5, gf5);
        output.println("pol5(x) = " + pol5);

        // prove irreducibility
        for (GaloisField.Element e : elems5) {
            output.print("pol5(" + e + ") = " + pol5.eval(e) + "  ");
        }
        output.println();

        ExtendedGaloisField gf25 = new ExtendedGaloisField(pol5);
        GaloisField.Element[] elems25 = gf25.element;

        // Complete
    }

    public static int gcd(int a, int b) {

        int res;

        while (!(b == 0)) {
            res = a % b;
            a = b;
            b = res;
        }

        return a;
    }

}
