//import { errorSub } from "./errorCode";
export function getCubicRoot(a3: number, a2: number, a1: number, A0: number, phase: string) {
// @customfunction
    try {

        /*'This work is adapted from work created by Tomas Co, Michigan Technological Univiersity
        Reference 6
        ' Computes the minimum real root of the cubic equation
        ' a3 x**3 + a2 x**2 + a1 x + a0 = 0 */
        
        let a: number = 0
        let b: number = 0
        let c: number = 0
        let d: number = 0
        let Z: number = 0
        let Q : number = 0
        let p : number = 0
        let h : number = 0
        let Y : number = 0
        let z1 : number = 0
        let z2 : number = 0
        let z3 : number = 0
        let c1 : number = 0
        let S1 : number = 0
        let m : number = 0
        let R : number = 0
        let Disc : number = 0
        let Theta : number = 0
        let fcnName: string = "getCubeRoot"
        let rootType: string = ""
        
        let myErrorMsg: string = "getCubicRoot"
        if(phase === "vapor".toLowerCase()){
            rootType = "largest"
        }

        if(phase === "liquid".toLowerCase()){
            rootType = "smallest"
        }

        if(phase === "") {
            return -500
        }

        if(a3 === 0) {
            myErrorMsg = "getCubicRoot error: The first term a3 equals zero!"
            throw new Error(`Function name: ${fcnName}, ${myErrorMsg}`)
        }
        
        a = a2 / a3
        b = a1 / a3
        c = A0 / a3
        p = ((-1)*a ** 2 / 3 + b) / 3
        
        if((9 * a * b - 2 * a ** 3 - 27 * c) === 0) {
            myErrorMsg = "getCubicRoot error: (9 * A * B - 2 * A ** 3 - 27 * C) equals zero!"
            throw new Error(`Function name: ${fcnName}, ${myErrorMsg}`)
        }
        
        Q = (9 * a * b - 2 * a ** 3 - 27 * c) / 54
        Disc = Q ** 2 + p ** 3
        
        if(Disc > 0) {
            h = Q + Disc ** (1 / 2)
            
            if(Math.abs(h)**(1 / 3) === 0) {
                myErrorMsg = "getCubicRoot error: (Abs(h)) ** (1 / 3) equals zero!"
                throw new Error(`Function name: ${fcnName}, ${myErrorMsg}`)
            }
            
            Y = (Math.abs(h)) ** (1 / 3)
            if(h < 0) {
                Y = -Y
            }
            Z = Y - p / Y - a / 3
            
        } else {
            Theta = Math.atan((-Disc) ** (1 / 2) / Q)
            c1 = Math.cos(Theta / 3)
            
            if(Q < 0) {
                S1 = Math.sin(Theta / 3)
                c1 = (c1 - S1 * 3 ** (1 / 2)) / 2
            }

            z1 = 2 * (-p) ** (1 / 2) * c1 - a / 3
            m = a + z1
            R = (m ** 2 - 4 * (b + m * z1)) ** (1 / 2)
            z2 = (-m + R) / 2
            z3 = (-m - R) / 2
            Z = z1

            if(rootType === "largest") {
                if(z2 > Z &&  z2 > 0) {
                    Z = z2
                }

                if(z3 > Z && z3 > 0) {
                    Z = z3
                }
            }

            if(rootType === "smallest") {
                if(z2 < Z &&  z2 > 0) {
                    Z = z2
                }

                if(z3 < Z && z3 > 0) {
                    Z = z3
                }
            }

        }
        
        return Z
        
     
    }

    catch(myErrorHandler) {
    return -500
    
    }
}