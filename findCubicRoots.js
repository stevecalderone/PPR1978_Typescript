"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
//import { errorSub } from "./errorCode";
function getCubicRoot(a3, a2, a1, A0, phase) {
    // @customfunction
    try {
        /*'This work is adapted from work created by Tomas Co, Michigan Technological Univiersity
        ' Computes the minimum real root of the cubic equation
        ' a3 x**3 + a2 x**2 + a1 x + a0 = 0 */
        var a = 0;
        var b = 0;
        var c = 0;
        var d = 0;
        var Z = 0;
        var Q = 0;
        var p = 0;
        var h = 0;
        var Y = 0;
        var z1 = 0;
        var z2 = 0;
        var z3 = 0;
        var c1 = 0;
        var S1 = 0;
        var m = 0;
        var R = 0;
        var Disc = 0;
        var Theta = 0;
        var fcnName = "getCubeRoot";
        var rootType = "";
        var myErrorMsg = "getCubicRoot";
        if (phase === "vapor".toLowerCase()) {
            rootType = "largest";
        }
        if (phase === "liquid".toLowerCase()) {
            rootType = "smallest";
        }
        if (phase === "") {
            return -500;
        }
        if (a3 === 0) {
            myErrorMsg = "getCubicRoot error: The first term a3 equals zero!";
            throw new Error("Function name: " + fcnName + ", " + myErrorMsg);
        }
        a = a2 / a3;
        b = a1 / a3;
        c = A0 / a3;
        p = ((-1) * Math.pow(a, 2) / 3 + b) / 3;
        if ((9 * a * b - 2 * Math.pow(a, 3) - 27 * c) === 0) {
            myErrorMsg = "getCubicRoot error: (9 * A * B - 2 * A ** 3 - 27 * C) equals zero!";
            throw new Error("Function name: " + fcnName + ", " + myErrorMsg);
        }
        Q = (9 * a * b - 2 * Math.pow(a, 3) - 27 * c) / 54;
        Disc = Math.pow(Q, 2) + Math.pow(p, 3);
        if (Disc > 0) {
            h = Q + Math.pow(Disc, (1 / 2));
            if (Math.pow(Math.abs(h), (1 / 3)) === 0) {
                myErrorMsg = "getCubicRoot error: (Abs(h)) ** (1 / 3) equals zero!";
                throw new Error("Function name: " + fcnName + ", " + myErrorMsg);
            }
            Y = Math.pow((Math.abs(h)), (1 / 3));
            if (h < 0) {
                Y = -Y;
            }
            Z = Y - p / Y - a / 3;
        }
        else {
            Theta = Math.atan(Math.pow((-Disc), (1 / 2)) / Q);
            c1 = Math.cos(Theta / 3);
            if (Q < 0) {
                S1 = Math.sin(Theta / 3);
                c1 = (c1 - S1 * Math.pow(3, (1 / 2))) / 2;
            }
            z1 = 2 * Math.pow((-p), (1 / 2)) * c1 - a / 3;
            m = a + z1;
            R = Math.pow((Math.pow(m, 2) - 4 * (b + m * z1)), (1 / 2));
            z2 = (-m + R) / 2;
            z3 = (-m - R) / 2;
            Z = z1;
            if (rootType === "largest") {
                if (z2 > Z && z2 > 0) {
                    Z = z2;
                }
                if (z3 > Z && z3 > 0) {
                    Z = z3;
                }
            }
            if (rootType === "smallest") {
                if (z2 < Z && z2 > 0) {
                    Z = z2;
                }
                if (z3 < Z && z3 > 0) {
                    Z = z3;
                }
            }
        }
        return Z;
    }
    catch (myErrorHandler) {
        return -500;
    }
}
exports.getCubicRoot = getCubicRoot;
//# sourceMappingURL=findCubicRoots.js.map